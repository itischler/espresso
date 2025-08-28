/*
 * Copyright (C) 2025 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "FFT_heffte_CPU.h"

#include <heffte.h>
#include <heffte_backends.h>
#include <heffte_geometry.h>

#include <utils/Vector.hpp>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>

namespace walberla {

template <typename FloatType> struct FFT_heffte_CPU<FloatType>::heffte_container {
  std::shared_ptr<heffte::box3d<>> m_box_in;
  std::shared_ptr<heffte::box3d<>> m_box_out;
  std::shared_ptr<heffte::fft3d<heffte::backend::fftw>> m_fft;
  std::shared_ptr<
      heffte::fft3d<heffte::backend::fftw>::buffer_container<ComplexType>>
      m_buffer;
};

template <typename FloatType>
FloatType greens_function(int x, int y, int z, Utils::Vector<int, 3> &dim) {
  if (x == 0u && y == 0u && z == 0u)
    return 0.;
  return -0.5 /
          (std::cos(2. * std::numbers::pi * FloatType(x) / FloatType(dim[0])) +
           std::cos(2. * std::numbers::pi * FloatType(y) / FloatType(dim[1])) +
           std::cos(2. * std::numbers::pi * FloatType(z) / FloatType(dim[2])) -
           3.) /
          FloatType(dim[0] * dim[1] * dim[2]);
}

template <typename T, std::size_t N>
auto to_array(Utils::Vector<T, N> const &vec) {
  std::array<T, N> res{};
  std::copy(vec.begin(), vec.end(), res.begin());
  return res;
}

inline int pos_to_linear_index(int x, int y, int z, Utils::Vector<int, 3> dim) {
  return (z * dim[1] + y) * dim[0] + x;
}

template <typename FloatType>
FFT_heffte_CPU<FloatType>::FFT_heffte_CPU(std::shared_ptr<LatticeWalberla> lattice,
                                          double permittivity)
    : m_lattice(std::move(lattice)), m_permittivity(permittivity){
  m_blocks = get_lattice().get_blocks();

  m_potential_field_with_ghosts_id = field::addToStorage<PotentialField>(
        get_lattice().get_blocks(), "potential field with ghosts", 0.0, field::fzyx,
        get_lattice().get_ghost_layers());

  heffte = std::make_shared<heffte_container>();
  auto dim = get_lattice().get_grid_dimensions();
  auto offset_vec = Utils::Vector3i({1, 1, 1});
  auto order = Utils::Vector3i({0, 1, 2});
  heffte->m_box_in = std::make_shared<heffte::box3d<>>(
      to_array(Utils::Vector3i({0, 0, 0})), to_array(dim - offset_vec),
      to_array(order));
  heffte->m_box_out = std::make_shared<heffte::box3d<>>(
      to_array(Utils::Vector3i({0, 0, 0})), to_array(dim - offset_vec),
      to_array(order));
  heffte->m_fft = std::make_shared<heffte::fft3d<heffte::backend::fftw>>(
      *(heffte->m_box_in), *(heffte->m_box_out), MPI_COMM_WORLD);
  heffte->m_buffer = std::make_shared<
      heffte::fft3d<heffte::backend::fftw>::buffer_container<ComplexType>>(
      heffte->m_fft->size_workspace());

  m_potential = std::vector<FloatType>(heffte->m_fft->size_inbox());
  m_greens = std::vector<FloatType>(heffte->m_fft->size_outbox());
  m_potential_fourier = std::vector<ComplexType>(heffte->m_fft->size_outbox());
  
  for (int x = 0; x < dim[0]; x++){
    for (int y = 0; y < dim[1]; y++){
      for (int z = 0; z < dim[2]; z++){
        m_greens[pos_to_linear_index(x,y,z,dim)] = greens_function<FloatType>(x, y, z, dim);
      }
    }
  }

  m_full_communication =
        std::make_shared<FullCommunicator>(get_lattice().get_blocks());
  m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PotentialField>>(
            m_potential_field_with_ghosts_id));
  reset_charge_field();
}

template <typename FloatType> void FFT_heffte_CPU<FloatType>::reset_charge_field() {
  auto dim = get_lattice().get_grid_dimensions();
  for (int x = 0; x < dim[0]; x++){
    for (int y = 0; y < dim[1]; y++){
      for (int z = 0; z < dim[2]; z++){
        m_potential[pos_to_linear_index(x,y,z,dim)] = FloatType(0.0);
      }
    }
  }
}

template <typename FloatType>
void FFT_heffte_CPU<FloatType>::add_charge_to_field(std::size_t id, double valency,
                                              bool is_double_precision) {
  auto dim = get_lattice().get_grid_dimensions();
  auto const factor = FloatType_c(valency) / FloatType_c(get_permittivity());
  const auto density_id = walberla::BlockDataID(id);
  for (auto &block : *get_lattice().get_blocks()) {
    auto density_field =
        block.template getData<PotentialField>(density_id);
    for (int x = 0; x < dim[0]; x++){
      for (int y = 0; y < dim[1]; y++){
        for (int z = 0; z < dim[2]; z++){
          m_potential[pos_to_linear_index(x,y,z,dim)] += factor * density_field->get(x,y,z);
        }
      }
    }
  }
}

template <typename FloatType> void FFT_heffte_CPU<FloatType>::solve() {
  auto dim = get_lattice().get_grid_dimensions();
  for (auto &block : *get_lattice().get_blocks()) {
    auto potential_with_ghosts =
        block.template getData<PotentialField>(m_potential_field_with_ghosts_id);
    heffte->m_fft->forward(m_potential.data(), m_potential_fourier.data(),
                           heffte->m_buffer->data());
    for (int i = 0; i < m_potential_fourier.size(); i++){
      m_potential_fourier[i] *= m_greens[i];
    }
    heffte->m_fft->backward(m_potential_fourier.data(), m_potential.data(),
                            heffte->m_buffer->data());
    for (int x = 0; x < dim[0]; x++){
      for (int y = 0; y < dim[1]; y++){
        for (int z = 0; z < dim[2]; z++){
          potential_with_ghosts->get(x,y,z) = m_potential[pos_to_linear_index(x,y,z,dim)];
        }
      }
    }
    ghost_communication();
  }
}

template class FFT_heffte_CPU<float>;
template class FFT_heffte_CPU<double>;

} // namespace walberla
