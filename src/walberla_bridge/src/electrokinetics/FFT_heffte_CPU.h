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

#pragma once

#include <blockforest/communication/UniformBufferedScheme.h>
#include <domain_decomposition/BlockDataID.h>
#include <field/AddToStorage.h>
#include <field/GhostLayerField.h>
#include <field/communication/PackInfo.h>

#include <stencil/D3Q27.h>

#include <walberla_bridge/LatticeWalberla.hpp>

#include "../src/electrokinetics/generated_kernels/EK_FieldAccessors_double_precision.h"
#include "../src/electrokinetics/generated_kernels/EK_FieldAccessors_single_precision.h"

#include <complex>
#include <cstddef>
#include <memory>
#include <type_traits>

namespace walberla {

template <typename FloatType> class FFT_heffte_CPU {
private:
  struct heffte_container;
  template <typename T> FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }

  using ComplexType = std::complex<FloatType>;
  using PotentialField = GhostLayerField<FloatType, 1>;

  std::shared_ptr<LatticeWalberla> m_lattice;
  double m_permittivity;

  walberla::BlockDataID m_potential_field_with_ghosts_id;
  std::vector<FloatType> m_greens;
  std::vector<FloatType> m_potential;
  std::vector<ComplexType> m_potential_fourier;

  std::shared_ptr<blockforest::StructuredBlockForest> m_blocks;
  std::shared_ptr<heffte_container> heffte;

  using FullCommunicator =
      blockforest::communication::UniformBufferedScheme<stencil::D3Q27>;
  std::shared_ptr<FullCommunicator> m_full_communication;

public:
  FFT_heffte_CPU(std::shared_ptr<LatticeWalberla> lattice, double permittivity);
  ~FFT_heffte_CPU() = default;

  void reset_charge_field();

  void add_charge_to_field(std::size_t id, double valency,
                           bool is_double_precision);

  [[nodiscard]] std::size_t get_potential_field_id() const noexcept {
    return static_cast<std::size_t>(m_potential_field_with_ghosts_id);
  }

  void solve();

  void set_permittivity(double permittivity) noexcept {
    m_permittivity = permittivity;
  }

  [[nodiscard]] double get_permittivity() const noexcept {
    return m_permittivity;
  }

  [[nodiscard]] auto const &get_lattice() const noexcept { return *m_lattice; }

private:
  void add_fields(PotentialField *field_out,
                  GhostLayerField<FloatType, 1> *field_add, FloatType factor);
  void ghost_communication() { (*m_full_communication)(); }
};

} // namespace walberla
