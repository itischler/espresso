/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
 * Header file for lb.cpp
 *
 * This is the header file for the Lattice Boltzmann implementation in lb.cpp
 */

#ifndef LB_H
#define LB_H

#include "config.hpp"

#include "lattice.hpp"

void mpi_set_lb_coupling_counter(int high, int low);

#ifdef LB

#include "errorhandling.hpp"

#include "utils.hpp"

#include <utils/Span.hpp>

/** \name Parameter fields for Lattice Boltzmann
 * The numbers are referenced in \ref mpi_bcast_lb_params
 * to determine what actions have to take place upon change
 * of the respective parameter. */
/*@{*/
#define LBPAR_DENSITY 0   /**< fluid density */
#define LBPAR_VISCOSITY 1 /**< fluid kinematic viscosity */
#define LBPAR_AGRID 2     /**< grid constant for fluid lattice */
#define LBPAR_TAU 3       /**< time step for fluid propagation */
#define LBPAR_FRICTION                                                         \
  4 /**< friction coefficient for viscous coupling between particles and fluid \
     */
#define LBPAR_EXTFORCE 5 /**< external force density acting on the fluid */
#define LBPAR_BULKVISC 6 /**< fluid bulk viscosity */

/** Note these are used for binary logic so should be powers of 2 */
#define LB_COUPLE_NULL 1
#define LB_COUPLE_TWO_POINT 2
#define LB_COUPLE_THREE_POINT 4

#ifdef ADDITIONAL_CHECKS
#endif
/*@}*/
/** Some general remarks:
 * This file implements the LB D3Q19 method to Espresso. The LB_Model
 * construction is preserved for historical reasons and might be removed
 * soon. It is constructed as a multi-relaxation time LB, thus all populations
 * are converted to modes, then collision is performed and transfered back
 * to population space, where the streaming is performed.
 *
 * For performance reasons it is clever to do streaming and collision at the
 * same time
 * because every fluid node has to be read and written only once. This increases
 * mainly cache efficiency.
 * Two alternatives are implemented: stream_collide and collide_stream.
 *
 * The hydrodynamic fields, corresponding to density, velocity and stress, are
 * stored in LB_FluidNodes in the array lbfields, the populations in lbfluid
 * which is constructed as 2 x (Nx x Ny x Nz) x 19 array.
 */

/** Description of the LB Model in terms of the unit vectors of the
 *  velocity sub-lattice and the corresponding coefficients
 *  of the pseudo-equilibrium distribution */
template <size_t N_vel = 19> struct LB_Model {
  /** number of velocities */
  static const constexpr int n_veloc = static_cast<int>(N_vel);

  /** unit vectors of the velocity sublattice */
  std::array<std::array<double, 3>, N_vel> c;

  /** coefficients in the pseudo-equilibrium distribution */
  std::array<std::array<double, 4>, N_vel> coeff;

  /** weights in the functional for the equilibrium distribution */
  std::array<double, N_vel> w;

  /** basis of moment space */
  std::array<std::array<double, N_vel>, N_vel> e_ki;

  /** normalization factors for the moment basis */
  std::array<double, N_vel> w_k;

  /** speed of sound squared */
  double c_sound_sq;
};

/** Data structure for fluid on a local lattice site */
struct LB_FluidNode {
#ifdef LB_BOUNDARIES
  /** flag indicating whether this site belongs to a boundary */
  int boundary;
  Vector3d slip_velocity = {};
#endif // LB_BOUNDARIES

  /** local force density */
  Vector3d force_density;
#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
  // For particle update, we need the force on the nodes in LBM
  // Yet, Espresso resets the force immediately after the LBM update
  // Therefore we save it here
  Vector3d force_density_buf;
#endif
};

/** Data structure holding the parameters for the Lattice Boltzmann system. */
typedef struct {
  /** number density (LJ units) */
  double rho;

  /** kinematic viscosity (LJ units) */
  double viscosity;

  /** bulk viscosity (LJ units) */
  double bulk_viscosity;

  /** lattice spacing (LJ units) */
  double agrid;

  /** time step for fluid propagation (LJ units)
   *  Note: Has to be larger than MD time step! */
  double tau;

  /** friction coefficient for viscous coupling (LJ units) */
  double friction;

  /** external force density applied to the fluid at each lattice site (MD
   * units) */
  Vector3d ext_force_density;
  double rho_lb_units;
  /** relaxation of the odd kinetic modes */
  double gamma_odd;
  /** relaxation of the even kinetic modes */
  double gamma_even;
  /** relaxation rate of shear modes */
  double gamma_shear;
  /** relaxation rate of bulk modes */
  double gamma_bulk;

  /** Flag determining whether lbpar.gamma_shear, gamma_odd, and gamma_even are
   * calculated
   *  from lbpar.gamma_shear in such a way to yield a TRT LB with minimized slip
   * at
   *  bounce-back boundaries */
  bool is_TRT;

  /** \name Derived parameters */
  /** Flag indicating whether fluctuations are present. */
  int fluct;
  /** amplitudes of the fluctuations of the modes */
  Vector<19, double> phi;
} LB_Parameters;

/** The DnQm model to be used. */
extern LB_Model<> lbmodel;

/** Struct holding the Lattice Boltzmann parameters */
extern LB_Parameters lbpar;

/** The underlying lattice */
extern Lattice lblattice;

/** Pointer to the velocity populations of the fluid.
 * lbfluid contains pre-collision populations, lbfluid_post
 * contains post-collision populations*/
using LB_Fluid = std::array<Utils::Span<double>, 19>;
extern LB_Fluid lbfluid;

/** Pointer to the hydrodynamic fields of the fluid */
extern std::vector<LB_FluidNode> lbfields;

/** Switch indicating momentum exchange between particles and fluid */
extern int transfer_momentum;

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Updates the Lattice Boltzmann system for one time step.
 * This function performs the collision step and the streaming step.
 * If external force densities are present, they are applied prior to the
 * collisions. If boundaries are present, it also applies the boundary
 * conditions.
 */
void lattice_boltzmann_update();

/** Performs a full initialization of
 *  the Lattice Boltzmann system. All derived parameters
 *  and the fluid are reset to their default values. */
void lb_init();

/** (Re-)initializes the derived parameters
 *  for the Lattice Boltzmann system.
 *  The current state of the fluid is unchanged. */
void lb_reinit_parameters();

/** (Re-)initializes the fluid. */
void lb_reinit_fluid();

/**
 * @brief Event handler for integration start.
 */
void lb_on_integration_start();

/** Calculates the equilibrium distributions.
    @param index Index of the local site
    @param rho local fluid density
    @param j local fluid speed
    @param pi local fluid pressure
*/
void lb_calc_n_from_rho_j_pi(const Lattice::index_t index, const double rho,
                             const std::array<double, 3> &j,
                             const std::array<double, 6> &pi);

/** Calculates the coupling of MD particles to the LB fluid.
 * This function  is called from \ref force_calc. The force is added
 * to the particle force and the corresponding momentum exchange is
 * applied to the fluid.
 * Note that this function changes the state of the fluid!
 */
void calc_particle_lattice_ia();

uint64_t lb_coupling_rng_state();
void lb_coupling_set_rng_state(uint64_t counter);

/** calculates the fluid velocity at a given position of the
 * lattice. Note that it can lead to undefined behaviour if the
 * position is not within the local lattice. */
Vector3d lb_lbfluid_get_interpolated_velocity(const Vector3d &p);

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
Vector3d lb_lbfluid_get_interpolated_force(const Vector3d &p);
#endif

void lb_calc_local_fields(Lattice::index_t index, double *rho, double *j,
                          double *pi);

/** Calculation of hydrodynamic modes.
 *
 *  @param index number of the node to calculate the modes for
 *  @param mode output pointer to a double[19]
 */
std::array<double, 19> lb_calc_modes(Lattice::index_t index);

/** Calculate the local fluid density.
 * The calculation is implemented explicitly for the special case of D3Q19.
 * @param index the local lattice site (Input).
 * @param rho   local fluid density
 */
inline void lb_calc_local_rho(Lattice::index_t index, double *rho) {
  // unit conversion: mass density
  if (!(lattice_switch & LATTICE_LB)) {
    runtimeErrorMsg() << "Error in lb_calc_local_rho in " << __FILE__
                      << __LINE__ << ": CPU LB not switched on.";
    *rho = 0;
    return;
  }

  double avg_rho = lbpar.rho * lbpar.agrid * lbpar.agrid * lbpar.agrid;

  *rho = avg_rho + lbfluid[0][index] + lbfluid[1][index] + lbfluid[2][index] +
         lbfluid[3][index] + lbfluid[4][index] + lbfluid[5][index] +
         lbfluid[6][index] + lbfluid[7][index] + lbfluid[8][index] +
         lbfluid[9][index] + lbfluid[10][index] + lbfluid[11][index] +
         lbfluid[12][index] + lbfluid[13][index] + lbfluid[14][index] +
         lbfluid[15][index] + lbfluid[16][index] + lbfluid[17][index] +
         lbfluid[18][index];
}

/** Calculate the local fluid momentum.
 * The calculation is implemented explicitly for the special case of D3Q19.
 * @param index The local lattice site (Input).
 * @param j local fluid speed
 */
inline void lb_calc_local_j(Lattice::index_t index, double *j) {
  if (!(lattice_switch & LATTICE_LB)) {
    runtimeErrorMsg() << "Error in lb_calc_local_j in " << __FILE__ << __LINE__
                      << ": CPU LB not switched on.";
    j[0] = j[1] = j[2] = 0;
    return;
  }

  j[0] = lbfluid[1][index] - lbfluid[2][index] + lbfluid[7][index] -
         lbfluid[8][index] + lbfluid[9][index] - lbfluid[10][index] +
         lbfluid[11][index] - lbfluid[12][index] + lbfluid[13][index] -
         lbfluid[14][index];
  j[1] = lbfluid[3][index] - lbfluid[4][index] + lbfluid[7][index] -
         lbfluid[8][index] - lbfluid[9][index] + lbfluid[10][index] +
         lbfluid[15][index] - lbfluid[16][index] + lbfluid[17][index] -
         lbfluid[18][index];
  j[2] = lbfluid[5][index] - lbfluid[6][index] + lbfluid[11][index] -
         lbfluid[12][index] - lbfluid[13][index] + lbfluid[14][index] +
         lbfluid[15][index] - lbfluid[16][index] - lbfluid[17][index] +
         lbfluid[18][index];
}

/** Calculate the local fluid fields.
 * The calculation is implemented explicitly for the special case of D3Q19.
 *
 * @param index   Index of the local lattice site.
 * @param rho     local fluid density
 * @param j       local fluid speed
 * @param pi      local fluid pressure
 */
void lb_calc_local_fields(Lattice::index_t index, double *rho, double *j,
                          double *pi);

#ifdef LB_BOUNDARIES
inline void lb_local_fields_get_boundary_flag(Lattice::index_t index,
                                              int *boundary) {

  if (!(lattice_switch & LATTICE_LB)) {
    runtimeErrorMsg() << "Error in lb_local_fields_get_boundary_flag in "
                      << __FILE__ << __LINE__ << ": CPU LB not switched on.";
    *boundary = 0;
    return;
  }

  *boundary = lbfields[index].boundary;
}
#endif

inline void lb_get_populations(Lattice::index_t index, double *pop) {
  for (int i = 0; i < lbmodel.n_veloc; ++i) {
    pop[i] = lbfluid[i][index] + lbmodel.coeff[i][0] * lbpar.rho;
  }
}

inline void lb_set_populations(Lattice::index_t index, double *pop) {
  for (int i = 0; i < lbmodel.n_veloc; ++i) {
    lbfluid[i][index] = pop[i] - lbmodel.coeff[i][0] * lbpar.rho;
  }
}
#endif

#if defined(LB) || defined(LB_GPU)
/* A C level interface to the LB fluid */
int lb_lbfluid_set_density(double *p_dens);
int lb_lbfluid_get_density(double *p_dens);
int lb_lbfluid_set_visc(double *p_visc);
int lb_lbfluid_set_bulk_visc(double *p_bulk_visc);
int lb_lbfluid_set_gamma_odd(double *p_gamma_odd);
int lb_lbfluid_set_gamma_even(double *p_gamma_even);
int lb_lbfluid_set_friction(double *p_friction);
int lb_lbfluid_set_couple_flag(int couple_flag);
int lb_lbfluid_set_agrid(double p_agrid);
int lb_lbfluid_set_ext_force_density(int component, double p_fx, double p_fy,
                                     double p_fz);
int lb_lbfluid_set_tau(double p_tau);
int lb_lbfluid_set_remove_momentum(void);
int lb_lbfluid_get_agrid(double *p_agrid);
int lb_lbfluid_get_tau(double *p_tau);
int lb_lbfluid_get_visc(double *p_visc);
int lb_lbfluid_get_bulk_visc(double *p_bulk_visc);
int lb_lbfluid_get_friction(double *p_friction);
int lb_lbfluid_get_couple_flag(int *couple_flag);
int lb_lbfluid_get_ext_force_density(double *p_f);
#ifdef SHANCHEN
int lb_lbfluid_set_shanchen_coupling(double *p_coupling);
int lb_lbfluid_set_mobility(double *p_mobility);
#endif
int lb_set_lattice_switch(int local_lattice_switch);

/* IO routines */
int lb_lbfluid_print_vtk_boundary(char *filename);
int lb_lbfluid_print_vtk_velocity(char *filename,
                                  std::vector<int> = {-1, -1, -1},
                                  std::vector<int> = {-1, -1, -1});
int lb_lbfluid_print_vtk_density(char **filename);
int lb_lbfluid_print_boundary(char *filename);
int lb_lbfluid_print_velocity(char *filename);

int lb_lbfluid_save_checkpoint(char *filename, int binary);
int lb_lbfluid_load_checkpoint(char *filename, int binary);

bool lb_lbnode_is_index_valid(const Vector3i &ind);
int lb_lbnode_get_rho(const Vector3i &ind, double *p_rho);
int lb_lbnode_get_u(const Vector3i &ind, double *u);
int lb_lbnode_get_pi(const Vector3i &ind, double *pi);
int lb_lbnode_get_pi_neq(const Vector3i &ind, double *pi_neq);
int lb_lbnode_get_boundary(const Vector3i &ind, int *p_boundary);
int lb_lbnode_get_pop(const Vector3i &ind, double *pop);

int lb_lbnode_set_rho(const Vector3i &ind, double *rho);
int lb_lbnode_set_u(const Vector3i &ind, double *u);
int lb_lbnode_set_pop(const Vector3i &ind, double *pop);

/** calculates the fluid velocity at a given position of the
 * lattice. Note that it can lead to undefined behaviour if the
 * position is not within the local lattice. This version of the function
 * can be called without the position needing to be on the local processor */
int lb_lbfluid_get_interpolated_velocity_global(Vector3d &p, double *v);

#endif

#endif /* _LB_H */
/*@}*/