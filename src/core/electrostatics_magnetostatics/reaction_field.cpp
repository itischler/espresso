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
 *
 *  Implementation of \ref reaction_field.hpp
 */
#include "reaction_field.hpp"

#ifdef ELECTROSTATICS
#include "communication.hpp"

Reaction_field_params rf_params{};

int rf_set_params(double kappa, double epsilon1, double epsilon2,
                  double r_cut) {
  rf_params.kappa = kappa;
  rf_params.epsilon1 = epsilon1;
  rf_params.epsilon2 = epsilon2;
  rf_params.r_cut = r_cut;
  rf_params.B = (2 * (epsilon1 - epsilon2) * (1 + kappa * r_cut) -
                 epsilon2 * kappa * kappa * r_cut * r_cut) /
                ((epsilon1 + 2 * epsilon2) * (1 + kappa * r_cut) +
                 epsilon2 * kappa * kappa * r_cut * r_cut);
  if (rf_params.epsilon1 < 0.0)
    return -1;

  if (rf_params.epsilon2 < 0.0)
    return -1;

  if (rf_params.r_cut < 0.0)
    return -2;

  mpi_bcast_coulomb_params();

  return 1;
}
#endif
