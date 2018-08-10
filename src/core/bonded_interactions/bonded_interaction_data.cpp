#include "bonded_interaction_data.hpp"
#include "bonded_interactions/thermalized_bond.hpp"

std::vector<Bonded_ia_parameters> bonded_ia_params;

/** calculate the maximal cutoff of bonded interactions, required to
    determine the cell size for communication. */
void recalc_maximal_cutoff_bonded() {
  int i;
  double max_cut_tmp;

  max_cut_bonded = 0.0;

  for (i = 0; i < bonded_ia_params.size(); i++) {
    switch (bonded_ia_params[i].type) {
    case BONDED_IA_FENE:
      max_cut_tmp =
          bonded_ia_params[i].p.fene.r0 + bonded_ia_params[i].p.fene.drmax;
      if (max_cut_bonded < max_cut_tmp)
        max_cut_bonded = max_cut_tmp;
      break;
    case BONDED_IA_HARMONIC:
      if ((bonded_ia_params[i].p.harmonic.r_cut > 0) &&
          (max_cut_bonded < bonded_ia_params[i].p.harmonic.r_cut))
        max_cut_bonded = bonded_ia_params[i].p.harmonic.r_cut;
      break;
    case BONDED_IA_THERMALIZED_DIST:
      if ((bonded_ia_params[i].p.thermalized_bond.r_cut > 0) && 
	  (max_cut_bonded < bonded_ia_params[i].p.thermalized_bond.r_cut))
    	max_cut_bonded = bonded_ia_params[i].p.thermalized_bond.r_cut;
      break;
    case BONDED_IA_RIGID_BOND:
      if (max_cut_bonded < sqrt(bonded_ia_params[i].p.rigid_bond.d2))
        max_cut_bonded = sqrt(bonded_ia_params[i].p.rigid_bond.d2);
      break;
#ifdef TABULATED
    case BONDED_IA_TABULATED:
      if (bonded_ia_params[i].p.tab.type == TAB_BOND_LENGTH &&
          max_cut_bonded < bonded_ia_params[i].p.tab.pot->cutoff())
        max_cut_bonded = bonded_ia_params[i].p.tab.pot->cutoff();
      break;
#endif
#ifdef IMMERSED_BOUNDARY
    case BONDED_IA_IBM_TRIEL:
      if (max_cut_bonded < bonded_ia_params[i].p.ibm_triel.maxDist)
        max_cut_bonded = bonded_ia_params[i].p.ibm_triel.maxDist;
      break;
#endif
    default:
      break;
    }
  }

  /* Bond angle and dihedral potentials do not contain a cutoff
     intrinsically. The cutoff for these potentials depends on the
     bond length potentials. For bond angle potentials nothing has to
     be done (it is assumed, that particles participating in a bond
     angle or dihedral potential are bound to each other by some bond
     length potential (FENE, Harmonic or tabulated)). For dihedral
     potentials (both normal and tabulated ones) it follows, that the
     cutoff is TWO TIMES the maximal cutoff! That's what the following
     lines assure. */
  max_cut_tmp = 2.0 * max_cut_bonded;
  for (i = 0; i < bonded_ia_params.size(); i++) {
    switch (bonded_ia_params[i].type) {
    case BONDED_IA_DIHEDRAL:
      max_cut_bonded = max_cut_tmp;
      break;
#ifdef TABULATED
    case BONDED_IA_TABULATED:
      if (bonded_ia_params[i].p.tab.type == TAB_BOND_DIHEDRAL)
        max_cut_bonded = max_cut_tmp;
      break;
#endif
    default:
      break;
    }
  }
}


void make_bond_type_exist(int type) {
  int i, ns = type + 1;
  const auto old_size = bonded_ia_params.size();
  if (ns <= bonded_ia_params.size()) {
    return;
  }
  /* else allocate new memory */
  bonded_ia_params.resize(ns);
  /* set bond types not used as undefined */
  for (i = old_size; i < ns; i++)
    bonded_ia_params[i].type = BONDED_IA_NONE;
}
