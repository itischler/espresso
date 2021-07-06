import sympy as sp
import pystencils as ps
import numpy as np
from pystencils.rng import random_symbol
from pystencils.field import Field
from pystencils_walberla import CodeGeneration, generate_sweep

#from lbmpy_walberla import generate_lattice_model, generate_boundary, generate_lb_pack_info

def grad(f):
    return sp.Matrix([ps.fd.diff(f, i) for i in range(3)])
def diffusion_equation(c_field, pot_field, D, z):
  return - D * grad(c_field) - D * z * c_field.center * grad(pot_field)
def create_advection_diffusion_method(c_field, v_field, pot_field, j_field, D, z):
    flux_eq = diffusion_equation(c_field, pot_field, D, z)
    fvm_eq = ps.fd.FVM1stOrder(c_field, flux=flux_eq)
    vof_adv = ps.fd.VOF(j_field, v_field, c_field)
    for adv, div in zip(vof_adv, fvm_eq.discrete_flux(j_field)):
        assert adv.lhs == div.lhs
        flux = ps.Assignment(adv.lhs, (adv.rhs + div.rhs))
    return ps.AssignmentCollection([flux])

def stencil_factor(stencil):
    factor = np.sqrt(1/(1+2*np.sqrt(2)))
    if(stencil == 27):
        factor = np.sqrt(1/(1+2*np.sqrt(2)+ 4.0/3.0 * np.sqrt(3)))
    if stencil == 6:
        factor = 1.0
    return factor

def add_fluctuations(flux, c_field, j_field, D, stencil_factor):
    rng_symbol_gen = random_symbol(flux.subexpressions, dim=3)
    for i in range(len(flux.main_assignments)):
        n = j_field.staggered_stencil[i]
        
        # calculate mean density
        dens = (c_field.neighbor_vector(n) + c_field.center_vector)[0]/2
        # multyply by smoothed haviside function so that fluctuation will not get bigger that the density
        dens *= sp.Max(0,sp.Min(1.0,c_field.neighbor_vector(n)[0]) * sp.Min(1.0,c_field.center_vector[0]))
        
        # lenght of the vector
        length = sp.sqrt(len(j_field.staggered_stencil[i]))
        
        # amplitude of the random fluctuations
        fluct = sp.sqrt(2*dens*D) * sp.sqrt(1/length) * stencil_factor
        # add fluctuations
        fluct *= 2 * (next(rng_symbol_gen)-0.5) * sp.sqrt(3)
        
        flux.main_assignments[i] = ps.Assignment(flux.main_assignments[i].lhs, flux.main_assignments[i].rhs + fluct)
    return flux
def add_ghostlayer_folding(flux, L):
    ''' Add the folding to the flux, so that the random numbers persist through the ghostlayers.'''
    fold = {ps.astnodes.LoopOverCoordinate.get_loop_counter_symbol(i):
            ps.astnodes.LoopOverCoordinate.get_loop_counter_symbol(i) % L[i] for i in range(3)}
    flux.subs(fold)
    return flux

def create_reaction_method(c_fields, r_flux_fields, r_rate_const, r_coefs, r_orders):
    reaction = r_rate_const
    for i in range(len(c_fields)):
        reaction *= sp.Pow(c_fields[i].center, r_orders[i])
    r_flux = []
    for i in range(len(c_fields)):
        r_flux.append(ps.Assignment(r_flux_fields[i].center, reaction * r_coefs[i]))
    return r_flux




with CodeGeneration() as ctx:
    stencil = 19
    kT = sp.symbols("kT")
    force_field = ps.fields("force(3): [3D]", layout='fzyx')
    c_field = ps.fields("c : float32[3D]", layout='fzyx')
    pot_field = ps.fields("Phi : float32[3D]", layout='fzyx')
    j_field = ps.fields(f"j({stencil//2}) : float32[3D]", layout='fzyx',
                              field_type=ps.FieldType.STAGGERED_FLUX)
    charge_field = ps.fields("q : float32[3D]", layout='fzyx')
    v_field = ps.fields("v(3) : float32[3D]", layout='fzyx')
    r_flux_field = ps.fields("r : float32[3D]", layout='fzyx')
    
    kT = sp.Symbol("kT")
    dt = sp.Symbol("dt")
    r_rate_const = sp.Symbol("gamma")
    D = sp.Symbol("D")
    z = sp.Symbol("z")
    L = sp.Matrix([ps.TypedSymbol(f'L_{i}', np.int) for i in range(3)])

    cpu_vectorize_info = {
        "instruction_set": "avx",
        "assume_inner_stride_one": True,
        "assume_aligned": True,
        "assume_sufficient_line_padding": True}


    # generate EK methods
    adv_diff_method = create_advection_diffusion_method(c_field, v_field, pot_field, j_field, D, z)
    
    generate_sweep(ctx,
                   'AdvectionDiffusion',
                   adv_diff_method,
                   ghost_layers_to_include=1,
                   staggered=True)

    adv_diff_fluct_method = add_fluctuations(adv_diff_method, c_field, j_field, D, stencil_factor(stencil))
    adv_diff_fluct_method = add_ghostlayer_folding(adv_diff_fluct_method, L)

    generate_sweep(ctx,
                   'AdvectionDiffusionFluctuation',
                   adv_diff_fluct_method,
                   ghost_layers_to_include=1,
                   staggered=True)

    fvm_eq = ps.fd.FVM1stOrder(c_field, flux=diffusion_equation(c_field, pot_field, D, z))
    continuity_method = fvm_eq.discrete_continuity(j_field)

    generate_sweep(ctx,
                   'ContinuityEquation',
                   continuity_method,
                   ghost_layers_to_include=1
                   staggered=False)

    gather_charges = ps.Assignment(charge_field.center, charge_field.center + z * c_field.center)

    generate_sweep(ctx,
                   'GatherCharges',
                   gather_charges,
                   staggered=False)

    # TODO FFT

    reaction_continuity = ps.Assignment(c_field.center, c_field.center + r_flux_field.center)

    generate_sweep(ctx,
                   'ReactionContinuity',
                   reaction_continuity,
                   staggered=False)

    c_fields = []
    r_flux_fields = []
    for i in range(1,6):
        c_fields.append(ps.fields(f"c_{i} : float32[3D]", layout='fzyx'))
        r_flux_fields.append(ps.fields(f"r_{i} : float32[3D]", layout='fzyx'))
        r_coefs = sp.Matrix([ps.TypedSymbol(f'n_{j}', np.float64) for j in range(i)])
        r_orders = sp.Matrix([ps.TypedSymbol(f'O_{j}', np.float64) for j in range(i)])

        reaction_flux = create_reaction_method(c_fields, r_flux_fields, r_rate_const, r_coefs, r_orders)

        generate_sweep(ctx,
                       f'ReactionFluxSpecies{i}',
                       reaction_flux,
                       staggered=False)
























