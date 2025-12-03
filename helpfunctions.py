import itertools
import copy
import numpy as np

def create_reactions_for_step_one(verbose=False):
    """
    Function to create all combinations of isotops with the following shape:

    OOPOO + HOH --> OOPOOO 
    
    OOPOO  = ["R1", "lhs", 16, 16, 16, 16, 0]
    OOPOOO = ["R1", "rhs", 16, 16, 16, 16, 18]

    """
    reactions = [] # [reactant, product, kdirect, kreverse]
    kd, kr = 1, 1 # ! HARDCODED

    # isomer O16 pure
    template = [16, 16, 16, 16]
    lhs_po4_16 = (tuple(template), (16,))
    lhs_po4_18 = (tuple(template), (18,))
    rhs_po4_16 = (tuple(template + [16]),)
    rhs_po4_18 = (tuple(template + [18]),)
    print(lhs_po4_16)
    reactions.append((lhs_po4_16, rhs_po4_16, kd, kr))
    reactions.append((lhs_po4_18, rhs_po4_18, kd, kr))

    # isomer O16 O18 mix
    template = [18, 16, 16, 16] # ! HARDCODED
    comb_po4 = list([list(o) for o in set(itertools.permutations(template))])
    for _lhs_po4 in comb_po4:
        lhs_po4_16 = (tuple(_lhs_po4), (16,))
        lhs_po4_18 = (tuple(_lhs_po4), (18,))
        rhs_po4_16 = (tuple(_lhs_po4 + [16]),)
        rhs_po4_18 = (tuple(_lhs_po4 + [18]),)
        reactions.append((lhs_po4_16, rhs_po4_16, kd, kr))
        reactions.append((lhs_po4_18, rhs_po4_18, kd, kr))

    if verbose: 
        for r in reactions:
            print(r)
        print(len(reactions))

    return reactions


def create_reactions_for_step_two(verbose=False):
    """
    Function to create all combinations of isotops with the following shape:

    OOPOOO --> OOPOO

    OOPOO  = [16, 16, 16, 16, 0]
    OOPOOO = [16, 16, 16, 16, 18]

    """
    
    reactions = [] # [reactant, product, kdirect, kreverse]
    kd, kr = 1, 1 # ! HARDCODED
    
    # isomer O16 pure
    template = [16, 16, 16, 16, 16] # ! HARDCODED 
    lhs_po5 = (tuple(template),)
    _rhs_po4_16 = template.copy()
    _rhs_po4_16[-1] = 0
    ## only one is created because if there are no O16, permutation are
    ## equivalent in terms of rate constants
    rhs_po4_16 = (tuple(_rhs_po4_16),)
    reactions.append((lhs_po5, rhs_po4_16, kd, kr))

    # isomer O16 O18 mix
    for template in [[18, 16, 16, 16, 16], [18, 18, 16, 16, 16]]: # ! HARDCODED
        comb_po5 = list([list(o) for o in set(itertools.permutations(template))])
        kd, kr = 1, 1 # ! HARDCODED
        for lhs_po5 in comb_po5:
           # the oxygen can only be lost from one specific position (arbitrarely [-1])
            _lhs_po5 = lhs_po5.copy()
            rhs_po4 = _lhs_po5
            rhs_po4[-1] = 0
            reactions.append(((tuple(lhs_po5),), (tuple(rhs_po4),), kd, kr))
    
    reactions_worep = list(set(reactions))
    if verbose:
        for r in reactions_worep:
            print(r)
        print(len(reactions_worep))
    return reactions

def convert_to_kinetx_notation(reactions, initial_conc):
    """
    Convert simplified form of isotopic reactions to a format that KiNetX can read
    """
    import scine_kinetx as kx

    # hash compounds to an index
    compounds = {}
    acc = 1
    for r in reactions:
        print("===========loop", r)
        lhs, rhs, _, _ = r
        for hs in [lhs, rhs]:
            for hsi in hs:
                if hsi in compounds.keys():
                    pass
                else:
                    compounds[hsi] = acc
                    acc += 1
    
    print(compounds)
    # create kinetx network object
    network_builder = kx.NetworkBuilder()
    n_compounds = len(compounds.keys()) + 1
    n_reactions = len(reactions) + 1 
    n_channels_per_reaction = 1
    network_builder.reserve(n_compounds, n_reactions, n_channels_per_reaction)
    
    
    concentrations = [0.9, 0.1] + [0 for _ in range(n_compounds-2)] # ! HARDCODED
    

    for i in range(n_compounds):
        network_builder.add_compound(1, str(i))
    for r in reactions:
        _lhs, _rhs, kd, kr = r
        lhs = [compounds[o] for o in _lhs]
        rhs = [compounds[o] for o in _rhs]
        if len(lhs) == 1: lhs = lhs + [0]
        if len(rhs) == 1: rhs = rhs + [0]
        print(_lhs, _rhs)
        print(lhs, rhs)
        network_builder.add_reaction([kd], [kr], [tuple(lhs)], [tuple(rhs)])
    network = network_builder.generate()
    solver = kx.Integrator.explicit_euler
    concentration_data, r_flux, rf_flux, rb_flux = kx.integrate(network, np.asarray(
        concentrations), 0.0, 1e+0, solver, 5000, 10000, 1e-10)

    return



