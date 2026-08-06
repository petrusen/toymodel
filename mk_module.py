import itertools
import copy
import numpy as np
import math

def _calculate_rate_constant(lhs, rhs, dkie):
    """
    Heuristic based rate constants. Reverse rate constant is considered 100
    times slower than direct rate constant.

    Pair subindex constants are in forward direction, otherwise they are in 
    backward direction.

    """
    
    k1p_18 = dkie["k1_16"] / dkie["kie1p"]
    k1s_18 = dkie["k1_16"] / dkie["kie1s"]

    k2p_18 = dkie["k2_16"] / dkie["kie2p"]
    k2s_18 = dkie["k2_16"] / dkie["kie2s"]

    k3p_18 = dkie["k3_16"] / dkie["kie3p"]
    k3s_18 = dkie["k3_16"] / dkie["kie3s"]

    k4p_18 = dkie["k4_16"] / dkie["kie4p"]
    k4s_18 = dkie["k4_16"] / dkie["kie4s"]

    penalty_k2 = dkie["penalty_k2"]
    penalty_k4 = dkie["penalty_k4"]


    if len(lhs) == 2:        # reaction 1
        r1, h2o = lhs        # organophosphate and water
        tsi = rhs[0]            # 'transition state' intermediate
        kf, kb = None, None
        assert len(h2o) == 1 # make sure its the h2o
        if 18 in h2o:
            kf = k1p_18      # kf depends only on h2o
            if tsi[-1] == 18: 
                kb = k2p_18 * penalty_k2
            elif tsi[-1] == 16:
                kb = k2s_18 * penalty_k2
        elif 16 in h2o: 
            if 18 in r1: 
                kf = k1s_18
                if tsi[-1] == 18:
                    kb = k2p_18 * penalty_k2
                elif tsi[-1] == 16:
                    kb = k2p_18 * penalty_k2
            else:
                kf = dkie["k1_16"]
                kb = dkie["k2_16"] * penalty_k2
            
    elif len(lhs) == 1: # reaction 2
        tsi = lhs[0]     # 'transition state' intermediate
        p1, ro = rhs  # phosphate and alcohol
        if 18 in ro:
            kf = k3p_18
            kb = k4p_18 * penalty_k4
        else: 
            if 18 in p1:
                kf = k3s_18
                kb = k4s_18 * penalty_k4
            else:
                kf = dkie["k3_16"]
                kb = dkie["k4_16"] * penalty_k4
    assert kf != None
    assert kb != None

    return kf, kb

def get_water_initial_concentrations(delta_O, initial_C=55.5, Rstd=0.002004):
    """
    Get the concentration of H2O(18) and H2O(16)
    """
    r0_h2o = (delta_O / 1000 + 1) * Rstd
    c16_0 = initial_C / (1 + r0_h2o)
    c18_0 = c16_0 * r0_h2o
    dinitialconc = {(16,): c16_0, (18,): c18_0}
    #assert sum([dinitialconc[k] for k in dinitialconc]) == initial_C
    assert math.isclose(sum(dinitialconc[k] for k in dinitialconc), initial_C, rel_tol=1e-9)
    return dinitialconc

def get_alcohol_initial_concentrations(delta_O, initial_C=55.5, Rstd=0.002004):
    """
    Get the concentration of H2O(18) and H2O(16)
    """
    r0_alc = (delta_O / 1000 + 1) * Rstd
    c16_0 = initial_C / (1 + r0_alc)
    c18_0 = c16_0 * r0_alc
    dinitialconc = {(12, 16): c16_0, (12, 18): c18_0}
    #assert sum([dinitialconc[k] for k in dinitialconc]) == initial_C
    assert math.isclose(sum(dinitialconc[k] for k in dinitialconc), initial_C, rel_tol=1e-9)
    return dinitialconc

def get_organophosphate_initial_concentrations(delta_O, initial_C=1e-3, Rstd=0.002004):
    """
    Get the concentration of organophosphates
    """
    r0_orgpo = (delta_O / 1000 + 1) * Rstd
    c12_16_4 = initial_C / (1 + (4 * r0_orgpo)/(1 - 3 * r0_orgpo))
    c12_16_3_18_1 = c12_16_4 * r0_orgpo / (1 - 3 * r0_orgpo)
    dinitialconc = {(12, 16, 16, 16, 16): c12_16_4, 
                    (12, 18, 16, 16, 16): c12_16_3_18_1,
                    (12, 16, 18, 16, 16): c12_16_3_18_1,
                    (12, 16, 16, 18, 16): c12_16_3_18_1,
                    (12, 16, 16, 16, 18): c12_16_3_18_1}
    #assert sum([dinitialconc[k] for k in dinitialconc]) == initial_C
    assert math.isclose(sum(dinitialconc[k] for k in dinitialconc), initial_C, rel_tol=1e-9)
    return dinitialconc

def get_phosphate_initial_concentrations(delta_O, initial_C=1e-3, Rstd=0.002004):
    """
    Get the concentration of organophosphates
    """
    r0_po = (delta_O / 1000 + 1) * Rstd
    c16_4 = initial_C / (1 + (4 * r0_po)/(1 - 3 * r0_po))
    c16_3_18_1 = c16_4 * r0_po / (1 - 3 * r0_po)
    dinitialconc = {(16, 16, 16, 16): c16_4,
                    (18, 16, 16, 16): c16_3_18_1,
                    (16, 18, 16, 16): c16_3_18_1,
                    (16, 16, 18, 16): c16_3_18_1,
                    (16, 16, 16, 18): c16_3_18_1}
    #assert sum([dinitialconc[k] for k in dinitialconc]) == initial_C
    assert math.isclose(sum(dinitialconc[k] for k in dinitialconc), initial_C, rel_tol=1e-9)
    return dinitialconc

def _get_filtered_scramble_reactions(reactions):
    """
    Filters out reactions where the O18 gets scrambled in the hypothetical
    transition state.
    """
    filtreactions = []
    for rxni in reactions:
        filtflag = True
        lhs, rhs, kd, kr = rxni
        lhsi, lhsj = lhs
        rhsi = rhs[0]
        assert len(lhsj) == 1
        assert len(rhsi) == 6
        # Search for consecutive 18s starting at any position except len(lst) - 2
        for i in range(len(rhsi) - 1):
            if rhsi[i] == 18 and rhsi[i+1] == 18:
                # Check if this pair is NOT in the last two positions
                if i != len(rhsi) - 2:
                    filtflag = False # Filter it out
        # Second O18 can only come from the reactive water
        if rhsi.count(18) == 2 and rhsi[-1] == 16:
            filtflag = False
        if filtflag:
            filtreactions.append(rxni)
    return filtreactions

def create_reactions_for_step_one(dkie, filtreactions=True, verbose=False):
    """
    Function to create all combinations of isotops with the following shape:

    OOPOO + HOH --> OOPOOO 
    
    OOPOO  = ["R1", "lhs", 16, 16, 16, 16, 0]
    OOPOOO = ["R1", "rhs", 16, 16, 16, 16, 18]

    """
    reactions = [] # [reactant, product, kdirect, kreverse]
    #kd, kr = 1, 0.001 # ! HARDCODED

    # isomer O16 O18 mix
    oxygens = [[16, 16, 16, 16, 16],[18, 16, 16, 16, 16], [18, 18, 16, 16, 16]]
    for template in oxygens: # only four out of five labile oxygens
        comb_po5 = list([list(o) for o in set(itertools.permutations(template))])
        for rhs_po5 in comb_po5:
           # the oxygen can only be lost from one specific position (arbitrarely [-1])
           _rhs_po5 = [12] + rhs_po5.copy()
           lhs_po4 = _rhs_po5.copy()
           h2o = lhs_po4[-1]
           c = lhs_po4[0]
           assert c == 12
           lhs_po4.pop()
           tmprhs = (tuple(_rhs_po5),)
           tmplhs = (tuple(lhs_po4), (h2o,))
           kd, kr = _calculate_rate_constant(tmplhs, tmprhs, dkie)
           reactions.append((tmplhs, tmprhs, kd, kr))
    if filtreactions:
        reactions = _get_filtered_scramble_reactions(reactions)
    if verbose: 
        for r in reactions:
            print(r)

    return reactions


def create_reactions_for_step_two(dkie, reactions1, verbose=False):
    """
    Function to create all combinations of isotops with the following shape:

    OOPOOO --> OOPOO

    OOPOO  = [16, 16, 16, 16, 0]
    OOPOOO = [16, 16, 16, 16, 18]

    """
    
    reactions = [] # [reactant, product, kdirect, kreverse]
    #kd, kr = 1, 0.001 # ! HARDCODED
    
    # isomer O16 O18 mix
    for rxni in reactions1:
        _lhs, _rhs, _kd, _kr = rxni
        lhs_po5 = list(_rhs[0][1:])
        #for lhs_po5 in comb_po5:
        _lhs_po5 = [12] + lhs_po5.copy()
        rhs_po4 = _lhs_po5.copy()
        o = rhs_po4[1]
        c = rhs_po4[0]
        assert c == 12
        rhs_po4.pop(0)
        rhs_po4.pop(0)
        tmplhs = (tuple(_lhs_po5),)
        tmprhs = (tuple(rhs_po4), (c, o))
        kd, kr = _calculate_rate_constant(tmplhs, tmprhs, dkie)
        reactions.append((tmplhs, tmprhs, kd, kr))
    
    if verbose:
        for r in reactions_worep:
            print(len(reactions_worep))
    return reactions

def calculate_concentrations_from_mk(reactions, initial_conc, timerange, verbose=True):
    """
    Convert simplified form of isotopic reactions to a format that KiNetX can read
    """
    import scine_kinetx as kx
    import networkx as nx
    # hash compounds to an index
    compounds = {}
    acc = 0
    for r in reactions:
        lhs, rhs, _, _ = r
        for hs in [lhs, rhs]:
            for hsi in hs:
                if hsi in compounds.keys():
                    pass
                else:
                    compounds[hsi] = acc
                    acc += 1

    if verbose: print("compounds", compounds)
    # create kinetx network object
    network_builder = kx.NetworkBuilder()
    n_compounds = len(compounds.keys())
    n_reactions = len(reactions) 
    n_channels_per_reaction = 1
    #network_builder.reserve(n_compounds, n_reactions, n_channels_per_reaction)
    
    initial_conc2 = {compounds[d]:initial_conc[d] for d in initial_conc}
    concentrations = []
    for ii in range(n_compounds):
        keys = initial_conc2.keys()
        if ii in keys:
            concentrations.append(initial_conc2[ii])
        else:
            concentrations.append(0)
    for i in range(n_compounds):
        idxcompound = i 
        network_builder.add_compound(1, str(idxcompound))
    edges = []
    for r in reactions:
        _lhs, _rhs, kd, kr = r
        lhs = [(compounds[o],1) for o in _lhs]
        rhs = [(compounds[o],1) for o in _rhs]
        lhs2 = [o for o in _lhs] 
        rhs2 = [o for o in _rhs]
        for e1 in lhs2:
            for e2 in rhs2:
                E1 = e1 #compounds[e1]
                E2 = e2 #compounds[e2]
                edges.append((E1, E2))
        network_builder.add_reaction([kd], [kr], lhs, rhs)
    #
    #import matplotlib.pyplot as plt
    #G = nx.Graph()
    #G.add_edges_from(edges)
    #plt.figure(figsize=(10,8))
    #pos = nx.spring_layout(G)
    #nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=600, arrows=True)
    #plt.savefig("test.png")
    
    concentration_data, time_data = [], []
    network = network_builder.generate()
    #solver = kx.Integrator.explicit_euler
    #solver = kx.Integrator.implicit_euler
    #solver = kx.Integrator.cash_karp_5
    solver = kx.Integrator.cvode_bdf
    tstart = 0
    batch_interval = 5000
    nbatches = 10000
    maxtime = 0
    convergence = 1e-10
    concentration_data.append(concentrations)
    time_data.append(timerange[0])
    for idx in range(len(timerange)-1):
        tstart = timerange[idx]
        maxtime = timerange[idx+1]
        dt = timerange[idx+1]
        try:
            concentration_tmp, r_flux, rf_flux, rb_flux = kx.integrate(network, np.asarray(
            concentrations),  tstart, dt, solver, batch_interval, nbatches, convergence, integrateByTime=True, maxTime=maxtime)
            col1, col2, col3 = zip(*concentration_tmp)
            concentrations = col1
            if all(x >= 0 for x in col1):
                concentration_data.append(col1)
                time_data.append(timerange[idx])
        except RuntimeError:
            break
            #concentration_data.append([None for _ in concentrations])


    return concentration_data, time_data, compounds


def _read_crn_data(file):
    """
    Not related to the toy model
    """
    import ast   
    data = []
    with open(file, "r") as infile:
        for line in infile:
            tmpline = line.split("   ")
            tmplist = []
            for t in tmpline:
                tt = ast.literal_eval(t)
                tmplist.append(tt)
            data.append(tmplist)
    return data


def calculate_mk_yields(nodesfile, edgesfile, labelsfile):
    """
    Not related to the toy model
    """
    import scine_kinetx as kx
    import numpy as np
    nodesdata = _read_crn_data(nodesfile)
    edgesdata = _read_crn_data(edgesfile)
    network_builder = kx.NetworkBuilder()
    concentrations_t0 = []
    for n in nodesdata:
        # its indexed from 0, instead of from 1, as in kinetx_legacy
        a, b, c, d = n
        idxcompound = a
        network_builder.add_compound(1, str(idxcompound-1))
        concentrations_t0.append(b)

    edges = []
    for e in edgesdata:
        r1, r2, p1, p2, kd, kr = e
        _lhs = (r1, r2)
        _rhs = (p1, p2)
        lhs = [(int(o)-1,1) for o in _lhs if o > 0]
        rhs = [(int(o)-1,1) for o in _rhs if o > 0]
        lhs2 = [o for o in _lhs]
        rhs2 = [o for o in _rhs]
        network_builder.add_reaction([kd], [kr], lhs, rhs)
    
    network = network_builder.generate()
    #solver = kx.Integrator.explicit_euler
    #solver = kx.Integrator.implicit_euler
    #solver = kx.Integrator.cash_karp_5
    solver = kx.Integrator.cvode_bdf
    tstart = 0
    batch_interval = 5000
    nbatches = 10000
    convergence = 1e-10
    dt = 1
    maxtime = 1000
    concentration_tmp, r_flux, rf_flux, rb_flux = kx.integrate(network, np.asarray(
    concentrations_t0),  tstart, dt, solver, batch_interval, nbatches, convergence)
    col1, col2, col3 = zip(*concentration_tmp)
    concentrations = col1


    import numpy as np
    import matplotlib.pyplot as plt
    from operator import itemgetter

    # Convert last row to floats
    fig, ax = plt.subplots(figsize=(4.8,7.4))
    labels = range(len(nodesdata))
    y_positions = np.arange(len(labels))
    tmp = sorted(zip(labels, concentrations))
    slabels = [a for a, _ in tmp]
    sconcentration_data = [c for _, c in tmp]
    sy_positions = np.arange(len(slabels))  # new positions in sorted order

    plt.barh(sy_positions, sconcentration_data)
    plt.yticks(sy_positions, slabels)
    plt.xlabel("Concentration")
    plt.tight_layout()
    plt.show()


    return concentrations


def _is_concentration_sparse(log_values, drop_threshold=60, verbose=True):
    """
    Detects numerically unstable downward jumps in concentration time series.
    Works in log space to handle extremely small values robustly.
    """
    log_values = np.asarray(log_values, dtype=float)

    diffs = np.diff(log_values)

    output = np.any(diffs < -drop_threshold)

    if verbose:
        print("Is there a concentration jump?", output)
        if output:
            print("Diffs:", diffs)
    return output


