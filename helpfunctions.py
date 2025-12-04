import itertools
import copy
import numpy as np

def _calculate_rate_constant(lhs, kie=1):
    """
    Heuristic based rate constants. Reverse rate constant is considered 100 
    times slower than direct rate constant.
    """
    #k18ps = kie * 1.0011 # O18 in reactive and non-reactive positions
    #k18p  = kie * 1.001  # O18 in reactive position (primary)
    #k18s  = kie * 1.01   # O18 in non-reactive position (secundary)
    #k16   = 1.0    # no O18
    k16   = 0.001    # no O18
    k18ps = k16 / kie * 1.0011 # O18 in reactive and non-reactive positions
    k18p  = k16 / kie * 1.001  # O18 in reactive position (primary)
    k18s  = k16 / kie * 1.01   # O18 in non-reactive position (secundary)

    if len(lhs) == 2: # reaction 1
        # write it up
        lhsi, lhsj = lhs
        assert len(lhsj) == 1 # make sure its the h2o 
        if 18 in lhsj:
            Ks = (k18p, k18p/100)
        elif 16 in lhsj:
            Ks = (k16, k16/100)
    elif len(lhs) == 1: # reaction 2
        lhs = lhs[0]
        if 18 not in lhs:
            Ks = (k16, k16/100)
        elif 18 in lhs:
            if lhs.index(18) == 5: # reactive site
                if lhs.count(18) == 1:
                    Ks = (k18p, k18p/100)
                elif lhs.count(18) == 2:
                    Ks = (k18ps, k18ps/100)
            elif lhs.index(18) != 5: # non reactive site
                Ks = (k18s, k18s/100)
    kd, kr = Ks
    return kd, kr

def create_reactions_for_step_one(verbose=False):
    """
    Function to create all combinations of isotops with the following shape:

    OOPOO + HOH --> OOPOOO 
    
    OOPOO  = ["R1", "lhs", 16, 16, 16, 16, 0]
    OOPOOO = ["R1", "rhs", 16, 16, 16, 16, 18]

    """
    reactions = [] # [reactant, product, kdirect, kreverse]
    #kd, kr = 1, 0.001 # ! HARDCODED

    # isomer O16 pure
    template = [16, 16, 16, 16]
    lhs_po4_16 = (tuple(template), (16,))
    lhs_po4_18 = (tuple(template), (18,))
    rhs_po4_16 = (tuple(template + [16]),)
    rhs_po4_18 = (tuple(template + [18]),)
    print(lhs_po4_16)
    kd, kr = _calculate_rate_constant(lhs_po4_16, kie=1)
    reactions.append((lhs_po4_16, rhs_po4_16, kd, kr))
    kd, kr = _calculate_rate_constant(lhs_po4_18, kie=1)
    reactions.append((lhs_po4_18, rhs_po4_18, kd, kr))

    # isomer O16 O18 mix
    template = [18, 16, 16, 16] # ! HARDCODED
    comb_po4 = list([list(o) for o in set(itertools.permutations(template))])
    for _lhs_po4 in comb_po4:
        lhs_po4_16 = (tuple(_lhs_po4), (16,))
        lhs_po4_18 = (tuple(_lhs_po4), (18,))
        rhs_po4_16 = (tuple(_lhs_po4 + [16]),)
        rhs_po4_18 = (tuple(_lhs_po4 + [18]),)

        kd, kr = _calculate_rate_constant(lhs_po4_16, kie=1)
        reactions.append((lhs_po4_16, rhs_po4_16, kd, kr))
        
        kd, kr = _calculate_rate_constant(lhs_po4_18, kie=1)
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
    #kd, kr = 1, 0.001 # ! HARDCODED
    
    # isomer O16 pure
    template = [16, 16, 16, 16, 16] # ! HARDCODED 
    lhs_po5 = (tuple(template),)
    _rhs_po4_16 = template.copy()
    _rhs_po4_16.pop()
    ## only one is created because if there are no O16, permutation are
    ## equivalent in terms of rate constants
    rhs_po4_16 = (tuple(_rhs_po4_16), (16,))
    kd, kr = _calculate_rate_constant(lhs_po5, kie=1)
    reactions.append((lhs_po5, rhs_po4_16, kd, kr))

    # isomer O16 O18 mix
    for template in [[18, 16, 16, 16, 16], [18, 18, 16, 16, 16]]: # ! HARDCODED
        comb_po5 = list([list(o) for o in set(itertools.permutations(template))])
        for lhs_po5 in comb_po5:
           # the oxygen can only be lost from one specific position (arbitrarely [-1])
            _lhs_po5 = lhs_po5.copy()
            kd, kr = _calculate_rate_constant(tuple((lhs_po5,)), kie=1)
            rhs_po4 = _lhs_po5
            h2o = rhs_po4[-1]
            rhs_po4.pop()
            reactions.append(((tuple(lhs_po5),), (tuple(rhs_po4), (h2o,)), kd, kr))
    
    #reactions_worep = list(set(reactions))
    if verbose:
        for r in reactions_worep:
            print(len(reactions_worep))
    return reactions

def convert_to_kinetx_notation(reactions, initial_conc, verbose=True):
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
        print("reactions", r)
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
    print("====", n_reactions)
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
        print("reactionss", _lhs, _rhs)
        print("reactions post", lhs, rhs)
        lhs2 = [o for o in _lhs] 
        rhs2 = [o for o in _rhs]
        for e1 in lhs2:
            for e2 in rhs2:
                E1 = e1 #compounds[e1]
                E2 = e2 #compounds[e2]
                edges.append((E1, E2))
        network_builder.add_reaction([kd], [kr], lhs, rhs)
    
    import matplotlib.pyplot as plt
    G = nx.Graph()
    G.add_edges_from(edges)
    plt.figure(figsize=(10,8))
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=600, arrows=True)
    plt.show()
    
    concentration_data = []
    network = network_builder.generate()
    solver = kx.Integrator.explicit_euler
    #solver = kx.Integrator.implicit_euler
    tstart = 0
    batch_interval = 5000
    nbatches = 10000
    convergence = 1e-10
    maxtime = 0
    timerange = np.logspace(-8, -2, num=100)
    concentration_data.append(concentrations)
    for idx in range(len(timerange)-1):
        tstart = timerange[idx]
        maxtime = timerange[idx+1]
        dt = timerange[idx+1]/10
        concentration_tmp, r_flux, rf_flux, rb_flux = kx.integrate(network, np.asarray(
        concentrations),  tstart, dt, solver, batch_interval, nbatches, convergence, integrateByTime=True, maxTime=maxtime)
        col1, col2, col3 = zip(*concentration_tmp)
        concentrations = col1
        concentration_data.append(col1)
    return concentration_data, timerange, compounds

def plot_ratio_vs_time(concentration_data, time_data, compounds, file, Rstd=0.002004):
    """
    Convert concentrations to substract ratio of O18/O16. Convention in the field.
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    Rsub = {"PO4": [], "PO5": []}
    for c in concentration_data:
        tmpo18, tmpo16 = [], []
        po4_18, po4_16 = [], []
        po5_18, po5_16 = [], []
        for d in compounds:
            if len(d) == 4: # PO4
                po4_18.append(d.count(18) * c[compounds[d]])
                po4_16.append(d.count(16) * c[compounds[d]])
            elif len(d) == 5:
                po5_18.append(d.count(18) * c[compounds[d]])
                po5_16.append(d.count(16) * c[compounds[d]])
        po4_sub = sum(po4_18) / sum(po4_16)
        Rsub["PO4"].append(po4_sub / Rstd - 1 * 1000)
        po5_sub = sum(po5_18) / sum(po5_16)
        Rsub["PO5"].append(po5_sub / Rstd - 1 * 1000)
    #plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.tab20.colors)
    fig, ax = plt.subplots()
    d_ind_stoich = {compounds[d]:d for d in compounds}
    plt.xscale("log")
    for d in Rsub:
        plt.plot(time_data, Rsub[d], label=d)
    plt.xlabel("Time (s)")
    plt.ylabel("Ratio substract O18/O16")
    plt.legend()
    plt.tight_layout()
    plt.savefig(file)


def plot_conc_vs_time(concentration_data, time_data, compounds, file):
    """
    Plot evolution of concentrations in time
    """
    import matplotlib.pyplot as plt
    #plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.tab20.colors)
    fig, ax = plt.subplots()
    concentration_dataT = np.array(concentration_data).T
    idx = 0
    d_ind_stoich = {compounds[d]:d for d in compounds}
    for c in concentration_dataT:
        if max(c) < 1e-8:
            plt.plot(time_data, c, color='tab:grey', linewidth=3)
        else:
            plt.plot(time_data, c, label=d_ind_stoich[idx])
        idx += 1
    for s, c in zip(compounds,concentration_data[-1]):
        print(s, c)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Time (s)")
    plt.ylabel("Concentration")
    plt.tight_layout()
    plt.legend(loc='best') #loc="upper left")
    plt.savefig(file)
    
    
