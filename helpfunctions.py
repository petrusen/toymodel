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
    #k16   = 0.001    # no O18
    
    #k16 = 1
    #k18ps = k16 / kie * 1.011 # O18 in reactive and non-reactive positions
    #k18p  = k16 / kie * 1.01  # O18 in reactive position (primary)
    #k18s  = k16 / kie * 0.999   # O18 in non-reactive position (secundary)

    k16 = 1
    k18p = 0.96
    k18s = 0.995
    k18sp = k18p * k18s

    rhspenalty = 100000
    if len(lhs) == 2: # reaction 1
        # write it up
        lhsi, lhsj = lhs
        assert len(lhsj) == 1 # make sure its the h2o 
        if 18 in lhsj:
            Ks = (k18p, k18p/rhspenalty)
        elif 16 in lhsj:
            Ks = (k16, k16/rhspenalty)
    elif len(lhs) == 1: # reaction 2
        lhs = lhs[0]
        if 18 not in lhs:
            Ks = (k16, k16/rhspenalty)
        elif 18 in lhs:
            if lhs[-1] == 18: # reactive site
                if lhs.count(18) == 1:
                    Ks = (k18p, k18p/rhspenalty)
                else:
                    Ks = (k18sp, k18sp/rhspenalty)
            elif lhs[-1] != 18: # non reactive site
                Ks = (k18s, k18s/rhspenalty)
    kd, kr = Ks
    #kd, kr = 1, 0.01
    return kd, kr

def create_reactions_for_step_one(kie, verbose=False):
    """
    Function to create all combinations of isotops with the following shape:

    OOPOO + HOH --> OOPOOO 
    
    OOPOO  = ["R1", "lhs", 16, 16, 16, 16, 0]
    OOPOOO = ["R1", "rhs", 16, 16, 16, 16, 18]

    """
    reactions = [] # [reactant, product, kdirect, kreverse]
    #kd, kr = 1, 0.001 # ! HARDCODED

    # isomer O16 pure
    template = [12, 16, 16, 16, 16]
    lhs_po4_16 = (tuple(template), (16,))
    lhs_po4_18 = (tuple(template), (18,))
    rhs_po4_16 = (tuple(template + [16]),)
    rhs_po4_18 = (tuple(template + [18]),)
    print(lhs_po4_16)
    kd, kr = _calculate_rate_constant(lhs_po4_16, kie)
    reactions.append((lhs_po4_16, rhs_po4_16, kd, kr))
    kd, kr = _calculate_rate_constant(lhs_po4_18, kie)
    reactions.append((lhs_po4_18, rhs_po4_18, kd, kr))

    # isomer O16 O18 mix
    template = [18, 16, 16, 16] # ! HARDCODED
    comb_po4 = list([list(o) for o in set(itertools.permutations(template))])
    for _lhs_po4 in comb_po4:
        lhs_po4_16 = (tuple([12] + _lhs_po4), (16,))
        lhs_po4_18 = (tuple([12] + _lhs_po4), (18,))
        rhs_po4_16 = (tuple([12] + _lhs_po4 + [16]),)
        rhs_po4_18 = (tuple([12] + _lhs_po4 + [18]),)

        kd, kr = _calculate_rate_constant(lhs_po4_16, kie)
        reactions.append((lhs_po4_16, rhs_po4_16, kd, kr))
        
        kd, kr = _calculate_rate_constant(lhs_po4_18, kie)
        reactions.append((lhs_po4_18, rhs_po4_18, kd, kr))

    if verbose: 
        for r in reactions:
            print(r)
        print(len(reactions))

    return reactions


def create_reactions_for_step_two(kie, verbose=False):
    """
    Function to create all combinations of isotops with the following shape:

    OOPOOO --> OOPOO

    OOPOO  = [16, 16, 16, 16, 0]
    OOPOOO = [16, 16, 16, 16, 18]

    """
    
    reactions = [] # [reactant, product, kdirect, kreverse]
    #kd, kr = 1, 0.001 # ! HARDCODED
    
    # isomer O16 pure
    template = [12, 16, 16, 16, 16, 16] # ! HARDCODED 
    lhs_po5 = (tuple(template),)
    _rhs_po4_16 = template.copy()
    _rhs_po4_16.pop()
    _rhs_po4_16.pop(0)
    ## only one is created because if there are no O16, permutation are
    ## equivalent in terms of rate constants
    rhs_po4_16 = (tuple(_rhs_po4_16), (12, 16))
    kd, kr = _calculate_rate_constant(lhs_po5, kie)
    reactions.append((lhs_po5, rhs_po4_16, kd, kr))

    # isomer O16 O18 mix
    for template in [[18, 16, 16, 16, 16], [18, 18, 16, 16, 16]]: # ! HARDCODED
        comb_po5 = list([list(o) for o in set(itertools.permutations(template))])
        for lhs_po5 in comb_po5:
           # the oxygen can only be lost from one specific position (arbitrarely [-1])
            _lhs_po5 = [12] + lhs_po5.copy()
            kd, kr = _calculate_rate_constant(tuple((_lhs_po5,)), kie)
            rhs_po4 = _lhs_po5.copy()
            h2o = rhs_po4[-1]
            c = rhs_po4[0]
            assert c == 12
            rhs_po4.pop()
            rhs_po4.pop(0)
            reactions.append(((tuple(_lhs_po5),), (tuple(rhs_po4), (c, h2o)), kd, kr))
    
    #reactions_worep = list(set(reactions))
    if verbose:
        for r in reactions_worep:
            print(len(reactions_worep))
    return reactions

def convert_to_kinetx_notation(reactions, initial_conc, concfile=False, timefile=False, verbose=True):
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
    
    #import matplotlib.pyplot as plt
    #G = nx.Graph()
    #G.add_edges_from(edges)
    #plt.figure(figsize=(10,8))
    #pos = nx.spring_layout(G)
    #nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=600, arrows=True)
    #plt.show()
    
    concentration_data = []
    network = network_builder.generate()
    #solver = kx.Integrator.explicit_euler
    #solver = kx.Integrator.implicit_euler
    #solver = kx.Integrator.cash_karp_5
    solver = kx.Integrator.cvode_bdf
    tstart = 0
    batch_interval = 5000
    nbatches = 10000
    maxtime = 0
    convergence = 1e-20
    #timerange = np.logspace(-8, -2, num=1000)
    timerange = np.logspace(-8, 7, num=1000)
    #timerange = np.linspace(0, 1e+4, num=1000)
    concentration_data.append(concentrations)
    for idx in range(len(timerange)-1):
        tstart = timerange[idx]
        maxtime = timerange[idx+1]
        dt = timerange[idx+1]/10
        try:
            concentration_tmp, r_flux, rf_flux, rb_flux = kx.integrate(network, np.asarray(
            concentrations),  tstart, dt, solver, batch_interval, nbatches, convergence, integrateByTime=True, maxTime=maxtime)
            col1, col2, col3 = zip(*concentration_tmp)
            concentrations = col1
            concentration_data.append(col1)
        except RuntimeError:
            pass #concentration_data.append([None for _ in concentrations])


    if (concfile is not False) and (timefile is not False):
        if verbose: print("Save concentration and time files")
        np.save(concfile, concentrations)
        np.save(timefile, timerange)

    return concentration_data, timerange, compounds

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

def plot_ratio_vs_time(concentration_data, time_data, compounds, file, Rstd=0.002004):
    """
    Convert concentrations to substract ratio of O18/O16. Convention in the field.
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    Rsub = {"CPO4": [], "CPO5": [], "PO4": []}
    for c in concentration_data[1:]:
        tmpo18, tmpo16 = [], []
        cpo4_18, cpo4_16 = [], []
        cpo5_18, cpo5_16 = [], []
        po4_18, po4_16 = [], []
        for d in compounds:
            if len(d) == 5:
                if 12 in d: # CPO4
                    cpo4_18.append(d.count(18) * c[compounds[d]])
                    cpo4_16.append(d.count(16) * c[compounds[d]])
                elif 12 not in d: # PO4
                    po4_18.append(d.count(18) * c[compounds[d]])
                    po4_16.append(d.count(16) * c[compounds[d]])
            elif len(d) == 6: #CPO5
                cpo5_18.append(d.count(18) * c[compounds[d]])
                cpo5_16.append(d.count(16) * c[compounds[d]])
        cpo4_sub = sum(cpo4_18) / sum(cpo4_16)
        Rsub["CPO4"].append(cpo4_sub / Rstd - 1 * 1000)
        po4_sub = sum(po4_18) / sum(po4_16)
        Rsub["PO4"].append(po4_sub / Rstd - 1 * 1000)
        cpo5_sub = sum(cpo5_18) / sum(cpo5_16)
        Rsub["CPO5"].append(cpo5_sub / Rstd - 1 * 1000)
    #plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.tab20.colors)
    fig, ax = plt.subplots()
    d_ind_stoich = {compounds[d]:d for d in compounds}
    plt.xscale("log")
    for d in Rsub:
        plt.plot(time_data[1:], Rsub[d], label=d)
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
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.tab20.colors)
    fig, ax = plt.subplots()
    concentration_dataT = np.array(concentration_data).T
    idx = 0
    d_ind_stoich = {compounds[d]:d for d in compounds}
    for d, c in zip(compounds, concentration_dataT):
        if idx < 20:
            plottype = '-'
            linewidth = 2
        elif 20 <= idx < 40:
            plottype = '--'
            linewidth = 3.5
        if len(d) == 6:
            plt.plot(time_data, c, plottype, label=d_ind_stoich[idx], linewidth=linewidth)
        else:
            plt.plot(time_data, c, plottype, label=d_ind_stoich[idx], linewidth=linewidth)
        idx += 1
    for s, c in zip(compounds,concentration_data[-1]):
        print(s, c)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Time (s)")
    plt.ylabel("Concentration")
    plt.tight_layout()
    plt.legend(loc="best", ncol=2)
    plt.savefig(file)
    plt.show()

def plot_concentrations(concentration_data, compounds):
    import numpy as np
    import matplotlib.pyplot as plt
   
    # Convert last row to floats
    fig, ax = plt.subplots(figsize=(4.8,7.4))
    labels = list(compounds.keys())
    y_positions = np.arange(len(labels))
   
    plt.barh(y_positions, concentration_data[-1])   #<-- horizontal bars
   
    plt.yticks(y_positions, labels)
    plt.xlabel("Concentration")
    plt.xscale("log")
    plt.tight_layout()
    plt.show() 
