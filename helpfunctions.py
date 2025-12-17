import itertools
import copy
import numpy as np

def _calculate_rate_constant(lhs, rhs, dkie, penalty):
    """
    Heuristic based rate constants. Reverse rate constant is considered 100
    times slower than direct rate constant.

    Pair subindex constants are in forward direction, otherwise they are in 
    backward direction.

    """
    #print("! penalty harcoded"); penalty = 1e-10  # hardcoded; can be changed   
    
    k1p_18 = dkie["k1_16"] / dkie["kie1p"]
    k1s_18 = dkie["k1_16"] / dkie["kie1s"]

    k2p_18 = dkie["k2_16"] / dkie["kie2p"]
    k2s_18 = dkie["k2_16"] / dkie["kie2s"]

    k3p_18 = dkie["k3_16"] / dkie["kie3p"]
    k3s_18 = dkie["k3_16"] / dkie["kie3s"]

    k4p_18 = k3p_18 * penalty
    k4s_18 = k3s_18 * penalty
    

    if len(lhs) == 2:        # reaction 1
        r1, h2o = lhs        # organophosphate and water
        tsi = rhs[0]            # 'transition state' intermediate
        kf, kb = None, None
        assert len(h2o) == 1 # make sure its the h2o
        if 18 in h2o:
            kf = k1p_18      # kf depends only on h2o
            if tsi[-1] == 18: 
                kb = k2p_18 * penalty
            elif tsi[-1] == 16:
                kb = k2s_18 * penalty
        elif 16 in h2o: 
            if 18 in r1: 
                kf = k1s_18
                if tsi[-1] == 18:
                    kb = k2p_18 * penalty
                elif tsi[-1] == 16:
                    kb = k2p_18 * penalty
            else:
                kf = dkie["k1_16"]
                kb = dkie["k2_16"] * penalty
            
    elif len(lhs) == 1: # reaction 2
        tsi = lhs[0]     # 'transition state' intermediate
        p1, ro = rhs  # phosphate and alcohol
        if 18 in ro:
            kf = k3p_18
            kb = k3p_18 * penalty
        else: 
            if 18 in p1:
                kf = k3s_18
                kb = k3s_18 * penalty
            else:
                kf = dkie["k3_16"]
                kb = dkie["k3_16"] * penalty
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
    return dinitialconc

def create_reactions_for_step_one(dkie, penalty, verbose=False):
    """
    Function to create all combinations of isotops with the following shape:

    OOPOO + HOH --> OOPOOO 
    
    OOPOO  = ["R1", "lhs", 16, 16, 16, 16, 0]
    OOPOOO = ["R1", "rhs", 16, 16, 16, 16, 18]

    """
    reactions = [] # [reactant, product, kdirect, kreverse]
    #kd, kr = 1, 0.001 # ! HARDCODED

    # isomer O16 O18 mix
    for template in [[16, 16, 16, 16, 16],[18, 16, 16, 16, 16], [18, 18, 16, 16, 16]]: # only four out of five labile oxygens
        comb_po5 = list([list(o) for o in set(itertools.permutations(template))])
        print("=======================================", comb_po5)
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
           kd, kr = _calculate_rate_constant(tmplhs, tmprhs, dkie, penalty)
           reactions.append((tmplhs, tmprhs, kd, kr))

    if verbose: 
        for r in reactions:
            print(r)

    return reactions


def create_reactions_for_step_two(dkie, penalty, verbose=False):
    """
    Function to create all combinations of isotops with the following shape:

    OOPOOO --> OOPOO

    OOPOO  = [16, 16, 16, 16, 0]
    OOPOOO = [16, 16, 16, 16, 18]

    """
    
    reactions = [] # [reactant, product, kdirect, kreverse]
    #kd, kr = 1, 0.001 # ! HARDCODED
    
    # isomer O16 O18 mix
    for template in [[16, 16, 16, 16, 16],[18, 16, 16, 16, 16], [18, 18, 16, 16, 16]]: # only four out of five labile oxygens
        comb_po5 = list([list(o) for o in set(itertools.permutations(template))])
        print("=======================================", comb_po5)
        for lhs_po5 in comb_po5:
                # the oxygen can only be lost from one specific position (arbitrarely [-1])
                _lhs_po5 = [12] + lhs_po5.copy()
                rhs_po4 = _lhs_po5.copy()
                h2o = rhs_po4[-1]
                c = rhs_po4[0]
                assert c == 12
                rhs_po4.pop()
                rhs_po4.pop(0)
                tmplhs = (tuple(_lhs_po5),)
                tmprhs = (tuple(rhs_po4), (c, h2o))
                kd, kr = _calculate_rate_constant(tmplhs, tmprhs, dkie, penalty)
                reactions.append((tmplhs, tmprhs, kd, kr))
    
    #reactions_worep = list(set(reactions))
    if verbose:
        for r in reactions_worep:
            print(len(reactions_worep))
    return reactions

def convert_to_kinetx_notation(reactions, initial_conc, timerange, verbose=True):
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
    import matplotlib.pyplot as plt
    G = nx.Graph()
    G.add_edges_from(edges)
    plt.figure(figsize=(10,8))
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=600, arrows=True)
    plt.savefig("test.png")
    
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
    #timerange = np.logspace(-8, -2, num=1000)
    #timerange = np.logspace(-8, 5, num=10000)
    #timerange = np.linspace(0, 1e+2, num=10)
    concentration_data.append(concentrations)
    time_data.append(0)
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
        print(str(idxcompound-1))
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
        print("reactions", [kd], [kr], lhs, rhs)
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

    print("CONC", concentrations)

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

def plot_ratio_vs_time(
    concentration_data,
    time_data,
    compounds,
    dkie,
    file=None,
    Rstd=0.002004,
    ax=None,
):
    """
    Convert concentrations to subtract ratio of O18/O16.
    Convention in the field.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Axes to plot into. If None, a new figure and axes are created.
    file : str, optional
        If provided, the figure is saved to this path.
    """
    import matplotlib.pyplot as plt

    # --- create axes if not provided ---
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots()
        created_fig = True
    else:
        fig = ax.figure

    # --- compute ratios ---
    Rsub = {"CPO4": [], "CPO5": [], "PO4": [], "O": [], "CO": []}

    for c in concentration_data:
        cpo4_18, cpo4_16 = [], []
        cpo5_18, cpo5_16 = [], []
        po4_18, po4_16 = [], []
        o_18, o_16 = [], []
        co_18, co_16 = [], []

        for d in compounds:
            idx = compounds[d]
            if len(d) == 1:      # O
                o_18.append(d.count(18) * c[idx])
                o_16.append(d.count(16) * c[idx])
            elif len(d) == 2:    # CO
                co_18.append(d.count(18) * c[idx])
                co_16.append(d.count(16) * c[idx])
            elif len(d) == 4:    # PO4
                po4_18.append(d.count(18) * c[idx])
                po4_16.append(d.count(16) * c[idx])
            elif len(d) == 5:    # CPO4
                cpo4_18.append(d.count(18) * c[idx])
                cpo4_16.append(d.count(16) * c[idx])
            elif len(d) == 6:    # CPO5
                cpo5_18.append(d.count(18) * c[idx])
                cpo5_16.append(d.count(16) * c[idx])
            else:
                raise ValueError("Unknown compound length")

        def _get_zero_safe_ratio(x18, x16):
            den = sum(x16)
            num = sum(x18)
            if den > 0 and num > 0:
                ratio = num / den
                delta = ((ratio / Rstd) - 1) * 1000
            else:
                delta = 0
            return delta

        for a, b, key in [
            (cpo4_18, cpo4_16, "CPO4"),
            (po4_18,  po4_16,  "PO4"),
            (cpo5_18, cpo5_16, "CPO5"),
            (o_18,    o_16,    "O"),
            (co_18,   co_16,   "CO"),
        ]:
            Rsub[key].append(_get_zero_safe_ratio(a, b))

    # --- plotting ---
    ax.set_xscale("log")
    for key in Rsub:
        ax.plot(time_data, Rsub[key], ".-", label=key)

    ax.set_xlabel("Time (s)")
    ax.set_ylabel(r"$\delta^{18}$O(S)")

    titlestr = _get_title(dkie)
    ax.set_title(titlestr)

    ax.legend()

    if created_fig:
        fig.tight_layout()
        if file is not None:
            fig.savefig(file)

    return ax



def _get_title(dkie):
    """
    Create title from the defined dkie. Done for keeping the record.
    """
    cnt = 1
    titlelist = []
    for d in dkie:
        if cnt == 4:
            cnt = 0
            titlelist.append(d+"="+str(dkie[d])+"\n")
        else:
            titlelist.append(d+"="+str(dkie[d]))
        cnt += 1
    titlestr = "  ".join(titlelist)
    return titlestr


def plot_conc_vs_time(concentration_data, time_data, compounds, dkie, file=None, ax=None):
    """
    Plot evolution of concentrations in time.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Axes to plot into. If None, a new figure and axes are created.
    file : str, optional
        If provided, the figure is saved to this path.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # Create axes if not provided
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.8, 8.4))
        created_fig = True
        ax.set_prop_cycle(color=plt.cm.tab20.colors)
    else:
        fig = ax.figure

    concentration_dataT = np.array(concentration_data).T
    idx = 0

    for d, c in zip(compounds, concentration_dataT):
        if max(c) < 1e-20:
            continue

        if idx < 20:
            plottype = '-'
            linewidth = 2
        elif 20 <= idx < 40:
            plottype = '--'
            linewidth = 3.5
        else:
            plottype = ':'
            linewidth = 2

        ax.plot(time_data, c, plottype, label=str(d), linewidth=linewidth)
        idx += 1

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Concentration")

    # Title
    titlestr = _get_title(dkie)
    ax.set_title(titlestr)

    ax.legend(loc="lower center", ncol=2)

    if created_fig:
        fig.tight_layout()
        if file is not None:
            fig.savefig(file)

    return ax


def plot_concentrations(concentration_data, compounds, dkie, outfile):
    """
    Plot final concentrations of all compounds
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from operator import itemgetter

    # Convert last row to floats
    fig, ax = plt.subplots(figsize=(4.8,7.4))
    labels = list(compounds.keys())
    y_positions = np.arange(len(labels))
    tmp = sorted(zip(labels, concentration_data[-1]), key=lambda x: len(x[0]), reverse=True)
    slabels = [a for a, _ in tmp]
    sconcentration_data = [c for _, c in tmp]
    sy_positions = np.arange(len(slabels))  # new positions in sorted order
    
    plt.barh(sy_positions, sconcentration_data)
    plt.yticks(sy_positions, slabels)
    plt.xlabel("Concentration")
    plt.xscale("log")

    # Get title
    titlestr = _get_title(dkie)
    plt.title(titlestr)

    plt.tight_layout()
    plt.savefig(outfile)
