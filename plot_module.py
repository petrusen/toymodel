import itertools
import copy
import numpy as np
import matplotlib.pyplot as plt

def _is_it_zero(x18, x16, Rstd=0.002004, tol=1e-9):
    """
    Returns True if either den or num is less than or equal to `tol`.
    """
    den = sum(x16)
    num = sum(x18)
    
    # Anything less than or equal to `tol` is considered zero
    return den <= tol or num <= tol

def _get_ratio(x18, x16, Rstd=0.002004):
    """
    Avoid dividing by zero.
    """
    den = sum(x16)
    num = sum(x18)
    ratio = num / den
    delta = ((ratio / Rstd) - 1) * 1000
    return delta

def plot_ratio_vs_time(
    concentration_data,
    time_data,
    compounds,
    dkie,
    file=None,
    Rstd=0.002004,
    ax=None,
    writecsv=False
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
    cnt = 0
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

        compounds_to_check = [
            (cpo4_18, cpo4_16, "CPO4"),
            (po4_18,  po4_16,  "PO4"),
            (cpo5_18, cpo5_16, "CPO5"),
            (o_18,    o_16,    "O"),
            (co_18,   co_16,   "CO"),
        ]
        
        all_valid = True
        
        for a, b, key in compounds_to_check:
            if _is_it_zero(a, b, Rstd=0.002004):
                all_valid = False
                break  # Stop checking early; at least one failed
        
        if all_valid:
            cnt += 1
            for a, b, key in compounds_to_check:
                Rsub[key].append(_get_ratio(a, b))

    # --- plotting ---
    ax.set_xscale("log")
    for key in Rsub:
        ax.plot(time_data[:cnt], Rsub[key], label=key)

    if writecsv is not False:
        safedata = []
        safedata.append(time_data[:cnt])
        for key in Rsub:
            safedata.append(Rsub[key])
        safedataT = np.array(safedata).T
        header = "Time (s),"+",".join(Rsub.keys())
        np.savetxt(writecsv.split(".")[0]+".csv", safedataT, delimiter=",", header=header)

    ax.set_xlabel("Time (s)")
    ax.set_ylabel(r"$\delta^{18}$O(S)")
    #ax.set_yscale("log")
    titlestr = _get_title(dkie)
    #ax.set_title(titlestr)

    ax.legend(ncols=5)

    if created_fig:
        fig.tight_layout()
        if file is not None:
            fig.savefig(file)

    return ax


def plot_ratio_vs_time_separate(
    concentration_data,
    time_data,
    compounds,
    dkie,
    file=None,
    Rstd=0.002004,
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

        for a, b, key in [
            (cpo4_18, cpo4_16, "CPO4"),
            (po4_18,  po4_16,  "PO4"),
            (cpo5_18, cpo5_16, "CPO5"),
            (o_18,    o_16,    "O"),
            (co_18,   co_16,   "CO"),
        ]:
            Rsub[key].append(_get_ratio(a, b))

    # --- plotting ---
    for key in Rsub:
        fig, ax = plt.subplots()
        #ax.ticklabel_format(axis='y', style='plain', useOffset=False)
        ax.set_xscale("log")
        ax.plot(time_data, Rsub[key], label=key)
        tmpname = file.split(".")[0] + "_O_vs_time" + "_" + key + ".png"
        ax.legend()
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(r"$\delta^{18}$O(S)")
        plt.tight_layout()
        plt.savefig(tmpname)

    return None


def plot_ratio_vs_reaction_progress(
    concentration_data,
    time_data,
    compounds,
    dkie,
    file=None,
    Rstd=0.002004,
    ax=None,
    reactionref=(12, 16, 16, 16, 16),
    writecsv=False
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
    Csumi = {"CPO4": [], "CPO5": [], "PO4": [], "O": [], "CO": []}
    cnt = 0
    for c in concentration_data:
        cpo4_18, cpo4_16 = [], []
        cpo5_18, cpo5_16 = [], []
        po4_18, po4_16 = [], []
        o_18, o_16 = [], []
        co_18, co_16 = [], []
        Csumj = {"CPO4": [], "CPO5": [], "PO4": [], "O": [], "CO": []}
        for d in compounds:
            idx = compounds[d]
            if len(d) == 1:      # O
                o_18.append(d.count(18) * c[idx])
                o_16.append(d.count(16) * c[idx])
                Csumj["O"].append(c[idx])
            elif len(d) == 2:    # CO
                co_18.append(d.count(18) * c[idx])
                co_16.append(d.count(16) * c[idx])
                Csumj["CO"].append(c[idx])
            elif len(d) == 4:    # PO4
                po4_18.append(d.count(18) * c[idx])
                po4_16.append(d.count(16) * c[idx])
                Csumj["PO4"].append(c[idx])
            elif len(d) == 5:    # CPO4
                cpo4_18.append(d.count(18) * c[idx])
                cpo4_16.append(d.count(16) * c[idx])
                Csumj["CPO4"].append(c[idx])
            elif len(d) == 6:    # CPO5
                cpo5_18.append(d.count(18) * c[idx])
                cpo5_16.append(d.count(16) * c[idx])
                Csumj["CPO5"].append(c[idx])
            else:
                raise ValueError("Unknown compound length")

        compounds_to_check = [
            (cpo4_18, cpo4_16, "CPO4"),
            (po4_18,  po4_16,  "PO4"),
            (cpo5_18, cpo5_16, "CPO5"),
            (o_18,    o_16,    "O"),
            (co_18,   co_16,   "CO"),
        ]

        all_valid = True

        for a, b, key in compounds_to_check:
            if _is_it_zero(a, b, Rstd=0.002004):
                all_valid = False
                break  # Stop checking early; at least one failed

        if all_valid:
            cnt += 1
            for key in Csumj.keys():
                Csumi[key].append(sum(Csumj[key]))
            for a, b, key in compounds_to_check:
                Rsub[key].append(_get_ratio(a, b))


    concentration_data_T = np.array(concentration_data).T
    
    # reaction progress reference selected by the user
    if isinstance(reactionref, tuple): 
        # --- plotting ---
        reacind = compounds[reactionref]
        reaction_progress = concentration_data_T[reacind][:cnt]
        reaction_progress_norm = [r/reaction_progress[0] for r in reaction_progress]
        for key in Rsub:
            ax.plot(reaction_progress_norm, Rsub[key], label=key)
        ax.set_xlabel("Reaction progress respect to {a}".format(a=str(reactionref)))
    elif isinstance(reactionref, str):
        # --- plotting ---
        reaction_progress_norm = [r/Csumi[reactionref][0] for r in Csumi[reactionref][:cnt]]
        for key in Rsub:
            ax.plot(reaction_progress_norm, Rsub[key], label=key)
        ax.set_xlabel("Reaction progress respect to {a}".format(a=reactionref))
    
    ax.set_ylabel(r"$\delta^{18}$O(S)")
    #ax.set_yscale("log")
    ax.invert_xaxis()
    titlestr = _get_title(dkie)
    #ax.set_title(titlestr)

    ax.legend(ncols=5)

    if created_fig:
        fig.tight_layout()
        if file is not None:
            fig.savefig(file)


    if writecsv is not False:
        safedata = []
        safedata.append(reaction_progress_norm)
        for key in Rsub:
            safedata.append(Rsub[key])
        safedataT = np.array(safedata).T
        header = "Reaction Progress {a}".format(a=reactionref)+","+",".join(Rsub.keys())
        np.savetxt(writecsv.split(".")[0]+".csv", safedataT, delimiter=",", header=header)

    return ax


def _get_title(dkie, cntthresh=6):
    """
    Create title from the defined dkie. Done for keeping the record.
    """
    cnt = 1
    titlelist = []
    for d in dkie:
        if cnt == cntthresh:
            cnt = 0
            titlelist.append(d+"="+str(dkie[d])+"\n")
        else:
            titlelist.append(d+"="+str(dkie[d]))
        cnt += 1
    titlestr = "  ".join(titlelist)
    return titlestr

def _create_extended_tab20():
    """
    Due to the large amount of compounds formed in the MKs simulations
    I have extended the color legend.
    """

    tab20 = plt.cm.tab20.colors  # returns 20 RGB tuples
    extended_colors = []

    for i in range(0, 20, 2):  # step by 2, since pairs are dark/light
        dark = np.array(tab20[i])
        light = np.array(tab20[i+1])

        # Create a medium color by blending dark and light (alpha blending)
        alpha = 0.5
        medium = alpha * dark + (1 - alpha) * light

        extended_colors.extend([dark, medium, light])

    return extended_colors

def plot_mk_information(dkie, handles, labels, ax):
    """
    Placeholder information plot
    """
    # Create axes if not provided
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.8, 8.4))
        created_fig = True
        colors = _create_extended_tab20() 
        ax.set_prop_cycle(color=colors)
        #ax.set_prop_cycle(color=plt.cm.tab20.colors)
    else:
        fig = ax.figure
    string = _get_title(dkie, cntthresh=1)
    ax.text(0.75, 0.475, string, transform=ax.transAxes, ha='left', va='center', fontsize=10)
    ax.axis('off')
    ax.legend(handles, labels, ncols=2, loc="center left")

    return ax


def plot_conc_vs_time(concentration_data, time_data, compounds, dkie, file=None, writecsv=False, ax=None):
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
    from cycler import cycler
    colors = _create_extended_tab20()
    ax.set_prop_cycle(color=colors)
    #ax.set_prop_cycle(cycler(color=plt.cm.tab20.colors))

    # Create axes if not provided
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.8, 8.4))
        created_fig = True
        colors = _create_extended_tab20()
        ax.set_prop_cycle(color=colors)
        #ax.set_prop_cycle(color=plt.cm.tab20.colors)
    else:
        fig = ax.figure

    concentration_dataT = np.array(concentration_data).T
    idx = 0
    dictcomp_conc = {}
    for d, c in zip(compounds, concentration_dataT):
        dictcomp_conc[d] = c
        if max(c) < 1e-20:
            continue
        if idx < 30:
            plottype = '-'
            linewidth = 2
        elif 30 <= idx < 60:
            plottype = '--'
            linewidth = 2.5
        else:
            plottype = ':'
            linewidth = 2

        ax.plot(time_data, c, plottype, label=str(d), linewidth=linewidth)
        idx += 1

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Concentration")

    if writecsv is not False:
        safedata = []
        safedata.append(time_data)
        for key in compounds:
            safedata.append(dictcomp_conc[key])
        safedataT = np.array(safedata).T
        np.savetxt(writecsv, safedataT, delimiter=",")


    if writecsv is not False:
        safedata = []
        safedata.append(time_data)
    
        for key in dictcomp_conc:
            safedata.append(dictcomp_conc[key])
    
        safedataT = np.array(safedata).T
    
        header = "Time (s)\t" + "\t".join([str(o) for o in dictcomp_conc.keys()])
    
        np.savetxt(
            writecsv,
            safedataT,
            delimiter="\t",
            header=header,
            comments=""   # prevents "#" before header
        )

    if created_fig:
        fig.tight_layout()
        if file is not None:
            fig.savefig(file)
    
    return ax




