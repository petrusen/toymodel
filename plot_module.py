import itertools
import copy
import numpy as np

def _get_zero_safe_ratio(x18, x16, Rstd=0.002004):
    """
    Avoid dividing by zero.
    """
    den = sum(x16)
    num = sum(x18)
    if den > 0 and num > 0:
        ratio = num / den
        delta = ((ratio / Rstd) - 1) * 1000
    else:
        delta = 0
    return delta

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
        ax.plot(time_data, Rsub[key], label=key)

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
            Rsub[key].append(_get_zero_safe_ratio(a, b))

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

        for a, b, key in [
            (cpo4_18, cpo4_16, "CPO4"),
            (po4_18,  po4_16,  "PO4"),
            (cpo5_18, cpo5_16, "CPO5"),
            (o_18,    o_16,    "O"),
            (co_18,   co_16,   "CO"),
        ]:
            Rsub[key].append(_get_zero_safe_ratio(a, b))
    
    concentration_data_T = np.array(concentration_data).T
    reacind = compounds[(12, 16, 16, 16, 16)]
    reaction_progress = concentration_data_T[reacind]
    reaction_progress_norm = [r/reaction_progress[0] for r in reaction_progress]
    # --- plotting ---
    #ax.set_xscale("log")
    for key in Rsub:
        ax.plot(reaction_progress_norm, Rsub[key], label=key)

    ax.set_xlabel("Reaction progress respect to (12, 16, 16, 16, 16)")
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

def plot_mk_information(dkie, handles, labels, ax):
    """
    Placeholder information plot
    """
    # Create axes if not provided
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.8, 8.4))
        created_fig = True
        ax.set_prop_cycle(color=plt.cm.tab20.colors)
    else:
        fig = ax.figure
    string = _get_title(dkie, cntthresh=1)
    string2 = "MK information:\n\n" + string
    ax.text(0.8, 0.5, string2, transform=ax.transAxes, ha='left', va='center', fontsize=12)
    ax.axis('off')
    ax.legend(handles, labels, ncols=2, loc="center left")

    return ax


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
    from cycler import cycler
    ax.set_prop_cycle(cycler(color=plt.cm.tab20.colors))

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
    #titlestr = _get_title(dkie)
    #ax.set_title(titlestr)

    ##ax.legend(loc="lower center", ncol=2)

    if created_fig:
        fig.tight_layout()
        if file is not None:
            fig.savefig(file)
    
    return ax




