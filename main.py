# Third party libraries
import matplotlib.pyplot as plt       
import numpy as np
import yaml

# Local libraries
from mk_module import (get_water_initial_concentrations, 
                       get_alcohol_initial_concentrations,
                       get_organophosphate_initial_concentrations,
                       get_phosphate_initial_concentrations,
                       create_reactions_for_step_one,
                       create_reactions_for_step_two,
                       calculate_concentrations_from_mk)
from plot_module import (plot_conc_vs_time, 
                         plot_ratio_vs_reaction_progress, 
                         plot_ratio_vs_time, 
                         plot_mk_information,
                         plot_ratio_vs_time_separate)

def main():

     # Load YAML
     with open("config.yaml", "r") as f:
         config = yaml.safe_load(f)
     
     # Extract values
     dkie = config["dkie"]
     dplot = config["plot"]
     
     # Recreate NumPy array
     tr = config["timerange"]
     timerange = np.logspace(tr["start"], tr["stop"], num=tr["num"])
     
     # Get initial concentrations
     c0_o = get_water_initial_concentrations(dkie["dh2o"], initial_C=dkie["conc_h2o"])
     c0_cpo = get_organophosphate_initial_concentrations(dkie["dcpo"], initial_C= dkie["conc_cpo"])
     c0_po = get_phosphate_initial_concentrations(dkie["dpo"], initial_C=dkie["conc_po"])
     c0_co = get_alcohol_initial_concentrations(dkie["dco"], initial_C=dkie["conc_co"])
     initial_conc = {}
     for dicti in [c0_o, c0_cpo, c0_po, c0_co]:
         initial_conc.update(dicti)
     
     # Get reactions for the two steps
     reactions1 = create_reactions_for_step_one(dkie, filtreactions=dkie["filtreactions"])
     reactions2 = create_reactions_for_step_two(dkie, reactions1)
     total_reac = reactions1 + reactions2
     print("TOTAL REACTIONS", len(total_reac), total_reac)
     # Set up and solve ODE equations
     simdata = calculate_concentrations_from_mk(total_reac, initial_conc, timerange)
     conc_data, t_data, compounds = simdata

     # Plot 5 separate figures for each dO vs time
     #plot_ratio_vs_time_separate(conc_data, t_data, compounds, dkie, file=outfile)
     
     # Plot four panel figure
     fig, axes = plt.subplots(2, 2, figsize=(14, 9))
     ax1, ax2, ax3, ax4 = axes[0, 0], axes[1, 0], axes[0, 1], axes[1, 1]
     plot_conc_vs_time(conc_data, t_data, compounds, dkie, writecsv=dplot["C_time"], ax=ax1)
     plot_ratio_vs_reaction_progress(conc_data, t_data, compounds, dkie, ax=ax4, 
                                     reactionref=dplot["rxnref"], writecsv=dplot["dO_rxnprogress"])
     plot_ratio_vs_time(conc_data, t_data, compounds, dkie, writecsv=dplot["dO_time"], ax=ax2)
     handles, labels = ax1.get_legend_handles_labels()
     plot_mk_information(dkie, handles, labels, ax=ax3)
     plt.tight_layout()
     plt.savefig(dplot["4plotlayout"])
     plt.show()

main()
