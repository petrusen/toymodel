# Third party libraries
import matplotlib.pyplot as plt       
import numpy as np
import yaml

# Local libraries
from mk_module import (get_water_initial_concentrations, 
                       get_organophosphate_initial_concentrations,
                       create_reactions_for_step_one,
                       create_reactions_for_step_two,
                       calculate_concentrations_from_mk)
from plot_module import (plot_conc_vs_time, 
                         plot_ratio_vs_reaction_progress, 
                         plot_ratio_vs_time, 
                         plot_mk_information)

def main():

     # Load YAML
     with open("config.yaml", "r") as f:
         config = yaml.safe_load(f)
     
     # Extract values
     dkie = config["dkie"]
     outfile = config["outfile"]
     
     # Recreate NumPy array
     tr = config["timerange"]
     timerange = np.logspace(tr["start"], tr["stop"], num=tr["num"])
     
     # Define initial parameters
     delta_O_h2o = dkie["dh2o"]
     delta_O_cpo = dkie["dcpo"]
     ch2o, ccpo = dkie["conc_h2o"], dkie["conc_cpo"]

     # Get initial concentrations
     initial_conc = get_water_initial_concentrations(delta_O_h2o, initial_C=ch2o)
     tmp = get_organophosphate_initial_concentrations(delta_O_cpo, initial_C=ccpo)
     initial_conc.update(tmp)

     # Get reactions for the two steps
     reactions1 = create_reactions_for_step_one(dkie, verbose=False)
     reactions2 = create_reactions_for_step_two(dkie, verbose=False)
     total_reac = reactions1 + reactions2
     
     # Set up and solve ODE equations
     simdata = calculate_concentrations_from_mk(total_reac, initial_conc, timerange)
     conc_data, t_data, compounds = simdata

     # Plot results
     fig, axes = plt.subplots(2, 2, figsize=(14, 9))
     ax1, ax2, ax3, ax4 = axes[0,0], [1,0], [0,1], [1,1]
     plot_conc_vs_time(conc_data, t_data, compounds, dkie, ax=ax1)
     plot_ratio_vs_reaction_progress(conc_data, t_data, compounds, dkie, ax=ax4)
     plot_ratio_vs_time(conc_data, t_data, compounds, dkie, ax=ax2)
     handles, labels = ax1.get_legend_handles_labels()
     plot_mk_information(dkie, handles, labels, ax=ax3)
     plt.tight_layout()
     plt.savefig(outfile)
     plt.show()

main()
