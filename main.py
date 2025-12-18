# third party libraries
import matplotlib.pyplot as plt       
import numpy as np
import yaml

# local libraries
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
     ############### INPUT DEFINITION ################
     # Load YAML
     with open("config.yaml", "r") as f:
         config = yaml.safe_load(f)
     
     # Extract values
     dkie = config["dkie"]
     outfile = config["outfile"]
     
     # Recreate NumPy array
     tr = config["timerange"]
     timerange = np.logspace(tr["start"], tr["stop"], num=tr["num"])
     #################################################
     
     delta_O_h2o = dkie["dh2o"]
     delta_O_cpo = dkie["dcpo"]
     ch2o, ccpo = dkie["conc_h2o"], dkie["conc_cpo"]
     initial_conc = get_water_initial_concentrations(delta_O_h2o, initial_C=ch2o)
     initial_conc.update(get_organophosphate_initial_concentrations(delta_O_cpo, initial_C=ccpo))
     reactions1 = create_reactions_for_step_one(dkie, verbose=False)
     reactions2 = create_reactions_for_step_two(dkie, verbose=False)
     total_reac = reactions1 + reactions2
     
     for r in total_reac:
         print(r)
 
     fig, axes = plt.subplots(
         2, 2,
         figsize=(14, 9),
     )

     concentration_data, time_data, compounds = calculate_concentrations_from_mk(total_reac, initial_conc, timerange)
     plot_conc_vs_time(concentration_data, time_data, compounds, dkie, ax=axes[0,0])
     plot_ratio_vs_reaction_progress(concentration_data, time_data, compounds, dkie, ax=axes[1,1])
     plot_ratio_vs_time(concentration_data, time_data, compounds, dkie, ax=axes[1,0])
     handles, labels = axes[0,0].get_legend_handles_labels()
     plot_mk_information(dkie, handles, labels, ax=axes[0,1])
     plt.tight_layout()
     plt.savefig(outfile)
     plt.show()

main()
