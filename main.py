from helpfunctions import *
import matplotlib.pyplot as plt

def main():
     ##### INPUT DEFINITION #####
     dkie = {"k1_16": 1,
             "k2_16": 1,
             "k3_16": 1,
             "kie1p": 1.03,
             "kie1s": 1,
             "kie2p": 1,
             "kie2s": 1,
             "kie3p": 1,
             "kie3s": 1,
             "penalty_k2": 1,
             "penalty_k4": 1,
             "dh2o": 20,
             "dcpo": 10,
             "conc_h2o": 55.5,
             "conc_cpo": 1e-3}
     outfile = "mk_scenario_1_2c.png"
     timerange = np.logspace(-8, -2, num=100)
     delta_O_h2o = dkie["dh2o"]
     delta_O_cpo = dkie["dcpo"]
     ############################
     ch2o, ccpo = dkie["conc_h2o"], dkie["conc_cpo"]
     initial_conc = get_water_initial_concentrations(delta_O_h2o, initial_C=ch2o)
     initial_conc.update(get_organophosphate_initial_concentrations(delta_O_cpo, initial_C=ccpo))
     print("INITIALCONC", initial_conc)
     reactions1 = create_reactions_for_step_one(dkie, verbose=False)
     reactions2 = create_reactions_for_step_two(dkie, verbose=False)
     total_reac = reactions1 + reactions2
     
     for r in total_reac:
         print(r)
 
    
     fig, axes = plt.subplots(
         2, 2,
         figsize=(14, 9),
     )

     concentration_data, time_data, compounds = convert_to_kinetx_notation(total_reac, initial_conc, timerange)
     plot_conc_vs_time(concentration_data, time_data, compounds, dkie, ax=axes[0,0])
     plot_ratio_vs_reaction_progress(concentration_data, time_data, compounds, dkie, ax=axes[1,1])
     plot_ratio_vs_time(concentration_data, time_data, compounds, dkie, ax=axes[1,0])
     handles, labels = axes[0,0].get_legend_handles_labels()
     plot_mk_information(dkie, handles, labels, ax=axes[0,1])
     #plot_concentrations(concentration_data, compounds, dkie, file3)
     plt.tight_layout()
     plt.savefig(outfile)
     plt.show()
main()
