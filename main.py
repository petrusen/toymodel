from helpfunctions import *

def main():
     ##### INPUT DEFINITION #####
     dkie = {"k1_16": 1,
             "k2_16": 1,
             "k3_16": 1,
             "kie1p": 1,
             "kie1s": 1,
             "kie2p": 1,
             "kie2s": 1,
             "kie3p": 1.005,
             "kie3s": 1.0005,
             }
     
     file1 = "ratio_vs_time.png"
     file2 = "conc_vs_time.png"
     file3 = "final_yields.png"
     timerange = np.logspace(-8, 5, num=100)
     penalty = 1e-10
     delta_O_h2o = 20
     delta_O_cpo = 10
     ############################


     initial_conc = get_water_initial_concentrations(delta_O_h2o)
     initial_conc.update(get_organophosphate_initial_concentrations(delta_O_cpo))
     reactions1 = create_reactions_for_step_one(dkie, penalty, verbose=False)
     reactions2 = create_reactions_for_step_two(dkie, penalty, verbose=False)
     total_reac = reactions1 + reactions2
     
     for r in total_reac:
         print(r)
     
     concentration_data, time_data, compounds = convert_to_kinetx_notation(total_reac, initial_conc, timerange)
     plot_ratio_vs_time(concentration_data, time_data, compounds, dkie, file1)
     plot_conc_vs_time(concentration_data, time_data, compounds, dkie, file2)
     plot_concentrations(concentration_data, compounds, dkie, file3)

main()
