from helpfunctions import *

file1 = "deleteme1.png"
file2 = "deleteme2.png"
kie = 1
reactions1 = create_reactions_for_step_one(kie, verbose=False)
reactions2 = create_reactions_for_step_two(kie, verbose=False)
concr1 = 1
concr2 = 2
p18 = 0.002
p16 = 0.998
initial_conc = {(12,16,16,16,16): concr1*p16, 
                (12,16,16,16,18): concr1*p18, 
                (16,): concr2*p16,
                (18,): concr2*p18}
total_reac = reactions1 + reactions2

for r in total_reac:
    print(r)

concentration_data, time_data, compounds = convert_to_kinetx_notation(total_reac, initial_conc)
#plot_ratio_vs_time(concentration_data, time_data, compounds, file1)
plot_conc_vs_time(concentration_data, time_data, compounds, file2)
plot_concentrations(concentration_data, compounds)
