from helpfunctions import *

file1 = "deleteme1.png"
file2 = "deleteme3.png"

reactions1 = create_reactions_for_step_one(verbose=False)
reactions2 = create_reactions_for_step_two(verbose=False)
concr1 = 1
concr2 = 0.8
p18 = 0.002
p16 = 0.9998
initial_conc = {(16,16,16,16): concr1*p16, 
                (16,16,16,18): concr1*p18, 
                (16,): concr2*p16,
                (18,): concr2*p18}
total_reac = reactions1

print(total_reac)
concentration_data, time_data, compounds = convert_to_kinetx_notation(total_reac, initial_conc)
#plot_ratio_vs_time(concentration_data, time_data, compounds, file1)
plot_conc_vs_time(concentration_data, time_data, compounds, file2)
