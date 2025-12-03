from helpfunctions import *

reactions1 = create_reactions_for_step_one(verbose=False)
reactions2 = create_reactions_for_step_two(verbose=False)
initial_conc = {(16,16,16,16,0): 0.9, (16,16,16,18,0): 0.1, (16,): 0.1, (18,): 0.9}
total_reac = reactions1 + reactions2 
convert_to_kinetx_notation(total_reac, initial_conc)
