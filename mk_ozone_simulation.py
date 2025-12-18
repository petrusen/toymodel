from helpfunctions import *

nodesfile = "/home/petrusen22/csv_files/nodes_dez2_buo_g250_vinfinfinf_bi.crn"
edgesfile = "/home/petrusen22/csv_files/edges_dez2_buo_g250_vinfinfinf_bi.crn"
labelfile = "/home/petrusen22/csv_files/labels_dez2_buo_g250_vinfinfinf_bi.crn"

#nodesfile = "/home/petrusen22/csv_files/nodes_test.crn"
#edgesfile = "/home/petrusen22/csv_files/edges_test.crn"
#labelfile = "/home/petrusen22/csv_files/labels_test.crn"

concentration = calculate_mk_yields(nodesfile, edgesfile, labelfile)

