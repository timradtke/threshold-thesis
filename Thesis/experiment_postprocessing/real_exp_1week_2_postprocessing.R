# Post processing of data for arm sequence visualization.
# This is for Experiment 2 on real data as used in thesis.

##########################################################

# load the algorithm results from Experiment 2

library(plyr)
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/"


