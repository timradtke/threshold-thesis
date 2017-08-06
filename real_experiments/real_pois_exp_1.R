#sessions_threshold <- 0
set.seed(35968399)
selected_products <- sample(1:196, 20)
niter <- 5000
nrounds <- 7580 # These two values give the maximum size data 
nahead <- 7580  # possible given the data set size we have
source("/Users/timradtke/Desktop/thesis_data/get_count_dataset.R")

head(pv_list[[1]])
colMeans(pv_list[[1]])
hist(pv_list[[1]][,8])
colMeans(pv_products_wide)
