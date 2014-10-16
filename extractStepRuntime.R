# extracting runtime

setwd("~/Dropbox/lsgm/scalability/")

corr3_projection = read.table("./runtime_corr3/projection_corr3.txt")
corr3_procrus = read.table("./runtime_corr3/procrusties_corr3.txt")
corr3_clust = read.table("./runtime_corr3/clustering_corr3.txt")
corr3_match = read.table("./runtime_corr3/matching_corr3.txt")

corr3_run_per_step = data.frame(cbind(corr3_projection, corr3_procrus, corr3_clust, corr3_match))
names(corr3_run_per_step) = c('embedding', 'procrusties', 'clustering', 'matching')

for (i in 1:4){
    ind = 1:nrow(corr3_clust)
    extract_ind = ind[seq(i, length(ind)+i-1, 4)]    
    print(colMeans(corr3_run_per_step[extract_ind,]))
    print(apply(corr3_run_per_step[extract_ind,], 2, sd))
}


projection = read.table("./runtime_corr3/projection_corr3.txt")
procrus = read.table("./runtime_corr3/procrusties_corr3.txt")
clust = read.table("./runtime_corr3/clustering_corr3.txt")
match = read.table("./runtime_corr3/matching_corr3.txt")

run_per_step = data.frame(cbind(projection, procrus, clust, match))
names(run_per_step) = c('Embedding', 'Procrusties', 'Clustering', 'Matching')

for (i in 1:4){
    ind = 1:nrow(clust)
    extract_ind = ind[seq(i, length(ind)+i-1, 4)]    
    print(colMeans(run_per_step[extract_ind,]))
    print(apply(run_per_step[extract_ind,], 2, sd))
}




