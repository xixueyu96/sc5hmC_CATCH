library(ClusterGVis)

# load data
data(exps)

# check
head(exps,3)

getClusters(exp = exps)

cm <- clusterData(exp = exps,
                  cluster.method = "mfuzz",
                  cluster.num = 8)
