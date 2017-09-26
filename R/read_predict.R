
topic_clus <- get(load("../utilities/moderns_Pinhasi/clus_2/model.rda"))

omega <- topic_clus$omega
theta <- topic_clus$theta

new_count <- matrix(0, ncol = dim(theta)[1], nrow = 5)
new_count[1, c(1, 100)] <- 1
new_count[2, c(1001)] <- 1
new_count[3, c(5000, 9000)] <- 1
new_count[4, c(1933)] <- 1

rownames(new_count) <- paste0("reads-", 1:5)
out <- read_memberships(topic_clus, new_count, method = "lik")
out <- read_memberships(topic_clus, new_count, method = "map")
out <- read_memberships(topic_clus, new_count, method = "independent")
out <- read_memberships(topic_clus, new_count, method = "independent-nostrand")

fit <- topic_clus
reads_data <- new_count

