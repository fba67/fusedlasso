set.seed(1); 
N <- 600
N1 <- 300
N2 <- 300

d <- matrix(NA, ncol= 40, nrow= N)

q <- list(); q[[1]] <- seq(10); q[[2]] <- seq(11,20); q[[3]] <- seq(21,30); q[[4]] <- seq(31,40)

if(F){
d[seq(N1), c(q[[1]], q[[2]])] <- rnorm(20 * N1, mean= 0, sd= .01)
d[seq(N1), c(q[[3]], q[[4]])] <- rnorm(20 * N1, mean= 10, sd= .01)

d[seq((N2 + 1), N), c(q[[1]], q[[2]])] <- rnorm(20 * N2, mean= 10, sd= .01)
d[seq((N2 + 1), N), c(q[[3]], q[[4]])] <- rnorm(20 * N2, mean= 0, sd= .01)
}


####################################################
####################################################
d[seq(N1), c(q[[1]], q[[2]], q[[3]])] <- 0
d[seq(N1), c(q[[4]])] <- rnorm(10 * N1, mean= 10, sd= .01)

d[seq((N2 + 1), N), c(q[[1]], q[[2]], q[[3]])] <- rnorm(30 * N2, mean= 10, sd= .01)
d[seq((N2 + 1), N), c(q[[4]])] <- 0
# y <- c(rnorm(N1, mean= 10, sd= .01), rnorm(N2, mean= 0, sd= .01))
y <- rowSums(d[, c(q[[3]], q[[4]])]) * 3 + 4
y[seq((N2 + 1), N)] <- -y[seq((N2 + 1), N)]
write.table(d, "fl_sim_features_1.txt", quote= F, row.names= F, col.names= F)
writeLines(as.character(y), "fl_sim_response_1.txt")

rm(list= c("d", "y"))
####################################################
####################################################
d <- matrix(NA, ncol= 40, nrow= N)
d[seq(N1), c(q[[1]], q[[2]], q[[3]])] <- 0
d[seq(N1), c(q[[1]])] <- 1
d[seq(N1), c(q[[4]])] <- rnorm(10 * N1, mean= 10, sd= .01)

d[seq((N2 + 1), N), c(q[[2]], q[[3]])] <- rnorm(20 * N2, mean= 10, sd= .01)
d[seq((N2 + 1), N), q[[1]]] <- 0
d[seq((N2 + 1), N), c(q[[4]])] <- 0
# y <- c(rnorm(N1, mean= 10, sd= .01), rnorm(N2, mean= 0, sd= .01))
y <- rowSums(d[, c(q[[3]], q[[4]])]) * 3 + 4
y[seq((N2 + 1), N)] <- -y[seq((N2 + 1), N)]
write.table(d, "fl_sim_features_2.txt", quote= F, row.names= F, col.names= F)
writeLines(as.character(y), "fl_sim_response_2.txt")

rm(list= c("d", "y"))
####################################################
####################################################
d <- matrix(NA, ncol= 40, nrow= N)
d[seq(N), c(q[[1]], q[[2]], q[[3]])] <- 0 + rnorm(10, 0, .01)
d[seq(N), c(q[[1]])] <- 1 + rnorm(10, 0, .01)
d[seq(N), c(q[[4]])] <- rnorm(10 * N, mean= 10, sd= .01)

y <- rowSums(d[, c(q[[3]], q[[4]])]) * 3 + 4

write.table(d, "fl_sim_features_3.txt", quote= F, row.names= F, col.names= F)
writeLines(as.character(y), "fl_sim_response_3.txt")
