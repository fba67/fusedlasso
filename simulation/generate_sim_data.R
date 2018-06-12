set.seed(1); 
N1 <- 300
N2 <- 300
N <- N1 + N2

if(F){
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
print("simulatioin 1")
d[seq(N1), c(q[[1]], q[[2]], q[[3]])] <- 0
d[seq(N1), c(q[[4]])] <- rnorm(10 * N1, mean= 10, sd= .01)

d[seq((N2 + 1), N), c(q[[1]], q[[2]], q[[3]])] <- rnorm(30 * N2, mean= 10, sd= .01)
d[seq((N2 + 1), N), c(q[[4]])] <- 0
# y <- c(rnorm(N1, mean= 10, sd= .01), rnorm(N2, mean= 0, sd= .01))
y <- rowSums(d[, q[[3]]]) * 2 + rowSums(d[, q[[4]]]) * 3 + 4
y[seq((N2 + 1), N)] <- -y[seq((N2 + 1), N)]
write.table(d, "fl_sim_features_1.txt", quote= F, row.names= F, col.names= F)
writeLines(as.character(y), "fl_sim_response_1.txt")

pdf("simulation_1.pdf")
par(mfrow= c(2,2))
plot(colMeans(d[seq(300), ]), main= "X", ylab= "avg. X[1:300, ]", pch= 20)
plot(y[seq(300)], main= "Y", pch= 20)
plot(colMeans(d[seq(301, 600), ]), main= "X", ylab= "avg. X[301:600, ]", pch= 20)
plot(y[seq(301, 600)], main= "Y", pch= 20)
dev.off()

rm(list= c("d", "y"))
####################################################
####################################################
print("simulatioin 2")
#set.seed(1); 
d <- matrix(NA, ncol= 40, nrow= N)
d[seq(N1), c(q[[1]], q[[2]], q[[3]])] <- 0
d[seq(N1), c(q[[1]])] <- 1
d[seq(N1), c(q[[4]])] <- rnorm(10 * N1, mean= 10, sd= .01)

d[seq((N2 + 1), N), c(q[[2]], q[[3]])] <- rnorm(20 * N2, mean= 10, sd= .01)
d[seq((N2 + 1), N), q[[1]]] <- 0
d[seq((N2 + 1), N), c(q[[4]])] <- 0
# y <- c(rnorm(N1, mean= 10, sd= .01), rnorm(N2, mean= 0, sd= .01))
y <- rowSums(d[, c(q[[1]], q[[4]])]) * 3 + 4
y[seq((N2 + 1), N)] <- -rowSums(d[seq((N2 + 1), N), c(q[[2]], q[[3]])]) * 3 + 4

write.table(d, "fl_sim_features_2.txt", quote= F, row.names= F, col.names= F)
writeLines(as.character(y), "fl_sim_response_2.txt")

pdf("simulation_2.pdf")
par(mfrow= c(2,2))
plot(colMeans(d[seq(300), ]), main= "X", ylab= "avg. X[1:300, ]", pch= 20)
plot(y[seq(300)], main= "Y", pch= 20)
plot(colMeans(d[seq(301, 600), ]), main= "X", ylab= "avg. X[301:600, ]", pch= 20)
plot(y[seq(301, 600)], main= "Y", pch= 20)
dev.off()

rm(list= c("d", "y"))
####################################################
####################################################
print("simulatioin 3")
#set.seed(1); 
d <- matrix(NA, ncol= 40, nrow= N)
d[seq(N), c(q[[1]], q[[2]], q[[3]])] <- 0 + rnorm(30 * N, mean= 0, sd= .01)
d[seq(N), c(q[[1]])] <- rnorm(10 * N, mean= 5, sd= .01)
#d[seq(N), c(q[[1]])] <- 1 + rnorm(10, 0, .01)
d[seq(N), c(q[[4]])] <- rnorm(10 * N, mean= 10, sd= .01)

y <- rowSums(d[, q[[1]][1:5]]) * 2 +  rowSums(d[, q[[4]]]) * 3 + 4 + rnorm(N, mean= 0, sd= .01)

write.table(d, "fl_sim_features_3.txt", quote= F, row.names= F, col.names= F)
writeLines(as.character(y), "fl_sim_response_3.txt")

pdf("simulation_3.pdf")
par(mfrow= c(2,2))
plot(colMeans(d[seq(300), ]), main= "X", ylab= "avg. X[1:300, ]", pch= 20)
plot(y[seq(300)], main= "Y", pch= 20)
plot(colMeans(d[seq(301, 600), ]), main= "X", ylab= "avg. X[301:600, ]", pch= 20)
plot(y[seq(301, 600)], main= "Y", pch= 20)
dev.off()

####################################################
####################################################
print("simulatioin 4")
#set.seed(1); 
d <- matrix(NA, ncol= 40, nrow= N)

segment.size <- 4
### Generate indices for segments of size segment.size
for(i in seq(10))
  q[[i]] <- seq((i-1) * segment.size + 1, i * segment.size)

### Assign ~0 to every other segment
d[, seq(20)] <- rnorm(N * 20, mean= 10, sd= 0.01)
d[, seq(21, 40)] <- rnorm(N * 20, mean= -10, sd= 0.01)
for(i in seq(1, 10, by= 2)){
  d[seq(N), q[[i]]] <- rnorm(segment.size * N, mean= 0, sd= .01)
}


y <- rowSums(d[, q[[1]]]) * 2 + rowSums(d[, q[[2]]]) * 2 + rowSums(d[, q[[4]]]) * 2 +
  rowSums(d[, q[[7]]]) * 3 + rowSums(d[q[[8]],]) * 3 + rowSums(d[q[[9]],]) * 3 + rnorm(N, mean= 0, sd= .01)

write.table(d, "fl_sim_features_4.txt", quote= F, row.names= F, col.names= F)
writeLines(as.character(y), "fl_sim_response_4.txt")

pdf("simulation_4.pdf")
par(mfrow= c(2,2))
plot(colMeans(d[seq(300), ]), main= "X", ylab= "avg. X[1:300, ]", pch= 20)
plot(y[seq(300)], main= "Y", pch= 20)
plot(colMeans(d[seq(301, 600), ]), main= "X", ylab= "avg. X[301:600, ]", pch= 20)
plot(y[seq(301, 600)], main= "Y", pch= 20)
dev.off()

####################################################
####################################################
print("simulatioin 5")
#set.seed(1); 
d <- matrix(NA, ncol= 40, nrow= N)

segment.size <- 4
### Generate indices for segments of size segment.size
for(i in seq(10))
  q[[i]] <- seq((i-1) * segment.size + 1, i * segment.size)

### Assign ~0 to every other segment
d[, seq(20)] <- rnorm(N * 20, mean= 10, sd= 0.01)
d[, seq(21, 40)] <- rnorm(N * 20, mean= -10, sd= 0.01)
for(i in seq(1, 10, by= 2)){
  d[seq(N), q[[i]]] <- rnorm(segment.size * N, mean= 0, sd= .01)
}


y <- rowSums(d[, q[[2]]]) * 2 + rowSums(d[, q[[4]]]) * 2 +
 rnorm(N, mean= 0, sd= .01)

write.table(d, "fl_sim_features_5.txt", quote= F, row.names= F, col.names= F)
writeLines(as.character(y), "fl_sim_response_5.txt")

pdf("simulation_5.pdf")
par(mfrow= c(2,2))
plot(colMeans(d[seq(300), ]), main= "X", ylab= "avg. X[1:300, ]", pch= 20)
plot(y[seq(300)], main= "Y", pch= 20)
plot(colMeans(d[seq(301, 600), ]), main= "X", ylab= "avg. X[301:600, ]", pch= 20)
plot(y[seq(301, 600)], main= "Y", pch= 20)
dev.off()

####################################################
####################################################
print("simulatioin 6")
#set.seed(1); 
d <- matrix(NA, ncol= 40, nrow= N)

segment.size <- 2
### Generate indices for segments of size segment.size
for(i in seq(20))
  q[[i]] <- seq((i-1) * segment.size + 1, i * segment.size)

### Assign ~0 to every other segment
d[, seq(20)] <- rnorm(N * 20, mean= 10, sd= 0.01)
d[, seq(21, 40)] <- rnorm(N * 20, mean= -5, sd= 0.01)

d[, c(seq(2, 13), seq(15, 20))] <- rnorm(N * length(c(seq(2, 13), seq(15, 20))), mean= 0, sd= .01)
d[, seq(22, 38)] <- rnorm(N * length(seq(22, 38)), mean= 0, sd= .01)
data.x <- d
print(dim(d))

y <- data.x[, 1] * 2 + (data.x[, 14]) * 2 + rowSums(data.x[, seq(21, 22)]) * 2 + rowSums(data.x[, seq(39, 40)]) * 2 + 
 rnorm(N, mean= 0, sd= .01)

write.table(d, "fl_sim_features_6.txt", quote= F, row.names= F, col.names= F)
writeLines(as.character(y), "fl_sim_response_6.txt")

pdf("simulation_6.pdf")
par(mfrow= c(2,2))
plot(colMeans(d[seq(300), ]), main= "X", ylab= "avg. X[1:300, ]", pch= 20)
plot(y[seq(300)], main= "Y", pch= 20)
plot(colMeans(d[seq(301, 600), ]), main= "X", ylab= "avg. X[301:600, ]", pch= 20)
plot(y[seq(301, 600)], main= "Y", pch= 20)
dev.off()
}
####################################################
####################################################
print("simulatioin 7")
load("../DEEP_allDS_divergenet_plus_newBPs/fl_01_HepG2_newBPs_plus_divergent_mRNA.RData")
data.mean <- colMeans(fl$cv.fl$bestobj$X[, c(seq(42, 121))]) # This subsets the HMs to only K4m3 and K27me3
#data.sd <- apply(fl$cv.fl$bestobj$X[, c(seq(42, 121))], 2, FUN= sd)
data.cov1 <- cov(fl$cv.fl$bestobj$X[, c(seq(42, 81))])
data.cov2 <- cov(fl$cv.fl$bestobj$X[, c(seq(82, 121))])
#set.seed(1); 
library(MASS)

### Generate indices for segments of size segment.size
d1 <- matrix(mvrnorm(N * 40, mu= data.mean[seq(40)], Sigma= data.cov1), ncol= 40, nrow= N, byrow=T)
d2 <- matrix(mvrnorm(N * 40, mu= data.mean[seq(41, 80)], Sigma= data.cov1), ncol= 40, nrow= N, byrow=T)
d <- cbind(d1, d2)
### Assign ~0 to every other segment
data.x <- d
print(dim(d))

y <- rowSums(data.x[, seq(15, 30)]) - rowSums(data.x[, seq(41, 45)]) - rowSums(data.x[, seq(50, 80)]) + 
 rnorm(N, mean= 0, sd= .01)

write.table(d, "fl_sim_features_7.txt", quote= F, row.names= F, col.names= F)
writeLines(as.character(y), "fl_sim_response_7.txt")

pdf("simulation_7.pdf")
par(mfrow= c(2,2))
plot(colMeans(d[seq(300), ]), main= "X", ylab= "avg. X[1:300, ]", pch= 20)
plot(y[seq(300)], main= "Y", pch= 20)
plot(colMeans(d[seq(301, 600), ]), main= "X", ylab= "avg. X[301:600, ]", pch= 20)
plot(y[seq(301, 600)], main= "Y", pch= 20)
dev.off()
