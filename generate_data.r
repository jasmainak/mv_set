source("distributions.R")
rpareto_inv <- function(n, theta, a) {
    u <- runif(n, min = 0, max = 1);
    return(theta / (u ^ (1 / a)));
}


rpareto_cond_inv <- function(x2, theta1, theta2, a) {
  u <- runif(length(x2), min = 0, max = 1);
  return(theta1 + theta1 / theta2 * x2 * (1 / (u ^ (1 / (a + 1))) - 1));
}

theta1 <- 5;   # Location parameter 1
theta2 <- 2;   # Location parameter 2

simulation.mc=function(n,r,rho,theta1,theta2,a){
x <- rmvexp(n, r, rho)
x2 <- rpareto_inv(n, theta = theta2, a = a);
x1 <- rpareto_cond_inv(x2, theta1, theta2, a);
df <- cbind.data.frame(x1 = x1, x2 = x2);
join=x-df

my_int<-n
nr<-as.integer(my_int)
mymat<-matrix(0,nr,2)
first=mymat[,1]
 for(i in 1:n){first[i+1]<-max(first[i]+join[i,1],0)}
first1=mymat[,2]
for(i in 1:n){first1[i+1]<-max(first1[i]+join[i,2],0)}

reg=cbind(first,first1)
w=sum(reg[,1] == 0.00000 & reg[,2] == 0.00000)
s=n/w
return(reg)}

r <- 0.75:0.8
rho <- matrix(c(1, 0.8, 0.8, 1), ncol=2)
a <- 3;        # Shape parameter
for (i in 1:5){
	n <- 10^i;     # Number of samples
	data = simulation.mc(n,r,rho,theta1,theta2,a)
	write.table(data, sprintf('regen1_%d.txt', n), sep='\t', row.names = FALSE, col.names = FALSE)
}

r <- 0.35:0.4
rho <- matrix(c(1, 0.8, 0.8, 1), ncol=2)
a <- 3.2;        # Shape parameter
for (i in 1:5){
	n <- 10^i;     # Number of samples
	data = simulation.mc(n,r,rho,theta1,theta2,a)
	write.table(data, sprintf('regen2_%d.txt', n), sep='\t', row.names = FALSE, col.names = FALSE)
}

r <- 0.25:0.4
rho <- matrix(c(1, 0.5, 0.5, 1), ncol=2)
a <- 3;        # Shape parameter
for (i in 1:5){
	n <- 10^i;     # Number of samples
	data = simulation.mc(n,r,rho,theta1,theta2,a)
	write.table(data, sprintf('regen3_%d.txt', n), sep='\t', row.names = FALSE, col.names = FALSE)
}


