#source('H:/Advising/Murphy_Shannon/Subgroup Detection/Kaizer_MEM_Weights.R')
source("C:/Users/mushanno/Desktop/Work/Dissertation1/CodeFromAlex/Kaizer_MEM_Weights.R")
library(pROC)


calc.weights_MEM(xvec=c(1,1.1), svec=c(1,1), nvec=c(100,100), prior='pi_e')
calc.weights_MEM(xvec=c(1,1.1), svec=c(1,1), nvec=c(100,100), prior='pi_n')

# sqrt(s1^2/n1 + s2^2/n2)

# Z=-2, setting 1
x1 <- 0
x2 <- 0.4
s1 <- s2 <- sqrt(2)
n1 <- n2 <- 100

z1 <- (x1-x2)/sqrt(s1^2/n1 + s2^2/n2); z1

e1 <- calc.weights_MEM(xvec=c(x1,x2), svec=c(s1,s2), nvec=c(n1,n2), prior='pi_e')
f1 <- calc.weights_MEM(xvec=c(x1,x2), svec=c(s1,s2), nvec=c(n1,n2), prior='pi_n')


# Z=-2, setting 2 (change x2, s1/s2)
x1 <- 0
x2 <- 4/5
s1 <- s2 <- sqrt(8)
n1 <- n2 <- 100

z2 <- (x1-x2)/sqrt(s1^2/n1 + s2^2/n2); z2

e2 <- calc.weights_MEM(xvec=c(x1,x2), svec=c(s1,s2), nvec=c(n1,n2), prior='pi_e')
f2 <- calc.weights_MEM(xvec=c(x1,x2), svec=c(s1,s2), nvec=c(n1,n2), prior='pi_n')


# Z=-2, setting 3 (change n1, n2 to larger; x2)
x1 <- 0
x2 <- 2*sqrt(2/7)/5
s1 <- s2 <- sqrt(2)
n1 <- n2 <- 350

z3 <- (x1-x2)/sqrt(s1^2/n1 + s2^2/n2); z3

e3 <- calc.weights_MEM(xvec=c(x1,x2), svec=c(s1,s2), nvec=c(n1,n2), prior='pi_e')
f3 <- calc.weights_MEM(xvec=c(x1,x2), svec=c(s1,s2), nvec=c(n1,n2), prior='pi_n')


# Z=-2, setting 4 (change n1, n2 to smaller; x2)
x1 <- 0
x2 <- 2*sqrt(2)/5
s1 <- s2 <- sqrt(2)
n1 <- n2 <- 50

z4 <- (x1-x2)/sqrt(s1^2/n1 + s2^2/n2); z4

e4 <- calc.weights_MEM(xvec=c(x1,x2), svec=c(s1,s2), nvec=c(n1,n2), prior='pi_e')
f4 <- calc.weights_MEM(xvec=c(x1,x2), svec=c(s1,s2), nvec=c(n1,n2), prior='pi_n')


# Z=-2, setting 5 (change n1, n2 to smaller; s1, s2 to bigger; x2)
x1 <- 0
x2 <- 4*sqrt(2)/5
s1 <- s2 <- sqrt(8)
n1 <- n2 <- 50

z5 <- (x1-x2)/sqrt(s1^2/n1 + s2^2/n2); z4

e5 <- calc.weights_MEM(xvec=c(x1,x2), svec=c(s1,s2), nvec=c(n1,n2), prior='pi_e')
f5 <- calc.weights_MEM(xvec=c(x1,x2), svec=c(s1,s2), nvec=c(n1,n2), prior='pi_n')


# Results indicate that pi_e changes with both x1-x2, s1, s2, n1, n2; pi_n is constant for given n1/n2
c(z1,z2,z3,z4,z5)
rbind(e1,e2,e3,e4,e5)
rbind(f1,f2,f3,f4,f5)

# Questions: can we find connection between pi_e and pi_n to use pi_e but apply correction
# Questions: or just use pi_n since it is constant for a given effect size
# Questions: could then incorporate n correction since large sample size does affect estimates
# Questions: impact of imbalance (n1/n2 versus n2/n1 matters?)
# Questions: use both pi_e and pi_n

# Formula to estimate needed x2 for any given n1, n2, s1, s2, z
x2 <- function(z,n1,n2,s1,s2){ z * (-sqrt((n1 * s2^2 + n2 * s1^2)/(n1 * n2))) }
z <- function(x1,x2,s1,s2,n1,n2){ (x1-x2) / sqrt(s1^2/n1 + s2^2/n2)  }

test <- NULL
for( i in seq(-4,4,by=0.1)){
	asd <- expand.grid( n1=c(50,100,150,200), n2=c(50,100,150,200), s1=c(0.5,1:5), s2=c(0.5,1:5))
	asd$x2 <- x2(z=i, n1=asd$n1, n2=asd$n1, s1=asd$s1, s2=asd$s2)
	asd$z <- round(z(x1=0, x2=asd$x2, n1=asd$n1, n2=asd$n1, s1=asd$s1, s2=asd$s2),1)
	test <- rbind(test, asd)
}

test$ratio <- (test$s1/test$n1)/(test$s2/test$n2)
test$ratio2 <- (test$s1^2/test$n1)/(test$s2^2/test$n2)

test$p1n <- NA
test$p1e <- NA

for( i in 1:nrow(test) ){ 
	test$p1n[i] <- calc.weights_MEM(xvec=c(0,test$x2[i]), svec=c(test$s1[i],test$s2[i]), nvec=c(test$n1[i], test$n2[i]), prior='pi_n')[1]
	test$p1e[i] <- calc.weights_MEM(xvec=c(0,test$x2[i]), svec=c(test$s1[i],test$s2[i]), nvec=c(test$n1[i], test$n2[i]), prior='pi_e')[1]
}	

par(mfrow=c(2,2))
plot(x=abs(test$z), y=test$p1e, ylim=c(0,1), xlab='|Z|', ylab='Pr(Not Exchangeable)', main='pi_e')

plot(x=abs(test$z), y=test$p1n, ylim=c(0,1), xlab='|Z|', ylab='Pr(Not Exchangeable)', main='pi_n')

plot(x=abs(test$z), y=test$p1e, ylim=c(0,1), xlab='|Z|', ylab='Pr(Not Exchangeable)', main='Both')
points(x=abs(test$z), y=test$p1n, col='gray65')
legend('bottomright', bty='n', pch=c(1,1), col=c('black','gray65'), legend=c('pi_e','pi_n'))

plot(x=abs(test$z), y=test$p1n, ylim=c(0,1), xlab='|Z|', ylab='Pr(Not Exchangeable)', main='Both', col='gray65')
points(x=abs(test$z), y=test$p1e)


# Regression fits
mod1 <- glm(p1n ~ n1 + n2 + s1 + s2 + abs(z), data=test)
summary(mod1)
hist(mod1$resid)

mod2 <- glm(p1n ~ ratio + abs(z), data=test)
summary(mod2)
hist(mod2$resid)

# Logistic regression to predict |Z|>=1
test$z_gte1 <- abs(test$z)>=1
test$perc_p1e <- test$p1e*100
test$perc_p1n <- test$p1n*100

test$col <- 'black'
test$col[which(test$z_gte1==T)] <- 'blue'
test$pch <- 3
test$pch[which(test$z_gte1==T)] <- 3

mod3 <- glm(z_gte1 ~ perc_p1e, data=test, family='binomial')
summary(mod3)
mean( test$z_gte1[which(abs(test$z) <= 2)] == (predict(mod3, type='response')>0.5) )


roc( test$z_gte1 ~ predict(mod3, type='response')) #0.7976455   
coords(roc( test$z_gte1 ~ predict(mod3, type='response')), best.method='youden', x='best')
table(zgte1=test$z_gte1, pred=predict(mod3, type='response')>0.7976455)
mean( test$z_gte1 == (predict(mod3, type='response')>0.7976455) )

plot(x=test$perc_p1e, y=predict(mod3, type='response'), col=test$col, pch=test$pch)
plot(x=test$perc_p1n, y=predict(mod3, type='response'), col=test$col, pch=test$pch)


#visualize and summarize results to look for other potential cutoffs
boxplot(perc_p1e ~ z_gte1, data=test[test$z < 2,], main = "|Z| < 2 Subset")
boxplot(perc_p1e ~ z_gte1, data=test, main = "Full Data Set")

max(test$perc_p1e[test$z == 0]) #71.48
quantile(test$perc_p1e[test$z == 0], 0.975) #66.129
max(test$perc_p1e[test$z_gte1 == FALSE]) #81.3376
quantile(test$perc_p1e[test$z_gte1 == FALSE], 0.975) #71.60678

# Strategy 1: simulation calibration, potentially using p*C_5 + (1-p)*C_80
# Strategy 2: general calibration above using range of values/combos:
#- Logistic regression with just pi_e
#- Logistic regression with pi_e + n1, n2, s1, s2, ratio, etc.
#- These are affected by the choice of threshold (e.g., >0.5 versus >Youden's J)
# Strategy 3: general calibration but focused moreso on range of values for trial
#- Could focus on power calculation assumptions
#- Could assume range of prevalence in each subgroup
# Strategy 4: consider different models for Z<0.5; 0.5 to 1.5; etc.




