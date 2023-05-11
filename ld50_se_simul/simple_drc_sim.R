## data:
Concentration <- as.numeric(c(0, 100, 200, 300, 400, 500, 600, 700, 800))
Survival <- as.integer(c(20, 20, 20, 18, 12, 8, 4, 0, 0))
Total <- as.integer(c(20, 20, 20, 20, 20, 20, 20, 20, 20))
Data <- data.frame(Concentration, Survival, Total)
Data$Dead <- Data$Total-Data$Survival
Data$Proport_Surv <- (Data$Total - Data$Dead) / Data$Total

mod <- glm(cbind(Survival, Dead) ~ Concentration, Data, family=binomial(link = "logit"))

summary(mod)

library(MASS)
dose.p(mod, p=0.5)

ec = 0.5
library(VGAM)
eta <- logit(ec) 
beta <- coef(mod)[1:2] 
ecx <- (eta - beta[1])/beta[2] 
pd <- -cbind(1, ecx)/beta[2] 
ff = as.matrix(vcov(mod)[1:2,1:2])
se <- sqrt(((pd %*% ff )*pd) %*% c(1, 1))
upper = (ecx+se*1.96)
lower = (ecx-se*1.96)
df1 = data.frame(ecx, lower, upper)

#############

# two examples, one for binomial and one for continuous data show that spread of data affects SE

## data:
Concentration <- as.numeric(c(0, 100, 200, 300, 400, 500, 600, 700, 800))
Survival <- as.integer(c(20, 20, 20, 10, 0, 0, 0, 0, 0))
Total <- as.integer(c(20, 20, 20, 20, 20, 20, 20, 20, 20))
Data <- data.frame(Concentration, Survival, Total)
Data$Dead <- Data$Total-Data$Survival
Data$Proport_Surv <- (Data$Total - Data$Dead) / Data$Total

mod <- glm(cbind(Survival, Dead) ~ Concentration, Data, family=binomial(link = "logit"))

summary(mod)

library(MASS)
dose.p(mod, p=0.5)

#############

## data:
Concentration <- as.numeric(c(0, 100, 200, 300, 400, 500, 600, 700, 800))
Survival <- as.integer(c(20, 20, 15, 10, 5, 0, 0, 0, 0))
Total <- as.integer(c(20, 20, 20, 20, 20, 20, 20, 20, 20))
Data <- data.frame(Concentration, Survival, Total)
Data$Dead <- Data$Total-Data$Survival
Data$Proport_Surv <- (Data$Total - Data$Dead) / Data$Total

mod <- glm(cbind(Survival, Dead) ~ Concentration, Data, family=binomial(link = "logit"))

summary(mod)

library(MASS)
dose.p(mod, p=0.5)

################

# in drc
library(drc)
Data.m1 <- drm(Survival/Total ~ Concentration,data=Data, weights = Total,fct=LL.3(), type = "binomial")

summary(Data.m1)

################
# does spread of the data in the y dimension affect SE of LD50 in the x dimension. Yes.

data("algae")

algae.m1 <- drm(vol~conc,data=algae,fct=LL.3())

summary(algae.m1)

conc <- sort(rep(c(0,1,10,100,1000),3))
growth <- c(90,100,110,90,100,110,40,50,60,4,5,6,1,2,3)
growth1 <- c(95,100,105,95,100,105,45,50,55,4.5,5,5.5,1.5,2,2.5)

test1 <- data.frame(conc,growth)

test2 <- data.frame(conc,growth1)

test1.m1 <- drm(growth~conc,data=test1,fct=LL.3())

summary(test1.m1)

test2.m1 <- drm(growth1~conc,data=test2,fct=LL.3())

summary(test2.m1)

####################################

# this example shows that number of individuals in a sample affects SE, use weights option with binomial

## data:
Concentration <- as.numeric(c(0, 100, 200, 300, 400, 500, 600, 700, 800))
Survival <- as.integer(c(20, 20, 15, 10, 5, 0, 0, 0, 0))
Total <- as.integer(c(20, 20, 20, 20, 20, 20, 20, 20, 20))
Data <- data.frame(Concentration, Survival, Total)
Data$Dead <- Data$Total-Data$Survival
Data$Proport_Surv <- (Data$Total - Data$Dead) / Data$Total

mod <- glm(cbind(Survival, Dead) ~ Concentration, Data, family=binomial(link = "logit"))

summary(mod)

library(MASS)
dose.p(mod, p=0.5)

#####################

# in drc
library(drc)
Data.m1 <- drm(Survival/Total ~ Concentration,data=Data, weights = Total,fct=LL.2(), type = "binomial")

summary(Data.m1)

# multiply sample size by 10

# smaller SE

## data:
Concentration <- as.numeric(c(0, 100, 200, 300, 400, 500, 600, 700, 800))
Survival <- as.integer(c(20, 20, 15, 10, 5, 0, 0, 0, 0)) * 10
Total <- as.integer(c(20, 20, 20, 20, 20, 20, 20, 20, 20)) * 10
Data <- data.frame(Concentration, Survival, Total)
Data$Dead <- Data$Total-Data$Survival
Data$Proport_Surv <- (Data$Total - Data$Dead) / Data$Total

mod <- glm(cbind(Survival, Dead) ~ Concentration, Data, family=binomial(link = "logit"))

summary(mod)

library(MASS)
dose.p(mod, p=0.5)

#####################

# in drc
library(drc)
Data.m1 <- drm(Survival/Total ~ Concentration,data=Data, weights = Total,fct=LL.2(), type = "binomial")

summary(Data.m1)