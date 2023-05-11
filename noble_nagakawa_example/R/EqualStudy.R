library(plyr)
library(MASS)
library(metafor)
library(robumeta)
library(Matrix)
library(nlme)
library(lme4)
library(MCMCglmm)
library(robumeta)
library(clubSandwich)

rm(list=ls())

#Set number of papers in the simulated meta-analysis data set#
no.paper = 20
#Upper and lower bound for the number of replicates in control and treatment for each study#
#Number of replicates is an integer draw uniformly between the upper and lower bound#
rep.low = 3
rep.up = 20
#Number of iterations of the simulation#
iteration = 5000

#standard deviation of the random study effect#
sd.e = c(0.1, 0.5, 1)
#correlation coefficient between studies within the same paper#
rho.e = c(0.1, 0.5, 0.9)
#Size parameter used to generate number of studies per paper from shifted binomial distribution#
#The three sizes correspond to means of 1.5, 5.5 and 15.5#
#Obtained from StudyNo_Distribution.R#
size = c(0.3339891,0.4522312,0.6287481)

#List experimental factors#
parm = expand.grid(size=size, sd.e=sd.e, rho.e=rho.e)
parm$mu = rep(c(0.5,4.5,14.5), times=9)
scenario = dim(parm)[1]

#Create array to store estimates from each method#
  parm.est = array(dim=c(scenario,iteration,11))
  parm.low = array(dim=c(scenario,iteration,11))
   parm.up = array(dim=c(scenario,iteration,11))
tau2.paper = array(dim=c(scenario,iteration,2))

set.seed(5)

for(k in 1:scenario){
	for(a in 1:iteration){
	#Generate number of studies for each of the 20 papers contained in the meta-analysis#
	no.study = 1 + rnbinom(n=no.paper, size=parm$size[k], mu=parm$mu[k])
	raw.data = data.frame()
	for(i in 1:no.paper){ #Loop through all papers#
		#The mean response variable in control and treatment is 10#
		#The value of mu does not matter as it cancels in the calculation of response ratio#
		mu = 10
		#alpha.t=1 means treatment has no effect#
		alpha.t = 1
		#correlation matrix for observed effect size in paper i#
		cor.e = matrix(data=parm$rho.e[k], nrow=no.study[i], ncol=no.study[i])
		diag(cor.e) = 1
		#Calculate covariance matrix from correlation matrix#
		cov.e = diag(parm$sd.e[k], nrow=no.study[i]) %*% cor.e %*% diag(parm$sd.e[k], nrow=no.study[i])
		log.e = mvrnorm(n=1, mu=rep(0,no.study[i]), Sigma=cov.e)
		for(j in 1:no.study[i]){ #Loop through each study#
			#alpha is the study-specific treatment effect#
			alpha = alpha.t * exp(log.e[j])
			no.rep = runif(n=1, rep.low, rep.up)
			sd.log.epsilon = runif(1,0.1,0.3)
			yc = mu * exp(rnorm(n=no.rep, mean=0, sd=sd.log.epsilon))
			yt = mu * alpha * exp(rnorm(n=no.rep, mean=0, sd=sd.log.epsilon))
			data.temp = data.frame(yc, yt, paper=rep(i,no.rep), study=rep(j, no.rep))
			raw.data = rbind(raw.data, data.temp)
		}
	}
	
	#Calculate response ratio and its variance for each study#
	effect.size = ddply(raw.data, c("paper","study"), summarise, log.ratio=log(mean(yt)/mean(yc)), 
													var.log.ratio=var(yt)/length(yt)/mean(yt)^2+var(yc)/length(yc)/mean(yc)^2)
	
	#Random effects model that assumes independence#
	mod1 = try(rma(yi=log.ratio, vi=var.log.ratio, mods=~1, method="REML", data=effect.size, test="knha"))
	parm.est[k,a,1] = as.numeric(try(mod1$b))
	parm.low[k,a,1] = as.numeric(try(mod1$ci.lb))
	parm.up[k,a,1] = as.numeric(try(mod1$ci.ub))

	#Calculate mean effect size from each paper#
	effect.size.paper = ddply(effect.size, c("paper"), summarise, log.ratio.p=rma(yi=log.ratio,vi=var.log.ratio, method="FE")$b, var.log.ratio.p=(rma(yi=log.ratio,vi=var.log.ratio, method="FE")$se)^2)
	mod2 = try(rma(yi=log.ratio.p, vi=var.log.ratio.p, mods=~1, method="REML", data=effect.size.paper, test="knha"))
	parm.est[k,a,2] = as.numeric(try(mod2$b))
	parm.low[k,a,2] = as.numeric(try(mod2$ci.lb))
	parm.up[k,a,2] = as.numeric(try(mod2$ci.ub))

	#Randomly choose one effect size from each paper#
	effect.size.random = ddply(effect.size, c("paper"), summarise, study.no=sample(study,size=1), log.ratio=log.ratio[study.no], var.log.ratio=var.log.ratio[study.no])
	mod3 = try(rma(yi=log.ratio, vi=var.log.ratio, mods=~1, method="REML", data=effect.size.random, test="knha"))
	parm.est[k,a,3] = as.numeric(try(mod3$b))
	parm.low[k,a,3] = as.numeric(try(mod3$ci.lb))
	parm.up[k,a,3] = as.numeric(try(mod3$ci.ub))
	
	#Model with an additive random paper effect#
	effect.size$study.unique = paste0(effect.size$paper,"-",effect.size$study)
	mod4 = try(rma.mv(yi=log.ratio, V=var.log.ratio, mods=~1, random=list(~1|paper,~1|study.unique), data=effect.size, test="t"))
	parm.est[k,a,4] = as.numeric(try(mod4$b))
	parm.low[k,a,4] = as.numeric(try(mod4$ci.lb))
	 parm.up[k,a,4] = as.numeric(try(mod4$ci.ub))
	#Store the variance component estimates, first element is paper, second element is study#
	tau2.paper[k,a,1] = as.numeric(try(mod4$sigma2[1]))	
	tau2.paper[k,a,2] = as.numeric(try(mod4$sigma2[2]))

	#Robust variance estimates method by Hedges#
	mod5 = try(robu(formula=log.ratio~1, data=effect.size, studynum=paper, var.eff.size=var.log.ratio))
	parm.est[k,a,5] = as.numeric(try(mod5$reg_table$b.r)) 
	parm.low[k,a,5] = as.numeric(try(mod5$reg_table$CI.L))
	parm.up[k,a,5] = as.numeric(try(mod5$reg_table$CI.U))

##################################################################
	## Added simulations
##################################################################

	# METHOD 6 # Use papers as DF when calculating CIs with metafor, method 8 #
		            papers <- mod4$s.nlevels [1]
		   parm.est[k,a,6] <- as.numeric(try(mod4$b))
		    parm.up[k,a,6] <- as.numeric(try(mod4$b + (mod4$se * qt(0.975, papers-1))))
		   parm.low[k,a,6] <- as.numeric(try(mod4$b - (mod4$se * qt(0.975, papers-1))))

	# METHOD 7 - clubsandwich vs robumeta should be the same 
		   mod6 <- coef_test(mod4, vcov="CR2", cluster = effect.size$paper)
		   SW_df <- mod6$df
			parm.est[k,a,7] <- as.numeric(try(mod6$beta))
		    parm.up[k,a,7] <- as.numeric(try(mod6$beta + (mod6$SE * qt(0.975, SW_df))))
		   parm.low[k,a,7] <- as.numeric(try(mod6$beta - (mod6$SE * qt(0.975, SW_df))))	

	# METHOD 8 - Saitterwaite DF
		   parm.est[k,a,8] <- as.numeric(try(mod4$b))
		    parm.up[k,a,8] <- as.numeric(try(mod4$b + (mod4$se * qt(0.975, SW_df))))
		   parm.low[k,a,8] <- as.numeric(try(mod4$b - (mod4$se * qt(0.975, SW_df))))

    # METHOD 9 - Robust with paper-1
		mod7 <- robust(mod4, cluster = effect.size$paper)
		   parm.est[k,a,9] <- as.numeric(try(mod7$b))
		    parm.up[k,a,9] <- as.numeric(try(mod7$b + (mod7$se * qt(0.975, papers-1))))
		   parm.low[k,a,9] <- as.numeric(try(mod7$b - (mod7$se * qt(0.975, papers-1))))

	# METHOD 10 - try to ditch DF correction.Method 8 and 7 should be identical
	mod8 = try(robu(formula=log.ratio~1, data=effect.size, studynum=paper, var.eff.size=var.log.ratio, small = FALSE))
	parm.est[k,a,10] = as.numeric(try(mod8$reg_table$b.r)) 
	parm.low[k,a,10] = as.numeric(try(mod8$reg_table$CI.L))
	 parm.up[k,a,10] = as.numeric(try(mod8$reg_table$CI.U))   

	# METHOD 11 = Use a Bayesian approach with MCMCglmm
		prior <- list(R = list(V = 1, nu = 0.002), 
              G = list(G1 = list(V = 1 , nu = 1, alpha.mu=0, alpha.V=25^2)))

	 mod9 <- MCMCglmm(log.ratio ~ 1, mev = effect.size$var.log.ratio, random = ~ paper, data = effect.size, prior = prior, verbose = FALSE)
		 parm.est[k,a,11] <- summary(mod9)$solutions[1]
		 parm.low[k,a,11] <- summary(mod9)$solutions[2]
		  parm.up[k,a,11] <- summary(mod9)$solutions[3]
	
	}
	#print when the kth scenario is done to track simulation progress#
	print(k)
}


save.image(file="./sim_check_rerun/EqualStudy.RData")


# Explore what the results are for the precision of the estimate across 27 sceneries; as defined in param. Take the SD of the estimate vector for all 10 methods
    # SD_result <- apply(parm.est, c(1,3), function(x) sd(x, na.rm = TRUE)) # take SD across the cols for each row and method
#full_SD_result <- cbind(parm, SD_result) # To make this easier to see, have a look at the specific param combinations for each simulation