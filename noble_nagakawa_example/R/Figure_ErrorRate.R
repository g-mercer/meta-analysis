rm(list = ls())

pacman::p_load(biwavelet, ggplot2, cowplot, gridGraphics, tidyverse)
#pdf(width = 7.365639 , height =  7.480176, file = "FigureS1.pdf")
#Define the layout of multiple panel figure using layout matrix#
#layout.mat = rbind(c(1:2),c(1:2),c(1:2),c(11:12),c(3:4),c(3:4),c(3:4),c(5:6),c(5:6),c(5:6),c(13:14),c(7:8),c(7:8),c(7:8),c(9:10),c(9:10),c(9:10))

method.col = c(rgb(0,0,0),rgb(0.9,0.6,0),rgb(0.8,0.4,0),rgb(0,0.45,0.7),rgb(0.35,0.7,0.9), rgb(0.35,0.5,0.5), rgb(0.35,0.6,0.8), rgb(0.40,0.5,0.8), rgb(0.40,0.2,0.3))

layout.mat = rbind(c(1:2),c(1:2),c(1:2),c(3:4),c(3:4),c(3:4),c(9:10),c(5:6),c(5:6),c(5:6),c(7:8),c(7:8),c(7:8))

#Set up plotting device and figure margins#

quartz(w=7.365639, h=7.480176)
layout(mat=layout.mat)
par(mar=c(0,0,0,0), oma=c(4,4,3,2))

#######################
##Equally correlated studies##
#######################
	load("EqualStudy.RData")
	#### Note. This is very dangerous because the objects are named the same and loaded, so if there are mis-haps, this will follow along the whole simulation. ####

	#Calculate error rate and the confidence interval of error rate#
	error_ES = array(dim=c(scenario, 9))
	error_ES.ci = array(dim=c(scenario, 9))

	# First, grab the relevant simulations for plotting
	relevant <- c(1:5, 6, 7, 8, 11)
	parm.low <- parm.low[,,relevant]
	 parm.up <-  parm.up[,,relevant]

	cover_ES = (parm.low<0) + (parm.up>0)

	for(i in 1:scenario){
		error_ES[i,] = 1-colSums(cover_ES[i,,]==2,na.rm=T)/iteration
		error_ES.ci[i,] = qnorm(0.975)*sqrt(error_ES[i,]*(1-error_ES[i,])/iteration)
	}

	#Find which sample corresponds to each paper size#
	size = array(dim=c(9,3))
	size[1,] = with(parm, which(mu==0.5 & rho.e==0.1))
	size[2,] = with(parm, which(mu==4.5 & rho.e==0.1))
	size[3,] = with(parm, which(mu==14.5 & rho.e==0.1))
	size[4,] = with(parm, which(mu==0.5 & rho.e==0.5))
	size[5,] = with(parm, which(mu==4.5 & rho.e==0.5))
	size[6,] = with(parm, which(mu==14.5 & rho.e==0.5))
	size[7,] = with(parm, which(mu==0.5 & rho.e==0.9))
	size[8,] = with(parm, which(mu==4.5 & rho.e==0.9))
	size[9,] = with(parm, which(mu==14.5 & rho.e==0.9))

	paper.mean = rep(c(1.5,5.5,15.5), times=3)
	paper.sd = rep(c(1.11,7.02,18.68), times=3)

	# This simply just can be updated to depict different numbers of methods in figures without having to update numbers throughout all the time. This equates to the number of columns (i.e., methods)
	method_num = 9
	for(i in c(1,3,7,9)){
		if(error_ES[size[i,1],1]>0.1){
			rate = error_ES[size[i,1],1]
			error_ES[size[i,1],1] = 999
		}
		plot(error_ES[size[i,1],]~c(1:method_num), xlim=c(0.5,29), ylim=c(0,0.13), pch=19, axes=F, col=method.col)
		arrows(x0=c(1:method_num),y0=error_ES[size[i,1],]-error_ES.ci[size[i,1],],y1=error_ES[size[i,1],]+error_ES.ci[size[i,1],], length=0, col=method.col)
		if(error_ES[size[i,1],1]>0.1){
			points(x=1, y=0.11, pch=19)
			text(x=1, y=0.1, labels=paste0(round(100*rate,0),"%"), cex=0.85)
			}
		box()
		abline(h=0.05, lty=2, col="grey")
		if(i %in% c(1,4,7)){axis(2, at=c(0,0.05,0.1), label=c("0%","5%","10%"), tck=-0.03)} #axis(2, at=c(0, 0.1, 0.20, 0.30, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), label=c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%"), tck=-0.03)
		for(j in 2:length(size[i,])){
			if(error_ES[size[i,j],1]>0.1){
			rate = error_ES[size[i,j],1]
			error_ES[size[i,j],1] = 0.5
			}
			points(error_ES[size[i,j],]~c(((j-1)*(method_num+1)+1):((j-1)*(method_num+1)+method_num)), pch=19, col=method.col)
			arrows(x0=c(((j-1)*(method_num+1)+1):((j-1)*(method_num+1)+method_num)), y0=error_ES[size[i,j],]-error_ES.ci[size[i,j],], y1=error_ES[size[i,j],]+error_ES.ci[size[i,j],], length=0, col=method.col)
			abline(v=(method_num+1)*(j-1), lty=2, col="grey")
			if(error_ES[size[i,j],1]>0.1){
				points(x=(j-1)*(method_num+1)+1, y=0.11, pch=19)
				text(x=(j-1)*(method_num+1)+1, y=0.1, labels=paste0(round(100*rate,0),"%"), cex=0.85)
				}
		}
		mtext(paste0("Mean:",paper.mean[i]),side=1,adj=0.95,line=-2.3, cex=0.75)
		mtext(paste0("SD:",paper.sd[i]),side=1, adj=0.95,line=-1.3, cex=0.75)
		mtext(expression(tau==0.1), side=3, adj=0.1, cex=0.8, line=-1.4)
		mtext(expression(tau==0.5), side=3, adj=0.5, cex=0.8, line=-1.4)
		mtext(expression(tau==1), side=3, adj=0.9, cex=0.8, line=-1.4)
		if(i %in% c(1,2,3)){mtext(expression(rho==0.1), side=1, adj=0.5, cex=0.8, line=-1.3)}
		if(i %in% c(4,5,6)){mtext(expression(rho==0.5), side=1, adj=0.5, cex=0.8, line=-1.3)}
		if(i %in% c(7,8,9)){mtext(expression(rho==0.9), side=1, adj=0.5, cex=0.8, line=-1.3)}
		#mtext(paste0("(",LETTERS[3+i],")"), side=1, adj=0.02, line=-1.35, cex=0.75)
	}


##################################
## Unequally correlated studies ##
##################################

	#Load results from simulations that allow tau to vary#
	#Because different simulation files using the same variable name, save results in new variable os that it can be plotted together with the unequal study with constant tau across paper scenario#
	load("UnequalStudy_VaryTau.RData")
	   error.varytau = array(dim=c(scenario,9))
	error.ci.varytau = array(dim=c(scenario,9))

	# First, grab the relevant simulations for plotting
	relevant <- c(1:5, 6, 7, 8, 11)
	parm.low <- parm.low[,,relevant]
	 parm.up <-  parm.up[,,relevant]

	cover.varytau = (parm.low<0) + (parm.up>0)

	for(i in 1:scenario){
		error.varytau[i,] = 1-colSums(cover.varytau[i,,]==2,na.rm=T)/iteration #Error rate#
		error.ci.varytau[i,] = qnorm(0.975)*sqrt(error.varytau[i,]*(1-error.varytau[i,])/iteration) #Lenght of CI of the error rate#  

		######### NOTE. There was an error here. THERE IS NO "error" OBJECT IN "UnequalStudy_VaryTau.RData", and the error object is called `error.varytau`. Since the whole workspace is loaded, this mean the error is coming from "EqualStudy.RData". They should have added in "error.varytau" object that they created. NOTE I AM CORRECTING WITH UNIQUE LABELLING OF COVER AND ERROR as a safeguard across the whole script. All simulations were run independently in three separate workspaces, so there should all be unique param.low and param.up objects, which are overwritten each load ################
	}

	load("UnequalStudy.RData")

	#Calculate error rate and the confidence interval of error rate#
	   error_US = array(dim=c(scenario, 9))
	error_US.ci = array(dim=c(scenario, 9))

	# First, grab the relevant simulations for plotting
	relevant <- c(1:5, 6, 7, 8, 11)
	parm.low <- parm.low[,,relevant]
	 parm.up <-  parm.up[,,relevant]

	cover_US = (parm.low<0) + (parm.up>0)
	
	for(i in 1:scenario){
		error_US[i,] = 1-colSums(cover_US[i,,]==2,na.rm=T)/iteration
		error_US.ci[i,] = qnorm(0.975)*sqrt(error_US[i,]*(1-error_US[i,])/iteration)
	}

	#Find which sample corresponds to each paper size#
	size = array(dim=c(6,3))
	size[1,] = with(parm, which(mu==0.5 & rho.e==0.1))
	size[2,] = with(parm, which(mu==4.5 & rho.e==0.1))
	size[3,] = with(parm, which(mu==14.5 & rho.e==0.1))
	size[4,] = with(parm, which(mu==0.5 & rho.e==0.6))
	size[5,] = with(parm, which(mu==4.5 & rho.e==0.6))
	size[6,] = with(parm, which(mu==14.5 & rho.e==0.6))

	paper.mean = rep(c(1.5,5.5,15.5), times=2)
	paper.sd = rep(c(1.11,7.02,18.68), times=2)

	# This simply just can be updated to depict different numbers of methods in figures without having to update numbers throughout all the time. This equates to the number of columns (i.e., methods)
	method_num = 9

	for(i in c(1,3,4,6)){
		if(error_US[size[i,1],1]>0.1){ #If error_US rate of method 1 too big to show, plot the value on the figure#
			rate = error_US[size[i,1],1]
			error_US[size[i,1],1] = 999
		} #ylim=c(0,0.13)
		plot(error_US[size[i,1],]~c(1:method_num), xlim=c(0.5,40), ylim=c(0,0.13), pch=19, axes=F, col=method.col)
		arrows(x0=c(1:method_num),y0=error_US[size[i,1],]-error_US.ci[size[i,1],],y1=error_US[size[i,1],]+error_US.ci[size[i,1],], length=0, col=method.col)
		box()
		abline(h=0.05, lty=2, col="grey")
		#Directly note value of error_US rate if too high to show on figure#
		if(error_US[size[i,1],1]>0.1){
			points(x=1, y=0.11, pch=19)
			text(x=1, y=0.1, labels=paste0(round(100*rate,0),"%"), cex=0.85)
			}
		if(i %in% c(1,4)){axis(2, at=c(0,0.05,0.1), label=c("0%","5%","10%"), tck=-0.03)}
		for(j in 2:length(size[i,])){
			if(error_US[size[i,j],1]>0.1){
			rate = error_US[size[i,j],1]
			error_US[size[i,j],1] = 0.5
			}
			points(error_US[size[i,j],]~c(((j-1)*(method_num+1)+1):((j-1)*(method_num+1)+method_num)), pch=19, col=method.col)
			arrows(x0=c(((j-1)*(method_num+1)+1):((j-1)*(method_num+1)+method_num)), y0=error_US[size[i,j],]-error_US.ci[size[i,j],], y1=error_US[size[i,j],]+error_US.ci[size[i,j],], length=0, col=method.col)
			abline(v=(method_num+1)*(j-1), lty=2, col="grey")
			if(error_US[size[i,j],1]>0.1){
				points(x=(j-1)*(method_num+1)+1, y=0.11, pch=19)
				text(x=(j-1)*(method_num+1)+1, y=0.1, labels=paste0(round(100*rate,0),"%"), cex=0.85)
				}
		}
		if(error.varytau[i,1]>0.1){
			rate = error.varytau[i,1]
			error.varytau[i,1] = 0.5
		}
		points(error.varytau[i,]~c(32:(31+method_num)), pch=19, col=method.col)
		arrows(x0=c(32:(31+method_num)), y0=error.varytau[i,]-error.ci.varytau[i,], y1=error.varytau[i,]+error.ci.varytau[i,], length=0, col=method.col)
		if(error.varytau[i,1]>0.1){
			points(x=32, y=0.11, pch=19)
			text(x=32, y=0.1, labels=paste0(round(100*rate,0),"%"), cex=0.85)
			}
		abline(v=31, lty=2, col="grey")
		mtext(paste0("Mean:",paper.mean[i]),side=1,adj=0.95,line=-2.3, cex=0.75)
		mtext(paste0("SD:",paper.sd[i]),side=1, adj=0.95,line=-1.3, cex=0.75)
		mtext(expression(tau==0.1), side=3, adj=0.06, cex=0.8, line=-1.4)
		mtext(expression(tau==0.5), side=3, adj=0.35, cex=0.8, line=-1.4)
		mtext(expression(tau==1), side=3, adj=0.64, cex=0.8, line=-1.4)
		mtext(expression(tau==0.1-1), side=3, adj=0.98, cex=0.77, line=-1.4)
		if(i %in% c(1,2,3)){mtext(expression(rho==0.1-0.4), side=1, adj=0.5, cex=0.8, line=-1.3)}
		if(i %in% c(4,5,6)){mtext(expression(rho==0.6-0.9), side=1, adj=0.5, cex=0.8, line=-1.3)}
		#mtext(paste0("(",LETTERS[12+i],")"), side=1, adj=0.02, line=-1.35, cex=0.75)
	}

#Label each section of figures with the type of non-independence#
#mtext("Experiment 1: independence", outer=T, line=0.5, font=2, adj=1)
mtext("Experiment 1: equal correlation", outer=T, line=0.5, font=2, adj=1)
mtext("Experiment 2: unequal correlation", outer=T, line=-26, font=2, adj=1)

#Label axes#
#mtext("Error rate %", side=2, outer=T, line=2.5, font=2)

#Legends#
#mtext("Methods:", outer=T, line=0.5, font=2, adj=0, cex=0.70)
#legend("topleft", legend=c("1","2","3","4","5", "6", "7", "8", "9"), pch=19, xpd=NA, inset = c(-0.8, -3.52), bty="n", horiz=T, cex=1.35, col=method.col[1:9])
#dev.off()
p2 <- recordPlot()

########################################
#Adding new figures
#######################################
# Lets just make a figure which shows the distribution of effects pooled across all scenarios. This gives less precedence to one particular scenario. What we really want to see is, how well do these perform, on average, across all the simulated scenarios. 

   combined_error <- as.data.frame(rbind(error_ES, error.varytau, error_US)[,-1]) # we don't care about 1, so we can ditch)

# Now, lets just plot across all the methods
   colnames(combined_error) <- c("AV", "One", "MLM", "RVE", "Papers_df", "CS", "SW_df", "Bayes")

combined_data <- tidyr::pivot_longer(combined_error, cols = c(1:8), names_to = "Method")
combined_data$Method <- factor(combined_data$Method, levels = c("AV", "One", "MLM", "RVE", "CS", "Papers_df", "SW_df", "Bayes"))

# Some summary stats
apply(combined_error, 2, range) #0.0344 for Bayes is lowest of all; 
apply(combined_error, 2, mean)
apply(combined_error, 2, median)

#pdf(width=7.365639, height = 5.207048, file = "Figure 1.pdf")
p1 <- ggplot(combined_data, aes(x=Method, y=value*100)) +
  geom_violin(aes(fill = Method), trim=FALSE) + 
  geom_jitter(shape=16, position=position_jitter(0.1), color = "black") +
  labs(y = "Error Rate %",
  		x = "Method") + 
  scale_fill_manual(breaks = c("AV", "One", "MLM", "RVE", "CS", "Papers_df", "SW_df", "Bayes"), 
                   values=c(rgb(0.9,0.6,0),rgb(0.8,0.4,0),rgb(0,0.45,0.7),rgb(0.35,0.7,0.9), rgb(0.35,0.6,0.8), rgb(0.35,0.5,0.5), rgb(0.40,0.5,0.8), rgb(0.40,0.2,0.3))) + 
  							#rgb(0.9,0.6,0),rgb(0.8,0.4,0),rgb(0,0.45,0.7),rgb(0.35,0.7,0.9), rgb(0.35,0.6,0.8), rgb(0.35,0.5,0.5), rgb(0.40,0.5,0.8), rgb(0.40,0.2,0.3)
  theme_classic() +         
  geom_hline(yintercept = 5, linetype = "dashed", colour = "darkgrey") +
  theme(axis.text.x = element_blank())
#dev.off()
#rgb(0.8,0.1,0.1), rgb(0.8,0.1,0.4)

############################################
# New combined figure
#plot_grid(p1, p2,
 #         nrow = 2, ncol = 1, labels = 'AUTO',
  #        hjust = 0, vjust = 1.5, label_size = 16, rel_heights = c(1,2))

plot_grid(p1, p2,
			 labels = 'AUTO',
             hjust = 0, vjust = 1.5, label_size = 16)


quartz.save(file = "Figure1_revised.pdf", type = "pdf")