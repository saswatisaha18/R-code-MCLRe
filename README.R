#########################################################################################
# README file accompanying paper "Testing effect of a drug using multiple nested models #
# for the dose-response".                                                               #
#                                                                                       #
# Functions refered to in this file can be found in file "DRevaluation.R                #
#                                                                                       #
# Below we will describe how to apply the method proposed in the paper to an example    #
# data set and how to simulate relevant distributions of the test statistics            #
#                                                                                       #
# Please make sure to have the following packages installed:                            #
# DoseFinding, MCPMod, stats4, plyr                                                     #
#                                                                                       #
# Made with R version 3.0.2                                                             #
#                                                                                       #
#########################################################################################

#To load the relevant functions and libraries, load the functions in Rcode_paper.R, for
#examply through:
source("DRevaluation.R")

###############################################
# Applying the method to an example data set  #
###############################################
#The example data discussed in the paper can be found in the DoseFinding package, or MCPMod
#package. Both packages were loaded with the previous command. Load the data as follows:

data(biom)

#For illustration, we can plot the average dose means:
plot(aggregate(biom$resp,list(biom$dose),mean),xlab="Dose",ylab="Response")

#apply new method from the paper to the data:
result <- Modsel(biom,estMED=TRUE)

#the selected model is stored under selmod and is in this case the power model:
result$selmod
#and the MED estimate is stored under med and here equals 0.23:
result$med

#we can add the estimated curve to the plot of the data (see above) using the estimated coefficients:
doses <- seq(0,1,0.01)  #sequence of doses to predict an outcome for
coefs <- result$coefs    #estimated coefficients of the power model
lines(doses,coefs[1]+coefs[2]*doses^coefs[3],col="red") #add the predicted curve to the plot

#note that, by default, the unrestricted model is included. To run the method without
#it, execute:
result_wu <- Modsel(biom,estMED=TRUE,critvals=c(3.979,2,2),addFull=FALSE)
#the results are the same

#to apply the MCP-Mod approach as done in Bretz et al. 2005, run:
models <- Mods(linear=NULL,linlog=NULL,quadratic=c(-0.854,-1),exponential=c(0.279,0.15),emax=0.2,doses=c(0,0.05,0.2,0.6,1),addArgs=list(off=1))
mm <- MCPMod(dose,resp,biom,models,Delta=0.4,alpha=0.05,selModel="maxT")

#For more details, execute:
summary(mm)

#we can also add the estimated emax model to the plot:
lines(doses,0.322+0.746*doses/(doses+0.142),col="blue")

#to do the same evaluation with MCP-Mod, but with model selection based on AIC, run:
mm_aic <- MCPMod(dose,resp,biom,models,Delta=0.4,alpha=0.05,selModel="AIC")

########################################################
# Simulating the distributions of the test statistics  #
########################################################
# Simulate theoretical distributions in GENERAL CASE (See Web Appendix A)

cl <- rchisq(1000,1)  #distribution of lratio comparing constant model and linear model
lp <- rchisq(1000,1)  #distribution of lratio comparing linear and power model
p4pl <- rchisq(1000,1) #distribution of lratio comparing power and 4pl model
fourplu <- rchisq(1000,1) #distribution of lratio comparing 4pl model and unrestricted model

#Distribution of the test statistics with the unrestricted model included:
T0u <- apply(cbind(cl,(cl+lp)/2,(cl+lp+p4pl)/3,(cl+lp+p4pl+fourplu)/4),1,max)
T1u <- apply(cbind(lp,(lp+p4pl)/2,(lp+p4pl+fourplu)/3),1,max)
T2u <- apply(cbind(p4pl,(p4pl+fourplu)/2),1,max)
T3u <- fourplu

#Distributions without the unrestricted model:
T0 <- apply(cbind(cl,(cl+lp)/2,(cl+lp+p4pl)/3),1,max)
T1 <- apply(cbind(lp,(lp+p4pl)/2),1,max)
T2 <- p4pl

#one can look at the quantiles of these distributions using the commands like:
quantile(T0u,0.95)


# Simulate theoretical DISTRIBUTIONS FOR THE CANDIDATE SET (See Web Appendix A)
#distributions for T_1, T_2 and T_3:
w21 <- rnorm(1000,0,1)
w32 <- rnorm(1000,0,1)
w43 <- rnorm(1000,0,1)

#Distributions with the unrestricted model included:
T1c <- apply(cbind(w21^2,(w21^2+w32^2*(w32>0))/2,(w21^2+w32^2+w43^2)/3),1,max)
T2c <- apply(cbind(w32^2*(w32>0),(w32^2+w43^2)/2),1,max)
T3c <- w43^2

#Distributions without the unrestricted model
T1c_wu <- apply(cbind(w21^2,(w21^2+w32^2*(w32>0))/2),1,max)
T2c_wu <- w32^2*(w32>0)

#distribution for T_0 for the candidate set:
#First we obtain the individual distributions of -2*log-likelihood ratio for comparing
#the constant model against the alternative models in the candidate set. The following
#command will return a data frame with these distributions:
res <- simdist(nsim=1000)

#to obtain the distribution of T_0, we take the maximum over the weighted distributions
#from above:

#with the unrestricted model included:
T0c <- apply(cbind(res$linear,res$power/2,res$fpl/3,res$full/4),1,max)
#without the unrestricted model:
T0c_wu <- apply(cbind(res$linear,res$power/2,res$fpl/3),1,max)

#now we can look at a quantile of interest of this distribution, say a 95%-quantile:
quantile(T0c,0.95)
quantile(T0c_wu,0.95)
