## This is the analysis of the biom data as in Bretz et al. (2005)
rm(list=ls())

set.seed(1)

library(chebsample)
library(mnll)
library(DoseFinding)

########################################################################
## fit models and calculate test statistics
pMCLRe<-function(data,doses,BLOCK_SIZE=100,MAX_BLOCKS=100,nSim=100)
{
  data <- data[order(data$dose),]
  
  ## calculate critical value
  getB <- function(n)
    t(svd(rep(1, n), n , n)$u[,-1])
  B <- getB(nrow(data))
  
  ## simulate standardized predictions
  ## emax model ("normal" interval)
  bndsEmax <- c(0.001,1.5)
  o <- list(x = data$dose, lambda = bndsEmax[1], upsilon = bndsEmax[2])
  q <- cheb_q(o, g_emax, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6)
  emax <- function(x, gamma) x / (x + gamma)
  tx_gamma <- function(x, gamma) stat(B, emax(x, gamma))
  sims_Emax <- sapply(esample(o, q, nSim), function(gamma) tx_gamma(x=o$x, gamma))
  
  ## emax model (larger interval)
  bndsEmax <- c(0.001,10)
  o <- list(x = data$dose, lambda = bndsEmax[1], upsilon = bndsEmax[2])
  q <- cheb_q(o, g_emax, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6)
  emax <- function(x, gamma) x / (x + gamma)
  tx_gamma <- function(x, gamma) stat(B, emax(x, gamma))
  sims_Emax2 <- sapply(esample(o, q, nSim), function(gamma) tx_gamma(x=o$x, gamma))
  
  ## exponential model
  bndsExpo <- c(0.1, 2)
  o <- list(x = data$dose, lambda = bndsExpo[1], upsilon = bndsExpo[2])
  q <- cheb_q(o, g_exponential, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6)
  expo <- function(x, gamma) exp(x / gamma)-1
  tx_gamma <- function(x, gamma) stat(B, expo(x, gamma))
  sims_Expo <- sapply(esample(o, q, nSim), function(gamma) tx_gamma(x=o$x, gamma))
  
  ## linear model (nothing to simulate, just one shape)
  sims_Lin <- stat(B, o$x)
  sims_LinLog <- stat(B, log(o$x+1))
  
  
  ## effect of including linear model
  sims_EmaxLin <- cbind(sims_Emax, sims_Lin)
  
  ## effect of including exponential model
  sims_EmaxLinExpo <- cbind(sims_Emax, sims_Lin, sims_Expo)
  
  
  ## effect of including linlog model
  sims_EmaxLinExpo <- cbind(sims_Emax, sims_Lin, sims_Expo, sims_LinLog)
  
mods <- vector("list", 4)
mods[[1]] <- fitMod(dose,resp, data=data, model="emax", bnds=bndsEmax)
mods[[2]] <- fitMod(dose,resp, data=data, model="linear")
mods[[3]] <- fitMod(dose,resp, data=data, model="exponential",bnds = bndsExpo)
mods[[4]] <- fitMod(dose,resp, data=data, model="linlog")

names(mods) <- lapply(mods, function(x) attr(x, "model"))
## function to calculate test statistic
getTstat <- function(mod, y)
{
  cf <- coef(mod)
  pred <- (predict(mod, predType = "l")-cf[1])/cf[2]
  B <- getB(length(pred))
  pr <- stat(B,pred)
  dat <- stat(B,y)
  sum(pr*dat)
}
tt <- lapply(mods, getTstat, y=data$resp)

## marginal p-value for different models
pvals <- numeric(4)
pvals[1] <- mnll_pvalue(tt$emax, sims_Emax, block_size = BLOCK_SIZE,
                        max_blocks = MAX_BLOCKS)
pvals[2] <- mnll_pvalue(tt$linear, sims_Lin,block_size = BLOCK_SIZE,
                        max_blocks = MAX_BLOCKS)
pvals[3] <- mnll_pvalue(tt$exponential, sims_Expo, block_size = BLOCK_SIZE,
                        max_blocks = MAX_BLOCKS)
pvals[4] <- mnll_pvalue(tt$exponential, sims_LinLog, block_size = BLOCK_SIZE,
                        max_blocks = MAX_BLOCKS)

## adjusted p-values for different models
adjpvals <- numeric(4)
sims <- cbind(sims_Emax, sims_Expo, sims_Lin, sims_LinLog)
adjpvals[1] <- mnll_pvalue(tt$emax, sims, block_size = BLOCK_SIZE,
                           max_blocks = MAX_BLOCKS)
adjpvals[2] <- mnll_pvalue(tt$linear, sims, block_size = BLOCK_SIZE,
                           max_blocks = MAX_BLOCKS)
adjpvals[3] <- mnll_pvalue(tt$exponential, sims, block_size = BLOCK_SIZE,
                           max_blocks = MAX_BLOCKS)
adjpvals[4] <- mnll_pvalue(tt$linlog, sims, block_size = BLOCK_SIZE,
                           max_blocks = MAX_BLOCKS)

## table
modNams <- c("Emax (gamma in [0.001,1.5])", "Linear", "Exponential (gamma in [0.1,2])", "Linlog")
parEst <- sapply(mods, function(x){
  cf <- coef(x)
  str <- sprintf("alpha=%.2f, beta=%.2f", cf[1], cf[2])
  if(length(cf) == 3)
    str <- paste(str, sprintf(", gamma=%.2f", cf[3]), sep="")
  str
})
out <- data.frame(model=modNams,parEst=parEst,teststatistic=round(as.numeric(tt),3), round(adjpvals,3), round(pvals,3))
out <- out[rev(order(out$teststatistic)),]
rownames(out) <- NULL
colnames(out) <- c("Model", "Parameter estimates", "test-statistic", "adj p-value", "p-value")
return(out)
}
