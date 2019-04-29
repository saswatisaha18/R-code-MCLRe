#Made with R version 3.0.2

#Required libraries
library(DoseFinding)
library(stats4)
library(plyr)

###################
#Fit power model  #
###################
#Function to fit the power model to the data, based on code from DoseFinding package
#(see also Pinheiro et al., Statist. Med. 2014, 33 1646-1661).

#Power function required for function fitPower
#INPUT
#dose: number specifying a dose
#e0, e1, pow, numbers giving value of model paramters
powermod <- function(dose, e0, e1, pow){
  e0 +  e1 * dose^pow
}

#INPUT:
#data: dataframe with a variable named "dose", containing the dose levels and a variable names "resp" containing the responses
#bnds: bounds for the non-linear parameter. If none are specified, bounds (0,10) are used

fitPower <- function(data,bnds=NULL){
  
  #define bounds
  if(is.null(bnds)){bnds<-c(0,10)}
  
  #data with doses and mean responses
  dataFit <- data.frame(dose=sort(unique(data$dose)),resp=as.numeric(tapply(data$resp,data$dose,mean)))
  
  #group sizes
  weights <- as.vector(table(data$dose))
  
  #bounds
  bnds <- matrix(bnds,nrow=1)
  
  #compute resXY
  dose <- dataFit$dose
  resp <- dataFit$resp
  
  m <- model.matrix(resp~1,dataFit)
  clinS <- diag(sqrt(weights))
  qrX <- qr(clinS %*% m)
  resXY <- as.numeric(qr.resid(qrX,sqrt(weights)*resp))
  
  #make a grid with starting values for pow
  grdnods <- (2*(1:30)-1)/(2*30)
  nodes <- matrix(grdnods*(bnds[2]-bnds[1])+bnds[1],ncol=1)
  
  #compute Zmat
  getPred <- function(vec, x) powermod(x,0,1,vec)
  Zmat <- apply(nodes, 1, getPred, x = dose)
  Zmat <- clinS%*%Zmat
  resZmat <- qr.resid(qrX,Zmat)
  
  #compute RSS
  colsms1 <- colSums(resZmat*resXY)
  colsms2 <- colSums(resZmat*resZmat)
  RSSvec <- sum(resXY*resXY)-(colsms1*colsms1)/colsms2
  
  #First optimization step:
  indMin <- which.min(RSSvec)
  strt <- nodes[indMin,]
  
  #starting values
  N<-30
  dif <- (bnds[2]-bnds[1])/N
  bnds[1] <- max(c(strt-1.1*dif),bnds[1])
  bnds[2] <- min(c(strt+1.1*dif),bnds[2])
  
  #local optimization
  optFunc <- function(par,x,qrX,resXY,clinS){
    Z <- powermod(x,0,1,par)
    Z <- clinS%*%Z
    resXZ <- try(qr.resid(qrX,Z))
    if (inherits(resXZ,"try-error")) {return(NA)}
    sumrsXYrsXZ <- sum(resXY*resXZ)
    sum(resXY*resXY)-sumrsXYrsXZ*sumrsXYrsXZ/sum(resXZ*resXZ)
  }
  
  optobj <- optimize(optFunc, c(bnds[1],bnds[2]), x=dose,qrX=qrX,resXY=resXY,clinS=clinS)
  coefs <- optobj$minimum
  
  #calculation of linear coefficients
  f0 <- powermod(dose,0,1,coefs)
  X <- cbind(1,f0)
  par0 <- qr.coef(qr(clinS%*%X),clinS%*%resp)
  par <- c(par0,coefs)
  names(par)<-c("int","fact","pow")
  df <- sum(weights)-length(par)
  
  #estimate linear model.
  m1 <- lm(data$resp~data$dose)
  coefs <- c(coef(m1),1)
  ll <- logLik(m1)
  w2<-1
  
  #compute log-likelihood
  llm2 <- function(int,fact,pow,s){
    -sum(log(1/sqrt(2*pi*s^2)*exp(-0.5*1/s^2*(data$resp-(int+fact*data$dose^pow))^2)))
  }
  
  RSS <- sum((data$resp-(par[1]+par[2]*data$dose^par[3]))^2)
  ll2 <- -llm2(par[1],par[2],par[3],sqrt(RSS/df))
  
  #if likelihood of power model is lower than likelihood of linear model, the model
  #has not converged and the power model estimate is put to the linear model estimate
  #since the linear model is a submodel of the power model
  if (ll2>ll){
    w2<-0
    ll <- ll2 
    coefs <- par
  }
  
  list("coefs"=coefs,"LogLik"=ll,"w2"=w2)
}

#OUTPUT
#Coefs: estimated model coefficients
#LogLik: log-likelihood value
#w2: indicates whether power model converged badly and a linear model was fitted instead (1=yes, 0=no)

###################
#fit 4PL model    #
###################
#function to fit the four-parameter logistic model.
#INPUT
#data: dataframe with a variable named "dose", containing the dose levels and a variable names "resp" containing the responses
#bnds: 2x2 matrix with bounds for the non-linear parameters theta_2 and theta_3 in the model. The first row corresponds to theta_2 and the second row to theta_3,
#in the parameterization theta_0+theta_1*dose^theta_3/(dose^theta_3+theta_2^theta_3)

fitSigEmax <- function(data, bnds=NULL){ 
  #estimate power model
  m2 <- fitPower(data)
  w3 <- 0
  
  #estimate model using DoseFinding package
  if(is.null(bnds)){
    bnds <- defBnds(1)$sigEmax
    bnds[2,2] <- 15
    bnds[2,1] <- 0
  }
  fit <- fitMod(dose,resp,data,"sigEmax",bnds=bnds)
  
  par <- coef(fit)
  names(par) <- c("int","fact","infl","pow")
  LL <- logLik(fit)
  
  #compare to estimate of power model. If log-likelihood is smaller than log-likelihood
  #of power model, power model is taken to be the best model estimate
  if(m2$LogLik>LL){
    LL <- m2$LogLik
    par <- m2$coefs
    w3 <- m2$w3
  }
  
  par <- par[c("int","fact","infl","pow")]
  
  list("coefs"=par,"LogLik"=LL,"w3"=w3)   
}

#OUTPUT
#coefs: estimated model coefficients
#LogLik: log-likelihood value
#w3: indicates whether four-parameter logistic model converged badly and the power model was fitted instead (1=yes, 0=no)

#######################################################################################
# Function which performs new model selection approach for the proposed candidate set #
#######################################################################################
#INPUT
#data: dataframe with a variable named "dose", containing the dose levels and a variable names "resp" containing the responses
#critvals: vector with critical values for the test statistics
#estMED: TRUE/FALSE, should an estimate of the MED be calculated?
#addFULL: should unrestricted (full) model be included?
Modsel <- function(data,critvals=c(4.020,2,2,2),estMED=FALSE, addFull=TRUE){
  resp <- data$resp
  dose <- data$dose
  
  #fit models
  m0 <- lm(resp~1)
  c0 <- coef(m0)
  m1 <- lm(resp~dose)
  c1 <- coef(m1)
  m2 <- fitPower(data)
  c2 <- m2$coefs
  m3 <- fitSigEmax(data)
  c3 <- m3$coefs
  m4 <- lm(resp~factor(dose))
  c4 <- coef(m4)
  
  ndose <- length(unique(dose))
  #the log likelihood values:
  ll <- c(logLik(m0),logLik(m1),m2$LogLik,m3$LogLik,logLik(m4))
  dgf <- c(2:4,ndose)
  
  if(addFull){
    #the test statistics when the unrestricted model is included:
    Tstat <- rep(0,4)
    for(i in 1:4){
      Tstat[i] <- max(2*(ll[(i+1):5]-ll[i])/(dgf[i:4]-i))
    }
    #selected model:
    selmod <- ifelse(Tstat[1]<critvals[1],"constant",ifelse(Tstat[2]<critvals[2],"linear",ifelse(Tstat[3]<critvals[3],"power",ifelse(Tstat[4]<critvals[4],"4PL","full"))))
  } else {
    #the test statistics when the unrestricted model is not included:
    Tstat <- rep(0,3)
    for(i in 1:3){
      Tstat[i] <- max(2*(ll[(i+1):4]-ll[i])/(dgf[i:3]-i))
    }
    selmod <- ifelse(Tstat[1]<critvals[1],"constant",ifelse(Tstat[2]<critvals[2],"linear",ifelse(Tstat[3]<critvals[3],"power","4PL")))
  }
  
  coefs <- if(selmod=="linear"){c1} else {if(selmod=="power"){c2} else {if(selmod=="4PL"){c3} else NULL}}
  coefs <- as.vector(coefs)
  
  med <- NA
  
  #estimate MED
  if(estMED){
    med <- findMED(data,model=selmod,coefs=coefs)
  }
  
  list("T"=Tstat, "selmod"=selmod, "med"=as.numeric(med), "coefs"=round(coefs,2))
}

#OUTPUT:
#Tstat: vector of test statistics
#selmod: the selected model
#med: MED estimate
#coefs: estimated model coefficients for selected model

####################
#Estimate MED      #
####################

#gradient of powermod, required for function findMED
powergrad <- function(dose, e0, e1, pow){
  cbind(1,dose^pow,e1*dose^pow*log(dose))
}


#estimates MED as defined in paper
#INPUT:
#data: dataframe with a variable named "dose", containing the dose levels and a variable names "resp" containing the responses
#model: string specifying for which model the MED should be estimated. Possible models are "constant", "linear", "power", "4PL", "full", "linlog", "quadratic".
#coefs: vector with model coefficients c(theta_0, theta_1, theta_2, theta_3) as given in Section 2.4 of paper, except for four-parameter logistic model,
#for which the following parameterization is used: theta_0+theta_1*dose^theta_3/(dose^theta_3+theta_2^theta_3)
#gamma: number, method will estimate 1-2*gamma confidence limit for MED estimation
#delta: number giving clinical relevance threshold for MED estimate

findMED <- function(data,model,coefs=NULL,gamma=0.1,delta=0.4){
  dose <- data$dose
  resp <- data$resp
  n <- nrow(data)
  mdose <- max(dose)
  
  #critical values of t-distribution
  crt<-qt(1-gamma,c(n-2,n-3,n-4))
  
  #constant model: no estimate can be given
  if(model=="constant"){
    med <- NA
  }
  
  #linear model:
  if(model=="linear"){
    med <- as.numeric(delta/coefs[2]) #estimate of dose with a clinical relevant effect of 0.4
    #check if lower confidence limit is above placebo effect
    if(med>mdose){med <- NA} else {  #med = NA if outside tested dose range
      #estimate variance-covariance matrix
      predLin <- coefs[1]+coefs[2]*dose
      sLin <- sum((resp-predLin)^2)/(n-2)
      fl <- function(d){
        c(1,d)
      }
      Fl <- t(sapply(dose,fl))
      cvmat <- try(solve(crossprod(Fl))*sLin,silent=TRUE)
      #if variance-covariance matrix cannot be estimated, try MCPMod function from MCPMod package, else return NA
      if (!inherits(cvmat, "matrix")) {
        models <- list(linear=NULL)
        mm <- MCPMod::MCPMod(data,models,clinRel=0.4,doseEst="MED2")
        med <- mm$tdose
        med <- ifelse(length(med)==0,NA,as.numeric(med))
      } else {
        #compute lower confidence interval
        testLin <- function(td){ 
          coefs[2]*td-crt[1]*sqrt(sLin)*sqrt((fl(td)%*%cvmat%*%t(t(fl(td)))))
        }
        #check if lower confidence interval is above placebo estimate, otherwise perform gridsearch for next lowest dose with high enough lower confidence level
        if(testLin(med)<0){
          gridLin <- seq(med,mdose,by=(mdose-med)/100)
          outLin <- sapply(gridLin,testLin)
          if (sum(outLin>0)==0){
            med <- NA #NA if med is too large or does not exist
          } else {
            indexLin <- min(which(outLin>0))
            med <- gridLin[indexLin]
          }
        }
      }
    }
  }
  
  #power model
  if(model=="power" | model=="powermod"){
    med <- as.numeric((delta/coefs[2])^(1/coefs[3]))
    if(med>mdose){med <- NA} else {
      predPow <- coefs[1]+coefs[2]*dose^coefs[3]
      sPow <- sum((resp-predPow)^2)/(n-3)
      fp <- function(d){
        if(d!=0){
          c(1,d^coefs[3],coefs[2]*d^coefs[3]*log(d))
        } else {c(1,0,0)}
      }
      Fp <- t(sapply(dose,fp))
      cvmat <- try(solve(crossprod(Fp))*sPow,silent=TRUE)
      if (!inherits(cvmat, "matrix")) {
        models <- list(powermod=c(coefs[1],coefs[2],coefs[3]))
        start <- list(powermod=c(e0=coefs[1],e1=coefs[2],pow=coefs[3]))
        mm <- MCPMod::MCPMod(data,models,clinRel=0.4,doseEst="MED2",start=start,uGrad=powergrad)
        med <- mm$tdose
        med <- ifelse(length(med)==0,NA,as.numeric(med))
      }
      else {
        testPow <- function(td){ 
          coefs[2]*td^coefs[3]-crt[2]*sqrt(sPow)*sqrt((fp(td)%*%cvmat%*%t(t(fp(td)))))
        }  
        if(testPow(med)<0){
          gridPow <- seq(med,mdose,by=(mdose-med)/100)
          outPow <- sapply(gridPow,testPow)
          if (sum(outPow>0)==0){
            med <- NA
          } else {
            indexPow <- min(which(outPow>0))
            med <- gridPow[indexPow]
          }
        }
      }
    }
  }
  
  #4PL model
  if(model=="4PL" | model=="sigEmax"){
    if(coefs[1]+coefs[2]<coefs[1]+delta){med <- NA} else {
      med <- as.numeric((delta*coefs[3]^coefs[4]/(coefs[2]-delta))^(1/coefs[4]))
      if(med>mdose){med <- NA} else {
        predFpl <- coefs[1]+coefs[2]*dose^coefs[4]/(dose^coefs[4]+coefs[3]^coefs[4])
        sFpl <- sum((resp-predFpl)^2)/(n-4)
        ff <- function(d){
          c(1,d^coefs[4]/(d^coefs[4]+coefs[3]^coefs[4]),-d^coefs[4]*coefs[4]*coefs[3]^(coefs[4]-1)*coefs[2]/(d^coefs[4]+coefs[3]^coefs[4])^2,ifelse(d==0,0,coefs[2]*d^coefs[4]*coefs[3]^coefs[4]*log(d/coefs[3])/(d^coefs[4]+coefs[3]^coefs[4])^2))
        }
        Ff <- t(sapply(dose,ff))
        cvmat <- try(solve(crossprod(Ff))*sFpl,silent=TRUE)
        if (!inherits(cvmat, "matrix")) {
          models <- list(sigEmax=c(coefs[3],coefs[4]))
          mm <- MCPMod::MCPMod(data,models,clinRel=0.4,doseEst="MED2")
          med <- mm$tdose
          med <- ifelse(length(med)==0,NA,as.numeric(med))
        }
        else {
          testFpl <- function(td){ 
            coefs[2]*td^coefs[4]/(td^coefs[4]+coefs[3]^coefs[4])-crt[3]*sqrt(sFpl)*sqrt((ff(td)%*%cvmat%*%t(t(ff(td)))))
          }
          if(testFpl(med)<0){
            gridFpl <- seq(med,mdose,by=(mdose-med)/100)
            outFpl <- sapply(gridFpl,testFpl)
            if (sum(outFpl>0)==0){
              med <- NA
            } else {
              indexFpl <- min(which(outFpl>0))
              med <- gridFpl[indexFpl]
            }
          }
        }
      }
    }
  }
  
  #linlog model
  if(model=="linlog"){
    med <- as.numeric(exp(delta/coefs[2])-1)
    if(med>mdose){med <- NA} else {
      predLinl <- coefs[1]+coefs[2]*log(dose+1)
      sLinl <- sum((resp-predLinl)^2)/(n-2)
      fll <- function(d){
        c(1,log(d+1))
      }
      Fll <- t(sapply(dose,fll))
      cvmat <- try(solve(crossprod(Fll))*sLinl,silent=TRUE)
      if (!inherits(cvmat, "matrix")) {
        models <- list(linlog=NULL)
        mm <- MCPMod::MCPMod(data,models,clinRel=0.4,doseEst="MED2")
        med <- mm$tdose
        med <- ifelse(length(med)==0,NA,as.numeric(med))
      } else {
        testLinl <- function(td){ 
          coefs[2]*td-crt[1]*sqrt(sLinl)*sqrt((fll(td)%*%cvmat%*%t(t(fll(td)))))
        }
        if(testLinl(med)<0){
          gridLinl <- seq(med,mdose,by=(mdose-med)/100)
          outLinl <- sapply(gridLinl,testLinl)
          if (sum(outLinl>0)==0){
            med <- NA
          } else {
            indexLinl <- min(which(outLinl>0))
            med <- gridLinl[indexLinl]
          }
        }
      }
    }
  }
  
  #quadratic model
  if(model=="quadratic"){
    med <- polyroot(c(-0.4,coefs[2],coefs[3]))
    med <- Re(med[abs(Im(med))<0.0000001])
    if (length(med)==0){med <- NA} else {
      med <- min(med[med>0])
      if(med>mdose){med<-NA} else {
        predq <- coefs[1]+coefs[2]*dose+coefs[3]*dose^2
        sq <- sum((resp-predq)^2)/(n-3)
        fq <- function(d){
          c(1,d,d^2)
        }
        Fq <- t(sapply(dose,fq))
        cvmat <- try(solve(crossprod(Fq))*sq,silent=TRUE)
        if (!inherits(cvmat, "matrix")) {
          models <- list(quadratic=c(-0.854))
          mm <- MCPMod::MCPMod(data,models,clinRel=0.4,doseEst="MED2")
          med <- mm$tdose
          med <- ifelse(length(med)==0,NA,as.numeric(med))
        }
        else {
          testq <- function(td){ 
            coefs[2]*td-crt[2]*sqrt(sq)*sqrt((fq(td)%*%cvmat%*%t(t(fq(td)))))
          }
          if(testq(med)<0){
            gridq <- seq(med,mdose,by=(mdose-med)/100)
            outq <- sapply(gridq,testq)
            if (sum(outq>0)==0){
              med <- NA
            } else {
              indexq <- min(which(outq>0))
              med <- gridq[indexq]
            }
          }
        }
      }
    }
  }
  
  #for the full model an MED estimate is not given
  if(model=="full"){
    med <- NA
  }
  
  list("med"=round(med,2))
}

#OUTPUT:
#med: MED estiamte

###################################################################
# Simulate theoretical critical values under non-identifiability  #
###################################################################
#function to calculate the supremum of one realisation of the process G^2 outlined in Web Appendix B for the candidate set
#INPUT:
#dose: vector with dose levels
#p: p_i given in Web Appendix B, needed in case of unequal group sizes.
G <- function(dose=c(0,0.05,0.2,0.6,1),p=NULL){
  k <- length(dose)
  eps <- rnorm(k)
  if(is.null(p)){
    p <- rep(1,k)
  }
  
  #Linear model
  g1lin <- 1/k*sum(p*dose)
  g2lin <- 1/k*sum(p*(dose-g1lin)^2)
  glin <- (1/sqrt(k)*sum((dose-g1lin)/sqrt(g2lin)*sqrt(p)*eps))^2
  outlin <- glin
  
  #Power model
  g1pow <- function(theta){
    1/k*sum(p*dose^theta)
  }
  g2pow <- function(theta){
    1/k*sum(p*(dose^theta-g1pow(theta))^2)
  }
  gpow <- function(theta){
    -(1/sqrt(k)*sum((dose^theta-g1pow(theta))/sqrt(g2pow(theta))*sqrt(p)*eps))^2
  }
  #use bounded optimization to find the supremum of process gpow, because it might
  #not be possible to evaluate the process at the boundaries of the parameter space:
  lowerpow <- 0.00001
  upperpow <- 100
  svalpow <- 0.1
  outpow <- -optim(svalpow,gpow,lower=lowerpow,upper=upperpow,method="L-BFGS-B")$value
  
  #4PL model
  g14pl <- function(theta){
    1/k*sum(p*dose^theta[1]/(dose^theta[1]+theta[2]^theta[1]))
  }
  g24pl <- function(theta){
    1/k*sum(p*(dose^theta[1]/(dose^theta[1]+theta[2]^theta[1])-g14pl(theta))^2)
  }
  g4pl <- function(theta){
    -(1/sqrt(k)*sum((dose^theta[1]/(dose^theta[1]+theta[2]^theta[1])-g14pl(theta))/sqrt(g24pl(theta))*sqrt(p)*eps))^2
  }
  lower4pl <- c(0.00001,0.00001)
  upper4pl <- c(100,100)
  sval4pl <- c(0.1,diff(range(dose))/2)
  out4pl <- -optim(sval4pl,g4pl,lower=lower4pl,upper=upper4pl,method="L-BFGS-B")$value
  
  #full model
  outfull <- sum(eps^2)-1/k*(sum(p*eps))^2
  
  list("linear"=outlin,"power"=outpow,"4pl"=out4pl,"full"=outfull)
}

#OUTPUT:
#realisation of sup(G) for linear, power, emax, 4pl and unrestricted (full) model

#simulate distribution of supremum of G^2
simdist <- function(nsim=1000,dose=c(0,0.05,0.2,0.6,1),p=NULL){
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  out <- matrix(0,nsim,4)
  for(i in 1:nrow(out)){
    out[i,]<-unlist(G(dose=dose,p=p))
    setTxtProgressBar(pb, i)
  }
  out <- as.data.frame(out)
  colnames(out)<-c("linear","power","fpl","full")
  return(out)
}

#OUTPUT
#dataframe containing correlated distributions of sup(G) for evaluating the constant
#model against the alternative models in the candidate set