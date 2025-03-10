######################################################################
#------------------Part1: Function for MSDB prior--------------------#
######################################################################
#The sub-functions that constitute the MSDB prior function------------
PS_score_label_func <- function(data.X.Z,cutgroup_n=5){
  data.X.Z$Z <- as.factor(data.X.Z$Z)
  data.X.Z$Z <- relevel(data.X.Z$Z, ref = "CTD")
  X.names <- grep("[Xx]", colnames(data.X.Z), value = TRUE, ignore.case = TRUE)
  PSformula <- paste(c("Z","~",paste(X.names,collapse = "+")),collapse = "") %>% as.formula()
  PSmodel <- multinom(PSformula,data = data.X.Z)
  PS_score <- fitted(PSmodel)
  data.X.Z$PS_score <- PS_score[,"CTD"]
  
  PS_label_func_all <- function(Z,PS_score,cutgroup_n=cutgroup_n){
    PS_label_func <- function(PS_score,cutpoint=NA,cutgroup_n=5){
      if(is.na(cutpoint[1])){
        cutpoint <- sort(PS_score)[seq(from=1,to=length(PS_score),length=(cutgroup_n+1)) %>% round(0)] %>% round(4)
      }
      PS_label <- cut(PS_score, cutpoint, labels = 1:cutgroup_n, include.lowest = TRUE)
      Ps_cutpoint <- cutpoint
      return(list(PS_label=PS_label,Ps_cutpoint=Ps_cutpoint))
    }
    c(PS_label_CTD,Ps_cutpoint) %<-% PS_label_func(PS_score[Z=="CTD"],cutgroup_n=cutgroup_n)
    c(PS_label_ETD,a) %<-% PS_label_func(PS_score[Z=="ETD"],cutpoint=Ps_cutpoint,cutgroup_n)
    c(PS_label_RWD,a) %<-% PS_label_func(PS_score[Z=="RWD"],cutpoint=Ps_cutpoint,cutgroup_n)
    return(c(PS_label_CTD,PS_label_ETD,PS_label_RWD))
  }
  data.X.Z$PS_label<- PS_label_func_all(Z=data.X.Z$Z,PS_score=data.X.Z$PS_score,cutgroup_n=cutgroup_n)
  return(data.X.Z)
}
K.label.func.3 <- function(data.T.d,data.X.Z,K=4){
  data.T.d.S.K <- merge(data.T.d,data.X.Z[,c("ID","PS_score","PS_label")],by = "ID")
  data.T.d.S.K <- data.T.d.S.K[is.na(data.T.d.S.K$PS_label)==F,]
  cut.time.func <- function(data.T.d,K=4){
    tdata.CTD <- data.T.d[data.T.d$Z=="CTD",]
    tdata.CTD <- tdata.CTD[order(tdata.CTD$Time),]
    time.max <- max(tdata.CTD$Time)
    tdata.CTD <- tdata.CTD[tdata.CTD$Time!=time.max,]
    n.death.all <- tdata.CTD$delta %>% sum()
    k.length <- round(n.death.all/K,digits = 0)
    count.point <- (1:K)*k.length
    cumsum.T <- cumsum(tdata.CTD$delta)
    cut.time <- vector()
    for (k in 1:(K-1)) {
      cut.time[k] <- (tdata.CTD$Time[cumsum.T==count.point[k]][1]+tdata.CTD$Time[cumsum.T==(count.point[k]+1)][1])/2 
    }
    cut.time <- c(0,cut.time,time.max)
    return(cut.time)
  }
  PWE.label.func <- function(data.T.d,cut.time,K=4){
    data.T.d$Time.label <- cut(data.T.d$Time, cut.time , labels = 1:K, include.lowest = TRUE)
    return(data.T.d)
  }
  cut.time <- cut.time.func(data.T.d,K=K)
  data.T.d.s.K <-  PWE.label.func(data.T.d = data.T.d.S.K,cut.time,K=K)
  names(cut.time) <- paste("cut.time",1:length(cut.time),sep = ".")
  return(list(data.T.d=data.T.d.s.K,cut.time=cut.time))
}
rE.func.for.Z.PS.K2 <- function(data.T.d,data.X.Z,cut.time,K=4){
  data.T.d.S.K <- data.T.d
  data.T.d.S.K <- data.T.d.S.K[is.na(data.T.d.S.K$Time.label)==F,]
  rE.func <- function(data.T.d,cut.time){
    r <- data.T.d %>% group_by(Time.label) %>% summarise(r=sum(delta))
    rr <- data.frame(Time.label=factor(1:K),r=rep(0,K)) 
    r <- bind_rows(r, rr %>% anti_join(r, by = "Time.label"))
    
    E_k <- vector()
    N <- vector()
    for (k in 1:K) {
      e1 <- sum(data.T.d$Time[as.numeric(data.T.d$Time.label)==k]-cut.time[k])
      e2 <- sum(as.numeric(data.T.d$Time.label)>k)*(cut.time[k+1]-cut.time[k])####杩欓噷鏄竴姝ラ噸瑕佺殑淇�??
      N[k] <- sum(as.numeric(data.T.d$Time.label)>=k)
      E_k[k] <- sum(c(e1,e2),na.rm = TRUE)
    }
    a <- data.frame(r=r,E=E_k,N=N,
                    Z=data.T.d$Z%>% unique(),
                    PS.label=data.T.d$PS_label%>% unique(),
                    group=data.T.d$group%>% unique())
    colnames(a) <- c("Time.label","r","E","N","Z","PS_label","group")
    return(a)
  }
  split_data.T.d.S.K_by.Z.PS.group <- split(data.T.d.S.K, list(data.T.d.S.K$Z, data.T.d.S.K$PS_label,data.T.d.S.K$group))
  rE_Z.PS.K <- list()
  for (i in 1:length(split_data.T.d.S.K_by.Z.PS.group)) {
    rE_Z.PS.K[[i]] <- rE.func(split_data.T.d.S.K_by.Z.PS.group[[i]],cut.time)
  }
  names(rE_Z.PS.K) <- names(split_data.T.d.S.K_by.Z.PS.group)
  b <- data.frame()
  for (i in 1:length(rE_Z.PS.K)) {
    b <- rbind(b,rE_Z.PS.K[[i]]) 
  }
  return(list(data.T.d=data.T.d.S.K,data.rE=b))
}
tauk.Z.K.S.func2 <- function(data.rE){
  split.data.rE <- split(data.rE,list(data.rE$PS_labe,data.rE$Time.label))
  edata <- list()
  for (i in 1:length(split.data.rE)) {
    datai=split.data.rE[[i]]
    rr=datai$r[datai$Z=="CTD"]
    EE=datai$E[datai$Z=="CTD"]
    if(datai$r[datai$Z=="ETD"]==0){P_value.ETD=NA}else{
      lambda.ETD <- mixgamma(c(1,datai$r[datai$Z=="ETD"],datai$E[datai$Z=="ETD"]))
      lambda.ETD.post <- postmix(lambda.ETD,m=rr/EE,n=EE)
      rmap_pred <- RBesT::preddist(lambda.ETD.post, n = EE)
      p_lower <- RBesT::pmix(rmap_pred, rr)
      P_value.ETD <- ifelse(p_lower < 0.5, 2 * p_lower, 
                            2 * (1 - p_lower))
    }
    if(datai$r[datai$Z=="RWD"]==0){P_value.RWD=NA}else{
      lambda.RWD <- mixgamma(c(1,datai$r[datai$Z=="RWD"],datai$E[datai$Z=="RWD"]))
      lambda.RWD.post <- postmix(lambda.RWD,m=rr/EE,n=EE)
      rmap_pred <- RBesT::preddist(lambda.RWD.post, n = EE)
      p_lower <- RBesT::pmix(rmap_pred, rr)
      P_value.RWD <- ifelse(p_lower < 0.5, 2 * p_lower, 
                            2 * (1 - p_lower))
    }
    edata[[i]] <- data.frame(e.ETD=P_value.ETD,e.RWD=P_value.RWD,
                             Time.label=datai$Time.label[1],PS_label=datai$PS_label[1])
  }
  tau <- data.frame(matrix(unlist(edata), ncol = 4, byrow = TRUE))
  colnames(tau) <- c("e.ETD","e.RWD","Time.label","PS_label")
  tau.median <- median(c(tau$e.RWD,tau$e.ETD),na.rm=T)
  tau.ETD <- data.frame(tau=c(tau$e.ETD/tau.median),Z="ETD",PS_label=tau$PS_label,Time.label=tau$Time.label)
  tau.RWD <- data.frame(tau=c(tau$e.RWD/tau.median),Z="RWD",PS_label=tau$PS_label,Time.label=tau$Time.label)
  tau=rbind(tau.ETD,tau.RWD)
  return(tau)
}
W.K.S.func2 <- function(data.rE,r0to=0.001){
  split.data.rE <- split(data.rE,list(data.rE$PS_labe,data.rE$Time.label))
  edata <- list()
  for (i in 1:length(split.data.rE)) {
    datai=split.data.rE[[i]]
    rr=datai$r[datai$Z=="CTD"]
    EE=datai$E[datai$Z=="CTD"]
    if(datai$r[datai$Z=="ETD"]==0){P_value.ETD=NA}else{
      lambda.ETD <- mixgamma(c(1,datai$r[datai$Z=="ETD"],datai$E[datai$Z=="ETD"]))
      lambda.ETD.post <- postmix(lambda.ETD,m=rr/EE,n=EE)
      rmap_pred <- RBesT::preddist(lambda.ETD.post, n = EE)
      p_lower <- RBesT::pmix(rmap_pred, rr)
      P_value.ETD <- ifelse(p_lower < 0.5, 2 * p_lower, 
                            2 * (1 - p_lower))
    }
    if(datai$r[datai$Z=="RWD"]==0){P_value.RWD=NA}else{
      lambda.RWD <- mixgamma(c(1,datai$r[datai$Z=="RWD"],datai$E[datai$Z=="RWD"]))
      lambda.RWD.post <- postmix(lambda.RWD,m=rr/EE,n=EE)
      rmap_pred <- RBesT::preddist(lambda.RWD.post, n = EE)
      p_lower <- RBesT::pmix(rmap_pred, rr)
      P_value.RWD <- ifelse(p_lower < 0.5, 2 * p_lower, 
                            2 * (1 - p_lower))
    }
    edata[[i]] <- data.frame(e.ETD=P_value.ETD,e.RWD=P_value.RWD,
                             Time.label=datai$Time.label[1],PS_label=datai$PS_label[1])
  }
  e.K.S.data <- data.frame(matrix(unlist(edata), ncol = 4, byrow = TRUE))
  colnames(e.K.S.data) <- c("e.ETD","e.RWD","Time.label","PS_label")
  
  data.N <- data.rE[,c("Time.label","PS_label","Z","N")]
  data.N <- pivot_wider(data.N,names_from = "Z", values_from = "N")
  data.m <- merge(e.K.S.data,data.N)
  w.func <- function(data.m){
    fenmu <- apply(data.frame(ETD=data.m$ETD*data.m$e.ETD,RWD=data.m$RWD*data.m$e.RWD), 1, sum,na.rm=TRUE)
    w.ETD <- (data.m$ETD*data.m$e.ETD)/(fenmu)
    w.RWD <- (data.m$RWD*data.m$e.RWD)/(fenmu)
    w.ETD[is.na(w.ETD)]=0
    w.RWD[is.na(w.RWD)]=0
    return(data.frame(Time.label=data.m$Time.label,PS_label=data.m$PS_label,
                      w.ETD=w.ETD,w.RWD=w.RWD))
  }
  w.Z.K.S<- w.func(data.m)
  return(w.Z.K.S)
}
lambda.prior.MCMC.func4 <- function(data.rE,data.R.Z.K.S,w.Z.K.S,w.for.ps=NA,seed=231115){
  library(cmdstanr)
  set_cmdstan_path("D:/cmdstan-2.35.0")
  split.data.R <- split(data.R.Z.K.S,list(data.R.Z.K.S$PS_label,data.R.Z.K.S$Time.label))
  split.data.rE <- split(data.rE,list(data.rE$PS_label,data.rE$Time.label))
  split.data.w <- split(w.Z.K.S,list(w.Z.K.S$PS_label,w.Z.K.S$Time.label))
  KK <- data.R.Z.K.S$Time.label %>% max()
  SS <- data.R.Z.K.S$PS_label %>% as.numeric() %>% max()
  modelcode_ETD <- {"
        data {
          real<lower=0> E_ETD;
          int<lower=0> r_ETD;
          real prior_prec_tau_ETD;}
        parameters {
          real logit_lambda_ETD;
          real mu;
          real<lower=0> tau_ETD;}
        model {
          // prior distributions
          tau_ETD ~ normal(0, prior_prec_tau_ETD)T[0,];
          mu ~ normal(0,10);
          // binomial likelihood
          logit_lambda_ETD ~ normal(mu,tau_ETD);
          r_ETD ~ poisson(exp(logit_lambda_ETD)*E_ETD);}"}
  modelcode_RWD <- {"
        data {
          real<lower=0> E_RWD;
          int<lower=0> r_RWD;
          real prior_prec_tau_RWD;}
        parameters {
          real logit_lambda_RWD;
          real mu;
          real<lower=0> tau_RWD;}
        model {
          // prior distributions
          tau_RWD ~ normal(0, prior_prec_tau_RWD)T[0,];
          mu ~ normal(0,10);
          // binomial likelihood
          logit_lambda_RWD ~ normal(mu,tau_RWD);
          r_RWD ~ poisson(exp(logit_lambda_RWD)*E_RWD);}"}
  mod_modelcode_ETD <- cmdstan_model(write_stan_file(modelcode_ETD))
  mod_modelcode_RWD <- cmdstan_model(write_stan_file(modelcode_RWD))
  
  lambda.SK.prior.MCMC <- list()
  forfunction <- function(i){
    set_cmdstan_path("D:/cmdstan-2.35.0")
    data.rE.for.model <- split.data.rE[[i]]
    data.tau.for.model <- split.data.R[[i]]
    data.w.for.model<- split.data.w[[i]]
    model_data_ETD <- list(
      E_ETD=data.rE.for.model[data.rE.for.model$Z=="ETD","E"],
      r_ETD=data.rE.for.model[data.rE.for.model$Z=="ETD","r"],
      prior_prec_tau_ETD=data.tau.for.model[data.tau.for.model$Z=="ETD","tau"]
    )
    model_data_RWD <- list(
      E_RWD=data.rE.for.model[data.rE.for.model$Z=="RWD","E"],
      r_RWD=data.rE.for.model[data.rE.for.model$Z=="RWD","r"],
      prior_prec_tau_RWD=data.tau.for.model[data.tau.for.model$Z=="RWD","tau"]
    )
    if(model_data_ETD$prior_prec_tau_ETD %>% is.na()){
      MCMC_lambda_ETD <- NA
    }else{
      if(model_data_ETD$E_ETD==0){MCMC_lambda_ETD <- NA}else{
        fit.stan_ETD <- mod_modelcode_ETD$sample(data = model_data_ETD,iter_sampling = 10000, chains = 4, parallel_chains = 4)
        MCMC_logit_lambda_ETD <- fit.stan_ETD$draws(variables = "logit_lambda_ETD") %>% as.vector()
        MCMC_lambda_ETD <- exp(MCMC_logit_lambda_ETD)}
    }
    if(model_data_RWD$prior_prec_tau_RWD %>% is.na()){
      MCMC_lambda_RWD <- NA
    }else{
      if(model_data_RWD$E_RWD==0){MCMC_lambda_RWD <- NA}else{
        fit.stan_RWD <- mod_modelcode_RWD$sample(data = model_data_RWD,iter_sampling = 10000, chains = 4, parallel_chains = 4)
        MCMC_logit_lambda_RWD <- fit.stan_RWD$draws(variables = "logit_lambda_RWD") %>% as.vector()
        MCMC_lambda_RWD <- exp(MCMC_logit_lambda_RWD)}
    }
    w.ETD <- data.w.for.model[,"w.ETD"] 
    w.RWD <- data.w.for.model[,"w.RWD"] 
    if(is.na(w.ETD)){MCMC_lambda_ETD <- NA}
    if(is.na(w.RWD)){MCMC_lambda_RWD <- NA}
    lambda.prior.MCMC <- NA
    if(is.na(MCMC_lambda_ETD[1])&is.na(MCMC_lambda_RWD[1])){lambda.prior.MCMC <- NA}else{
      lambda.prior.MCMC <- apply(data.frame(MCMC_ETD=w.ETD*MCMC_lambda_ETD,
                                            MCMC_RWD=w.RWD*MCMC_lambda_RWD), 1, sum,na.rm=TRUE)}
    return(lambda.prior.MCMC)
  }
  lambda.SK.prior.MCMC <- foreach(
    i=names(split.data.rE),   
    .combine=cbind,  
    .packages = c("MASS", "RBesT","rstan","dplyr","VGAM","tidyr","cmdstanr") 
  ) %do%  forfunction(i)
  
  data <- data.rE[data.rE$Z=="CTD",]
  if(is.list(w.for.ps)==FALSE){
    split.data <- split(data,list(data$Time.label))
    for (i in 1:length(split.data)) {
      split.data[[i]]["N"] <- split.data[[i]]["N"]/sum(split.data[[i]]["N"])
      colnames(split.data[[i]])[which(colnames(split.data[[i]])=="N")] <- "w.for.PS"
    }}else{
      split.data <- w.for.ps
    }
  lambda.K.prior.MCMC <- list()
  for (i in 1:KK) {
    p <- (1:SS)+SS*(i-1)
    lambda.K.prior.MCMC[[i]] <- apply(lambda.SK.prior.MCMC[,p], 1, function(row){ a <-  sum(row * split.data[[i]][["w.for.PS"]],na.rm = TRUE);return(a)})
  }
  lambda.K.prior.MCMC <- do.call(data.frame, lambda.K.prior.MCMC)
  colnames(lambda.K.prior.MCMC) <- 1:ncol(lambda.K.prior.MCMC)
  return(list(lambda.K.prior.MCMC=lambda.K.prior.MCMC,w.for.ps=split.data))
}
logHR_distr_CTD_func <- function(data.T.d,data.rE){
  data.T.d_CTD <- data.T.d[data.T.d$Z=="CTD",]
  data.rE_CTD <- data.rE[data.rE$Z=="CTD",] %>% group_by(Time.label,group) %>% summarise(rsum=sum(r),Esum=sum(E),Nsum=sum(N))
  logHR_distr_CTD <- list()
  mean_logHR_K_CTD <- vector()
  sigma_logHR_K_CTD <- vector()
  for(i in 1:K){
    lambda.T.mcmc<- rmix(mixgamma(c(1,data.rE_CTD$rsum[data.rE_CTD$Time.label==i&data.rE_CTD$group==1],
                                    data.rE_CTD$Esum[data.rE_CTD$Time.label==i&data.rE_CTD$group==1])),100000) 
    lambda.C.mcmc<- rmix(mixgamma(c(1,data.rE_CTD$rsum[data.rE_CTD$Time.label==i&data.rE_CTD$group==0],
                                    data.rE_CTD$Esum[data.rE_CTD$Time.label==i&data.rE_CTD$group==0])),100000) 
    logHR_distr_CTD[[i]]<- mixfit(log(lambda.T.mcmc/lambda.C.mcmc),type = "norm",Nc=1)
    mean_logHR_K_CTD[i] <- logHR_distr_CTD[[i]]["m",]
    sigma_logHR_K_CTD[i] <- logHR_distr_CTD[[i]]["s",]
  }
  mean_bind <- sum(mean_logHR_K_CTD*(1/sigma_logHR_K_CTD^2))/sum(1/sigma_logHR_K_CTD^2)
  sigma_bind <- sqrt(1/sum(1/sigma_logHR_K_CTD^2))
  logHR_distr_CTD <- mixnorm(c(1,mean_bind,sigma_bind))
  return(logHR_distr_CTD)
}
lambda.MCMC.to.distr2 <- function(lambda.K.prior.MCMC,N.mix=3){
  mixgamma.K <- list()
  mixgamma.distr.K <- list()
  for (k in 1:ncol(lambda.K.prior.MCMC)) {
    mixgamma <- mixfit(lambda.K.prior.MCMC[,k],type = "gamma",Nc=N.mix)
    mixgamma.K[[k]] <- mixgamma
  }
  return(mixgamma.K=mixgamma.K)
}
lambda.unimform.func <- function(data.T.d,cut.time){
  data.T.d.CTD <- data.T.d[data.T.d$Z=="CTD",]
  fit <- survreg(Surv(Time, delta) ~ 1, data = data.T.d.CTD, dist = "weibull")
  shape <- 1/fit$scale   
  scale <- exp(fit$coefficients[1]) 
  rate <- 1/scale 
  lambda.unimform <- vector()
  cut.time.1 <- cut.time
  for (i in 2:length(cut.time.1)) {
    lambda.unimform[i-1] <- ((cut.time.1[i]^shape-cut.time.1[i-1]^shape)/((scale^shape)*(cut.time.1[i]-cut.time.1[i-1]))) %>% as.numeric()
  }
  return(lambda.unimform)
}
Wvfunc.lambda.k.SAM <- function(mixgamma.K,data.rE.group,K=4,SAM_delta=0.1){
  data <- data.rE.group[data.rE.group$Z=="CTD",] %>% group_by(Time.label) %>% summarise(rsum=sum(r),Esum=sum(E))
  yangb.distr <- list()
  WV.myMAP.k1 <- list()
  for (k in 1:K) {
    lambda.CTD.K.distr <- mixgamma(c(1,data$rsum[k],data$Esum[k]))
    yangb.distr[[k]] <- lambda.CTD.K.distr
    mixgamma.KK <- mixgamma.K[[k]]
    uninforma.distr <- mixgamma(c(1,(data$rsum[k]/data$Esum[k]),1),param = "mn")
    WV.myMAP.k1[[k]] <- 
      foreach(SAM_deltai=SAM_delta,.combine = rbind,.packages = "SAMprior"
      )%do%{ wv<- SAM_weight(if.prior = mixgamma.KK,delta = SAM_deltai,data = rbind(data$rsum[k],data$Esum[k]))
      }
  }
  WV.matrix <- matrix(unlist(WV.myMAP.k1),nrow=K,byrow = T)
  result_func <- function(i){
    lambda.prior.MY.MAP1 <- list()
    lambda.post.MY.MAP1 <- list()
    WV.myMAP1 <- vector()
    mean.MyMAP.post1 <- vector()
    ess_group1=vector()
    WV_iK <- as.vector(WV.matrix[,i])
    for (k in 1:K) {
      WV_i <- WV_iK[k]
      mixgamma.KK <- mixgamma.K[[k]]
      uninforma.distr <- mixgamma(c(1,(data$rsum[k]/data$Esum[k]),1),param = "mn")
      if(WV_i==1){ lambda.prior.k.MY.MAP1=uninforma.distr
      }else{if(WV_i==0){lambda.prior.k.MY.MAP1=mixgamma.KK
      }else{lambda.prior.k.MY.MAP1 <- mixcombine(mixgamma.KK,uninforma.distr,weight=c(1-WV_i, WV_i))}}
      lambda.prior.MY.MAP1[[k]] <- lambda.prior.k.MY.MAP1
      likelihood(lambda.prior.k.MY.MAP1) <- "poisson"
      lambda.post.k.MY.MAP1 <- postmix(lambda.prior.k.MY.MAP1,m=(data$rsum[k]/data$Esum[k]),n=data$Esum[k])
      lambda.post.MY.MAP1[[k]] <- lambda.post.k.MY.MAP1
      WV.myMAP1[k] <- WV_i
      mean.MyMAP.post1[k] <- summary(lambda.post.k.MY.MAP1)["mean"]
      likelihood(lambda.post.k.MY.MAP1)="exp"
      ess_group1[k]=ess(lambda.post.k.MY.MAP1,method = "elir") 
    }
    ESS_borrow1=sum(ess_group1)-sum(data$rsum)
    return(list(yangb.distr=yangb.distr,
                lambda.prior.MY.MAP1=lambda.prior.MY.MAP1,WV.myMAP1=WV.myMAP1,
                lambda.post.MY.MAP1=lambda.post.MY.MAP1,mean.MyMAP.post1=mean.MyMAP.post1,
                ESS_borrow1=ESS_borrow1))}
  result_all<- foreach(i=1:ncol(WV.matrix),.packages = c("RBesT"))%do%result_func(i)
  return(result_all)
}
Wvfunc.lambda.k.BOX <- function(mixgamma.K,data.rE.group,K=4,P_cutpoint=0.7){
  data <- data.rE.group[data.rE.group$Z=="CTD",] %>% group_by(Time.label) %>% summarise(rsum=sum(r),Esum=sum(E))
  yangb.distr <- list()
  WV.myMAP.k1 <- list()

  W_byBOX_Pvalue<- function (map_obj, vague_prior, E, r,P_cutpoint) {
    library(dplyr)
    library(RBesT)
    library(foreach)
    Wvfunction <- function(W){
      rmap <-  mixcombine(map_obj,vague_prior,weight=c(1-W, W))
      likelihood(rmap) <- "poisson"
      rmap_pred <- RBesT::preddist(rmap, n = E)
      x <- RBesT::rmix(rmap_pred,30000)
      d_value <- dmix(rmap_pred, x)
      d_obs <- dmix(rmap_pred,r)
      P_value <- mean(d_value<d_obs)
      return(P_value)
    }
    P_value_vec <- foreach(W=seq(0,1,0.01),.combine = rbind,.packages = "RBesT")%do%Wvfunction(W)
    pvalue_func <- function(P_cutpoint){
      Wv=seq(from=0,to=1,by=0.01)
      if(length(Wv[P_value_vec>P_cutpoint])>0){
        WV.myMAP.k1<- Wv[P_value_vec>P_cutpoint] %>% min
      }else{WV.myMAP.k1<- 1}
      return(WV.myMAP.k1)
    }
    WV.myMAP.k1 <- foreach(P_cutpointi=P_cutpoint,.combine = cbind)%do%pvalue_func(P_cutpointi)
    return(WV.myMAP.k1)
  }
  for (k in 1:K) {
    lambda.CTD.K.distr <- mixgamma(c(1,data$rsum[k],data$Esum[k]))
    yangb.distr[[k]] <- lambda.CTD.K.distr
    mixgamma.KK <- mixgamma.K[[k]]
    uninforma.distr <- mixgamma(c(1,(data$rsum[k]/data$Esum[k]),1),param = "mn")
    WV.myMAP.k1[[k]] <- W_byBOX_Pvalue(mixgamma.KK, uninforma.distr, E=data$Esum[k], r=data$rsum[k],P_cutpoint)
  }
  WV.matrix <- matrix(unlist(WV.myMAP.k1),nrow=K,byrow = T)
  result_func <- function(i){
    lambda.prior.MY.MAP1 <- list()
    lambda.post.MY.MAP1 <- list()
    WV.myMAP1 <- vector()
    mean.MyMAP.post1 <- vector()
    ess_group1=vector()
    WV_iK <- as.vector(WV.matrix[,i])
    for (k in 1:K) {
      WV_i <- WV_iK[k]
      mixgamma.KK <- mixgamma.K[[k]]
      uninforma.distr <- mixgamma(c(1,(data$rsum[k]/data$Esum[k]),1),param = "mn")
      if(WV_i==1){ lambda.prior.k.MY.MAP1=uninforma.distr
      }else{if(WV_i==0){lambda.prior.k.MY.MAP1=mixgamma.KK
      }else{lambda.prior.k.MY.MAP1 <- mixcombine(mixgamma.KK,uninforma.distr,weight=c(1-WV_i, WV_i))}}
      lambda.prior.MY.MAP1[[k]] <- lambda.prior.k.MY.MAP1
      likelihood(lambda.prior.k.MY.MAP1) <- "poisson"
      lambda.post.k.MY.MAP1 <- postmix(lambda.prior.k.MY.MAP1,m=(data$rsum[k]/data$Esum[k]),n=data$Esum[k])
      lambda.post.MY.MAP1[[k]] <- lambda.post.k.MY.MAP1
      WV.myMAP1[k] <- WV_i
      mean.MyMAP.post1[k] <- summary(lambda.post.k.MY.MAP1)["mean"]
      likelihood(lambda.post.k.MY.MAP1)="exp"
      ess_group1[k]=ess(lambda.post.k.MY.MAP1,method = "elir") 
    }
    ESS_borrow1=sum(ess_group1)-sum(data$rsum)
    return(list(yangb.distr=yangb.distr,
                lambda.prior.MY.MAP1=lambda.prior.MY.MAP1,WV.myMAP1=WV.myMAP1,
                lambda.post.MY.MAP1=lambda.post.MY.MAP1,mean.MyMAP.post1=mean.MyMAP.post1,
                ESS_borrow1=ESS_borrow1))}
  result_all<- foreach(i=1:ncol(WV.matrix),.packages = c("RBesT"))%do%result_func(i)
  return(result_all)
}
Wvfunc.lambda.k.EBrmap <- function(mixgamma.K,data.rE.group,K=4,P_cutpoint=0.7){
  data <- data.rE.group[data.rE.group$Z=="CTD",] %>% group_by(Time.label) %>% summarise(rsum=sum(r),Esum=sum(E))
  
  yangb.distr <- list()
  WV.myMAP.k1 <- list()

  EB_rMAP_pallar<- function (map_obj, vague_prior, E, r,P_cutpoint) {
    library(dplyr)
    library(RBesT)
    library(foreach)
    Wvfunction <- function(W){
      rmap <-  mixcombine(map_obj,vague_prior,weight=c(1-W, W))
      likelihood(rmap) <- "poisson"
      rmap_pred <- RBesT::preddist(rmap, n = E)
      p_lower <- RBesT::pmix(rmap_pred, r)
      P_value <- ifelse(p_lower < 0.5, 2 * p_lower, 
                        2 * (1 - p_lower))
      return(P_value)
    }
    P_value_vec <- foreach(W=seq(0,1,0.01),.combine = rbind,.packages = "RBesT")%do%Wvfunction(W)
    pvalue_func <- function(P_cutpoint){
      Wv=seq(from=0,to=1,by=0.01)
      if(length(Wv[P_value_vec>P_cutpoint])>0){
        WV.myMAP.k1<- Wv[P_value_vec>P_cutpoint] %>% min
      }else{WV.myMAP.k1<- 1}
      return(WV.myMAP.k1)
    }
    WV.myMAP.k1 <- foreach(P_cutpointi=P_cutpoint,.combine = cbind)%do%pvalue_func(P_cutpointi)
    return(WV.myMAP.k1)
  }
  for (k in 1:K) {
    lambda.CTD.K.distr <- mixgamma(c(1,data$rsum[k],data$Esum[k]))
    yangb.distr[[k]] <- lambda.CTD.K.distr
    mixgamma.KK <- mixgamma.K[[k]]
    uninforma.distr <- mixgamma(c(1,(data$rsum[k]/data$Esum[k]),1),param = "mn")
    WV.myMAP.k1[[k]] <- EB_rMAP_pallar(mixgamma.KK, uninforma.distr, E=data$Esum[k], r=data$rsum[k],P_cutpoint)
  }
  WV.matrix <- matrix(unlist(WV.myMAP.k1),nrow=K,byrow = T)
  result_func <- function(i){
    lambda.prior.MY.MAP1 <- list()
    lambda.post.MY.MAP1 <- list()
    WV.myMAP1 <- vector()
    mean.MyMAP.post1 <- vector()
    ess_group1=vector()
    WV_iK <- as.vector(WV.matrix[,i])
    for (k in 1:K) {
      WV_i <- WV_iK[k]
      mixgamma.KK <- mixgamma.K[[k]]
      uninforma.distr <- mixgamma(c(1,(data$rsum[k]/data$Esum[k]),1),param = "mn")
      if(WV_i==1){ lambda.prior.k.MY.MAP1=uninforma.distr
      }else{if(WV_i==0){lambda.prior.k.MY.MAP1=mixgamma.KK
      }else{lambda.prior.k.MY.MAP1 <- mixcombine(mixgamma.KK,uninforma.distr,weight=c(1-WV_i, WV_i))}}
      lambda.prior.MY.MAP1[[k]] <- lambda.prior.k.MY.MAP1
      likelihood(lambda.prior.k.MY.MAP1) <- "poisson"
      lambda.post.k.MY.MAP1 <- postmix(lambda.prior.k.MY.MAP1,m=(data$rsum[k]/data$Esum[k]),n=data$Esum[k])
      lambda.post.MY.MAP1[[k]] <- lambda.post.k.MY.MAP1
      WV.myMAP1[k] <- WV_i
      mean.MyMAP.post1[k] <- summary(lambda.post.k.MY.MAP1)["mean"]
      likelihood(lambda.post.k.MY.MAP1)="exp"
      ess_group1[k]=ess(lambda.post.k.MY.MAP1,method = "elir") 
    }
    ESS_borrow1=sum(ess_group1)-sum(data$rsum)
    return(list(yangb.distr=yangb.distr,
                lambda.prior.MY.MAP1=lambda.prior.MY.MAP1,WV.myMAP1=WV.myMAP1,
                lambda.post.MY.MAP1=lambda.post.MY.MAP1,mean.MyMAP.post1=mean.MyMAP.post1,
                ESS_borrow1=ESS_borrow1))}
  result_all<- foreach(i=1:ncol(WV.matrix),.packages = c("RBesT"))%do%result_func(i)
  return(result_all)
}
Wvfunc.lambda.k.myPvalue2 <- function(mixgamma.K,data.rE.group,K=4,P_cutpoint=0.7){
  data <- data.rE.group[data.rE.group$Z=="CTD",] %>% group_by(Time.label) %>% summarise(rsum=sum(r),Esum=sum(E))
  
  yangb.distr <- list()
  WV.myMAP.k1 <- list()
  
  W_by_myPvalue <- function (map_obj, vague_prior, E, r,P_cutpoint) {
    library(dplyr)
    library(RBesT)
    library(foreach)
    Wvfunction <- function(W){
      rmap <-  mixcombine(map_obj,vague_prior,weight=c(1-W, W))
      likelihood(rmap) <- "poisson"
      likelihood(map_obj) <- "poisson"
      likelihood(vague_prior) <- "poisson"
      rmappost <- postmix(rmap,m=r/E,n=E)
      rmap_pred <- RBesT::preddist(rmappost, n = E)
      p_lower <- RBesT::pmix(rmap_pred, r)
      P_value <- ifelse(p_lower < 0.5, 2 * p_lower, 
                        2 * (1 - p_lower))
      if(W==0){rmappost=postmix(map_obj,m=r/E,n=E)}
      if(W==1){rmappost=postmix(vague_prior,m=r/E,n=E)}
      likelihood(rmappost) <- "exp"
      ess <- ess(rmappost,"elir")
      return(c(P_value,ess))
    }
    P_value_vec <- foreach(W=seq(0,1,0.01),.combine = rbind,.packages = "RBesT")%do%Wvfunction(W)
    pvalue_func <- function(P_cutpoint){
      Wv=seq(from=0,to=1,by=0.01)
      if(length(Wv[P_value_vec[,1]>P_cutpoint])>0){
        lambdaGdata <- data.frame(Wv,P_value_vec)
        colnames(lambdaGdata) <- c("Wv","P","ess")
        lambdaGdata <- lambdaGdata %>% filter(P>P_cutpoint)
        WV.myMAP.k1<- lambdaGdata$Wv[which.max(lambdaGdata$ess)] 
      }else{WV.myMAP.k1<- 1}
      return(WV.myMAP.k1)
    }
    WV.myMAP.k1 <- foreach(P_cutpointi=P_cutpoint,.combine = cbind)%do%pvalue_func(P_cutpointi)
    return(WV.myMAP.k1)
  }
  
  for (k in 1:K) {
    lambda.CTD.K.distr <- mixgamma(c(1,data$rsum[k],data$Esum[k]))
    yangb.distr[[k]] <- lambda.CTD.K.distr
    mixgamma.KK <- mixgamma.K[[k]]
    uninforma.distr <- mixgamma(c(1,(data$rsum[k]/data$Esum[k]),1),param = "mn")
    WV.myMAP.k1[[k]] <- W_by_myPvalue(mixgamma.KK, uninforma.distr, E=data$Esum[k], r=data$rsum[k],P_cutpoint)
  }
  WV.matrix <- matrix(unlist(WV.myMAP.k1),nrow=K,byrow = T)
  result_func <- function(i){
    lambda.prior.MY.MAP1 <- list()
    lambda.post.MY.MAP1 <- list()
    WV.myMAP1 <- vector()
    mean.MyMAP.post1 <- vector()
    ess_group1=vector()
    WV_iK <- as.vector(WV.matrix[,i])
    for (k in 1:K) {
      WV_i <- WV_iK[k]
      mixgamma.KK <- mixgamma.K[[k]]
      uninforma.distr <- mixgamma(c(1,(data$rsum[k]/data$Esum[k]),1),param = "mn")
      if(WV_i==1){ lambda.prior.k.MY.MAP1=uninforma.distr
      }else{if(WV_i==0){lambda.prior.k.MY.MAP1=mixgamma.KK
      }else{lambda.prior.k.MY.MAP1 <- mixcombine(mixgamma.KK,uninforma.distr,weight=c(1-WV_i, WV_i))}}
      lambda.prior.MY.MAP1[[k]] <- lambda.prior.k.MY.MAP1
      likelihood(lambda.prior.k.MY.MAP1) <- "poisson"
      lambda.post.k.MY.MAP1 <- postmix(lambda.prior.k.MY.MAP1,m=(data$rsum[k]/data$Esum[k]),n=data$Esum[k])
      lambda.post.MY.MAP1[[k]] <- lambda.post.k.MY.MAP1
      WV.myMAP1[k] <- WV_i
      mean.MyMAP.post1[k] <- summary(lambda.post.k.MY.MAP1)["mean"]
      likelihood(lambda.post.k.MY.MAP1)="exp"
      ess_group1[k]=ess(lambda.post.k.MY.MAP1,method = "elir") 
    }
    ESS_borrow1=sum(ess_group1)-sum(data$rsum)
    return(list(yangb.distr=yangb.distr,
                lambda.prior.MY.MAP1=lambda.prior.MY.MAP1,WV.myMAP1=WV.myMAP1,
                lambda.post.MY.MAP1=lambda.post.MY.MAP1,mean.MyMAP.post1=mean.MyMAP.post1,
                ESS_borrow1=ESS_borrow1))}
  result_all<- foreach(i=1:ncol(WV.matrix),.packages = c("RBesT"))%do%result_func(i)
  return(result_all)
}
find_median_survival_time <- function(lambdas, break_points) {
  PWE_survival_function <- function(time, lambdas, break_points) {
    surv <- 1
    for (i in 1:length(lambdas)) {
      t_min <- break_points[i]
      t_max <- break_points[i + 1]
      if (time <= t_max) {
        surv <- surv * exp(-lambdas[i] * (min(time, t_max) - t_min))
        break
      } else {
        surv <- surv * exp(-lambdas[i] * (t_max - t_min))}}
    return(surv)
  }
  optimize(function(time) abs(PWE_survival_function(time, lambdas, break_points) - 0.5), 
           interval = c(break_points[[1]], break_points[[length(break_points)]]), tol = 1e-10)$minimum
}
sample_of_medianTHR <- function(robust.MY.MAP.dist.K,cut.time){
  K=length(robust.MY.MAP.dist.K)
  smi_N <- 10000
  robust.MY.MAP.sample <- matrix(nrow = smi_N,ncol = K)
  for (k in 1:K) {
    robust.MY.MAP.sample[,k] <- rmix(robust.MY.MAP.dist.K[[k]],smi_N)
  }
  forfunction <- function(i) {
    median_survival_time <- find_median_survival_time(lambdas=robust.MY.MAP.sample[i,], break_points=cut.time)
    return(median_survival_time)
  }
  median_survival_time <- foreach(
    i=1:smi_N,   
    .combine=rbind,  
    .packages = c("RBesT"),
    .export="find_median_survival_time"
  ) %do% forfunction(i)
  median_survival_time=as.vector(median_survival_time)
  return(median_survival_time)
} 
#MSDB prior function-------------------------------------------------
#Convert the data into a format that meets the requirements
#data.X.Z,data.T.d are examples of data and the detail requirement are: 
#The input data needs to use the same variable name as in the examples
data.X.Z <- read.csv(file = "D:/MSDB_prior/data_X_example.csv",row.names = 1)
data.T.d <- read.csv(file = "D:/MSDB_prior/data_T_example.csv",row.names = 1)
View(data.X.Z);View(data.T.d)
#Parameter setting: number of time interval; number of propensity scores stratification; 
#Consistency parameter for SAM method, BOX method, EB-rMAP method, and MSDB method; Seed number
K=5;S=4;SAM_delta=0.05;P_cutpoint=0.85;P_cutpoint_my=0.85;seed=666
MSDB_prior <- function(data.X.Z,data.T.d,K,S,SAM_delta,P_cutpoint,P_cutpoint_my,seed){
  data.X.Z_C <- PS_score_label_func(data.X.Z[data.X.Z$group==0,],cutgroup_n=S)
  data.X.Z_T <- PS_score_label_func(data.X.Z[data.X.Z$group==1,],cutgroup_n=S)
  c(data.T.d_C, cut.time_C) %<-% K.label.func.3(data.T.d = data.T.d[data.T.d$group==0,],data.X.Z_C,K=K) 
  c(data.T.d_T, cut.time_T) %<-% K.label.func.3(data.T.d = data.T.d[data.T.d$group==1,],data.X.Z_T,K=K) 
  c(data.T.d_C,data.rE_C) %<-% rE.func.for.Z.PS.K2(data.T.d_C,data.X.Z_C,cut.time_C,K=K)
  c(data.T.d_T,data.rE_T) %<-% rE.func.for.Z.PS.K2(data.T.d_T,data.X.Z_T,cut.time_T,K=K)
  data.R.Z.K.S_C <- tauk.Z.K.S.func2(data.rE_C)
  data.R.Z.K.S_T <- tauk.Z.K.S.func2(data.rE_T)
  w.Z.K.S_C <- W.K.S.func2(data.rE_C)
  w.Z.K.S_T <- W.K.S.func2(data.rE_T)

  c(lambda.K.prior.MCMC_T,w.for.ps) %<-% lambda.prior.MCMC.func4(data.rE_T,
                                                                 data.R.Z.K.S_T,
                                                                 w.Z.K.S_T,w.for.ps=NA,seed=seed)
  c(lambda.K.prior.MCMC_C,w.for.ps) %<-% lambda.prior.MCMC.func4(data.rE_C,
                                                                 data.R.Z.K.S_C,
                                                                 w.Z.K.S_C,w.for.ps=w.for.ps,seed=seed)

  mixgamma.K_T %<-% lambda.MCMC.to.distr2(lambda.K.prior.MCMC_T)
  mixgamma.K_C %<-% lambda.MCMC.to.distr2(lambda.K.prior.MCMC_C)
  lambda.unimform_T <- lambda.unimform.func(data.T.d[data.T.d$group==1,],cut.time_T)
  lambda.unimform_C <- lambda.unimform.func(data.T.d[data.T.d$group==0,],cut.time_C)
  
  result_SAM_T <- Wvfunc.lambda.k.SAM(mixgamma.K_T,data.rE_T,K=K,SAM_delta=SAM_delta)
  result_SAM_C <- Wvfunc.lambda.k.SAM(mixgamma.K_C,data.rE_C,K=K,SAM_delta=SAM_delta)

  result_EBrmap_T<- Wvfunc.lambda.k.EBrmap(mixgamma.K_T,data.rE_T,K=K,P_cutpoint=P_cutpoint)
  result_EBrmap_C<- Wvfunc.lambda.k.EBrmap(mixgamma.K_C,data.rE_C,K=K,P_cutpoint=P_cutpoint)

  result_BOX_T<- Wvfunc.lambda.k.BOX(mixgamma.K_T,data.rE_T,K=K,P_cutpoint=P_cutpoint)
  result_BOX_C<- Wvfunc.lambda.k.BOX(mixgamma.K_C,data.rE_C,K=K,P_cutpoint=P_cutpoint)

  result_mypvalue_T<- Wvfunc.lambda.k.myPvalue2(mixgamma.K_T,data.rE_T,K=K,P_cutpoint=P_cutpoint_my)
  result_mypvalue_C<- Wvfunc.lambda.k.myPvalue2(mixgamma.K_C,data.rE_C,K=K,P_cutpoint=P_cutpoint_my)

  result_func <- function(result_gamma_T,result_gamma_C){
    c(yangb.distr_T,lambda.prior.MY.MAP1_T,WV.myMAP1_T,lambda.post.MY.MAP1_T,mean.MyMAP.post1_T,ESS_borrow1_T) %<-% result_gamma_T
    c(yangb.distr_C,lambda.prior.MY.MAP1_C,WV.myMAP1_C,lambda.post.MY.MAP1_C,mean.MyMAP.post1_C,ESS_borrow1_C) %<-% result_gamma_C
    
    ###result1 for pppcut1
    medianT_T.sample=sample_of_medianTHR(lambda.post.MY.MAP1_T,cut.time_T)
    medianT_C.sample=sample_of_medianTHR(lambda.post.MY.MAP1_C,cut.time_C)
    mius_medianT.sample=medianT_T.sample-medianT_C.sample
    if(quantile(mius_medianT.sample,0.05)<0){H.power=0}else{H.power=1}
    CI95_low <- quantile(mius_medianT.sample,0.025)
    CI95_up <- quantile(mius_medianT.sample,0.975)
    
    H_sim1 <- H.power;
    medianTD_postmean_sim1  <- mean(mius_medianT.sample)
    medianT.T1 <- mean(medianT_T.sample)
    medianT.C1 <- mean(medianT_C.sample)
    ESS_T1=ESS_borrow1_T
    ESS_C1=ESS_borrow1_C

    return(c(CI95_low,CI95_up,H_sim1,medianTD_postmean_sim1,medianT.T1,medianT.C1,ESS_T1,ESS_C1))
  }
  result_SAM_all <- foreach(result_gamma_T=result_SAM_T,
                            result_gamma_C=result_SAM_C,
                            .packages = c("RBesT","overlapping","dplyr" ,"survival" )
  )%do%result_func(result_gamma_T,result_gamma_C)
  result_EBrmap_all <- foreach(result_gamma_T=result_EBrmap_T,
                               result_gamma_C=result_EBrmap_C,
                               .packages = c("RBesT","overlapping","dplyr" ,"survival" )
  )%do%result_func(result_gamma_T,result_gamma_C)
  result_BOX_all <- foreach(result_gamma_T=result_BOX_T,
                            result_gamma_C=result_BOX_C,
                            .packages = c("RBesT","overlapping","dplyr" ,"survival" )
  )%do%result_func(result_gamma_T,result_gamma_C)
  result_myPvalue_all <- foreach(result_gamma_T=result_mypvalue_T,
                                 result_gamma_C=result_mypvalue_C,
                                 .packages = c("RBesT","overlapping","dplyr" ,"survival" )
  )%do%result_func(result_gamma_T,result_gamma_C)
  
  result_SAM_all <- unlist(result_SAM_all)
  result_EBrmap_all <- unlist(result_EBrmap_all)
  result_BOX_all <- unlist(result_BOX_all)
  result_myPvalue_all <- unlist(result_myPvalue_all)
  
  names(result_SAM_all) <- paste(c("CI95_low","CI95_up","H_sim1","medianTD_postmean_sim1","medianT.T1",
                                   "medianT.C1","ESS_T1","ESS_C1") %>% rep(length(SAM_delta)),"_","SAM","_", SAM_delta %>% rep(each=8),sep = "")
  names(result_EBrmap_all) <- paste(c("CI95_low","CI95_up","H_sim1","medianTD_postmean_sim1","medianT.T1",
                                      "medianT.C1","ESS_T1","ESS_C1") %>% rep(length(P_cutpoint)),"_","EBrmap","_", P_cutpoint %>% rep(each=8),sep = "")
  names(result_BOX_all) <- paste(c("CI95_low","CI95_up","H_sim1","medianTD_postmean_sim1","medianT.T1",
                                   "medianT.C1","ESS_T1","ESS_C1") %>% rep(length(P_cutpoint)),"_","BOX","_", P_cutpoint %>% rep(each=8),sep = "")
  names(result_myPvalue_all) <- paste(c("CI95_low","CI95_up","H_sim1","medianTD_postmean_sim1","medianT.T1",
                                        "medianT.C1","ESS_T1","ESS_C1") %>% rep(length(P_cutpoint_my)),"_","myPvalue","_", P_cutpoint_my %>% rep(each=8),sep = "")
  result_all <- c(result_SAM_all,result_EBrmap_all,result_BOX_all,result_myPvalue_all)
  return(result_all)
}

######################################################################
#------------------Part2: Application for MSDB prior-----------------#
######################################################################
#Function for application data generation based on literature
simu_ETD_data_f <- function(n){
  rnorm_f <- function(n,median_value,Q1,Q3,lim=1){
    std_dev <- (Q3 - Q1) / (2 * 0.6745)
    rn <- vector()
    while(length(rn)<n){
      rni <- rnorm(1,median_value,std_dev)
      if(lim[1]==1){rn <- c(rn,rni);next}
      if(rni>lim[2]|rni<lim[1]){next}
      rn <- c(rn,rni)
    }
    return(rn)
  }
  age_T_ETD <- rnorm_f(n,68,60,74,lim=c(18,100))  
  age_C_ETD <- rnorm_f(n,66,59,71,lim=c(18,100))
  sex_T_ETD <- rbinom(n,1,0.58)
  sex_C_ETD <- rbinom(n,1,0.46)
  
  COPD_T_ETD  <- rbinom(n,1,0.1)
  COPD_C_ETD  <- rbinom(n,1,0.11)
  eGFR_T_ETD <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(23,39,38))
  eGFR_C_ETD <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(33,34,33))
  
  ASCT_T_ETD <- rbinom(n,1,0.54)
  ASCT_C_ETD <- rbinom(n,1,0.59)
  T.initial_T_ETD  <- rnorm_f(n,4.46,2.6,7.2)
  T.initial_C_ETD  <- rnorm_f(n,4.09,2.9,7.0)
  Type.mye_T_ETD <- sample(c("IgA","IgG","Light.chain","other"), size = n, replace = TRUE, prob = c(0.22,0.66,0.1,0.02))
  Type.mye_C_ETD <- sample(c("IgA","IgG","Light.chain","other"), size = n, replace = TRUE, prob = c(0.27,0.65,0.07,0.01))
  ISS.stage1_T_ETD <- sample(c(1,2,3,0), size = n, replace = TRUE, prob = c(0.23,0.32,0.27,0.18))
  ISS.stage1_C_ETD <- sample(c(1,2,3,0), size = n, replace = TRUE, prob = c(0.27,0.31,0.29,0.13))
  ISS.stage2_T_ETD <- sample(c(1,2,3), size = n, replace = TRUE, prob = c(0.42,0.34,0.22))
  ISS.stage2_C_ETD <- sample(c(1,2,3), size = n, replace = TRUE, prob = c(0.33,0.37,0.28))
  Cyt.risk_T_ETD <- sample(c("high","standard","missing"), size = n, replace = TRUE, prob = c(0.16,0.67,0.18))
  Cyt.risk_C_ETD <- sample(c("high","standard","missing"), size = n, replace = TRUE, prob = c(0.24,0.51,0.26))
  pre.line_T_ETD  <- rnorm_f(n,3,2,4) %>% round()%>% abs
  pre.line_C_ETD  <- rnorm_f(n,3,2,4) %>% round()%>% abs
  Pre.al_T_ETD <- rbinom(n,1,0.90)
  Pre.al_C_ETD <- rbinom(n,1,0.97)
  Pre.pro_T_ETD <- rbinom(n,1,1)
  Pre.pro_C_ETD <- rbinom(n,1,1)
  Pre.le_T_ETD <- rbinom(n,1,1)
  Pre.le_C_ETD <- rbinom(n,1,1)
  
  #Patients refractory to previous therapy
  Last.line_T_ETD <- rbinom(n,1,0.97)
  Last.line_C_ETD <- rbinom(n,1,0.99)
  Imm.drug_T_ETD <- rbinom(n,1,0.96)
  Imm.drug_C_ETD <- rbinom(n,1,0.94)
  Len_T_ETD <- rbinom(n,1,0.94)
  Len_C_ETD <- rbinom(n,1,0.92)
  Pro.in_T_ETD <- rbinom(n,1,0.77)
  Pro.in_C_ETD <- rbinom(n,1,0.75)
  Lena.pro_T_ETD <- rbinom(n,1,0.72)
  Lena.pro_C_ETD <- rbinom(n,1,0.70)
  Lenal.last_T_ETD <- rbinom(n,1,0.60)
  Lenal.last_C_ETD <- rbinom(n,1,0.58)
  
  ETD_T_data <- {data.frame(age=age_T_ETD,
                            sex=sex_T_ETD,
                            COPD=COPD_T_ETD,
                            eGFR=eGFR_T_ETD,
                            ASCT=ASCT_T_ETD,
                            T.initial=T.initial_T_ETD,
                            Type.mye=Type.mye_T_ETD,
                            ISS.stage1=ISS.stage1_T_ETD,
                            ISS.stage2=ISS.stage2_T_ETD,
                            Cyt.risk=Cyt.risk_T_ETD,
                            pre.line=pre.line_T_ETD,
                            Pre.al=Pre.al_T_ETD,
                            Pre.pro=Pre.pro_T_ETD,
                            Pre.le=Pre.le_T_ETD,
                            #Patients refractory to previous therapy
                            Last.line=Last.line_T_ETD,
                            Imm.drug=Imm.drug_T_ETD,
                            Len=Len_T_ETD,
                            Pro.in=Pro.in_T_ETD,
                            Lena.pro=Lena.pro_T_ETD,
                            Lenal.last=Lenal.last_T_ETD,
                            group=1)}
  ETD_C_data <-{ data.frame(age=age_C_ETD,
                            sex=sex_C_ETD,
                            COPD=COPD_C_ETD,
                            eGFR=eGFR_C_ETD,
                            ASCT=ASCT_C_ETD,
                            T.initial=T.initial_C_ETD,
                            Type.mye=Type.mye_C_ETD,
                            ISS.stage1=ISS.stage1_C_ETD,
                            ISS.stage2=ISS.stage2_C_ETD,
                            Cyt.risk=Cyt.risk_C_ETD,
                            pre.line=pre.line_C_ETD,
                            Pre.al=Pre.al_C_ETD,
                            Pre.pro=Pre.pro_C_ETD,
                            Pre.le=Pre.le_C_ETD,
                            #Patients refractory to previous therapy
                            Last.line=Last.line_C_ETD,
                            Imm.drug=Imm.drug_C_ETD,
                            Len=Len_C_ETD,
                            Pro.in=Pro.in_C_ETD,
                            Lena.pro=Lena.pro_C_ETD,
                            Lenal.last=Lenal.last_C_ETD,
                            group=0)}
  
  data_ETD <- rbind(ETD_T_data,ETD_C_data)
  #Previous autoloqous stem-cell transplantation=Pre.aut.ste.cel.t
  #Years since initial diagnosis=T.initial
  #Type of myeloma at diagnosis=Type.mye
  #ISS stage at diagnosis=ISS.stage1
  #ISS stage at study entry=ISS.stage2
  #Cytogenetic risk at study entry=Cyt.risk
  #Previous lines of therapy=pre.line
  #Previous alkylating drug=Pre.al
  #Previous proteasome inhibitors=Pre.pro
  #Previous lenalidomide=Pre.le
  #Patients refractory to previous therapy
  ##Last line of therapy=Last.line
  ##Immunomodulatory drug=Imm.drug
  ##Lenalidomide=Len
  ##Protease inhibitor=Pro.in
  ##Lenalidomide and protease inhibitor=Lena.pro
  ##Lenalidomide last line=Lenal.last
  return(data_ETD)
}
simu_CTD_data_f <- function(n){
  rnorm_f <- function(n,median_value,Q1,Q3,lim=1){
    std_dev <- (Q3 - Q1) / (2 * 0.6745)
    rn <- vector()
    while(length(rn)<n){
      rni <- rnorm(1,median_value,std_dev)
      if(lim[1]==1){rn <- c(rn,rni);next}
      if(rni>lim[2]|rni<lim[1]){next}
      rn <- c(rn,rni)
    }
    return(rn)
  }
  datagroup <- function(n){
    age_C_ETD <<- rnorm_f(n,61,44,77,lim=c(18,100))
    sex_C_ETD <<- rbinom(n,1,0.595)
    # COPD_C_ETD  <<- rbinom(n,1,0.11)
    # eGFR_C_ETD <<- sample(c(0,1,2), size = n, replace = TRUE, prob = c(33,34,33))
    ASCT_C_ETD <<- rbinom(n,1,0.108)
    T.initial_C_ETD  <<- rnorm_f(n,3.4,2.2,6.2)
    #Type.mye_C_ETD <<- sample(c("IgA","IgG","Light.chain","other"), size = n, replace = TRUE, prob = c(0.27,0.65,0.07,0.01))
    #ISS.stage1_C_ETD <<- sample(c(1,2,3,0), size = n, replace = TRUE, prob = c(0.27,0.31,0.29,0.13))
    ISS.stage2_C_ETD <<- sample(c(1,2,3), size = n, replace = TRUE, prob = c(0.257,0.567,0.176))
    Cyt.risk_C_ETD <<- sample(c("high","standard","missing"), size = n, replace = TRUE, prob = c(0.662,0.338,0))
    pre.line_C_ETD  <<- rnorm_f(n,3,1,9) %>% round() %>% abs
    # Pre.al_C_ETD <<- rbinom(n,1,0.97)
    # Pre.pro_C_ETD <<- rbinom(n,1,1)
    # Pre.le_C_ETD <<- rbinom(n,1,1)
    
    #Patients refractory to previous therapy
    # Last.line_C_ETD <<- rbinom(n,1,0.99)
    # Imm.drug_C_ETD <<- rbinom(n,1,0.94)
    Len_C_ETD <<- rbinom(n,1,0.432)
    # Pro.in_C_ETD <<- rbinom(n,1,0.75)
    # Lena.pro_C_ETD <<- rbinom(n,1,0.70)
    # Lenal.last_C_ETD <<- rbinom(n,1,0.58)
  }
  datagroup(n)
  ETD_C_data <-{ data.frame(age=age_C_ETD,
                            sex=sex_C_ETD,
                            # COPD=COPD_C_ETD,
                            # eGFR=eGFR_C_ETD,
                            ASCT=ASCT_C_ETD,
                            T.initial=T.initial_C_ETD,
                            # Type.mye=Type.mye_C_ETD,
                            # ISS.stage1=ISS.stage1_C_ETD,
                            ISS.stage2=ISS.stage2_C_ETD,
                            Cyt.risk=Cyt.risk_C_ETD,
                            pre.line=pre.line_C_ETD,
                            # Pre.al=Pre.al_C_ETD,
                            # Pre.pro=Pre.pro_C_ETD,
                            # Pre.le=Pre.le_C_ETD,
                            # #Patients refractory to previous therapy
                            # Last.line=Last.line_C_ETD,
                            # Imm.drug=Imm.drug_C_ETD,
                            Len=Len_C_ETD,
                            # Pro.in=Pro.in_C_ETD,
                            # Lena.pro=Lena.pro_C_ETD,
                            # Lenal.last=Lenal.last_C_ETD,
                            group=0)}
  datagroup(n)
  ETD_T_data <- { data.frame(age=age_C_ETD,
                             sex=sex_C_ETD,
                             # COPD=COPD_C_ETD,
                             # eGFR=eGFR_C_ETD,
                             ASCT=ASCT_C_ETD,
                             T.initial=T.initial_C_ETD,
                             # Type.mye=Type.mye_C_ETD,
                             # ISS.stage1=ISS.stage1_C_ETD,
                             ISS.stage2=ISS.stage2_C_ETD,
                             Cyt.risk=Cyt.risk_C_ETD,
                             pre.line=pre.line_C_ETD,
                             # Pre.al=Pre.al_C_ETD,
                             # Pre.pro=Pre.pro_C_ETD,
                             # Pre.le=Pre.le_C_ETD,
                             # #Patients refractory to previous therapy
                             # Last.line=Last.line_C_ETD,
                             # Imm.drug=Imm.drug_C_ETD,
                             Len=Len_C_ETD,
                             # Pro.in=Pro.in_C_ETD,
                             # Lena.pro=Lena.pro_C_ETD,
                             # Lenal.last=Lenal.last_C_ETD,
                             group=1)}
  
  data_ETD <- rbind(ETD_T_data,ETD_C_data)
  #Previous autoloqous stem-cell transplantation=ASCT
  #Years since initial diagnosis=T.initial
  #Type of myeloma at diagnosis=Type.mye
  #ISS stage at diagnosis=ISS.stage1
  #ISS stage at study entry=ISS.stage2
  #Cytogenetic risk at study entry=Cyt.risk
  #Previous lines of therapy=pre.line
  #Previous alkylating drug=Pre.al
  #Previous proteasome inhibitors=Pre.pro
  #Previous lenalidomide=Pre.le
  #Patients refractory to previous therapy
  ##Last line of therapy=Last.line
  ##Immunomodulatory drug=Imm.drug
  ##Lenalidomide=Len
  ##Protease inhibitor=Pro.in
  ##Lenalidomide and protease inhibitor=Lena.pro
  ##Lenalidomide last line=Lenal.last
  return(data_ETD)
}
simu_RWD_data_f <- function(n){
  rnorm_f <- function(n,median_value,Q1,Q3,lim=1){
    std_dev <- (Q3 - Q1) / (2 * 0.6745)
    rn <- vector()
    while(length(rn)<n){
      rni <- rnorm(1,median_value,std_dev)
      if(lim[1]==1){rn <- c(rn,rni);next}
      if(rni>lim[2]|rni<lim[1]){next}
      rn <- c(rn,rni)
    }
    return(rn)
  }
  age_T_ETD <- rnorm_f(n,63.5,49,79,lim=c(18,100))
  age_C_ETD <- rnorm_f(n,63.5,40,90,lim=c(18,100))
  sex_T_ETD <- rbinom(n,1,0.563)
  sex_C_ETD <- rbinom(n,1,0.632)
  
  # COPD_T_ETD  <- rbinom(n,1,0.1)
  # COPD_C_ETD  <- rbinom(n,1,0.11)
  # eGFR_T_ETD <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(23,39,38))
  # eGFR_C_ETD <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(33,34,33))
  
  ASCT_T_ETD <- rbinom(n,1,0.188)
  ASCT_C_ETD <- rbinom(n,1,0.162)
  # T.initial_T_ETD  <- rnorm_f(n,4.46,2.6,7.2)
  # T.initial_C_ETD  <- rnorm_f(n,4.09,2.9,7.0)
  T.initial_T_ETD  <- rbinom(n,1,0.521)
  T.initial_C_ETD  <- rbinom(n,1,0.353)
  
  Type.mye_T_ETD <- sample(c("IgA","IgG","Light.chain","other"), size = n, replace = TRUE, prob = c(0.25,0.417,0.208,0.125))
  Type.mye_C_ETD <- sample(c("IgA","IgG","Light.chain","other"), size = n, replace = TRUE, prob = c(0.162,0.515,0.191,0.132))
  # ISS.stage1_T_ETD <- sample(c(1,2,3,0), size = n, replace = TRUE, prob = c(0.23,0.32,0.27,0.18))
  # ISS.stage1_C_ETD <- sample(c(1,2,3,0), size = n, replace = TRUE, prob = c(0.27,0.31,0.29,0.13))
  ISS.stage2_T_ETD <- sample(c(1,2,3), size = n, replace = TRUE, prob = c(0.4,445,0.155))
  ISS.stage2_C_ETD <- sample(c(1,2,3), size = n, replace = TRUE, prob = c(0.414,0.414,0.172))
  Cyt.risk_T_ETD <- sample(c("high","standard","missing"), size = n, replace = TRUE, prob = c(0.5,0.4,0.1))
  Cyt.risk_C_ETD <- sample(c("high","standard","missing"), size = n, replace = TRUE, prob = c(0.42,0.38,0.2))
  pre.line_T_ETD  <- rnorm_f(n,2,1,5) %>% round()%>% abs
  pre.line_C_ETD  <- rnorm_f(n,2,1,4) %>% round()%>% abs
  # Pre.al_T_ETD <- rbinom(n,1,0.90)
  # Pre.al_C_ETD <- rbinom(n,1,0.97)
  # Pre.pro_T_ETD <- rbinom(n,1,1)
  # Pre.pro_C_ETD <- rbinom(n,1,1)
  Pre.le_T_ETD <- rbinom(n,1,0.896)
  Pre.le_C_ETD <- rbinom(n,1,0.838)
  
  #Patients refractory to previous therapy
  # Last.line_T_ETD <- rbinom(n,1,0.97)
  # Last.line_C_ETD <- rbinom(n,1,0.99)
  # Imm.drug_T_ETD <- rbinom(n,1,0.96)
  # Imm.drug_C_ETD <- rbinom(n,1,0.94)
  Len_T_ETD <- rbinom(n,1,0.771)
  Len_C_ETD <- rbinom(n,1,0.662)
  Pro.in_T_ETD <- rbinom(n,1,0.625)
  Pro.in_C_ETD <- rbinom(n,1,0.676)
  Lena.pro_T_ETD <- rbinom(n,1,0.521)
  Lena.pro_C_ETD <- rbinom(n,1,0.441)
  # Lenal.last_T_ETD <- rbinom(n,1,0.60)
  # Lenal.last_C_ETD <- rbinom(n,1,0.58)
  
  ETD_T_data <- {data.frame(age=age_T_ETD,
                            sex=sex_T_ETD,
                            # COPD=COPD_T_ETD,
                            # eGFR=eGFR_T_ETD,
                            ASCT=ASCT_T_ETD,
                            # T.initial=T.initial_T_ETD,
                            Type.mye=Type.mye_T_ETD,
                            # ISS.stage1=ISS.stage1_T_ETD,
                            ISS.stage2=ISS.stage2_T_ETD,
                            Cyt.risk=Cyt.risk_T_ETD,
                            pre.line=pre.line_T_ETD,
                            # Pre.al=Pre.al_T_ETD,
                            # Pre.pro=Pre.pro_T_ETD,
                            Pre.le=Pre.le_T_ETD,
                            #Patients refractory to previous therapy
                            # Last.line=Last.line_T_ETD,
                            # Imm.drug=Imm.drug_T_ETD,
                            Len=Len_T_ETD,
                            Pro.in=Pro.in_T_ETD,
                            Lena.pro=Lena.pro_T_ETD,
                            # Lenal.last=Lenal.last_T_ETD,
                            group=1)}
  ETD_C_data <-{ data.frame(age=age_C_ETD,
                            sex=sex_C_ETD,
                            # COPD=COPD_C_ETD,
                            # eGFR=eGFR_C_ETD,
                            ASCT=ASCT_C_ETD,
                            # T.initial=T.initial_C_ETD,
                            Type.mye=Type.mye_C_ETD,
                            # ISS.stage1=ISS.stage1_C_ETD,
                            ISS.stage2=ISS.stage2_C_ETD,
                            Cyt.risk=Cyt.risk_C_ETD,
                            pre.line=pre.line_C_ETD,
                            # Pre.al=Pre.al_C_ETD,
                            # Pre.pro=Pre.pro_C_ETD,
                            Pre.le=Pre.le_C_ETD,
                            #Patients refractory to previous therapy
                            # Last.line=Last.line_C_ETD,
                            # Imm.drug=Imm.drug_C_ETD,
                            Len=Len_C_ETD,
                            Pro.in=Pro.in_C_ETD,
                            Lena.pro=Lena.pro_C_ETD,
                            # Lenal.last=Lenal.last_C_ETD,
                            group=0)}
  
  data_ETD <- rbind(ETD_T_data,ETD_C_data)
  #Previous autoloqous stem-cell transplantation=Pre.aut.ste.cel.t
  #Years since initial diagnosis=T.initial
  #Type of myeloma at diagnosis=Type.mye
  #ISS stage at diagnosis=ISS.stage1
  #ISS stage at study entry=ISS.stage2
  #Cytogenetic risk at study entry=Cyt.risk
  #Previous lines of therapy=pre.line
  #Previous alkylating drug=Pre.al
  #Previous proteasome inhibitors=Pre.pro
  #Previous lenalidomide=Pre.le
  #Patients refractory to previous therapy
  ##Last line of therapy=Last.line
  ##Immunomodulatory drug=Imm.drug
  ##Lenalidomide=Len
  ##Protease inhibitor=Pro.in
  ##Lenalidomide and protease inhibitor=Lena.pro
  ##Lenalidomide last line=Lenal.last
  return(data_ETD)
}
#Parameter setting
K=2;S=3 
P_cutpoint_my=c(0.8,0.85,0.9)
N_CTD=80;N_ETD=150;N_RWD=100
canshu <- c(0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44)
seed_simu=1:500
simu_time=1:500
#smulation function
forfunction_result_ALL <- function(i){
  seed=seed_simu[i]
  set.seed(seed)
  {
    data_ETD_X <- simu_ETD_data_f(N_ETD)
    data_CTD_X <- simu_CTD_data_f(N_CTD)
    data_RWD_X <- simu_RWD_data_f(N_RWD)
    inter_name <- intersect(intersect(colnames(data_ETD_X),colnames(data_CTD_X)),colnames(data_RWD_X))
    
    data_ETD_X <- data_ETD_X %>% select(inter_name)
    data_CTD_X <- data_CTD_X %>% select(inter_name)
    data_RWD_X <- data_RWD_X %>% select(inter_name)
    
    data.simu.T_example <- function(data.X.Z,Ture_effect,alpha_ETD,
                                    beta,shape=1.5,T.shuif=5,seed=666){
      set.seed(seed)
      library(survival)
      library(dplyr)
      library(dummies)#build Dummy Variables
      X <- data.X.Z %>% select(-group) 
      X <- dummy.data.frame(X, names = c("Type.mye","Cyt.risk")) %>% as.matrix()
      
      group <- data.X.Z$group 
      beta1 <- rep(beta,ncol(X))
      
      t_median=X%*%beta1+alpha_ETD+Ture_effect*group
      scale=t_median/(log(2)^(1/shape))
      Time<- rweibull(nrow(X), shape = shape, scale = scale)
      delta=1
      time_data <- data.frame(S.time=Time,delta=delta)
      # time_data <- data.frame(Time=Time,delta=delta)
      
      time_data <- data.frame(time_data,data.X.Z)
      time_data$C.time <- T.shuif
      
      time_data$Time <- pmin(time_data$S.time,T.shuif)
      time_data$delta <- ifelse(time_data$S.time<time_data$C.time,1,0)
      data.T.d <- time_data[,c("Time","delta")]
      data.T.d <- cbind(data.T.d,data.X.Z)
      
      return(data.T.d)
    }
    data_ETD_T <- data.simu.T_example(data_ETD_X,Ture_effect=0.40,alpha_ETD=-0.55,beta=0.015)
    data_RWD_T <- data.simu.T_example(data_RWD_X,Ture_effect=ture_RWD,alpha_ETD=-0.47,beta=0.015)
    data_CTD_T <- data.simu.T_example(data_CTD_X,Ture_effect=ture_CTD,alpha_ETD=-0.47,beta=0.015)
    
    data.X.Z <- rbind(data.frame(data_CTD_X,Z1=0,Z2=0,Z="CTD"),
                      data.frame(data_RWD_X,Z1=0,Z2=1,Z="RWD"),
                      data.frame(data_ETD_X,Z1=1,Z2=0,Z="ETD")) %>% cbind(data.frame(ID=1:(nrow(data_CTD_X)+nrow(data_ETD_X)+nrow(data_RWD_X))))
    data.T.d <- rbind(data.frame(data_CTD_T,Z1=0,Z2=0,Z="CTD"),
                      data.frame(data_RWD_T,Z1=0,Z2=1,Z="RWD"),
                      data.frame(data_ETD_T,Z1=1,Z2=0,Z="ETD")) %>% cbind(data.frame(ID=1:(nrow(data_CTD_X)+nrow(data_ETD_X)+nrow(data_RWD_X))))
    data.T.d <- data.T.d[which(is.nan(data.T.d$Time)==F),]}
  data.X.Z_C <- PS_score_label_func_my(data.X.Z[data.X.Z$group==0,],cutgroup_n=S)
  data.X.Z_T <- PS_score_label_func_my(data.X.Z[data.X.Z$group==1,],cutgroup_n=S)
  c(data.T.d_C, cut.time_C) %<-% K.label.func.3_my(data.T.d = data.T.d[data.T.d$group==0,],data.X.Z_C,K=K) 
  c(data.T.d_T, cut.time_T) %<-% K.label.func.3_my(data.T.d = data.T.d[data.T.d$group==1,],data.X.Z_T,K=K) 
  c(data.T.d_C,data.rE_C) %<-% rE.func.for.Z.PS.K2_my(data.T.d_C,data.X.Z_C,cut.time_C,K=K)
  c(data.T.d_T,data.rE_T) %<-% rE.func.for.Z.PS.K2_my(data.T.d_T,data.X.Z_T,cut.time_T,K=K)
  data.R.Z.K.S_C <- tauk.Z.K.S.func2_my(data.rE_C)
  data.R.Z.K.S_T <- tauk.Z.K.S.func2_my(data.rE_T)
  w.Z.K.S_C <- W.K.S.func2_my(data.rE_C)
  w.Z.K.S_T <- W.K.S.func2_my(data.rE_T)
  
  c(lambda.K.prior.MCMC_T,w.for.ps) %<-% lambda.prior.MCMC.func4_my(data.rE_T,
                                                                    data.R.Z.K.S_T,
                                                                    w.Z.K.S_T,w.for.ps=NA,seed=seed)
  c(lambda.K.prior.MCMC_C,w.for.ps) %<-% lambda.prior.MCMC.func4_my(data.rE_C,
                                                                    data.R.Z.K.S_C,
                                                                    w.Z.K.S_C,w.for.ps=w.for.ps,seed=seed)
  
  mixgamma.K_T %<-% lambda.MCMC.to.distr2_my(lambda.K.prior.MCMC_T)
  mixgamma.K_C %<-% lambda.MCMC.to.distr2_my(lambda.K.prior.MCMC_C)
  lambda.unimform_T <- lambda.unimform.func_my(data.T.d[data.T.d$group==1,],cut.time_T)
  lambda.unimform_C <- lambda.unimform.func_my(data.T.d[data.T.d$group==0,],cut.time_C)
  
  result_mypvalue_T<- Wvfunc.lambda.k.myPvalue2_my(mixgamma.K_T,data.rE_T,K=K,P_cutpoint=P_cutpoint_my)
  result_mypvalue_C<- Wvfunc.lambda.k.myPvalue2_my(mixgamma.K_C,data.rE_C,K=K,P_cutpoint=P_cutpoint_my)
  
  result_func <- function(result_gamma_T,result_gamma_C){
    c(yangb.distr_T,lambda.prior.MY.MAP1_T,WV.myMAP1_T,lambda.post.MY.MAP1_T,mean.MyMAP.post1_T,ESS_borrow1_T) %<-% result_gamma_T
    c(yangb.distr_C,lambda.prior.MY.MAP1_C,WV.myMAP1_C,lambda.post.MY.MAP1_C,mean.MyMAP.post1_C,ESS_borrow1_C) %<-% result_gamma_C
    
    ###result1 for pppcut1
    medianT_T.sample=sample_of_medianTHR(lambda.post.MY.MAP1_T,cut.time_T)
    medianT_C.sample=sample_of_medianTHR(lambda.post.MY.MAP1_C,cut.time_C)
    mius_medianT.sample=medianT_T.sample-medianT_C.sample
    if(quantile(mius_medianT.sample,0.21)<0){H.power=0}else{H.power=1}
    
    H_sim1 <- H.power
    medianTD_postmean_sim1  <- mean(mius_medianT.sample)
    medianT.T1 <- mean(medianT_T.sample)
    medianT.C1 <- mean(medianT_C.sample)
    ESS_T1=ESS_borrow1_T
    ESS_C1=ESS_borrow1_C
    gc()
    return(c(H_sim1,medianTD_postmean_sim1,medianT.T1,medianT.C1,ESS_T1,ESS_C1))
  }
  
  result_myPvalue_all <- foreach(result_gamma_T=result_mypvalue_T,
                                 result_gamma_C=result_mypvalue_C,
                                 .packages = c("RBesT","overlapping","dplyr" ,"survival" )
  )%do%result_func(result_gamma_T,result_gamma_C)
  
  result_myPvalue_all <- unlist(result_myPvalue_all)
  
  names(result_myPvalue_all) <- paste(c("H_sim1","medianTD_postmean_sim1","medianT.T1",
                                        "medianT.C1","ESS_T1","ESS_C1") %>% rep(length(P_cutpoint_my)),"_","myPvalue","_", P_cutpoint_my %>% rep(each=6),sep = "")
  return(result_myPvalue_all)
}

for(ii in 1:length(canshu)){
  ture_RWD=ture_CTD=canshu[ii]
  
  library(doParallel)
  library(foreach)
  cl<- makeCluster(10)
  registerDoParallel(cl)
  error_messages <- NULL
  result_simu1 <- foreach(
    i=simu_time,
    .combine=cbind,  
    .packages = c("SAMprior","base","tidyr","MASS","VGAM","RBesT","dplyr","keras",
                  "overlapping","nnet" ,"knitr","survival","distr","tidyr",
                  "zeallot","EBrmap","doParallel","foreach","BH","cmdstanr")
  ) %dopar% tryCatch({     forfunction_result_ALL(i)   },error = function(e) {error_messages<<-rbind(error_messages,c(e$message,i));print(e$message)})   
  stopCluster(cl)
  
  result_name <- paste("D:/MSDB_prior/","example","_MYMAP_",ture_CTD,".CSV",sep = "")
  write.csv(result_simu1,file = result_name)
}

######################################################################
#------------------Part3: Simulation for MSDB prior------------------#
######################################################################
#simulation 
data.simu.X.Z.HR2 <- function(n.for.model_now,n.for.model_histr,beta_RWD=0.1,beta_ETD=0.2,seed){#鎹㈡�????????濮嬬殑鏁版嵁妯℃嫙鎬濊矾??????杩囧钩绉婚噺鏉ヤ骇鐢熷熀绾夸笉鍧囪�??鐨勬�?????
  set.seed(seed)
  mean.X.norm <- 1
  sigma.X.norm <- 0.2
  prob=0.3
  X.CTD.norm <- rnorm(n.for.model_now*2, mean.X.norm, sigma.X.norm)
  X.RWD.norm <- rnorm(n.for.model_histr*2, mean.X.norm+beta_RWD, sigma.X.norm)
  X.ETD.norm <- rnorm(n.for.model_histr*2, mean.X.norm+beta_ETD, sigma.X.norm)
  X.CTD.bin <- rbinom(n.for.model_now*2,1,prob)
  X.RWD.bin <- rbinom(n.for.model_histr*2,1,prob+beta_RWD)
  X.ETD.bin <- rbinom(n.for.model_histr*2,1,prob+beta_ETD)
  group1 <- rep(c(0,1),each=n.for.model_now)
  group2 <- rep(c(0,1),each=n.for.model_histr)
  data.CTD <- data.frame(X.norm=X.CTD.norm,X.bin=X.CTD.bin,Z1=0,Z2=0,Z="CTD",group=group1)
  data.RWD <- data.frame(X.norm=X.RWD.norm,X.bin=X.RWD.bin,Z1=0,Z2=1,Z="RWD",group=group2)
  data.ETD <- data.frame(X.norm=X.ETD.norm,X.bin=X.ETD.bin,Z1=1,Z2=0,Z="ETD",group=group2)
  data.all <- rbind(data.CTD,data.ETD,data.RWD)
  data.all$Z <- factor(data.all$Z)
  data.all$ID <- 1:nrow(data.all)
  return(data.all)
}
data.simu.T.delta.NHR3 <- function(data.X.Z,alpha_CTD,alpha_ETD,alpha_RWD,
                                   beta,censor.parm.list,c.rate,shape=2,T.shuif=5,seed){
  set.seed(seed)
  library(survival)
  X <- as.matrix(data.X.Z[,grep("X", colnames(data.X.Z))])
  Z12 <- as.matrix(data.X.Z[,grep("Z", colnames(data.X.Z))][,1:2]) 
  group <- data.X.Z$group 
  beta1 <- rep(beta,ncol(X))
  
  t_median=1+X%*%beta1+alpha_ETD*Z12[,1]*group+alpha_RWD*Z12[,2]*group+alpha_CTD*group
  scale=t_median/(log(2)^(1/shape))
  Time<- rweibull(nrow(X), shape = shape, scale = scale)
  delta=1
  time_data <- data.frame(S.time=Time,delta=delta)

  time_data <- data.frame(time_data,data.X.Z)
  c(c.para.CTD.C,c.para.ETD.C,c.para.RWD.C,c.para.CTD.T,c.para.ETD.T,c.para.RWD.T) %<-% censor.parm.list[censor.parm.list$c.rate==c.rate,1:6]
  n.for.model_CTD <- length(intersect(which(time_data$group==0),which(time_data$Z=="CTD")))
  n.for.model_histr <- length(intersect(which(time_data$group==0),which(time_data$Z=="ETD")))
  time_data$C.time <- T.shuif
  time_data$C.time[intersect(which(time_data$group==0),which(time_data$Z=="CTD"))
  ] <- runif(n.for.model_CTD,0.01,c.para.CTD.C)
  time_data$C.time[intersect(which(time_data$group==0),which(time_data$Z=="ETD"))
  ] <- runif(n.for.model_histr,0.01,c.para.ETD.C)
  time_data$C.time[intersect(which(time_data$group==0),which(time_data$Z=="RWD"))
  ] <- runif(n.for.model_histr,0.01,c.para.RWD.C)
  time_data$C.time[intersect(which(time_data$group==1),which(time_data$Z=="CTD"))
  ] <- runif(n.for.model_CTD,0.01,c.para.CTD.T)
  time_data$C.time[intersect(which(time_data$group==1),which(time_data$Z=="ETD"))
  ] <- runif(n.for.model_histr,0.01,c.para.ETD.T)
  time_data$C.time[intersect(which(time_data$group==1),which(time_data$Z=="RWD"))
  ] <- runif(n.for.model_histr,0.01,c.para.RWD.T)
  
  time_data$Time <- pmin(time_data$S.time,time_data$C.time,T.shuif)
  time_data$delta <- ifelse(time_data$S.time<apply(data.frame(C.Time=time_data$C.time,T.shuif=T.shuif),1,min),1,0)
  data.T.d <- time_data[,c("Time","delta")]
  data.T.d <- cbind(data.T.d,data.X.Z)
  
  return(data.T.d)
}
censor.parm.for.censor.rate.weibull.NHR3 <- function(beta_RWD,beta_ETD,alpha_CTD,alpha_ETD,alpha_RWD,beta,
                                                     c.rate.list=c(0.2,0.3,0.5),shape=2,T.shuif=5,seed=666){
  set.seed(seed)
  data.X.Z <- data.simu.X.Z.HR2(n.for.model_now = 10000,n.for.model_histr = 10000,beta_RWD=beta_RWD,beta_ETD=beta_ETD,seed)
  func.t.d <- function(data.X.Z,n.for.model=10000,alpha_CTD,alpha_ETD,alpha_RWD,beta,censor.rate,shape,scale_base){
    library(survival)
    X <- as.matrix(data.X.Z[,grep("X", colnames(data.X.Z))])
    Z12 <- as.matrix(data.X.Z[,grep("Z", colnames(data.X.Z))][,1:2]) 
    group <- data.X.Z$group 
    beta1 <- rep(beta,ncol(X))
    
    t_median=1+X%*%beta1+alpha_ETD*Z12[,1]*group+alpha_RWD*Z12[,2]*group+alpha_CTD*group
    scale=t_median/(log(2)^(1/shape))
    Time<- rweibull(nrow(X), shape = shape, scale = scale)
    delta <- rep(1, length(Time))
    time_data <- data.frame(S.time=Time,delta=delta)
    
    time_data$C.time <- runif(n.for.model,0.01,censor.rate)
    time_data$delta <- ifelse(time_data$S.time<apply(data.frame(C.Time=time_data$C.time,T.shuif=T.shuif),1,min),1,0)
    
    time_data <- data.frame(time_data,data.X.Z)
    censor.rate.CTD_C <- 1-sum(time_data$delta[intersect(which(time_data$Z=="CTD"),which(time_data$group==0))])/n.for.model
    censor.rate.ETD_C <- 1-sum(time_data$delta[intersect(which(time_data$Z=="ETD"),which(time_data$group==0))])/n.for.model
    censor.rate.RWD_C <- 1-sum(time_data$delta[intersect(which(time_data$Z=="RWD"),which(time_data$group==0))])/n.for.model
    censor.rate.CTD_T <- 1-sum(time_data$delta[intersect(which(time_data$Z=="CTD"),which(time_data$group==1))])/n.for.model
    censor.rate.ETD_T <- 1-sum(time_data$delta[intersect(which(time_data$Z=="ETD"),which(time_data$group==1))])/n.for.model
    censor.rate.RWD_T <- 1-sum(time_data$delta[intersect(which(time_data$Z=="RWD"),which(time_data$group==1))])/n.for.model
    
    return(censor.rate=list(censor.rate.CTD_C=censor.rate.CTD_C,
                            censor.rate.ETD_C=censor.rate.ETD_C,
                            censor.rate.RWD_C=censor.rate.RWD_C,
                            censor.rate.CTD_T=censor.rate.CTD_T,
                            censor.rate.ETD_T=censor.rate.ETD_T,
                            censor.rate.RWD_T=censor.rate.RWD_T))
  }
  parm <- seq(from=0.01,to=30,by=0.01)
  forfunction <- function(i) {
    censor.rate <- func.t.d(data.X.Z,n.for.model=10000,alpha_CTD,alpha_ETD,alpha_RWD,beta,censor.rate=parm[i],shape,scale_base)
    censor.parm.list <-c(censor.rate$censor.rate.CTD_C,
                         censor.rate$censor.rate.ETD_C,
                         censor.rate$censor.rate.RWD_C,
                         censor.rate$censor.rate.CTD_T,
                         censor.rate$censor.rate.ETD_T,
                         censor.rate$censor.rate.RWD_T)
    names(censor.parm.list) <- c("censor.rate.CTD_C",
                                 "censor.rate.ETD_C",
                                 "censor.rate.RWD_C",
                                 "censor.rate.CTD_T",
                                 "censor.rate.ETD_T",
                                 "censor.rate.RWD_T")
    return(censor.parm.list)
  }
  library(doParallel)
  library(foreach)
  cl<- makeCluster(34)
  registerDoParallel(cl)       
  censor.parm.list <- foreach(
    i=1:length(parm),   
    .combine=rbind,  
    .packages = c("base","tidyr") 
    
  ) %dopar% forfunction(i)
  stopCluster(cl)
  censor.parm.list <- data.frame(censor.parm.list,parm=parm)

  censor.parm.list_last <- data.frame()
  for (i in 1:length(c.rate.list)) {
    parm_calcu<- function(col,parm,c.rate){
      max(parm[col>=c.rate])
    }
    censor.parm.list_last <- 
      rbind(censor.parm.list_last,
            apply(censor.parm.list[,1:6], 2, parm_calcu,parm=censor.parm.list$parm,c.rate=c.rate.list[i]))
  }
  censor.parm.list_last <- cbind(censor.parm.list_last,c.rate.list)
  colnames(censor.parm.list_last) <- c(colnames(censor.parm.list)[1:6],"c.rate")
  return(censor.parm.list_last)
}
forfunction_MCMC <- function(i){
  seed=seed_simu[i]
  data.X.Z <- data.simu.X.Z.HR2(n.for.model_now = n.for.model_now,n.for.model_histr = n.for.model_histr,beta_RWD = beta_RWD,beta_ETD = beta_ETD,seed)
  data.T.d=data.simu.T.delta.NHR3(data.X.Z,alpha_CTD,alpha_ETD,alpha_RWD,beta,censor.parm.list,c.rate,shape=2,T.shuif=5,seed)
  data.T.d0=data.simu.T.delta.NHR3(data.X.Z,alpha_CTD = 0,alpha_ETD,alpha_RWD,beta,censor.parm.list,c.rate,shape=2,T.shuif=5,seed)
  data.X.Z_C <- PS_score_label_func(data.X.Z[data.X.Z$group==0,],cutgroup_n=S)
  data.X.Z_T <- PS_score_label_func(data.X.Z[data.X.Z$group==1,],cutgroup_n=S)
  c(data.T.d_C, cut.time_C) %<-% K.label.func.3(data.T.d = data.T.d[data.T.d$group==0,],data.X.Z_C,K=K) 
  c(data.T.d_T, cut.time_T) %<-% K.label.func.3(data.T.d = data.T.d[data.T.d$group==1,],data.X.Z_T,K=K) 
  c(data.T.d_T0, cut.time_T0) %<-% K.label.func.3(data.T.d = data.T.d0[data.T.d0$group==1,],data.X.Z_T,K=K) 
  c(data.T.d_C,data.rE_C) %<-% rE.func.for.Z.PS.K2(data.T.d_C,data.X.Z_C,cut.time_C,K=K)
  c(data.T.d_T,data.rE_T) %<-% rE.func.for.Z.PS.K2(data.T.d_T,data.X.Z_T,cut.time_T,K=K)
  c(data.T.d_T0,data.rE_T0) %<-% rE.func.for.Z.PS.K2(data.T.d_T0,data.X.Z_T,cut.time_T0,K=K)
  data.R.Z.K.S_C <- tauk.Z.K.S.func2(data.rE_C)
  data.R.Z.K.S_T <- tauk.Z.K.S.func2(data.rE_T)
  data.R.Z.K.S_T0 <- tauk.Z.K.S.func2(data.rE_T0)
  w.Z.K.S_C <- W.K.S.func2(data.rE_C)
  w.Z.K.S_T <- W.K.S.func2(data.rE_T)
  w.Z.K.S_T0 <- W.K.S.func2(data.rE_T0)
  
  c(lambda.K.prior.MCMC_T,w.for.ps) %<-% lambda.prior.MCMC.func4(data.rE_T,
                                                                 data.R.Z.K.S_T,
                                                                 w.Z.K.S_T,w.for.ps=NA,seed=seed)
  c(lambda.K.prior.MCMC_C,w.for.ps) %<-% lambda.prior.MCMC.func4(data.rE_C,
                                                                 data.R.Z.K.S_C,
                                                                 w.Z.K.S_C,w.for.ps=w.for.ps,seed=seed)
  c(lambda.K.prior.MCMC_T0,w.for.ps0) %<-% lambda.prior.MCMC.func4(data.rE_T0,
                                                                   data.R.Z.K.S_T0,
                                                                   w.Z.K.S_T0,w.for.ps=w.for.ps,seed=seed)
  return(cbind(lambda.K.prior.MCMC_T,lambda.K.prior.MCMC_C,lambda.K.prior.MCMC_T0))
}
forfunction_result_ALL <- function(i){
  seed=seed_simu[i]
  data.X.Z <- data.simu.X.Z.HR2(n.for.model_now = n.for.model_now,n.for.model_histr = n.for.model_histr,beta_RWD = beta_RWD,beta_ETD = beta_ETD,seed)
  data.T.d=data.simu.T.delta.NHR3(data.X.Z,alpha_CTD,alpha_ETD,alpha_RWD,beta,censor.parm.list,c.rate,shape=2,T.shuif=5,seed)
  data.T.d0=data.simu.T.delta.NHR3(data.X.Z,alpha_CTD = 0,alpha_ETD,alpha_RWD,beta,censor.parm.list,c.rate,shape=2,T.shuif=5,seed)
  data.X.Z_C <- PS_score_label_func(data.X.Z[data.X.Z$group==0,],cutgroup_n=S)
  data.X.Z_T <- PS_score_label_func(data.X.Z[data.X.Z$group==1,],cutgroup_n=S)
  c(data.T.d_C, cut.time_C) %<-% K.label.func.3(data.T.d = data.T.d[data.T.d$group==0,],data.X.Z_C,K=K) 
  c(data.T.d_T, cut.time_T) %<-% K.label.func.3(data.T.d = data.T.d[data.T.d$group==1,],data.X.Z_T,K=K) 
  c(data.T.d_T0, cut.time_T0) %<-% K.label.func.3(data.T.d = data.T.d0[data.T.d0$group==1,],data.X.Z_T,K=K) 
  c(data.T.d_C,data.rE_C) %<-% rE.func.for.Z.PS.K2(data.T.d_C,data.X.Z_C,cut.time_C,K=K)
  c(data.T.d_T,data.rE_T) %<-% rE.func.for.Z.PS.K2(data.T.d_T,data.X.Z_T,cut.time_T,K=K)
  c(data.T.d_T0,data.rE_T0) %<-% rE.func.for.Z.PS.K2(data.T.d_T0,data.X.Z_T,cut.time_T0,K=K)
  data.R.Z.K.S_C <- tauk.Z.K.S.func2(data.rE_C)
  data.R.Z.K.S_T <- tauk.Z.K.S.func2(data.rE_T)
  data.R.Z.K.S_T0 <- tauk.Z.K.S.func2(data.rE_T0)
  w.Z.K.S_C <- W.K.S.func2(data.rE_C)
  w.Z.K.S_T <- W.K.S.func2(data.rE_T)
  w.Z.K.S_T0 <- W.K.S.func2(data.rE_T0)
  
  lambda.K.prior.MCMC_T <-list_MCMC[,1:K]
  lambda.K.prior.MCMC_C <- list_MCMC[,(K+1):(2*K)]
  lambda.K.prior.MCMC_T0 <- list_MCMC[,(2*K+1):(3*K)]
  mixgamma.K_T %<-% lambda.MCMC.to.distr2(lambda.K.prior.MCMC_T)
  mixgamma.K_C %<-% lambda.MCMC.to.distr2(lambda.K.prior.MCMC_C)
  mixgamma.K_T0 %<-% lambda.MCMC.to.distr2(lambda.K.prior.MCMC_T0)
  lambda.unimform_T <- lambda.unimform.func(data.T.d[data.T.d$group==1,],cut.time_T)
  lambda.unimform_T0 <- lambda.unimform.func(data.T.d0[data.T.d0$group==1,],cut.time_T0)
  lambda.unimform_C <- lambda.unimform.func(data.T.d[data.T.d$group==0,],cut.time_C)
  
  result_SAM_T <- Wvfunc.lambda.k.SAM(mixgamma.K_T,data.rE_T,K=K,SAM_delta=SAM_delta)
  result_SAM_C <- Wvfunc.lambda.k.SAM(mixgamma.K_C,data.rE_C,K=K,SAM_delta=SAM_delta)
  result_SAM_T0 <- Wvfunc.lambda.k.SAM(mixgamma.K_T0,data.rE_T0,K=K,SAM_delta=SAM_delta)
  
  result_EBrmap_T<- Wvfunc.lambda.k.EBrmap(mixgamma.K_T,data.rE_T,K=K,P_cutpoint=P_cutpoint)
  result_EBrmap_C<- Wvfunc.lambda.k.EBrmap(mixgamma.K_C,data.rE_C,K=K,P_cutpoint=P_cutpoint)
  result_EBrmap_T0<- Wvfunc.lambda.k.EBrmap(mixgamma.K_T0,data.rE_T0,K=K,P_cutpoint=P_cutpoint)
  
  result_BOX_T<- Wvfunc.lambda.k.BOX(mixgamma.K_T,data.rE_T,K=K,P_cutpoint=P_cutpoint)
  result_BOX_C<- Wvfunc.lambda.k.BOX(mixgamma.K_C,data.rE_C,K=K,P_cutpoint=P_cutpoint)
  result_BOX_T0<- Wvfunc.lambda.k.BOX(mixgamma.K_T0,data.rE_T0,K=K,P_cutpoint=P_cutpoint)
  
  result_mypvalue_T<- Wvfunc.lambda.k.myPvalue2(mixgamma.K_T,data.rE_T,K=K,P_cutpoint=P_cutpoint_my)
  result_mypvalue_C<- Wvfunc.lambda.k.myPvalue2(mixgamma.K_C,data.rE_C,K=K,P_cutpoint=P_cutpoint_my)
  result_mypvalue_T0<- Wvfunc.lambda.k.myPvalue2(mixgamma.K_T0,data.rE_T0,K=K,P_cutpoint=P_cutpoint_my)
  
  result_func <- function(result_gamma_T,result_gamma_T0,result_gamma_C){
    c(yangb.distr_T,lambda.prior.MY.MAP1_T,WV.myMAP1_T,lambda.post.MY.MAP1_T,mean.MyMAP.post1_T,ESS_borrow1_T) %<-% result_gamma_T
    c(yangb.distr_T0,lambda.prior.MY.MAP1_T0,WV.myMAP1_T0,lambda.post.MY.MAP1_T0,mean.MyMAP.post1_T0,ESS_borrow1_T0) %<-% result_gamma_T0
    c(yangb.distr_C,lambda.prior.MY.MAP1_C,WV.myMAP1_C,lambda.post.MY.MAP1_C,mean.MyMAP.post1_C,ESS_borrow1_C) %<-% result_gamma_C
    
    medianT_T.sample=sample_of_medianTHR(lambda.post.MY.MAP1_T,cut.time_T)
    medianT_T0.sample=sample_of_medianTHR(lambda.post.MY.MAP1_T0,cut.time_T0)
    medianT_C.sample=sample_of_medianTHR(lambda.post.MY.MAP1_C,cut.time_C)
    mius_medianT.sample=medianT_T.sample-medianT_C.sample
    mius_medianT0.sample=medianT_T0.sample-medianT_C.sample
    if(quantile(mius_medianT.sample,0.05)<0){H.power=0}else{H.power=1}
    if(quantile(mius_medianT0.sample,0.05)<0){H.T1E=0}else{H.T1E=1}
    if(quantile(mius_medianT.sample,0.05)<alpha_CTD & alpha_CTD<quantile(mius_medianT.sample,0.975)){CI.cover=1}else{CI.cover=0}
    CI95_low <- quantile(mius_medianT.sample,0.025)
    CI95_up <- quantile(mius_medianT.sample,0.975)
    
    H_sim1 <- H.power;H_sim0 <- H.T1E;
    medianTD_postmean_sim1  <- mean(mius_medianT.sample)
    medianT.T1 <- mean(medianT_T.sample)
    medianT.C1 <- mean(medianT_C.sample)
    ESS_T1=ESS_borrow1_T
    ESS_C1=ESS_borrow1_C
    
    biasT1=medianT.T1-find_median_survival_time(lambda.unimform_T, cut.time_T)
    biasC1=medianT.C1-find_median_survival_time(lambda.unimform_C, cut.time_C)
    gc()
    return(c(CI95_low,CI95_up,CI.cover,H_sim1,H_sim0,medianTD_postmean_sim1,medianT.T1,medianT.C1,ESS_T1,ESS_C1,biasT1,biasC1))
  }

  result_SAM_all <- foreach(result_gamma_T=result_SAM_T,
                            result_gamma_T0=result_SAM_T0,
                            result_gamma_C=result_SAM_C,
                            #.combine = rbind,
                            .packages = c("RBesT","overlapping","dplyr" ,"survival" )
  )%do%result_func(result_gamma_T,result_gamma_T0,result_gamma_C)
  result_EBrmap_all <- foreach(result_gamma_T=result_EBrmap_T,
                               result_gamma_T0=result_EBrmap_T0,
                               result_gamma_C=result_EBrmap_C,
                               #.combine = rbind,
                               .packages = c("RBesT","overlapping","dplyr" ,"survival" )
  )%do%result_func(result_gamma_T,result_gamma_T0,result_gamma_C)
  result_BOX_all <- foreach(result_gamma_T=result_BOX_T,
                            result_gamma_T0=result_BOX_T0,
                            result_gamma_C=result_BOX_C,
                            #.combine = rbind,
                            .packages = c("RBesT","overlapping","dplyr" ,"survival" )
  )%do%result_func(result_gamma_T,result_gamma_T0,result_gamma_C)
  result_myPvalue_all <- foreach(result_gamma_T=result_mypvalue_T,
                                 result_gamma_T0=result_mypvalue_T0,
                                 result_gamma_C=result_mypvalue_C,
                                 #.combine = rbind,
                                 .packages = c("RBesT","overlapping","dplyr" ,"survival" )
  )%do%result_func(result_gamma_T,result_gamma_T0,result_gamma_C)
  
  
  result_SAM_all <- unlist(result_SAM_all)
  result_EBrmap_all <- unlist(result_EBrmap_all)
  result_BOX_all <- unlist(result_BOX_all)
  result_myPvalue_all <- unlist(result_myPvalue_all)
  
  names(result_SAM_all) <- paste(c("CI95_low","CI95_up","CI.cover","H_sim1","H_sim0","medianTD_postmean_sim1","medianT.T1",
                                   "medianT.C1","ESS_T1","ESS_C1","biasT1","biasC1") %>% rep(length(SAM_delta)),"_","SAM","_", SAM_delta %>% rep(each=12),sep = "")
  names(result_EBrmap_all) <- paste(c("CI95_low","CI95_up","CI.cover","H_sim1","H_sim0","medianTD_postmean_sim1","medianT.T1",
                                      "medianT.C1","ESS_T1","ESS_C1","biasT1","biasC1") %>% rep(length(P_cutpoint)),"_","EBrmap","_", P_cutpoint %>% rep(each=12),sep = "")
  names(result_BOX_all) <- paste(c("CI95_low","CI95_up","CI.cover","H_sim1","H_sim0","medianTD_postmean_sim1","medianT.T1",
                                   "medianT.C1","ESS_T1","ESS_C1","biasT1","biasC1") %>% rep(length(P_cutpoint)),"_","BOX","_", P_cutpoint %>% rep(each=12),sep = "")
  names(result_myPvalue_all) <- paste(c("CI95_low","CI95_up","CI.cover","H_sim1","H_sim0","medianTD_postmean_sim1","medianT.T1",
                                        "medianT.C1","ESS_T1","ESS_C1","biasT1","biasC1") %>% rep(length(P_cutpoint_my)),"_","myPvalue","_", P_cutpoint_my %>% rep(each=12),sep = "")
  result_all <- c(result_SAM_all,result_EBrmap_all,result_BOX_all,result_myPvalue_all)
  return(result_all)
}
#Parameter setting 1:fixed value parameter
K=5;S=4;beta=0.5;alpha_CTD=0.25;c.rate=0.3;T.shuif=5
n.for.model_now=300;n.for.model_histr=300
SAM_delta=c(0.03,0.05,0.1)
P_cutpoint=c(0.7,0.75,0.8,0.85,0.9)
P_cutpoint_my=c(0.7,0.75,0.8,0.85,0.9)
seed_simu=1:500
simu_time=1:500
#Parameter setting 2:parameter settings representing different simulation scenarios
canshu=rbind(data.frame(a=0.05,b=0.1,expand.grid(sort(c(seq(from=-0.2,to=0.2,by=0.05),-0.001)),c(0.05,0.1,0.15))),
             data.frame(a=0.15,b=0.2,expand.grid(sort(c(seq(from=-0.2,to=0.2,by=0.05),-0.001)),c(0.05,0.1,0.15))))
colnames(canshu)=c("beta_RWD","beta_ETD","alpha_RWD","diff")
library(foreach);library(zeallot)
foreach(i=1:nrow(canshu),.combine=rbind)%do%{if(canshu[i,3]<0)canshu[i,c(1,2,4)]=canshu[i,c(1,2,4)]*-1}
#500 simulation cycles with different seeds in every scenario
for(ii in 1:nrow(canshu)){
  c(beta_RWD,beta_ETD,alpha_RWD,diff)%<-%canshu[ii,]
  alpha_ETD=alpha_RWD+diff
  BLD_level <- ifelse(abs(beta_RWD)==0.05,"LOW","HIGH")
  Sign <- ifelse(alpha_RWD>=0,"P","N")
  seed=666
  censor.parm.list=censor.parm.for.censor.rate.weibull.NHR3(beta_RWD,beta_ETD,alpha_CTD,alpha_ETD,alpha_RWD,beta,
                                                            c.rate.list=c(0.2,0.3,0.5),shape=2,T.shuif=5,seed=seed)
  library(doParallel)
  library(foreach)
  cl<- makeCluster(34)
  registerDoParallel(cl)
  error_messages <- NULL
  lambda.SK.prior.MCMC <- foreach(
    i=simu_time,   
    .packages = c("SAMprior","base","tidyr","MASS","VGAM","RBesT","dplyr","keras",
                  "overlapping","nnet" ,"knitr","survival","distr","tidyr",
                  "zeallot","EBrmap","doParallel","foreach","BH","cmdstanr")
  ) %dopar% tryCatch({
    forfunction_MCMC(i)
  },error = function(e) {error_messages<<-rbind(error_messages,c(e$message,i));print(e$message)})
  stopCluster(cl)
  list_lambda.SK.prior.MCMC <- lambda.SK.prior.MCMC[which(lapply(lambda.SK.prior.MCMC, class)=="data.frame")]
  
  library(doParallel)
  library(foreach)
  cl<- makeCluster(34)
  registerDoParallel(cl)
  result_simu_i500 <- foreach(
    i=which(lapply(lambda.SK.prior.MCMC, class)=="data.frame"),   
    list_MCMC=list_lambda.SK.prior.MCMC,
    .combine=cbind,  
    .packages = c("SAMprior","base","tidyr","MASS","VGAM","RBesT","dplyr","keras",
                  "overlapping","nnet" ,"knitr","survival","distr","tidyr",
                  "zeallot","EBrmap","doParallel","foreach","BH","cmdstanr")
  ) %dopar% tryCatch({     forfunction_result_ALL(i)   },error = function(e) {error_messages<<-rbind(error_messages,c(e$message,i));print(e$message)})   
   
  stopCluster(cl)
  filename1 <- paste("D:/MSDB_prior/",BLD_level,"_MYMAP.ALL_",alpha_RWD,"-",abs(diff),"_",Sign,".CSV",sep = "")
  write.csv(result_simu_i500,file = filename1 )
}


