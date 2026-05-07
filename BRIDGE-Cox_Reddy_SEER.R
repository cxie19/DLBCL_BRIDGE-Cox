suppressMessages(library(dplyr)) 
suppressMessages(library(reticulate))
suppressMessages(library(survival))
suppressMessages(library(TransCox))
suppressMessages(library(doParallel))
suppressMessages(library(tdROC))
suppressMessages(library(glmnet))
suppressMessages(library(survC1))
suppressMessages(library(fastDummies))

use_condaenv("trans_cox", required = TRUE) 
source_python(system.file("python", "TransCoxFunction.py", package = "TransCox"))

#### Read in the datasets ####
# the datasets could be downloaded from the online supplementary files
save.path <- "~" # specify the directory
reddy_base <- read.csv(paste0(save.path,"Reddy dataset.csv"))
seer_base <- read.csv(paste0(save.path,"SEER dataset.csv"))
nci <- read.csv(paste0(save.path,"NCI validation cohort.csv"))

# remove row names
reddy_base <- reddy_base[,-1]
seer_base <- seer_base[,-1]
nci <- nci[,-1]
pData = reddy_base
aData = seer_base

aData$status <- ifelse(aData$status==1,2,1)
pData$status.tuning <- pData$status
pData$status <- ifelse(pData$status==1,2,1)

colnames(nci)[colnames(nci)=="os_years"] <- "time"
colnames(nci)[colnames(nci)=="os_status_0censor_1event"] <- "status"
colnames(nci)[colnames(nci)=="stage_IPI_0no_1yes"] <- "stage_IPI"  
colnames(nci)[colnames(nci)=="ECOG_PS_IPI_0no_1yes"] <-'ECOG.IPI'
colnames(nci)[colnames(nci)=="LDH_IPI_0no_1yes"] <-'LDH.IPI'
colnames(nci)[colnames(nci)=="extranodal_IPI_0no_1yes"] <-'multiple.extranodal.IPI'

X_var <- c("age_cat_le40","age_cat_41_60","age_cat_gt75","stage_IPI")
Z_var <- c( "LDH.IPI","ECOG.IPI","multiple.extranodal.IPI")

#### Model fitting ####
#### stage 1 ####
Cout <- GetAuxSurv(aData, cov = X_var)
Pout <- GetPrimaryParam(pData, q = Cout$q, estR = Cout$estR)

# criterion is auc
learning_rate_vec =c(seq(0.0001, 0.001,by=0.0001),0.002,0.003,0.004,0.005)
nsteps_vec = seq(400,1000,by=100)
SelLR_By_AUC <- matrix(NA,nrow=length(learning_rate_vec)*length(nsteps_vec),ncol=3)
index <- 0
for (i in learning_rate_vec){
  for (j in nsteps_vec){
    index <- index +1
    Tres <- runTransCox_one(Pout, l1 = 0.1, l2 = 0.1,
                            learning_rate = i, nsteps = j,
                            cov = X_var)
    
    if(sum(Tres$new_IntH<0)==0){
      dis_time <- Tres$time
      
      # Baseline cumulative hazards estimation
      newh <- data.frame(time=Tres$time,new_h=Tres$new_IntH)
      trans_base_H_time_event <- data.frame(cum_hazard=cumsum(newh$new_h),time=newh$time)
    
      trans_base_S <- exp(-trans_base_H_time_event$cum_hazard)
      # base_S at one time point
      trans_base_S <- ifelse(sum(dis_time<=max(pData$time))!=0,trans_base_S[sum(dis_time<=max(pData$time))],0)
      trans_predict_surv <- sapply(seq(nrow(pData)),function(x)
        trans_base_S^exp(Tres$new_beta%*%t(pData[x,X_var])))
      trans_predict_surv <- t(trans_predict_surv)
      
      SelLR_By_AUC[index,] <- try({
        trans_auc <- tdROC( X = 1-trans_predict_surv, Y = pData$time, delta =  pData$status.tuning, tau = max(pData$time), span = 0.01, nboot = 0)
        
        SelLR_By_AUC[index,] <- c(i,j,trans_auc$main_res$AUC.empirical)
      },
      silent = TRUE)
    }else{
      SelLR_By_AUC[index,] <- c(i,j,NA)
    }
  }
}
SelLR_By_AUC <- SelLR_By_AUC[!is.na(SelLR_By_AUC[,3]),]

lambda1_vec = c(0.1, 0.5, seq(1, 10, by = 0.5))
lambda2_vec = c(0.1, 0.5, seq(1, 10, by = 0.5))
SelParam_By_AUC <- matrix(NA,nrow=length(lambda1_vec)*length(lambda2_vec),ncol=3)
index <- 0
for(i in lambda1_vec){
  
  for (j in lambda2_vec){
    index <- index+1
    Tres <- runTransCox_one(Pout, l1 = i,
                            l2 = j,
                            learning_rate = SelLR_By_AUC[SelLR_By_AUC[,3]==max(SelLR_By_AUC[,3]),1][1],
                            nsteps = SelLR_By_AUC[SelLR_By_AUC[,3]==max(SelLR_By_AUC[,3]),2][1],
                            cov = X_var)
    if(sum(Tres$new_IntH<0)==0){
      dis_time <- Tres$time
      
      newh <- data.frame(time=Tres$time,new_h=Tres$new_IntH)
      trans_base_H_time_event <- data.frame(cum_hazard=cumsum(Tres$new_IntH),time=Tres$time)
      
      trans_base_S <- exp(-trans_base_H_time_event$cum_hazard)
      # base_S at one time point
      trans_base_S <-  ifelse(sum(dis_time<=max(pData$time))!=0,trans_base_S[sum(dis_time<=max(pData$time))],0)
      trans_predict_surv <- sapply(seq(nrow(pData)),function(x)
        trans_base_S^exp(Tres$new_beta%*%t(pData[x,X_var])))
      trans_predict_surv <- t(trans_predict_surv)
      
      # compute auc
      SelParam_By_AUC[index,] <- try({
        trans_auc <- tdROC( X = 1-trans_predict_surv, Y = pData$time, delta =  pData$status.tuning, tau = max(pData$time), span = 0.01, nboot = 0)
        SelParam_By_AUC[index,] <- c(i,j,trans_auc$main_res$AUC.empirical)
      },
      silent = TRUE)
    }else{
      SelParam_By_AUC[index,] <- c(i,j,NA)
    }
  }
}
SelParam_By_AUC <- SelParam_By_AUC[!is.na(SelParam_By_AUC[,3]),]

est_result <- runTransCox_one(Pout, l1 = SelParam_By_AUC[SelParam_By_AUC[,3]==max(SelParam_By_AUC[,3]),1][1],
                              l2 = SelParam_By_AUC[SelParam_By_AUC[,3]==max(SelParam_By_AUC[,3]),2][1],
                              learning_rate = SelLR_By_AUC[SelLR_By_AUC[,3]==max(SelLR_By_AUC[,3]),1][1],
                              nsteps = SelLR_By_AUC[SelLR_By_AUC[,3]==max(SelLR_By_AUC[,3]),2][1],
                              cov = X_var)
beta <- est_result$new_beta
pData$offset <-  sapply(1:nrow(pData),function(x) t(beta)%*%t(pData[x,X_var]))

#### stage 2 ####
reddy.cox <- coxph(Surv(time,status.tuning) ~ ECOG.IPI+LDH.IPI+multiple.extranodal.IPI+offset(offset),data = pData)
summary(reddy.cox)

#### stage 3 for subgroup ####
# offset 
pData$lp_reddy <- predict(
  reddy.cox,
  newdata = pData,
  type   = "lp"    # linear predictor = X β̂  (including your original offset term)
)

# COO
pData_subgroup <- pData[!is.na(pData$ABC.GCB..RNAseq.),]
table(pData_subgroup$ABC.GCB..RNAseq.)
# create dummy variables
pData_subgroup <- dummy_cols(pData_subgroup, 
                             select_columns =c("ABC.GCB..RNAseq."), 
                             remove_selected_columns = FALSE, 
                             remove_most_frequent_dummy  = FALSE)
reddy.subgroupcoo.cox <- coxph(
  Surv(time, status) ~ ABC.GCB..RNAseq._GCB + ABC.GCB..RNAseq._Unclassified + offset(lp_reddy),
  data = pData_subgroup
)
summary(reddy.subgroupcoo.cox)

# MYD88 
table(pData$MYD88,useNA = "always")
reddy.subgroupmyd88.cox <- coxph(
  Surv(time, status) ~ MYD88 + offset(lp_reddy),
  data = pData
)
summary(reddy.subgroupmyd88.cox)

#### External validation #####
# add offset values for individuals
nci$offset <- sapply(1:nrow(nci),function(x) t(beta)%*%t(nci[x,X_var]))
nci$lp_reddy <- predict(
  reddy.cox,
  newdata = nci,
  type   = "lp"    # linear predictor = X β̂  (including your original offset term)
)

# All patients
predict_whole_test <- survfit(reddy.cox,nci)
# C-index for overall survival
nci$risk_score <- predict(reddy.cox, newdata=nci,type="lp") 
C=Inf.Cval(nci[,c("time","status","risk_score")],tau=max(pData$time[pData$status.tuning==1]), itr=nrow(nci)) 
round(c(C$Dhat,C$low95, C$upp95), digits=3)

# GCB subgroup 
nci_gcb <- nci[nci$Gene_expression_subgroup=="GCB",]
nci_gcb$ABC.GCB..RNAseq._GCB <- 1
nci_gcb$ABC.GCB..RNAseq._Unclassified <- 0
predict_whole_test <- survfit(reddy.subgroupcoo.cox,nci_gcb)
C=Inf.Cval(nci_gcb[,c("time","status","risk_score")],tau=max(nci_gcb$time[nci_gcb$status==1]), itr=nrow(nci_gcb)) 
round(c(C$Dhat,C$low95, C$upp95), digits=3)

# ABC subgroup
nci_abc <- nci[nci$Gene_expression_subgroup=="ABC",]
nci_abc$ABC.GCB..RNAseq._GCB <- 0
nci_abc$ABC.GCB..RNAseq._Unclassified <- 0
predict_whole_test <- survfit(reddy.subgroupcoo.cox,nci_abc)
C=Inf.Cval(nci_abc[,c("time","status","risk_score")],tau=max(nci_abc$time[nci_abc$status==1]), itr=nrow(nci_abc)) 
round(c(C$Dhat,C$low95, C$upp95), digits=3)

# MYD88 subgroup
nci_myd88 <- nci[nci$MYD88==1,]
predict_whole_test <- survfit(reddy.subgroupmyd88.cox,nci_myd88)
C=Inf.Cval(nci[,c("time","status","risk_score")],tau=max(nci_myd88$time[nci_myd88$status==1]), itr=nrow(nci_myd88)) 
round(c(C$Dhat,C$low95, C$upp95), digits=3)

#### AUC and Brier score: same code for all patients and subgroups based on object predict_whole_test ####
# AUC and Brier score for 2-, 3-, 5- year survival probabilities
each_time_metric <- foreach(metric_time=c(2,3,5),.combine = "rbind")%do%{
  
  predict_surv_whole_point <- predict_whole_test[["surv"]][sum(predict_whole_test[["time"]]<=metric_time),]
  auc_whole <- tdROC( X = 1-predict_surv_whole_point, Y = nci$time, delta = nci$status, tau = metric_time, span = 0.01, nboot = 1000)
  metrics <- c(auc_whole$main_res$AUC.empirical,auc_whole$boot_res$bAUC.empirical$CIlow,auc_whole$boot_res$bAUC.empirical$CIhigh,
               auc_whole$calibration_res[1],auc_whole$boot_res$bBS$CIlow,auc_whole$boot_res$bBS$CIhigh) #
  return(c(metric_time,round(metrics,3)))
}

colnames(each_time_metric) <- c("Time","AUC","AUC_CIlow","AUC_CIhigh","BS","BS_CIlow","BS_CIhigh")
each_time_metric

