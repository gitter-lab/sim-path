library(survival)
library(R.matlab)
library(impute)
library(pamr)
library(mclust)


## Load data

# Data was kindly provided by Győrffy et al. 
# For more information please see: http://www.kmplot.com
library_data <- '~/projects/cancer_strat/data/nbs_release_v0.2'

OV_expression <- read.table(file=paste(library_data,'/data/OV_expression/ovar_1394_mas5_normalized.txt',sep=""),header=1,sep="\t",row.names=1)
OV_exp_pheno <- read.table(file=paste(library_data,'/data/OV_expression/expdesc_ovar_111101_final_batch.txt',sep=""),header=1,sep="\t")
OV_HM90_K4 <- read.table(file=paste(library_data,'/data/NBS_precalculated_stratification/NBS_HM90_netnmf_K4_label.tsv',sep=""),head=FALSE,sep="\t",row.names=1)


## Orgainze data
OV_exp_pheno$AffyID <- gsub('-','.',OV_exp_pheno$AffyID)
OV_exp_pheno$AffyID <- sub('55000240','X55000240',OV_exp_pheno$AffyID)
OV_exp_pheno$AffyID <- sub('HG_U133A_96_','HG.U133A_96.',OV_exp_pheno$AffyID)
rownames(OV_exp_pheno) <- OV_exp_pheno$AffyID

join_id_list <- intersect(colnames(OV_expression),OV_exp_pheno[,1])

OV_exp_pheno <- OV_exp_pheno[join_id_list,]
OV_expression <- OV_expression[,join_id_list]

colnames(OV_expression) <- OV_exp_pheno$Patient.ID
rownames(OV_exp_pheno) <- OV_exp_pheno$Patient.ID

## 5-Year survival

zselect60 <- which(OV_exp_pheno[,4] > 60)
OV_exp_pheno[zselect60,4] <- 60
OV_exp_pheno[zselect60,3] <- 0

## Impute missing values

OV_expression_i <- impute.knn(as.matrix(OV_expression))

## Setup data sets 
TCGA_names <- rownames(OV_HM90_K4)
TCGA_names <- TCGA_names[TCGA_names %in% rownames(OV_exp_pheno)]

OV_TCGA_sm <- OV_expression_i$data[,TCGA_names]
OV_TCGA_sm_labels <- OV_HM90_K4[TCGA_names,]
OV_TCGA_sm_pheno <- OV_exp_pheno[TCGA_names,]

OV_meta_full <- OV_expression_i$data[,setdiff(colnames(OV_expression_i$data),TCGA_names)]
OV_meta_full_pheno <- OV_exp_pheno[setdiff(rownames(OV_exp_pheno),TCGA_names),]

## Train model 

pamr.model.NBS.TCGA <- pamr.train(list(x=OV_TCGA_sm,y=factor(OV_TCGA_sm_labels)))
pamr.cvres.NBS.TCGA <- pamr.cv(pamr.model.NBS_TCGA,list(x=OV_TCGA_sm,y=factor(OV_TCGA_sm_labels)))

OV_TCGA_sm_labels.predict <- pamr.predict(pamr.model.NBS.TCGA,OV_TCGA_sm,threshold=0.075,type="class")
adjustedRandIndex(OV_TCGA_sm_labels.predict,OV_TCGA_sm_labels)

## Plot survival trained model

fitOrg <- survfit(Surv(OV_TCGA_sm_pheno$OS_time, OV_TCGA_sm_pheno$OS_event..1.death.) ~ factor(OV_HM90_K4[TCGA_names,]) )
plot(fitOrg,col=c("blue","green","magenta","cyan"),xlim=c(0,60),xlab='Time(months)',ylab="Probablity of survival",lwd=3)
llgened.label <- paste('Cluster ',unique(OV_TCGA_sm_labels.predict)," (",table(OV_TCGA_sm_labels.predict),')',sep="")
legend(x=5,y=0.3,llgened.label,col=c("blue","green","magenta","cyan"),lty=1,lwd=3)

## Apply to meta anlaysis 

OV_meta_full_labels.predict <- pamr.predict(pamr.model.NBS.TCGA,OV_meta_full,threshold=0.075,type="class")
cphm <- coxph(Surv(OV_meta_full_pheno$OS_time, OV_meta_full_pheno$OS_event..1.death.) ~ factor(OV_meta_full_labels.predict))
summary(cphm)

fitMeta <- survfit(Surv(OV_meta_full_pheno$OS_time, OV_meta_full_pheno$OS_event..1.death.) ~ factor(OV_meta_full_labels.predict) )
plot(fitMeta,col=c("blue","green","magenta","cyan"),xlim=c(0,60),xlab='Time(months)',ylab="Probablity of survival",lwd=3)
llgened.label <- paste('Cluster ',unique(OV_meta_full_labels.predict )," (",table(OV_meta_full_labels.predict ),')',sep="")
legend(x=5,y=0.3,llgened.label,col=c("blue","green","magenta","cyan"),lty=1,lwd=3)

## Tothill

beffect <- 4
pbatch <- which(OV_exp_pheno[,11] == beffect & OV_exp_pheno[,7] == 1)

OV_tothil <- OV_expression_i$data[,pbatch]
OV_tothil_pheno <- OV_exp_pheno[pbatch,]

pp_pamr <- pamr.predict(pamr.model.NBS.TCGA,OV_tothil,threshold=0.075,type="class")

fcov.stage <- factor(OV_tothil_pheno[,8])
fcov.grade <- factor(OV_tothil_pheno[,9])
fcov.debulk <- factor(OV_tothil_pheno[,10])
fcov.debulk[is.na(fcov.debulk)] <- 1

tothil_data = na.omit(data.frame(time=OV_tothil_pheno[,4],event=OV_tothil_pheno[,3],debulk=fcov.debulk,stage=fcov.stage,grade=fcov.grade,clustFact=factor(pp_pamr)))
cphm_full <- coxph(Surv(time,event) ~ clustFact + stage + grade + debulk,tothil_data)
cphm_base <- coxph(Surv(time,event) ~ stage + grade + debulk,tothil_data)
summary(cphm_full)
anova(cphm_full,cphm_base,test='Chisq')

fit1 <- survfit(Surv(OV_tothil_pheno[,4], OV_tothil_pheno[,3]) ~ factor(pp_pamr))

plot(fit1,col=c("blue","green","magenta","cyan"),xlim=c(0,60),xlab='Time(months)',ylab="Probablity of survival",lwd=3)
llgened.label <- paste('Cluster ',levels(pp_pamr)," (",table(pp_pamr),')',sep="")
legend(x=5,y=0.3,llgened.label,col=c("blue","green","magenta","cyan"),lty=1,lwd=3)

## Bonome 

beffect <- 6
pbatch <- which(OV_exp_pheno[,11] == beffect & !is.na(OV_exp_pheno[,10]) )

OV_bonome <- OV_expression_i$data[,pbatch]
OV_bonome_pheno <- OV_exp_pheno[pbatch,]

pp_pamr <- pamr.predict(pamr.model.NBS.TCGA,OV_bonome,threshold=0.075,type="class")

#fcov.stage <- factor(OV_bonome_pheno[,8])
#fcov.grade <- factor(OV_bonome_pheno[,9])
fcov.debulk <- factor(OV_bonome_pheno[,10])
fcov.debulk[is.na(fcov.debulk)] <- 1

bonome_data = na.omit(data.frame(time=OV_bonome_pheno[,4],event=OV_bonome_pheno[,3],debulk=fcov.debulk,clustFact=factor(pp_pamr)))
cphm_full <- coxph(Surv(time,event) ~ clustFact + debulk,bonome_data)
cphm_base <- coxph(Surv(time,event) ~ debulk,bonome_data)
summary(cphm_full)
anova(cphm_full,cphm_base,test='Chisq')

fit1 <- survfit(Surv(OV_bonome_pheno[,4], OV_bonome_pheno[,3]) ~ factor(pp_pamr))

plot(fit1,col=c("blue","green","magenta","cyan"),xlim=c(0,60),xlab='Time(months)',ylab="Probablity of survival",lwd=3)
llgened.label <- paste('Cluster ',levels(pp_pamr)," (",table(pp_pamr),')',sep="")
legend(x=5,y=0.3,llgened.label,col=c("blue","green","magenta","cyan"),lty=1,lwd=3)

