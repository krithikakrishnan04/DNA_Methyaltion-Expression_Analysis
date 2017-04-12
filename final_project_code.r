#############################################
# Download and prepare the data
#
# Data description:
#   PAAD.met - PAAD methylation data - 485577 by 28 (25 samples)
#   PAAD.exp - PAAD RNA-seq data - 20531 by 183 (full samples)
#   clinical_paad_data - clinical data - 24 by 78 (24 samples or full samples)
#############################################

#############################################
# NOTE:
#   We only need to match the samples by the first 3 digits
#   of barcode (eg. TCGA-F2-7276)
#############################################

rm(list=ls())
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("TCGAbiolinks","SummarizedExperiment"))

library(TCGAbiolinks)
library(SummarizedExperiment)

#---------------------------------------------------------
# 1. Download DNA methylation data with TCGAbiolinks
#---------------------------------------------------------
DownloadFlag=F  # flip this flag to T to download data from TCGA

if(DownloadFlag) {
  path <- "paad"
  batchsample <- "TCGA-2J-AAB1-01"
  
  query <- TCGAquery(tumor = "PAAD", level = 3, platform = "HumanMethylation450", sample = batchsample)
  
  # How many sample we have?
  length(unlist(strsplit(query$barcode,",")))
  
  # Download the TCGA data
  TCGAdownload(query, path =path)
  
  df=T
  if(df){
    PAAD.met <- TCGAprepare(query, dir = path,
                            save = TRUE,
                            filename = "metPAAD.rda",summarizedExperiment = F,
                            add.subtype = TRUE)
  } else {
    # summarizedExperiment version (default)
    PAAD.met <- TCGAprepare(query, dir = path,
                            save = TRUE,
                            filename = "metPAAD.rda",summarizedExperiment = T,
                            add.subtype = TRUE)
  }
  
  #---------------------------------------------------------
  # 2. download RNA expression data
  #---------------------------------------------------------
  query.rna <- TCGAquery(tumor="PAAD",level=3, platform="IlluminaHiSeq_RNASeqV2", sample = batchsample)
  
  # How many sample we have?
  length(unlist(strsplit(query.rna$barcode,",")))
  
  TCGAdownload(query.rna,path=path,type = "rsem.genes.normalized_results")
  
  df=T
  if(df) {
    PAAD.exp <- TCGAprepare(query.rna, dir=path, save = TRUE,
                            type = "rsem.genes.normalized_results",
                            filename = "expPAAD.rda", summarizedExperiment = F,
                            add.subtype = TRUE)
  } else {
    PAAD.exp <- TCGAprepare(query.rna, dir=path, save = TRUE,
                            type = "rsem.genes.normalized_results",
                            filename = "expPAAD.rda", summarizedExperiment = T,
                            add.subtype = TRUE)
  }
  
  #
  # Add a CrossReact column to the the data frame
  # also add a SNP column
  
  # load nonSpecific
  nonSpecificTable=read.table("Final project/RData/nonspecific probes Illumina450k.table",sep="\t",header = T)
  
  # load SNP
  load("Final project/RData/SNPinfo.RData")
  
  # Identify which records to take from SNPinfo
  SNPflag=(rowSums(!is.na(SNPinfo[,3:6]))!=0)
  badSNPs=rownames(SNPinfo)[SNPflag]
  matchDX=match(badSNPs,rownames(PAAD.met))
  
  PAAD.met$SNPflag=F
  PAAD.met$SNPflag[matchDX]=T
  
  if(F){ # double check our work
    CheckBad=rownames(PAAD.met)[which(PAAD.met$SNPflag)]
    CheckMatch=match(CheckBad,rownames(SNPinfo))
    SNPinfo[CheckMatch,]
  }
  
  # now, create a CrossReact flag
  matchDX=match(nonSpecificTable[,1],rownames(PAAD.met))
  PAAD.met$CrossReactflag=F
  PAAD.met$CrossReactflag[matchDX]=T
  
  # type I/II probe designation
  require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  class(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  MyAnn=getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
  matchDX=match(MyAnn$Name,rownames(PAAD.met))
  
  PAAD.met$Type=NA
  PAAD.met$Type[matchDX]=MyAnn$Type
  
  #-----------------------------
  # We can obtain clinical data
  # by TCGAquery_clinic
  #-----------------------------
  
  clinical_paad_data_all <- TCGAquery_clinic("paad","clinical_patient")
  clinical_paad_data <- TCGAquery_clinic("paad","clinical_patient",sample=unlist(strsplit(query$barcode,",")))
  
  save.image(file="Final project/RData/LoadData.RData")
  
}

if(!DownloadFlag) load("Final project/RData/LoadData.RData")

################################
### Pre-Processing Data
################################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
rm(list=ls())

load("Final project/RData/LoadData.RData")

################################################
# Look at our data set first, and manipulation
# these can be used in our data decription
#------------------------------------------------
# ## Note: we have 3 dataset:
#         1. Methylation data
#         2. RNA-seq data
#         3. Clinical data
################################################

#--------------
# Methylation
#--------------
# See how many type I we have in met data
table(PAAD.met$Type)

# select Type I only (met)
paad.met <- subset(PAAD.met,Type=="I")

#--------------
# RNA-seq
#--------------
# remove unknown gene (exp)
J=dim(PAAD.exp)[1]
genenames=rep("",J)
for(j in 1:J) genenames[j]=unlist(strsplit(rownames(PAAD.exp)[j],"\\|"))[1]
rows.to.keep <- which(genenames!="?")
paad.exp <- PAAD.exp[rows.to.keep,]

#--------------
# Clinical
#--------------

# look at the data
# matched to met data
attach(clinical_paad_data)
table(age_at_initial_pathologic_diagnosis)
table(tobacco_smoking_history)
table(pathologic_stage)
table(histological_type)
table(gender)

median(as.numeric(age_at_initial_pathologic_diagnosis))

detach(clinical_paad_data)

# unmatched (all samples)
attach(clinical_paad_data_all)
table(age_at_initial_pathologic_diagnosis)
table(tobacco_smoking_history)
table(pathologic_stage)
table(histological_type)
table(gender)

median(as.numeric(age_at_initial_pathologic_diagnosis))

detach(clinical_paad_data_all)

# Recode the covariates by making cut-off
# Age group
clinical_paad_data$agegrp <- ifelse(clinical_paad_data$age_at_initial_pathologic_diagnosis > 65, 1,0)
clinical_paad_data_all$agegrp <- ifelse(clinical_paad_data_all$age_at_initial_pathologic_diagnosis > 65, 1,0)

# Smoking group
clinical_paad_data$smoking <- ifelse(as.numeric(clinical_paad_data$tobacco_smoking_history) <= 1, 0, 1)
clinical_paad_data_all$smoking <- ifelse(as.numeric(clinical_paad_data_all$tobacco_smoking_history) <= 1, 0, 1)

# pathologic_stage
clinical_paad_data$path_stage[clinical_paad_data$pathologic_stage=="Stage I"|
                                clinical_paad_data$pathologic_stage=="Stage IA"|
                                clinical_paad_data$pathologic_stage=="Stage IB"] <- "Stage I"

clinical_paad_data$path_stage[clinical_paad_data$pathologic_stage=="Stage II"|
                                clinical_paad_data$pathologic_stage=="Stage IIA"|
                                clinical_paad_data$pathologic_stage=="Stage IIB"] <- "Stage II"

clinical_paad_data$path_stage[clinical_paad_data$pathologic_stage=="Stage III"|
                                clinical_paad_data$pathologic_stage=="Stage IV"] <- "Stage III & IV"

clinical_paad_data_all$path_stage[clinical_paad_data_all$pathologic_stage=="Stage I"|
                                    clinical_paad_data_all$pathologic_stage=="Stage IA"|
                                    clinical_paad_data_all$pathologic_stage=="Stage IB"] <- "Stage I"

clinical_paad_data_all$path_stage[clinical_paad_data_all$pathologic_stage=="Stage II"|
                                    clinical_paad_data_all$pathologic_stage=="Stage IIA"|
                                    clinical_paad_data_all$pathologic_stage=="Stage IIB"] <- "Stage II"

clinical_paad_data_all$path_stage[clinical_paad_data_all$pathologic_stage=="Stage III"|
                                    clinical_paad_data_all$pathologic_stage=="Stage IV"] <- "Stage III & IV"

# histological_type
clinical_paad_data$hist_type <- ifelse(clinical_paad_data$histological_type =="Pancreas-Adenocarcinoma Ductal Type", 1,0)
clinical_paad_data_all$hist_type <- ifelse(clinical_paad_data_all$histological_type =="Pancreas-Adenocarcinoma Ductal Type", 1,0)

## Note: we need recode other covariates..

#------------------------------------
# Manipulate 3 data set for analysis
#------------------------------------
# match clinical sample to exp data (183/185)
sample.exp <-  substr(colnames(paad.exp),1,12)

sample.cli <- clinical_paad_data_all$bcr_patient_barcode

length(sample.exp)  # 183
length(unique(sample.exp)) # 178
length(sample.cli) # 185
length(unique(sample.cli)) # 185

## note: there are 183-178=5 replicated patients (maybe multiple samples per patient)
## we need to remove them

sample.select <- unique(intersect(sample.exp,sample.cli))

paad.cli <- clinical_paad_data_all[which(clinical_paad_data_all$bcr_patient_barcode %in% sample.select),]
paad.rna <- paad.exp[,-which(duplicated(sample.exp))]

# then sort by barcode
paad.cli <- paad.cli[order(paad.cli$bcr_patient_barcode),]

colnames(paad.rna) <- substr(colnames(paad.rna),1,12)

paad.rna <- paad.rna[, order(names(paad.rna))]

# Note: so far we have 3 dataset, include info in our data description:
#         1. Methylation data (paad.met) 38/40
#         2. RNA-seq data (paad.rna) 178
#         3. Clinical data (paad.cli) 178 

# save.image("Final project/RData/ReadyToAnalysis.RData")

load("Final project/RData/ReadyToAnalysis.RData")

#------------------------------------------------
# Methylation
#   - Filter by previous RNA-seq results
#   - Build a linear model for all 5 covariates
#------------------------------------------------

#-----------------------------------
# Load p-values and rank them
#-----------------------------------

load("Final project/RData/pval_lm_rna.RData")
load("Final project/RData/pval_rna.RData")

# adjust by FDR first
myPvals_Lm <- p.adjust(myPvals_Lm, method="fdr", n=length(myPvals_Lm))

# rank p-values
bestK=order(myPvals_Lm)[1:100]

# find the most significant gene
bestK.gene=rep("",100)
for(j in 1:100) bestK.gene[j]=unlist(strsplit(rownames(paad.rna[bestK,])[j],"\\|"))[1]

# map to methylation data
paad.met.best <- subset(paad.met, Gene_Symbol %in% bestK.gene)

#-------------------------------------------------
# use the "best" set of methylation data to test
#-------------------------------------------------
nrow(na.omit(paad.met.best)) # 849
paad.met.best <- na.omit(paad.met.best)

# not run
paad.met.best.b <- na.omit(paad.met.best)[,-c(1:3,44:46)]

# match clinical sample to exp data (38/40)
sample.met <-  substr(colnames(paad.met.best.b),1,12)
sample.flag <- unique(intersect(sample.met,sample.cli))
paad.cli.match <- paad.cli[which(paad.cli$bcr_patient_barcode %in% sample.flag),]

# remove duplicate patitent
paad.met.best.b<- paad.met.best.b[,-which(duplicated(sample.met))]
paad.met.best.match <- paad.met.best.b[,which(substr(colnames(paad.met.best.b),1,12) %in% sample.flag)]

# then sort by barcode
paad.cli.match <- paad.cli.match[order(paad.cli.match$bcr_patient_barcode),]
colnames(paad.met.best.match) <- substr(colnames(paad.met.best.match),1,12)
paad.met.best.match <- paad.met.best.match[, order(names(paad.met.best.match))]

# save and load 2nd manipulated data
# save(paad.met.best,paad.met.best.match,paad.cli.match,file="Final project/RData/ReadyToAnalysis_lm.RData")

dim(paad.met.best.match)
dim(paad.cli.match)

load("Final project/RData/ReadyToAnalysis_lm.RData")

# adopt previous code
age.grp <- as.factor(paad.cli.match$agegrp)
smoking <- as.factor(paad.cli.match$smoking)
path.stage <- as.factor(paad.cli.match$path_stage)
hist.type <- as.factor(paad.cli.match$hist_type)
gender <- as.factor(paad.cli.match$gender)
#------------------------------------------------
# Do a lm scan
#------------------------------------------------
X=as.matrix(paad.met.best.match)
myLmp=function(i){
  cat(i,"...",fill=F)
  lm <- lm(X[i,]~age.grp + smoking + path.stage + hist.type)
  f <- summary(lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
myPvals_Lm_met=mapply(myLmp,1:(dim(X)[1]))
myPvals_Lm_met <- p.adjust(myPvals_Lm_met, method="fdr", n=length(myPvals_Lm_met))
DMR <- data.frame(Gene_Symbol=paad.met.best$Gene_Symbol, pval=myPvals_Lm_met)
GeneList.met <- unique(DMR$Gene_Symbol)

Gpval.met <- c()
for (i in 1:length(GeneList.met)) {
  idx = which(DMR$Gene_Symbol %in% GeneList.met[i])
  Gpval.met[i] <- min(DMR[idx,]$pval)
}

DMR.g <- data.frame(Gene_Symbol=GeneList.met, pval=Gpval.met)

#------------------------------------------------
# Match and compare met vs. rna
#------------------------------------------------
match.idx <- match(GeneList.met,bestK.gene)
bestK.gene[match.idx]
Gpval.rna <- myPvals_Lm[bestK][match.idx]
df <- data.frame(Gene_Symbol=DMR.g$Gene_Symbol, pval.met=DMR.g$pval, pval.rna=Gpval.rna)
p.lm <- ggplot(data = df, aes(x = log10(df$pval.met), y = log10(df$pval.rna))) +
  geom_point() +
  geom_hline(aes(yintercept = log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = log10(0.05)), linetype = "dashed") +
  xlab("DNA methylation \nLog10(FDR corrected P values)") +
  ylab("RNA-seq \nLog10(FDR corrected P values)") + 
  theme_minimal()
x11()
print(p.lm)

pdf("p_lm.pdf")
print(p.lm)
dev.off()

#------------------------------------------------
# RNA-seq
# Do a quick kruskal-wallis scan
# Marginal test for each covariates
#------------------------------------------------
load("Final project/RData/ReadyToAnalysis.RData")
X=as.matrix(paad.rna)
#----
age.grp <- as.factor(paad.cli$agegrp)
myKrusk=function(i){
  cat(i,"...",fill=F)
  kruskal.test(x=X[i,],g=age.grp)$p.value
}
myPvals_age=mapply(myKrusk,1:(dim(X)[1]))

#----
smoking <- as.factor(paad.cli$smoking)
myKrusk=function(i){
  cat(i,"...",fill=F)
  kruskal.test(x=X[i,],g=smoking)$p.value
}

myPvals_smoke=mapply(myKrusk,1:(dim(X)[1]))

#----
path.stage <- as.factor(paad.cli$path_stage)
myKrusk=function(i){
  cat(i,"...",fill=F)
  kruskal.test(x=X[i,],g=path.stage)$p.value
}

myPvals_path=mapply(myKrusk,1:(dim(X)[1]))

#----
hist.type <- as.factor(paad.cli$hist_type)
myKrusk=function(i){
  cat(i,"...",fill=F)
  kruskal.test(x=X[i,],g=hist.type)$p.value
}
myPvals_hist=mapply(myKrusk,1:(dim(X)[1]))

#----
gender <- as.factor(paad.cli$gender)
myKrusk=function(i){
  cat(i,"...",fill=F)
  kruskal.test(x=X[i,],g=gender)$p.value
}

myPvals_gender=mapply(myKrusk,1:(dim(X)[1]))

# save(myPvals_age,myPvals_smoke,myPvals_path,myPvals_hist,myPvals_gender,file="Final project/RData/pval_rna.RData")

#------------------------------------------------
# Build a linear model for all 5 covariates
#------------------------------------------------
myLmp=function(i){
  cat(i,"...",fill=F)
  lm <- lm(X[i,]~age.grp + smoking + path.stage + hist.type)
  f <- summary(lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
myPvals_Lm=mapply(myLmp,1:(dim(X)[1]))

# save(myPvals_Lm, file="Final project/RData/pval_lm_rna.RData")

#------------------------------------------------
# Methylation
#   - Filter by previous RNA-seq results
#   - Do a quick kruskal-wallis scan
#     Marginal test for each covariates
#------------------------------------------------
load("Final project/RData/pval_lm_rna.RData")
load("Final project/RData/pval_rna.RData")

##################
# age
##################
# adjust by FDR first
myPvals_age <- p.adjust(myPvals_age, method="fdr", n=length(myPvals_age))

# rank p-values
bestK=order(myPvals_age)[1:100]

# find the most significant gene
bestK.gene=rep("",100)
for(j in 1:100) bestK.gene[j]=unlist(strsplit(rownames(paad.rna[bestK,])[j],"\\|"))[1]

# map to methylation data
paad.met.best <- subset(paad.met, Gene_Symbol %in% bestK.gene)
#-------------------------------------------------
# use the "best" set of methylation data to test
#-------------------------------------------------
nrow(na.omit(paad.met.best)) 
paad.met.best <- na.omit(paad.met.best)
# not run
paad.met.best.b <- na.omit(paad.met.best)[,-c(1:3,44:46)]

# match clinical sample to exp data (38/40)
sample.met <-  substr(colnames(paad.met.best.b),1,12)
sample.flag <- unique(intersect(sample.met,sample.cli))
paad.cli.match <- paad.cli[which(paad.cli$bcr_patient_barcode %in% sample.flag),]

# remove duplicate patitent
paad.met.best.b<- paad.met.best.b[,-which(duplicated(sample.met))]
paad.met.best.match <- paad.met.best.b[,which(substr(colnames(paad.met.best.b),1,12) %in% sample.flag)]

# then sort by barcode
paad.cli.match <- paad.cli.match[order(paad.cli.match$bcr_patient_barcode),]

colnames(paad.met.best.match) <- substr(colnames(paad.met.best.match),1,12)
paad.met.best.match <- paad.met.best.match[, order(names(paad.met.best.match))]

# save and load 2nd manipulated data
# save(paad.met.best,paad.met.best.match,paad.cli.match,file="Final project/RData/ReadyToAnalysis_age.RData")
dim(paad.met.best.match)
dim(paad.cli.match)

load("Final project/RData/ReadyToAnalysis_age.RData")

# adopt previous code
age.grp <- as.factor(paad.cli.match$agegrp)
smoking <- as.factor(paad.cli.match$smoking)
path.stage <- as.factor(paad.cli.match$path_stage)
hist.type <- as.factor(paad.cli.match$hist_type)
gender <- as.factor(paad.cli.match$gender)
#------------------------------------------------
# Do a age scan
#------------------------------------------------
X=as.matrix(paad.met.best.match)
myKrusk=function(i){
  cat(i,"...",fill=F)
  kruskal.test(x=X[i,],g=age.grp)$p.value
}
myPvals_age_met=mapply(myKrusk,1:(dim(X)[1]))
myPvals_age_met <- p.adjust(myPvals_age_met, method="fdr", n=length(myPvals_age_met))
DMR <- data.frame(Gene_Symbol=paad.met.best$Gene_Symbol, pval=myPvals_age_met)
GeneList.met <- unique(DMR$Gene_Symbol)
Gpval.met <- c()
for (i in 1:length(GeneList.met)) {
  idx = which(DMR$Gene_Symbol %in% GeneList.met[i])
  Gpval.met[i] <- min(DMR[idx,]$pval)
}
DMR.g <- data.frame(Gene_Symbol=GeneList.met, pval=Gpval.met)
#------------------------------------------------
# Match and compare met vs. rna
#------------------------------------------------
match.idx <- match(GeneList.met,bestK.gene)
bestK.gene[match.idx]
Gpval.rna <- myPvals_age[bestK][match.idx]
df <- data.frame(Gene_Symbol=DMR.g$Gene_Symbol, pval.met=DMR.g$pval, pval.rna=Gpval.rna)
p.age <- ggplot(data = df, aes(x = log10(df$pval.met), y = log10(df$pval.rna))) +
  geom_point() +
  geom_hline(aes(yintercept = log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = log10(0.05)), linetype = "dashed") +
  xlab("DNA methylation \nLog10(FDR corrected P values)") +
  ylab("RNA-seq \nLog10(FDR corrected P values)") + 
  ggtitle("For age group") +
  theme_minimal()
x11()
print(p.age)

pdf("p_age.pdf")
print(p.age)
dev.off()

load("Final project/RData/pval_lm_rna.RData")
load("Final project/RData/pval_rna.RData")

#--- omit another covariates---#


############################################
# For full data set
# use TCGAbiolink for analysis
#############################################
rm(list=ls())
load("RData/LoadData.RData")
load("RData/LoadData_cli.RData")
load("RData/LoadData_met_df.RData")

# attach clinical infomation to both data set
length(colData(paad.exp)$patient)
length(colData(paad.met)$patient)
length(paad.cli$bcr_patient_barcode)

names(paad.cli)[79:82]

matchDX= match(colData(paad.exp)$patient, paad.cli$bcr_patient_barcode) 
matchDX_met= match(colData(paad.met)$patient, paad.cli$bcr_patient_barcode) 

colData(paad.exp)$agegrp <- paad.cli$agegrp[matchDX]
colData(paad.exp)$smoking <- paad.cli$smoking[matchDX]
colData(paad.exp)$path_stage <- paad.cli$path_stage[matchDX]
colData(paad.exp)$hist_type <- paad.cli$hist_type[matchDX]
colData(paad.exp)$gender <- paad.cli$gender[matchDX]

colData(paad.met)$agegrp <- paad.cli$agegrp[matchDX_met]
colData(paad.met)$smoking <- paad.cli$smoking[matchDX_met]
colData(paad.met)$path_stage <- paad.cli$path_stage[matchDX_met]
colData(paad.met)$hist_type <- paad.cli$hist_type[matchDX_met]
colData(paad.met)$gender <- paad.cli$gender[matchDX_met]

#--------------------------------
# load nonSpecific
nonSpecificTable=read.table("RData/nonspecific probes Illumina450k.table",sep="\t",header = T)

# load SNP
load("RData/SNPinfo.RData")

# Identify which records to take from SNPinfo
SNPflag=(rowSums(!is.na(SNPinfo[,3:6]))!=0)
badSNPs=rownames(SNPinfo)[SNPflag]
matchDX=match(badSNPs,rownames(PAAD.met))

PAAD.met$SNPflag=F
PAAD.met$SNPflag[matchDX]=T

if(F){ # double check our work
  CheckBad=rownames(PAAD.met)[which(PAAD.met$SNPflag)]
  CheckMatch=match(CheckBad,rownames(SNPinfo))
  SNPinfo[CheckMatch,]
}

# now, create a CrossReact flag
matchDX=match(nonSpecificTable[,1],rownames(PAAD.met))
PAAD.met$CrossReactflag=F
PAAD.met$CrossReactflag[matchDX]=T

# type I/II probe designation
require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

class(IlluminaHumanMethylation450kanno.ilmn12.hg19)

MyAnn=getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
matchDX=match(MyAnn$Name,rownames(PAAD.met))

PAAD.met$Type=NA
PAAD.met$Type[matchDX]=MyAnn$Type

# select Type I only (met)
# PAAD.met <- subset(PAAD.met,Type=="I")

# filter Type I
TypeI.DX <- which(PAAD.met$Type=="I")
PAAD.met <- subset(PAAD.met,Type=="I")
paad.met <- paad.met[TypeI.DX,]

# filter SNPflag
SNP.DX <- which(PAAD.met$SNPflag==F)
paad.met <- paad.met[SNP.DX,]
save(paad.cli,paad.met,paad.exp,query.met,paad.query.exp,file="RData/ManipulatedData.RData")

#-----------------------------
# Age group
#-----------------------------
rm(list=ls())
load("RData/ManipulatedData.RData")

##########################################################
# removing the samples without subtype classification
paad.aux <- subset(paad.met,select = colData(paad.met)$agegrp %in% c("65AndAbove","LessThan65"))

##########################################################
TCGAvisualize_meanMethylation(paad.aux,
                              groupCol = "agegrp",
                              subgroupCol = "gender",
                              group.legend  = "Groups",
                              subgroup.legend = "Gender",
                              filename = "paad_mean_agegrp.png")

# na.omit
paad.aux <- subset(paad.aux,subset = (rowSums(is.na(assay(paad.aux))) == 0))

# Volcano plot
paad.aux <- TCGAanalyze_DMR(paad.aux, groupCol = "agegrp",
                            group1 = "65AndAbove",
                            group2="LessThan65",
                            p.cut = 0.05,
                            diffmean.cut = 0.2,
                            legend = "State",
                            plot.filename = "paad_agegrp_met_volcano.png")

#----------------------
paad.exp.aux <- subset(paad.exp, select = colData(paad.exp)$agegrp %in% c("65AndAbove","LessThan65"))

idx <- colData(paad.exp.aux)$agegrp %in% c("65AndAbove")
idx2 <- colData(paad.exp.aux)$agegrp %in% c("LessThan65")

dataPrep <- TCGAanalyze_Preprocessing(object = paad.exp.aux,
                                      cor.cut = 0.0)

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  qnt.cut = 0.25,
                                  method='quantile')

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,idx],
                            mat2 = dataFilt[,idx2],
                            Cond1type = "65AndAbove",
                            Cond2type = "LessThan65",
                            method = "glmLRT")

TCGAVisualize_volcano(dataDEGs$logFC,dataDEGs$FDR,
                      filename = "paad_agegrp_exp_volcano.png",
                      x.cut = 3,
                      y.cut = 0.05,
                      names = rownames(dataDEGs),
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      title = "Volcano plot (LessThan65 vs 65AndAbove)")

#---------------------
starburst_agegrp <- TCGAvisualize_starburst(paad.aux, dataDEGs,"65AndAbove","LessThan65",
                                            filename = "paad_agegrp_starburst.png",
                                            met.p.cut = 0.05,
                                            exp.p.cut = 0.05,
                                            diffmean.cut = 0.25,
                                            logFC.cut = 3,
                                            names = TRUE)
