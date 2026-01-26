###############
## Libraries ##
###############

library(survival)
library(survminer)
library(tidyverse)
library(dplyr)
library(lubridate)
library(grid)
library(gridExtra)
library(readr)

##########################
## Read expression data ##
##########################


expr <- read.csv2("TMM conversion/TMM_TCGA_HNSC_counts_NoNormal_log2_filtered.csv", sep = ";", as.is = T, check.names = F)

#check if there is any gene duplicates 
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- tibble::rownames_to_column(expr, "genes")
dup<-expr[duplicated(expr$genes),] #check if there is any gene duplicates 

#check if there is any patient duplicates 
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.data.frame(t(expr))
expr <- tibble::rownames_to_column(expr, "id")
dup<-expr[duplicated(expr$id),]

trims <- expr %>% select(c("id", "ENSG00000119401", "ENSG00000136997"))

trims <- trims[!grepl("-06",trims$id),]
trims$id <- gsub("-01A", "",trims$id)
trims$id <- gsub("-01B", "",trims$id)

dup <- as.data.frame(duplicated(trims$id))


###################
## Clinical info ##
###################

info2 <- read.csv2("raw data/tcga_hnsc_updated_locations_03112023.csv")
info2 <- info2[,c(3,5)]
colnames(info2)[1] <- "Patient ID"

trims <- merge(trims, info2, by.x = "id",by.y = "Patient ID")

info <- read.csv2("raw data/Clinical_info_HNSC.csv", sep = ";", as.is = T, check.names = F)
head(info)


dss <- info[,-c(15,16, 41,42, 34, 35)]

pfi <- info[,-c(15:18,34,35)]
#Merging data
trim_dss <- merge(trims, dss, by.x = "id", by.y = "Patient ID")

trim_pfi <- merge(trims, pfi, by.x = "id", by.y = "Patient ID")


trim_dss$TRIM32_expression <- ifelse(trim_dss$ENSG00000119401 >= median(trim_dss$ENSG00000119401), 'High', "Low")

trim_pfi$TRIM32_expression <- ifelse(trim_pfi$ENSG00000119401 >= median(trim_pfi$ENSG00000119401), 'High', "Low")


trim_dss<- subset(trim_dss, trim_dss$`Months of disease-specific survival` > 0)

trim_pfi<- subset(trim_pfi, trim_pfi$`Progress Free Survival (Months)` > 0)

trim_dss <- trim_dss[!is.na(trim_dss$`Disease-specific Survival status`),]
trim_pfi <- trim_pfi[!is.na(trim_pfi$`Progression Free Status`),]

trim_dss$`Disease-specific Survival status` <- gsub("\\:.*","",trim_dss$`Disease-specific Survival status`)

trim_pfi$`Progression Free Status` <- gsub("\\:.*","",trim_pfi$`Progression Free Status`)


trim_dss$`Months of disease-specific survival` <- as.numeric(trim_dss$`Months of disease-specific survival`)

trim_pfi$`Progress Free Survival (Months)`<- as.numeric(trim_pfi$`Progress Free Survival (Months)`)


trim_dss$years<- trim_dss$`Months of disease-specific survival`/12 #Make DSS.months into years

trim_pfi$years<- trim_pfi$`Progress Free Survival (Months)`/12


trim_dss$`Disease-specific Survival status`<- as.integer(trim_dss$`Disease-specific Survival status`)

trim_pfi$`Progression Free Status` <- as.integer(trim_pfi$`Progression Free Status`)

#Cox regression: DSS 

head(trim_dss)

#order reference
head(trim_dss)
str(trim_dss)

table(trim_dss$`Disease-specific Survival status`,trim_dss$`Tumor Stage Code`)

trim_dss$TRIM32_expression <- factor(trim_dss$TRIM32_expression,
                                     levels = c("Low", "High"))
trim_dss$`Neoplasm Disease Stage American Joint Committee on Cancer Code` <- factor(trim_dss$`Neoplasm Disease Stage American Joint Committee on Cancer Code`, 
         levels = c("STAGE I", "STAGE II", "STAGE III", "STAGE IVA", "STAGE IVB", "STAGE IVC"),
         labels = c("STAGE I", "STAGE II", "STAGE III", "STAGE IV", "STAGE IV", "STAGE IV"))

trim_dss$`Neoplasm Histologic Grade` <- factor(trim_dss$`Neoplasm Histologic Grade`,
                                               levels = c("GX", "G1", "G2", "G3"))
trim_dss$`Neoplasm Histologic Grade`[trim_dss$`Neoplasm Histologic Grade`== "G4"] <- NA

trim_dss$`American Joint Committee on Cancer Tumor Stage Code`[trim_dss$`American Joint Committee on Cancer Tumor Stage Code` == "T0"] <- NA

trim_dss$`American Joint Committee on Cancer Tumor Stage Code` <- factor(trim_dss$`American Joint Committee on Cancer Tumor Stage Code`,
                                                                         levels = c("T1", "T2", "T3", "T4", "T4A", "T4B", "TX"),
                                                                         labels = c("T1", "T2", "T3", "T4", "T4", "T4", "TX"))


trim_dss$`Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code` <- factor(trim_dss$`Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`,
                                                                         levels = c("N0", "N1", "N2", "N2A", "N2B", "N2C", "N3", "NX"),
                                                                         labels = c("N0", "N1", "N2", "N2", "N2", "N2", "N3", "NX"))


trim_dss$site_of_resection_or_biopsy <- factor(trim_dss$site_of_resection_or_biopsy,
                                               levels = c("Tongue, NOS", "Pharynx, NOS", "Overlapping lesion of lip, oral cavity and pharynx",
                                                          "Mouth, NOS", "Larynx, NOS", "Gum, NOS", "Floor of mouth, NOS", "Cheek mucosa"))

trim_dss$`American Joint Committee on Cancer Tumor Stage Code` <- relevel(trim_dss$`American Joint Committee on Cancer Tumor Stage Code`, ref = "T1")
levels(trim_dss$`American Joint Committee on Cancer Tumor Stage Code`)



colnames(trim_dss)[colnames(trim_dss) == 'American Joint Committee on Cancer Tumor Stage Code'] <- "Tumor Stage Code"
colnames(trim_dss)[colnames(trim_dss) == 'Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code'] <- "Lymph Node Code"
colnames(trim_dss)[colnames(trim_dss) == 'American Joint Committee on Cancer Metastasis Stage Code'] <- "Metastasis Stage Code"

cox.model <- coxph(Surv(time = trim_dss$years, event = trim_dss$`Disease-specific Survival status`)~
                     `Sex`+
                     `Diagnosis Age`+
                     site_of_resection_or_biopsy+
                     `Subtype`+
                     `Lymph Node Code`+
                     `Metastasis Stage Code`+
                     `Aneuploidy Score` +
                     `Tumor Stage Code`+
                     `Ragnum Hypoxia Score`+`Winter Hypoxia Score` +
                     `Buffa Hypoxia Score` +
                     TRIM32_expression,
                   data = trim_dss)
cox.model

p <- ggforest(cox.model, data = trim_dss, fontsize = 1)
p

png("multivariate_cox_trim32_hnsc_dss.png", res= 200, height = 2500, width = 3000)
print(p)
dev.off()

pdf("multivariate_cox_trim32_hnsc_dss_2.pdf", height = 20, width = 15)
print(p)
dev.off()

#Cox regression:PFI 
head(trim_pfi)

table(trim_pfi$`Progression Free Status`, trim_pfi$`Tumor Stage Code`)

trim_pfi$TRIM32_expression <- factor(trim_pfi$TRIM32_expression,
                                     levels = c("Low", "High"))

trim_pfi$site_of_resection_or_biopsy <- factor(trim_pfi$site_of_resection_or_biopsy,
                                               levels = c("Tongue, NOS", "Pharynx, NOS", "Overlapping lesion of lip, oral cavity and pharynx",
                                                          "Mouth, NOS", "Larynx, NOS", "Gum, NOS", "Floor of mouth, NOS", "Cheek mucosa"))


colnames(trim_pfi)[colnames(trim_pfi) == 'American Joint Committee on Cancer Tumor Stage Code'] <- "Tumor Stage Code"
colnames(trim_pfi)[colnames(trim_pfi) == 'Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code'] <- "Lymph Node Code"
colnames(trim_pfi)[colnames(trim_pfi) == 'American Joint Committee on Cancer Metastasis Stage Code'] <- "Metastasis Stage Code"

trim_pfi$`Tumor Stage Code`[trim_pfi$`Tumor Stage Code` == "T0"] <- NA
trim_pfi$`Tumor Stage Code`<- factor(trim_pfi$`Tumor Stage Code`,levels = c("T1", "T2", "T3", "T4", "T4A", "T4B", "TX"),
                                                                 labels = c("T1", "T2", "T3", "T4", "T4", "T4", "TX"))


trim_pfi$`Lymph Node Code` <- factor(trim_pfi$`Lymph Node Code`,levels = c("N0", "N1", "N2", "N2A", "N2B", "N2C", "N3", "NX"),
                                                                labels = c("N0", "N1", "N2", "N2", "N2", "N2", "N3", "NX"))

cox.model <- coxph(Surv(time = trim_pfi$years, event = trim_pfi$`Progression Free Status`)~
                     `Sex`+
                     `Diagnosis Age`+
                     site_of_resection_or_biopsy+
                     `Subtype`+
                     `Lymph Node Code`+
                     `Metastasis Stage Code`+
                     `Aneuploidy Score` +
                     `Tumor Stage Code`+
                     `Ragnum Hypoxia Score`+`Winter Hypoxia Score` +
                     `Buffa Hypoxia Score` +
                     TRIM32_expression,
                   data = trim_pfi)

cox.model

p <- ggforest(cox.model, data = trim_dss, fontsize = 1)
p

png("multivariate_cox_trim32_hnsc_pfi.png", res= 200, height = 2500, width = 3000)
print(p)
dev.off()

pdf("multivariate_cox_trim32_hnsc_pfi_2.pdf", height = 20, width = 15)
print(p)
dev.off()
