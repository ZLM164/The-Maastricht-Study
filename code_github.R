setwd('D:/document/ALL project/polygneic risk score/AV_635_Zhewen_Ren')

install.packages("remotes")
install.packages("readxl")
install.packages("dplyr")
install.packages(c("survival", "survminer"))
install.packages('mlogit')
install.packages('quantreg')

library(remotes)
library(readxl)
library(dplyr)
library(foreign)
library("survival")
library("survminer")
library(mlogit)
library(quantreg)

####################################################################################Interaction of IHL and PRS for TG clearance on TG
### read dataset
TG <- readxl::read_xlsx("TG_genes.xlsx")
TG_yes <- subset(TG,TG$`Role in TG metabolism`==1)
TG_no <- subset(TG,TG$`Role in TG metabolism`==2)
data <- read.csv("AV635.csv", na.strings = c("", " ",'NA'))

### exclude people with missing data
TG_clearance <- data[,c(1,767:769,774:783,785:787,791,803,827:832,841:846,852:858,which(substring(names(data),4) %in% TG_yes$SNP))]
TG_clearance <- subset(TG_clearance,is.na(TG_clearance$liverfat_adj)==FALSE) # 4007 without liver MRI
TG_clearance <- subset(TG_clearance, TG_clearance$liverfat_adj != -999) # 159 invalid liver fat data
TG_clearance <- subset(TG_clearance,is.na(TG_clearance$batch)== FALSE) #411 people without genotyping
TG_clearance <- TG_clearance[complete.cases(TG_clearance[,39:56]),] #593 people without full genotyping data
TG_clearance <- subset(TG_clearance,TG_clearance$NIT_alcoholtot != -999) # 129 people with without alcohol intake data
TG_clearance <- subset(TG_clearance,TG_clearance$NIT_alcoholtot != -888) # 78 people with implausible alcohol intake 

### generate allele dosage
for (i in seq(from = 39, to = 55, by = 2)) {
  for (j in 1:nrow(TG_clearance)) {
    if ((TG_clearance[j,i] == TG_clearance[j,i+1]) & (TG_clearance[j,i] == TG$`Effect allele`[which(TG$SNP == substring(colnames(TG_clearance)[i],4))])){
      TG_clearance[j,56+(i-37)/2] <- 2
    } else if ((TG_clearance[j,i] == TG_clearance[j,i+1]) & (TG_clearance[j,i] == TG$`Non-efect allele`[which(TG$SNP == substring(colnames(TG_clearance)[i],4))])){
      TG_clearance[j,56+(i-37)/2] <- 0
    } else {
      TG_clearance[j,56+(i-37)/2] <- 1
    }
  }
}

for (i in seq(from = 39, to = 55, by = 2)) {
  names(TG_clearance)[56+(i-37)/2] <-  substring(colnames(TG_clearance)[i],4)
}


### generate PRS
for (i in 57:ncol(TG_clearance)){
  TG_clearance[,i+9] <- (TG_yes$Beta[which(TG_yes$SNP == colnames(TG_clearance)[i])]) * TG_clearance[,i] 
}
TG_clearance$PRS <- rowSums(TG_clearance[,c(66:74)])


### generate new variable
#T2D
T2D <- NULL
for (i in 1:nrow(TG_clearance)){
  if (TG_clearance$N_Diabetes_WHO_2[i]==2){
    T2D <- c(T2D,1)
  }
  else{T2D <- c(T2D,0)
  }
} 
TG_clearance["T2D"]<- T2D 
TG_clearance$T2D=as.factor(TG_clearance$T2D)

# steatosis
steotosis <- NULL
for (i in 1:nrow(TG_clearance)){
  if (TG_clearance$liverfat_adj[i] < 5.89){
    steotosis <- c(steotosis,0)
  }else { 
    steotosis <- c(steotosis,1)
  }
}
TG_clearance['steotosis'] <- steotosis
TG_clearance$steotosis = as.factor(TG_clearance$steotosis)

# hypertriglyceridemia
htg <- NULL
for (i in 1:nrow(TG_clearance)){
  if (TG_clearance$Triglyc[i] >= 2.3){
    htg <- c(htg, 1)
  }else {
    htg <- c(htg,0)
  }
}
TG_clearance["hypertrigly"] <- htg
TG_clearance$hypertrigly = as.factor(TG_clearance$hypertrigly)

# PRS and TG in tertile
quantile(TG_clearance$PRS,c((1/3),(2/3)))  
PRS_tertile <- cut(TG_clearance$PRS,breaks = c(-Inf,0.4480,0.5415,Inf),labels = c(1,2,3),include.lowest = T)
TG_clearance['PRS_tertile'] <- PRS_tertile

quantile(TG_clearance$Triglyc,c((1/3),(2/3)))  
TG_tertile <- cut(TG_clearance$Triglyc,breaks = c(-Inf,0.95, 1.44, Inf),labels = c(1,2,3),include.lowest = T)
TG_clearance['TG_tertile'] <- TG_tertile

### association analysis
# IHL-TG
#Model 1
lm1 <- lm(formula = log10(Triglyc) ~ liverfat_adj,data = TG_clearance)
hist(residuals(lm1))
sum1 <- summary(lm1)
coef1 <- data.frame(sum1$coefficients)
confint1 <- data.frame(confint(lm1))
paste(round(exp(coef1$Estimate[2]),3),"(", round(exp(confint1$X2.5..[2]),3),",",round(exp(confint1$X97.5..[2]),3),")")

#Model 2
lm1 <- lm(formula = log10(Triglyc) ~ liverfat_adj + Age + SEX + T2D + MRI_lagtime,data = TG_clearance)
hist(residuals(lm1))
sum1 <- summary(lm1)
coef1 <- data.frame(sum1$coefficients)
confint1 <- data.frame(confint(lm1))
paste(round(exp(coef1$Estimate[2]),3),"(", round(exp(confint1$X2.5..[2]),3),",",round(exp(confint1$X97.5..[2]),3),")")

#Model 3   
lm1 <- lm(formula = log10(Triglyc) ~ liverfat_adj + Age + SEX + T2D + MRI_lagtime + liverfat_adj * MRI_lagtime + med_LP + NIT_alcoholtot,data = TG_clearance)
hist(residuals(lm1))
sum1 <- summary(lm1)
print(sum1)
coef1 <- data.frame(sum1$coefficients)
confint1 <- data.frame(confint(lm1))
paste(round(exp(coef1$Estimate[2]),3),"(", round(exp(confint1$X2.5..[2]),3),",",round(exp(confint1$X97.5..[2]),3),")")

# PRS-TG
# Model 1 
lm1 <- lm(formula = log10(Triglyc) ~ PRS,data = TG_clearance)
hist(residuals(lm1))
sum1 <- summary(lm1)
coef1 <- data.frame(sum1$coefficients)
confint1 <- data.frame(confint(lm1))
paste(round(exp(coef1$Estimate[2]),3),"(", round(exp(confint1$X2.5..[2]),3),",",round(exp(confint1$X97.5..[2]),3),")")

# Model 2
lm1 <- lm(formula = log10(Triglyc) ~ PRS + Age + SEX + T2D + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3
          +C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3,data = TG_clearance)
hist(residuals(lm1))
sum1 <- summary(lm1)
coef1 <- data.frame(sum1$coefficients)
confint1 <- data.frame(confint(lm1))
paste(round(exp(coef1$Estimate[2]),3),"(", round(exp(confint1$X2.5..[2]),3),",",round(exp(confint1$X97.5..[2]),3),")")

# Model 3
lm1 <- lm(formula = log10(Triglyc) ~ PRS + Age + SEX + T2D + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3
          +C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot,data = TG_clearance)
hist(residuals(lm1))
sum1 <- summary(lm1)
print(sum1)
coef1 <- data.frame(sum1$coefficients)
confint1 <- data.frame(confint(lm1))
paste(round(exp(coef1$Estimate[2]),3),"(", round(exp(confint1$X2.5..[2]),3),",",round(exp(confint1$X97.5..[2]),3),")")

### Interaction analysis-- IHL and PRS are continuous
# Model 1
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj*PRS, data = TG_clearance)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)]),3),")",sep = '')

# Model 2
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj*PRS + Age + SEX + T2D +MRI_lagtime + liverfat_adj * MRI_lagtime +
                       +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3, data = TG_clearance)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

# Model 3
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + SEX + T2D + MRI_lagtime + liverfat_adj * MRI_lagtime +
                       +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot, data = TG_clearance)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
print(sum_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),5),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

# additionally adjusted for bmi
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + SEX + T2D + MRI_lagtime + liverfat_adj * MRI_lagtime +
                       +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot + bmi, data = TG_clearance)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
print(sum_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),5),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

### interaction analysis-steatosis and PRS in tertile
group <- NULL
for (i in 1:nrow(TG_clearance)){
  if (TG_clearance$PRS_tertile[i] ==1 & TG_clearance$steotosis[i] == 0) {
    group <- c(group,1)
  } else if  (TG_clearance$PRS_tertile[i] ==1 & TG_clearance$steotosis[i] == 1) {
    group <- c(group,2)
  } else if  (TG_clearance$PRS_tertile[i] ==2 & TG_clearance$steotosis[i] == 0) {
    group <- c(group,3)
  } else if  (TG_clearance$PRS_tertile[i] ==2 & TG_clearance$steotosis[i] == 1) {
    group <- c(group,4)
  } else if  (TG_clearance$PRS_tertile[i] ==3 & TG_clearance$steotosis[i] == 0) {
    group <- c(group,5)
  } else if  (TG_clearance$PRS_tertile[i] ==3 & TG_clearance$steotosis[i] == 1) {
    group <- c(group,6)
  }
}
TG_clearance$group <- group
TG_clearance$group = as.factor(TG_clearance$group)


# model 3
logistic_model <- glm(hypertrigly ~ group + Age + SEX + T2D + MRI_lagtime + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 + steotosis * MRI_lagtime +
                        C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot,data = TG_clearance,family = "binomial"(link='logit'))
summary(logistic_model)
OR <- round(exp(coef(logistic_model)),3)
CI <- round(exp(confint(logistic_model)),3)

### Interaction analysis-- IHL and PRS are continuous-stratification by sex
MALE <- subset(TG_clearance,SEX == 1)
FEMALE <- subset(TG_clearance,SEX == 2)

# for male
# Model 1
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj*PRS, data = MALE)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)]),3),")",sep = '')

# Model 2
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + T2D + MRI_lagtime + liverfat_adj * MRI_lagtime
                     +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3, data = MALE)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

# Model 3
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + T2D + MRI_lagtime
                     +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot + liverfat_adj * MRI_lagtime, data = MALE)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
print(sum_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

# for female
# Model 1
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj*PRS, data = FEMALE)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)]),3),")",sep = '')

# Model 2
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + T2D + MRI_lagtime + liverfat_adj * MRI_lagtime
                     +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3, data = FEMALE)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

# Model 3
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + T2D + MRI_lagtime
                     +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot + liverfat_adj * MRI_lagtime, data = FEMALE)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
print(sum_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')


### Interaction analysis-- IHL and PRS are continuous-stratification by T2D status
T2D <- subset(TG_clearance,T2D == 1)
non_T2D <- subset(TG_clearance,T2D == 0)

#for T2D
# Model 1
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj*PRS, data = T2D)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)]),3),")",sep = '')

# Model 2
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + SEX + MRI_lagtime
                     +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3+ liverfat_adj * MRI_lagtime, data = T2D)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

# Model 3
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + SEX + MRI_lagtime
                     +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot + liverfat_adj * MRI_lagtime, data = T2D)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
print(sum_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

#for non-T2D
# Model 1
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj*PRS, data = non_T2D)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)]),3),")",sep = '')

# Model 2
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + SEX + MRI_lagtime
                     +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3+ liverfat_adj * MRI_lagtime, data = non_T2D)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

# Model 3
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + SEX + MRI_lagtime
                     +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot + liverfat_adj * MRI_lagtime, data = non_T2D)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
print(sum_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

### Interaction analysis-- IHL and PRS are continuous-stratification by lipid modification
lipid_modify <- subset(TG_clearance, med_LP == 1)
non_lipid_modify <- subset(TG_clearance, med_LP == 0)

# for individuals with lipid modification
# Model 1
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj*PRS, data = lipid_modify)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)]),3),")",sep = '')

# Model 2
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + T2D + SEX + MRI_lagtime
                     +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3+ liverfat_adj * MRI_lagtime, data = lipid_modify)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

# Model 3
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + T2D + SEX + MRI_lagtime
                     +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3 + NIT_alcoholtot + liverfat_adj * MRI_lagtime, data =lipid_modify)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
print(sum_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

# for individuals without lipid modification
# Model 1
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj*PRS, data = non_lipid_modify)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)]),3),")",sep = '')

# Model 2
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + T2D + SEX + MRI_lagtime
                     +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3+ liverfat_adj * MRI_lagtime, data = non_lipid_modify)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

# Model 3
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + T2D + SEX + MRI_lagtime
                     +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3 + NIT_alcoholtot + liverfat_adj * MRI_lagtime, data = non_lipid_modify)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
print(sum_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

### interation analysis-steatosis and PRS intertile-stratified by sex

# for male model 3
logistic_model <- glm(hypertrigly ~ group + Age + T2D + MRI_lagtime + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 + 
                        C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot + steotosis * MRI_lagtime,data = MALE,family = "binomial"(link='logit'))
summary(logistic_model)
OR <- round(exp(coef(logistic_model)),1)
CI <- round(exp(confint(logistic_model)),1)

# for female model 3
logistic_model <- glm(hypertrigly ~ group + Age + T2D + MRI_lagtime + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 + 
                        C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot + steotosis * MRI_lagtime,data = FEMALE,family = "binomial"(link='logit'))
summary(logistic_model)
OR <- round(exp(coef(logistic_model)),1)
CI <- round(exp(confint(logistic_model)),1)

### interation analysis-steatosis and PRS intertile-stratified by T2D status
# for T2D model 3
logistic_model <- glm(hypertrigly ~ group + Age + SEX + MRI_lagtime + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 + 
                        C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot + steotosis * MRI_lagtime,data = T2D,family = "binomial"(link='logit'))
summary(logistic_model)
OR <- round(exp(coef(logistic_model)),2)
CI <- round(exp(confint(logistic_model)),2)

# for non-T2D model 3
logistic_model <- glm(hypertrigly ~ group + Age + SEX + MRI_lagtime + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 + 
                        C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot + steotosis * MRI_lagtime,data = non_T2D,family = "binomial"(link='logit'))
summary(logistic_model)
OR <- round(exp(coef(logistic_model)),2)
CI <- round(exp(confint(logistic_model)),2)


### interation analysis-steatosis and PRS intertile-stratified by lipid modification
# for lipid-modification model 3
logistic_model <- glm(hypertrigly ~ group + Age + SEX + T2D + MRI_lagtime + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 + 
                        C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3 +NIT_alcoholtot + steotosis * MRI_lagtime,data = lipid_modify,family = "binomial"(link='logit'))
summary(logistic_model)
OR <- round(exp(coef(logistic_model)),1)
CI <- round(exp(confint(logistic_model)),1)


# for non-lipid-modification model 3
logistic_model <- glm(hypertrigly ~ group + Age + SEX + T2D + MRI_lagtime + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 + 
                        C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3 +NIT_alcoholtot + steotosis * MRI_lagtime,data = non_lipid_modify,family = "binomial"(link='logit'))
summary(logistic_model)
OR <- round(exp(coef(logistic_model)),1)
CI <- round(exp(confint(logistic_model)),1)

##############################################################################Interaction of IHL and PRS for TG production on TG
### exclude people with missing data
non_TG_clearance <- data[,c(1,767:769,774:783,785:787,791,803,827:832,841:846,852:858,which(substring(names(data),4) %in% TG_no$SNP))]
non_TG_clearance <- subset(non_TG_clearance,is.na(non_TG_clearance$liverfat_adj)==FALSE) # 4007 without liver MRI
non_TG_clearance <- subset(non_TG_clearance, non_TG_clearance$liverfat_adj != -999) # 159 invalid liver fat data
non_TG_clearance <- subset(non_TG_clearance,is.na(non_TG_clearance$batch)== FALSE) #411 people without genotyping
non_TG_clearance <- non_TG_clearance[complete.cases(non_TG_clearance[,39:62]),] #642 people without full genotyping data
non_TG_clearance <- subset(non_TG_clearance,non_TG_clearance$NIT_alcoholtot != -999) # 126 people with without alcohol intake data
non_TG_clearance <- subset(non_TG_clearance,non_TG_clearance$NIT_alcoholtot != -888) # 73 people with implausible alcohol intake 

### generate allele dosage
for (i in seq(from = 39, to = 61, by = 2)) {
  for (j in 1:nrow(non_TG_clearance)) {
    if ((non_TG_clearance[j,i] == non_TG_clearance[j,i+1]) & (non_TG_clearance[j,i] == TG_no$`Effect allele`[which(TG_no$SNP == substring(colnames(non_TG_clearance)[i],4))])){
      non_TG_clearance[j,62+(i-37)/2] <- 2
    } else if ((non_TG_clearance[j,i] == non_TG_clearance[j,i+1]) & (non_TG_clearance[j,i] == TG_no$`Non-efect allele`[which(as.vector(TG_no$SNP) == substring(colnames(non_TG_clearance)[i],4))])){
      non_TG_clearance[j,62+(i-37)/2] <- 0
    } else {
      non_TG_clearance[j,62+(i-37)/2] <- 1
    }
  }
}

for (i in seq(from = 39, to = 61, by = 2)) {
  names(non_TG_clearance)[62+(i-37)/2] <-  substring(colnames(non_TG_clearance)[i],4)
}

### generate PRS
for (i in 63:ncol(non_TG_clearance)){
  non_TG_clearance[,i+12] <- (TG_no$Beta[which(TG_no$SNP == colnames(non_TG_clearance)[i])]) * non_TG_clearance[,i] 
}
non_TG_clearance$PRS <- rowSums(non_TG_clearance[,c(75:86)])


### generate new variable
#T2D
T2D <- NULL
for (i in 1:nrow(non_TG_clearance)){
  if (non_TG_clearance$N_Diabetes_WHO_2[i]==2){
    T2D <- c(T2D,1)
  }
  else{T2D <- c(T2D,0)
  }
} 
non_TG_clearance["T2D"]<- T2D 
non_TG_clearance$T2D=as.factor(non_TG_clearance$T2D)

#### association analysis:PRS-TG
# Model 1 
lm1 <- lm(formula = log10(Triglyc) ~ PRS,data = non_TG_clearance)
hist(residuals(lm1))
sum1 <- summary(lm1)
coef1 <- data.frame(sum1$coefficients)
confint1 <- data.frame(confint(lm1))
paste(round(exp(coef1$Estimate[2]),3),"(", round(exp(confint1$X2.5..[2]),3),",",round(exp(confint1$X97.5..[2]),3),")")

# Model 2
lm1 <- lm(formula = log10(Triglyc) ~ PRS + Age + SEX + T2D + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3
          +C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3,data = non_TG_clearance)
hist(residuals(lm1))
sum1 <- summary(lm1)
coef1 <- data.frame(sum1$coefficients)
confint1 <- data.frame(confint(lm1))
paste(round(exp(coef1$Estimate[2]),3),"(", round(exp(confint1$X2.5..[2]),3),",",round(exp(confint1$X97.5..[2]),3),")")

# Model 3
lm1 <- lm(formula = log10(Triglyc) ~ PRS + Age + SEX + T2D + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3
          +C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot,data = non_TG_clearance)
hist(residuals(lm1))
sum1 <- summary(lm1)
print(sum1)
coef1 <- data.frame(sum1$coefficients)
confint1 <- data.frame(confint(lm1))
paste(round(exp(coef1$Estimate[2]),3),"(", round(exp(confint1$X2.5..[2]),3),",",round(exp(confint1$X97.5..[2]),3),")")


### Interaction analysis-- IHL and PRS are continuous
# Model 1
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj*PRS, data = non_TG_clearance)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)]),3),")",sep = '')

# Model 2
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj*PRS + Age + SEX + T2D + MRI_lagtime + liverfat_adj * MRI_lagtime +
                       +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3, data = non_TG_clearance)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

# Model 3
lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + SEX + T2D + MRI_lagtime+ liverfat_adj * MRI_lagtime
                     + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                       C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot , data = non_TG_clearance)
hist(residuals(lm_interaction))
sum_interaction <- summary(lm_interaction)
print(sum_interaction)
coef_interaction <- data.frame(sum_interaction$coefficients)
confint_interaction <- data.frame(confint(lm_interaction))
paste(round(exp(coef_interaction$Estimate[2]),3)," ", "(", round(exp(confint_interaction$X2.5..[2]),3),","," ", round(exp(confint_interaction$X97.5..[2]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[3]),3)," ", "(", round(exp(confint_interaction$X2.5..[3]),3),",", round(exp(confint_interaction$X97.5..[3]),3),")",sep = '')
paste(round(exp(coef_interaction$Estimate[nrow(coef_interaction)-1]),3)," ", "(", round(exp(confint_interaction$X2.5..[nrow(coef_interaction)-1]),3),",", round(exp(confint_interaction$X97.5..[nrow(coef_interaction)-1]),3),")",sep = '')

#################################################################Interaction of IHL and PRS for 9 TG-associated genes (from all 142 SNPs) at random
### select 9 genes at random for 500 times
set.seed(123)
TG_sample <- replicate(500, {
  sample(TG$SNP, 9, replace = FALSE)
})

# generate new variable
T2D <- NULL
for (i in 1:nrow(data)){
  if (data$N_Diabetes_WHO_2[i]==2){
    T2D <- c(T2D,1)
  }
  else{T2D <- c(T2D,0)
  }
} 
data["T2D"]<- T2D 
data$T2D=as.factor(data$T2D)

### repeat the main analysis for 1000 times with 1000 sets of TG-associated genes at random
IHL_TG1 <- NULL
IHL_TG2 <- NULL
IHL_TG3 <- NULL
PRS_TG1 <- NULL
PRS_TG2 <- NULL
PRS_TG3 <- NULL
interaction1 <- NULL
interaction2 <- NULL
interaction3 <- NULL
association <- NULL

for (i in 1:500) {
  TG_new <- data[,c(1,767:769,774:783,785:787,791,803,827:832,841:846,852:858,876,which(substring(names(data),4) %in% TG_sample[,i]))]
  TG_new <- subset(TG_new,is.na(TG_new$liverfat_adj)==FALSE) 
  TG_new <- subset(TG_new, TG_new$liverfat_adj != -999) 
  TG_new <- subset(TG_new,is.na(TG_new$batch)== FALSE) 
  TG_new <- TG_new[complete.cases(TG_new[,40:57]),] 
  TG_new <- subset(TG_new,TG_new$NIT_alcoholtot != -999) 
  TG_new <- subset(TG_new,TG_new$NIT_alcoholtot != -888) 
  
  for (a in seq(from = 40, to = 56, by = 2)) {
    for (b in 1:nrow(TG_new)) {
      if ((TG_new[b,a] == TG_new[b,a+1]) & (TG_new[b,a] == TG$`Effect allele`[which(TG$SNP == substring(colnames(TG_new)[a],4))])){
        TG_new[b,57+(a-38)/2] <- 2
      } else if ((TG_new[b,a] == TG_new[b,a+1]) & (TG_new[b,a] == TG$`Non-efect allele`[which(TG$SNP == substring(colnames(TG_new)[a],4))])){
        TG_new[b,57+(a-38)/2] <- 0
      } else {
        TG_new[b,57+(a-38)/2] <- 1
      }
    }
  }
  
  for (d in seq(from = 40, to = 56, by = 2)) {
    names(TG_new)[57+(d-38)/2] <-  substring(colnames(TG_new)[d],4)
  }
  
  for (c in 58:ncol(TG_new)){
    TG_new[,c+9] <- (TG$Beta[which(TG$SNP == colnames(TG_new)[c])]) * TG_new[,c] 
  }
  PRS <- rowSums(TG_new[,c(67:75)])
  
  # association between PRS and TG
  lm1 <- lm(formula = log10(Triglyc) ~ PRS + Age + SEX + T2D + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3
            +C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot,data = TG_new)
  sum1 <- summary(lm1)
  coef1 <- data.frame(sum1$coefficients)
  association <- c(association, coef1$Estimate[2])
  
  # interaction analysis
  #model 1
  lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS, data = TG_new)
  sum_interaction <- summary(lm_interaction)
  coef_interaction <- data.frame(sum_interaction$coefficients)
  IHL_TG1 <- c(IHL_TG1, coef_interaction$Estimate[2])
  PRS_TG1 <- c(PRS_TG1,coef_interaction$Estimate[3])
  interaction1 <- c(interaction1,coef_interaction$Estimate[nrow(coef_interaction)])
  
  #model 2
  lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + SEX + T2D + MRI_lagtime + liverfat_adj * MRI_lagtime +
                         +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                         C9_1KGP3 + C10_1KGP3, data = TG_new)
  sum_interaction <- summary(lm_interaction)
  coef_interaction <- data.frame(sum_interaction$coefficients)
  IHL_TG2 <- c(IHL_TG2, coef_interaction$Estimate[2])
  PRS_TG2 <- c(PRS_TG2,coef_interaction$Estimate[3])
  interaction2 <- c(interaction2,coef_interaction$Estimate[nrow(coef_interaction)-1])
  
  
  # model 3
  lm_interaction <- lm(formula = log10(Triglyc) ~ liverfat_adj + PRS + liverfat_adj * PRS + Age + SEX + T2D + MRI_lagtime + liverfat_adj * MRI_lagtime +
                         +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                         C9_1KGP3 + C10_1KGP3 + med_LP + NIT_alcoholtot, data = TG_new)
  sum_interaction <- summary(lm_interaction)
  coef_interaction <- data.frame(sum_interaction$coefficients)
  confint_interaction <- data.frame(confint(lm_interaction))
  IHL_TG3 <- c(IHL_TG3, coef_interaction$Estimate[2])
  PRS_TG3 <- c(PRS_TG3,coef_interaction$Estimate[3])
  interaction3 <- c(interaction3,coef_interaction$Estimate[nrow(coef_interaction)-1])

  print(i)
}

### get mean value of association and corresponding CI
n_bootstrap <- 5000  
bootstrap_means <- numeric(n_bootstrap)  

# Bootstrap sampling
set.seed(123)  
for (i in 1:n_bootstrap) {
  sample_data <- sample(association, size = length(association), replace = TRUE)
  bootstrap_means[i] <- mean(sample_data)
}

# mean association
mean_value <- mean(bootstrap_means)
exp(mean_value) # 1.211

# 95% CI
ci_lower <- quantile(bootstrap_means, 0.025)  
ci_upper <- quantile(bootstrap_means, 0.975)  
exp(ci_lower) # 1.195
exp(ci_upper) # 1.227


### get mean value of IHL_TG1 and corresponding CI
# Bootstrap sampling
set.seed(123)  
for (i in 1:n_bootstrap) {
  sample_data <- sample(IHL_TG1, size = length(IHL_TG1), replace = TRUE)
  bootstrap_means[i] <- mean(sample_data)
}

# mean association
mean_value <- mean(bootstrap_means)
exp(mean_value) 

# 95% CI
ci_lower <- quantile(bootstrap_means, 0.025)  
ci_upper <- quantile(bootstrap_means, 0.975)  
exp(ci_lower) 
exp(ci_upper)

### get mean value of PRS_TG1 and corresponding CI
# Bootstrap sampling
set.seed(123)  
for (i in 1:n_bootstrap) {
  sample_data <- sample(PRS_TG1, size = length(PRS_TG1), replace = TRUE)
  bootstrap_means[i] <- mean(sample_data)
}

# mean association
mean_value <- mean(bootstrap_means)
exp(mean_value) 

# 95% CI
ci_lower <- quantile(bootstrap_means, 0.025)  
ci_upper <- quantile(bootstrap_means, 0.975)  
exp(ci_lower) 
exp(ci_upper) 

### get mean value of the interaction term1 and corresponding CI
# Bootstrap sampling
set.seed(123)  
for (i in 1:n_bootstrap) {
  sample_data <- sample(interaction1, size = length(interaction1), replace = TRUE)
  bootstrap_means[i] <- mean(sample_data)
}

# mean association
mean_value <- mean(bootstrap_means)
exp(mean_value) 

# 95% CI
ci_lower <- quantile(bootstrap_means, 0.025)  
ci_upper <- quantile(bootstrap_means, 0.975)  
exp(ci_lower) 
exp(ci_upper)

# p value
alpha <- 0.05 

z_alpha_half <- qnorm(1 - alpha / 2)  
SE <- (ci_upper - ci_lower) / (2 * z_alpha_half)

# Z statistic
z_stat <- mean_value / SE

# p value
p_value <- 2 * (1 - pnorm(abs(z_stat)))  
p_value 

### get mean value of IHL_TG2 and corresponding CI
# Bootstrap sampling
set.seed(123)  
for (i in 1:n_bootstrap) {
  sample_data <- sample(IHL_TG2, size = length(IHL_TG2), replace = TRUE)
  bootstrap_means[i] <- mean(sample_data)
}

# mean association
mean_value <- mean(bootstrap_means)
exp(mean_value) 

# 95% CI
ci_lower <- quantile(bootstrap_means, 0.025)  
ci_upper <- quantile(bootstrap_means, 0.975)  
exp(ci_lower) 
exp(ci_upper)

### get mean value of PRS_TG2 and corresponding CI
# Bootstrap sampling
set.seed(123)  
for (i in 1:n_bootstrap) {
  sample_data <- sample(PRS_TG2, size = length(PRS_TG2), replace = TRUE)
  bootstrap_means[i] <- mean(sample_data)
}

# mean association
mean_value <- mean(bootstrap_means)
exp(mean_value) 

# 95% CI
ci_lower <- quantile(bootstrap_means, 0.025)  
ci_upper <- quantile(bootstrap_means, 0.975)  
exp(ci_lower) 
exp(ci_upper) 

### get mean value of the interaction term2 and corresponding CI
# Bootstrap sampling
set.seed(123)  
for (i in 1:n_bootstrap) {
  sample_data <- sample(interaction2, size = length(interaction2), replace = TRUE)
  bootstrap_means[i] <- mean(sample_data)
}

# mean association
mean_value <- mean(bootstrap_means)
exp(mean_value) 

# 95% CI
ci_lower <- quantile(bootstrap_means, 0.025)  
ci_upper <- quantile(bootstrap_means, 0.975)  
exp(ci_lower) 
exp(ci_upper)

# p value
alpha <- 0.05 

z_alpha_half <- qnorm(1 - alpha / 2)  
SE <- (ci_upper - ci_lower) / (2 * z_alpha_half)

# Z statistic
z_stat <- mean_value / SE

# p value
p_value <- 2 * (1 - pnorm(abs(z_stat)))  
p_value 

### get mean value of IHL_TG3 and corresponding CI
# Bootstrap sampling
set.seed(123)  
for (i in 1:n_bootstrap) {
  sample_data <- sample(IHL_TG3, size = length(IHL_TG3), replace = TRUE)
  bootstrap_means[i] <- mean(sample_data)
}

# mean association
mean_value <- mean(bootstrap_means)
exp(mean_value) 

# 95% CI
ci_lower <- quantile(bootstrap_means, 0.025)  
ci_upper <- quantile(bootstrap_means, 0.975)  
exp(ci_lower) 
exp(ci_upper)

### get mean value of PRS_TG3 and corresponding CI
# Bootstrap sampling
set.seed(123)  
for (i in 1:n_bootstrap) {
  sample_data <- sample(PRS_TG3, size = length(PRS_TG3), replace = TRUE)
  bootstrap_means[i] <- mean(sample_data)
}

# mean association
mean_value <- mean(bootstrap_means)
exp(mean_value) 

# 95% CI
ci_lower <- quantile(bootstrap_means, 0.025)  
ci_upper <- quantile(bootstrap_means, 0.975)  
exp(ci_lower) 
exp(ci_upper) 

### get mean value of the interaction term3 and corresponding CI
# Bootstrap sampling
set.seed(123)  
for (i in 1:n_bootstrap) {
  sample_data <- sample(interaction3, size = length(interaction3), replace = TRUE)
  bootstrap_means[i] <- mean(sample_data)
}

# mean association
mean_value <- mean(bootstrap_means)
exp(mean_value) 

# 95% CI
ci_lower <- quantile(bootstrap_means, 0.025)  
ci_upper <- quantile(bootstrap_means, 0.975)  
exp(ci_lower) 
exp(ci_upper)

# p value
alpha <- 0.05 

z_alpha_half <- qnorm(1 - alpha / 2)  
SE <- (ci_upper - ci_lower) / (2 * z_alpha_half)

# Z statistic
z_stat <- mean_value / SE

# p value
p_value <- 2 * (1 - pnorm(abs(z_stat)))  
p_value 


############################################################# Interaction of IHL and PRS for TG clearance on incident CVD
# exclude ineligible people
CVD <- merge(TG_clearance,data[,c(1,788,804:823,859:875)], by = 'randomID')
CVD <- subset(CVD, is.na(CVD$LD_CVD_event) == FALSE) # 52 with missing on  follow-up data
CVD <- subset(CVD,CVD$N_CVD == 0) # 509 with CVD at baseline
CVD_corrected <- subset(CVD,(CVD$LD_CVD_event==1 & (CVD$LD_CVD_Surv_time_standard_default> CVD$MRI_lagtime))|CVD$LD_CVD_event==0) # 30 with CVD happening before MRI measurement
CVD_corrected <- CVD_corrected[complete.cases(CVD_corrected[,c(22,27, 30:31,36)]),] # 2 with missing data on OSBP
table(CVD_corrected$LD_CVD_event) # 318/3217
CVD_corrected$LD_CVD_event = as.integer(CVD_corrected$LD_CVD_event)
CVD$LD_CVD_event = as.numeric(CVD$LD_CVD_event)

### IHL and incident CVD
res.cox <- coxph(Surv(LD_CVD_Surv_time_standard_default,LD_CVD_event) ~ liverfat_adj + 
                   Age + SEX + T2D + MRI_lagtime + liverfat_adj * MRI_lagtime + med_LP + NIT_alcoholtot + med_HT + smoking_3cat + OSBP + LDL, data = CVD_corrected)
summary(res.cox)

### PRS and incident CVD
res.cox <- coxph(Surv(LD_CVD_Surv_time_standard_default,LD_CVD_event) ~  PRS + 
                   Age + SEX + T2D + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                   C9_1KGP3 + C10_1KGP3 +  med_LP + NIT_alcoholtot + med_HT + smoking_3cat + OSBP + LDL, data = CVD_corrected)
summary(res.cox)

# person-year
colSums(data.frame(CVD_corrected$LD_CVD_Surv_time_standard_default[1:3217]))

### interaction analysis on incident CVD- continuous IHL and PRS 
# model 1
res.cox <- coxph(Surv(LD_CVD_Surv_time_standard_default,LD_CVD_event) ~ liverfat_adj + PRS + liverfat_adj * PRS, data = CVD_corrected)
summary(res.cox)

# model 2
res.cox <- coxph(Surv(LD_CVD_Surv_time_standard_default,LD_CVD_event) ~ liverfat_adj + PRS + liverfat_adj * PRS+
                   Age + SEX + T2D + MRI_lagtime + liverfat_adj * MRI_lagtime + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                   C9_1KGP3 + C10_1KGP3, data = CVD_corrected)
summary(res.cox)

# model 3
res.cox <- coxph(Surv(LD_CVD_Surv_time_standard_default,LD_CVD_event) ~ liverfat_adj + PRS + liverfat_adj * PRS+
                   Age + SEX + T2D + MRI_lagtime + liverfat_adj * MRI_lagtime + C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + 
                   C9_1KGP3 + C10_1KGP3 +  med_LP + NIT_alcoholtot + med_HT + smoking_3cat + OSBP + LDL, data = CVD_corrected)
summary(res.cox)


### interation analysis on incident CVD -steatosis and PRS intertile
res.cox <- coxph(Surv(LD_CVD_Surv_time_standard_default, LD_CVD_event) ~ group + Age + SEX + T2D + MRI_lagtime+ steotosis * MRI_lagtime +
                   +C1_1KGP3 + C2_1KGP3 + C3_1KGP3 + C4_1KGP3 + C5_1KGP3 + C6_1KGP3 +C7_1KGP3 + C8_1KGP3 + C9_1KGP3 + C10_1KGP3 +  
                   med_LP + NIT_alcoholtot + med_HT + smoking_3cat + OSBP + LDL, data = CVD_corrected)
summary(res.cox)
