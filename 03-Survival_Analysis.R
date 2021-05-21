#最佳評分閥值，區分高低分組
##
library(maxstat)
library(survminer)

#確定最佳切點
#StromalScore
TNBC_surv.cut<-surv_cutpoint(tmp3,
                             time = "times",
                             event = "patient.vital_status",
                             variables = c("StromalScore"))
summary(TNBC_surv.cut)
#-------------cutpoint statistic
#StromalScore 916.1294  2.165501
#繪製

plot(TNBC_surv.cut,"StromalScore",palette="npg")
#ImmuneScore
TNBC_surv2.cut<-surv_cutpoint(tmp3,
                             time = "times",
                             event = "patient.vital_status",
                             variables = c("ImmuneScore"))
summary(TNBC_surv2.cut)
#------------cutpoint statistic
#ImmuneScore 1642.413  1.367714

plot(TNBC_surv2.cut,"ImmuneScore",palette="npg")
#
TNBC_surv3.cut<-surv_cutpoint(tmp3,
                              time = "times",
                              event = "patient.vital_status",
                              variables = c("StromalScore","ImmuneScore"))
summary(TNBC_surv3.cut)

###分高低組
TNBC_surv.cat <- surv_categorize(TNBC_surv.cut)
TNBC_surv.cat<-cbind(tmp3$bcr_patient_barcode,TNBC_surv.cat,tmp3$StromalScore)
TNBC_surv2.cat <- surv_categorize(TNBC_surv2.cut)
TNBC_surv2.cat<-cbind(tmp3$bcr_patient_barcode,TNBC_surv2.cat,tmp3$ImmuneScore)
TNBC_surv3.cat <- surv_categorize(TNBC_surv3.cut)

write.csv(TNBC_surv.cat,file = "TNBC_surv.cat.csv")
write.csv(TNBC_surv2.cat,file = "TNBC_surv2.cat.csv")
table(TNBC_surv.cat$StromalScore)
#high low 
# 12  103
table(TNBC_surv2.cat$ImmuneScore)
#high low 
# 12  103

library(survival)
#StromalScore
fit <- survfit(Surv(times, patient.vital_status) ~ StromalScore,
               data = TNBC_surv.cat)

fit$time=(fit$time)/30
ggsurvplot(fit, 
           risk.table = TRUE,
           pval = TRUE,
           conf.int = F,
           # point estimaes of survival curves
           palette = "aaas",
           xlim = c(0,140), 
           # survival estimates
           break.time.by = 10,
           xlab="Months",
           legend = c(0.85,0.75),
           legend.labs = c("StromalScore=high", "StromalScore=low"),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE)
#ImmuneScore
fit2 <- survfit(Surv(times, patient.vital_status) ~ ImmuneScore,
               data = TNBC_surv2.cat)
fit2$time=(fit2$time)/30
ggsurvplot(fit2, 
           risk.table = TRUE,
           pval = TRUE,
           conf.int = F,
           # point estimaes of survival curves
           palette = "aaas",
           xlim = c(0,140), 
           # survival estimates
           break.time.by = 10,
           xlab="Months",
           legend = c(0.85,0.75),
           legend.labs = c("ImmuneScore=high", "ImmuneScore=low"),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE)
#race
fit3 <- survfit(Surv(times, patient.vital_status) ~ race,
                data = tmp3text)

ggsurvplot(fit3, 
           risk.table = TRUE,
           pval = TRUE,
           conf.int = TRUE,
           # point estimaes of survival curves
           xlim = c(0,4000), 
           # survival estimates
           break.time.by = 500,
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE )

#age
fit4 <- survfit(Surv(times, patient.vital_status) ~ age_status,
                data = tmp3text)
ggsurvplot(fit4, 
           risk.table = TRUE,
           pval = TRUE,
           conf.int = TRUE,
           # point estimaes of survival curves
           xlim = c(0,4000), 
           # survival estimates
           break.time.by = 500,
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE )

#tumor_stage
library(survival)
library(survminer)
table<-summary(fit5)$table
write.csv(table,file = "table.csv")
fit5 <- survfit(Surv(times, patient.vital_status) ~ tumor_stage,
                data = tmp3text)
fit5$time=(fit5$time)/30
ggsurvplot(fit5, 
           risk.table = TRUE,
           pval = TRUE,
           conf.int = F,
           # point estimaes of survival curves
           palette = "lancet",
           # survival estimates
           break.time.by = 10,
           legend.labs = c("Tumor_stage=not reported", "Stage I","Stage II","Stage III","Stage IV"),
           xlab="Months",
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE)
#StromalScore+ImmuneScore
fit11<- survfit(Surv(times, patient.vital_status) ~ StromalScore+ImmuneScore,
               data = TNBC_surv3.cat)
fit11$time=(fit11$time)/30
ggsurvplot(fit11, 
           risk.table = TRUE,
           pval = TRUE,
           conf.int = F,
           # point estimaes of survival curves
           palette = "aaas",
           # survival estimates
           break.time.by = 10,
           xlab="Months",
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE)