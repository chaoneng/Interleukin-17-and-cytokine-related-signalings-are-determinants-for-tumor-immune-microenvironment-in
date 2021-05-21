#cibersoft
install.packages("remotes")
library(remotes)
install.packages("matrixStats")
remotes::install_github("icbi-lab/immunedeconv")
library(immunedeconv)
source("CIBERSORT.R")
LM22.file <- "LM22.txt"
exp.file <- "cibersoftHL.txt"
TME.results = CIBERSORT(LM22.file, exp.file, perm = 1000, QN = TRUE)
write.csv(TME.results,file = "TME.results.csv")
library(ggplot2)
TME.results<-as.data.frame(TME.results)
exp.file2 <-"cibersoftLH.txt"
TME.results2 = CIBERSORT(LM22.file, exp.file2, perm = 1000, QN = TRUE)
write.csv(TME.results2,file = "TME.results2.csv")
#圖
library(dplyr)
library(tidyr)
library(tibble)
TME_box<-TME.results%>%
  as.data.frame()%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols=2:23,
               names_to="celltype",
               values_to="Proportion")
library(ggplot2)
ggplot(TME_box,aes(sample,Proportion,fill=celltype))+
  geom_bar(position = "stack",stat = "identity")+
  theme_bw()+
  guides(fill=guide_legend(ncol = 1))
TME_box2<-TME.results2%>%
  as.data.frame()%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols=2:23,
               names_to="celltype",
               values_to="Proportion")
library(ggplot2)
ggplot(TME_box2,aes(sample,Proportion,fill=celltype))+
  geom_bar(position = "stack",stat = "identity")+
  theme_bw()+
  guides(fill=guide_legend(ncol = 1))