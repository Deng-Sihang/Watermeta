rm(list=ls())
library(lme4)
library(car)
library(plyr)
library(vegan)
library(dplyr)
library(picante)
library(ieggr)
library(ggcor)
library(psych)
library(corrplot)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(reshape2)
library(tidyverse)
library(scales)
library(GGally)
library(cowplot)
library(betapart)
library(nlme)
library(emmeans)
library(tidyr) 
library(patchwork)
library(lmerTest)
library(metafor)
save.wd <- iwd(choose.dir())

###### Publication bias -------------------------------
dat <- read.table("Datasets-Github.csv",sep=",",header = TRUE)
dat$Replication_weighted <- (dat$N_C*dat$N_T)/(dat$N_C+dat$N_T)
dat <- dat%>%mutate(Site=paste(Latitude,Longitude,sep="_"))
#length(unique(dat$Site))

### Egger’s test + Rosenthal’s fail safe number
result.f <- c()
for (i in 1:13){
    dat1 <- data.frame(Site=dat$Site,Title1=dat$Title1,Title2=dat$Title2,Title3=dat$Title3,
                       Depth=dat$Depth,Type=dat$Type,N_C=dat$N_C,N_T=dat$N_T,
                       Replication_weighted=dat$Replication_weighted,
                       RR=dat[,(i*2+37)],Variance_weighted=dat[,(i*2+38)]) 
    dat1 <- na.omit(dat1)
    dat1 <- dat1[!duplicated(dat1[,10:11]),]
    if (i %in% 3:11){dat1 <- dat1[dat1$Depth!="C",]}
    
    for (type in unique(dat1$Type)){
      dat2 <- dat1[dat1$Type==type,]
      if ((nrow(dat2)>1)){
        result <- try(result <- regtest(RR, Variance_weighted, model = "rma", predictor = 'sei', data=dat2), silent = TRUE)
        if (inherits(result, "try-error")) {result <- regtest(RR, 1/Replication_weighted, model = "rma", predictor = 'sei', data=dat2)}
      
        z <- result$zval
        p <- result$pval
        n <- as.numeric(nrow(dat2)) 
        
        result <- fsn(RR, Variance_weighted, data = dat2)
        Nfs <- result$fsnum
        Type <- ifelse((p>0.05),"NO",
                           ifelse(((Nfs>(5*n+10))&(p<0.05)),"NO","YES"))
        
        output <- data.frame(index=colnames(dat)[(i*2+37)], type = type,n=n,Egger.z=z,Egger.p=p,
                             Nfs=Nfs,Type=Type)
        result.f <- rbind(result.f,output)
      }else{result.f <- result.f}
  }
}

write.csv(result.f, file = "1-Publication bias.csv")

###### LnRR calculation -------------------------------
library(lme4)
library(lmerTest)
library(multcomp)

result.f <- c()
### 1 ### Combined all microbial taxa
for (i in (1:13)){
  dat1 <- data.frame(Site=dat$Site,Title1=dat$Title1,Title2=dat$Title2,Title3=dat$Title3,
                     Depth=dat$Depth,Type=dat$Type,Magnitude=dat$Magnitude/100,
                     N_C=dat$N_C,N_T=dat$N_T,
                     Replication_weighted=dat$Replication_weighted,
                     RR=dat[,(i*2+37)],Variance_weighted=dat[,(i*2+38)]) 

  ### Response ratio
  dat1 <- na.omit(dat1)
  dat1 <- dat1[!duplicated(dat1[,(ncol(dat1)-1):ncol(dat1)]),]
  if (i %in% 3:11){dat1 <- dat1[dat1$Depth!="C",]}
  
  result <- c()
  for (Type in unique(dat1$Type)){
    dat2 <- dat1[(dat1$Type==Type),]
    test <- rma.mv(RR,Variance_weighted,data=dat2,
                   method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                   mods=~1,random=list(~1|Title1/Site))
    
    result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Response ratio",Micro="All",
                          Type=Type,Estimate=test$b[1],SE=test$se[1],
                         Z=test$zval,p=test$pval,CL=test$ci.lb,CU=test$ci.ub)
    
    result1$n <- nrow(dat2)
    result1$s <- length(unique(dat2$Title1))
    result <- rbind(result,result1)
  }
  
  test <- rma.mv(RR,Variance_weighted,data=dat1,
                 method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                 mods=~Type,random=list(~1|Title1/Site))
  result$QM <- QM <- test$QM
  result$QMp <- QMp <- test$QMp
  
  ### Sensitivity
  dat1$RR <- dat1$RR/abs(dat1$Magnitude)
  dat1$Variance_weighted <- dat1$Variance_weighted/(abs(dat1$Magnitude)*abs(dat1$Magnitude))
  
  results <- c()
  for (Type in unique(dat1$Type)){
    dat2 <- dat1[(dat1$Type==Type),]
    test <- rma.mv(RR,Variance_weighted,data=dat2,
                   method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                   mods=~1,random=list(~1|Title1/Site))
    
    result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Sensitivity",Micro="All",
                          Type=Type,Estimate=test$b[1],SE=test$se[1],
                          Z=test$zval,p=test$pval,CL=test$ci.lb,CU=test$ci.ub)
    
    result1$n <- nrow(dat2)
    result1$s <- length(unique(dat2$Title1))
    results <- rbind(results,result1)
  }
  
  test <- rma.mv(RR,Variance_weighted,data=dat1,
                 method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                 mods=~Type,random=list(~1|Title1/Site))
  results$QM <- QM <- test$QM
  results$QMp <- QMp <- test$QMp
  
  result <- rbind(result,results)
  ### Percent change
  result$percent <- (exp(result$Estimate) - 1) * 100  
  result$percent_lower <- (exp(result$CL) - 1) * 100 
  result$percent_upper <- (exp(result$CU) - 1) * 100  
  result.f <- rbind(result.f,result)
}

write.csv(result.f, file = "2-Total effect-Combined.csv")

result.f <- c()
### 2 ### Single microbial response
### For diversity (1:2), homogenization and differentiation (12:13)
for (i in (1:2)){
  dat1 <- data.frame(Site=dat$Site,Title1=dat$Title1,Title2=dat$Title2,Title3=dat$Title3,
                     Type=dat$Type,Microbical_type=dat$Microbical_types,Magnitude=dat$Magnitude/100,
                     N_C=dat$N_C,N_T=dat$N_T,
                     Replication_weighted=dat$Replication_weighted,
                     RR=dat[,(i*2+37)],Variance_weighted=dat[,(i*2+38)]) 
  
  ### Response ratio
  dat1 <- na.omit(dat1)
  dat1 <- dat1[!duplicated(dat1[,(ncol(dat1)-1):ncol(dat1)]),]

  for (Micro in unique(dat1$Microbical_type)){
    dat2 <- dat1[(dat1$Microbical_type==Micro),]
    
    result <- c()
    for (Type in unique(dat2$Type)){
      dat3 <- dat2[(dat2$Type==Type),]
      if (nrow(dat3) > 4){
        test <- rma.mv(RR,Variance_weighted,data=dat3,
                         method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                         mods=~1,random=list(~1|Title1/Site))
          
          result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Response ratio",Micro=Micro,
                                Type=Type,Estimate=test$b[1],SE=test$se[1],
                                Z=test$zval,p=test$pval,CL=test$ci.lb,CU=test$ci.ub)
          
          result1$n <- nrow(dat3)
          result1$s <- length(unique(dat3$Title1))
          result <- rbind(result,result1)
      }
    }
        
        ### Sensitivity
        dat2$RR <- dat2$RR/abs(dat2$Magnitude)
        dat2$Variance_weighted <- dat2$Variance_weighted/(abs(dat2$Magnitude)*abs(dat2$Magnitude))
        
        results <- c()
        for (Type in unique(dat2$Type)){
          dat3 <- dat2[(dat2$Type==Type),]
          if (nrow(dat3) > 4){
            test <- rma.mv(RR,Variance_weighted,data=dat3,
                           method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                           mods=~1,random=list(~1|Title1/Site))
            
            result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Sensitivity",Micro=Micro,
                                  Type=Type,Estimate=test$b[1],SE=test$se[1],
                                  Z=test$zval,p=test$pval,CL=test$ci.lb,CU=test$ci.ub)
            
            result1$n <- nrow(dat3)
            result1$s <- length(unique(dat3$Title1))
            results <- rbind(results,result1)
          }
        }
        
        result <- rbind(result,results)
        ### Percent change
        result$percent <- (exp(result$Estimate) - 1) * 100  
        result$percent_lower <- (exp(result$CL) - 1) * 100 
        result$percent_upper <- (exp(result$CU) - 1) * 100  
        result.f <- rbind(result.f,result)
  }
}

write.csv(result.f, file = "2-Total effect-Single.csv")

result.f <- c()
### 3 ### Omnibus test
### For diversity (1:2), homogenization and differentiation (12:13)
for (i in (1:2)){
  dat1 <- data.frame(Site=dat$Site,Title1=dat$Title1,Title2=dat$Title2,Title3=dat$Title3,
                     Type=dat$Type,Microbical_type=dat$Microbical_types,Magnitude=dat$Magnitude/100,
                     N_C=dat$N_C,N_T=dat$N_T,
                     Replication_weighted=dat$Replication_weighted,
                     RR=dat[,(i*2+37)],Variance_weighted=dat[,(i*2+38)]) 
  
  ### Response ratio
  dat1 <- na.omit(dat1)
  dat1 <- dat1[!duplicated(dat1[,(ncol(dat1)-1):ncol(dat1)]),]
  
  result <- c()
  ### Calculate P-value for between-group heterogeneity (across treatment)
  test <- rma.mv(RR,Variance_weighted,data=dat1,
                 method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                 mods=~Type,random=list(~1|Title1/Site))
  QM <- test$QM
  QMp <- test$QMp
  
  result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Response ratio",Treat="All",
                        Micro="All",Type="across treatment",QM=QM,QMp=QMp)
  result <- rbind(result,result1) 
  
  dat1.asymmetry <- dat1
  dat1.asymmetry[dat1.asymmetry$Type=="P-",]$RR <- dat1.asymmetry[dat1.asymmetry$Type=="P-",]$RR*(-1)
  ### Calculate P-value for between-group heterogeneity (asymmetry index)
  test <- rma.mv(RR,Variance_weighted,data=dat1.asymmetry,
                 method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                 mods=~Type,random=list(~1|Title1/Site))
  QM <- test$QM
  QMp <- test$QMp
  
  result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Response ratio",Treat="All",
                        Micro="All",Type="asymmetry index",QM=QM,QMp=QMp)
  result <- rbind(result,result1) 
  
  ### Different microbe
  for (Micro in unique(dat1$Microbical_type)){
    dat2 <- dat1[(dat1$Microbical_type==Micro),]
    ### Calculate P-value for between-group heterogeneity (across treatment)
    test <- rma.mv(RR,Variance_weighted,data=dat2,
                   method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                   mods=~Type,random=list(~1|Title1/Site))
    QM <- test$QM
    QMp <- test$QMp
    
    result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Response ratio",Treat="All",
                          Micro=Micro,Type="across treatment",QM=QM,QMp=QMp)
    result <- rbind(result,result1) 
    
    dat2.asymmetry <- dat2
    dat2.asymmetry[dat2.asymmetry$Type=="P-",]$RR <- dat2.asymmetry[dat2.asymmetry$Type=="P-",]$RR*(-1)
    ### Calculate P-value for between-group heterogeneity (asymmetry index)
    test <- rma.mv(RR,Variance_weighted,data=dat2.asymmetry,
                   method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                   mods=~Type,random=list(~1|Title1/Site))
    QM <- test$QM
    QMp <- test$QMp
    
    result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Response ratio",Treat="All",
                          Micro=Micro,Type="asymmetry index",QM=QM,QMp=QMp)
    result <- rbind(result,result1) 
  }
  
  ### Different treatments
  for (Type in unique(dat1$Type)){
    dat2 <- dat1[(dat1$Type==Type),]
    ### Calculate P-value for between-group heterogeneity (across microbe)
    test <- rma.mv(RR,Variance_weighted,data=dat2,
                   method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                   mods=~Microbical_type,random=list(~1|Title1/Site))
    QM <- test$QM
    QMp <- test$QMp
    
    result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Response ratio",Treat=Type,
                          Micro="All",Type="across microbe",QM=QM,QMp=QMp)
    result <- rbind(result,result1) 
  }
    
  ### Sensitivity
  dat1$RR <- dat1$RR/abs(dat1$Magnitude)
  dat1$Variance_weighted <- dat1$Variance_weighted/(abs(dat1$Magnitude)*abs(dat1$Magnitude))
    
  ### Calculate P-value for between-group heterogeneity (across treatment)
  test <- rma.mv(RR,Variance_weighted,data=dat1,
                 method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                 mods=~Type,random=list(~1|Title1/Site))
  QM <- test$QM
  QMp <- test$QMp
  
  result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Sensitivity",Treat="All",
                        Micro="All",Type="across treatment",QM=QM,QMp=QMp)
  result <- rbind(result,result1) 
  
  dat1.asymmetry <- dat1
  dat1.asymmetry[dat1.asymmetry$Type=="P-",]$RR <- dat1.asymmetry[dat1.asymmetry$Type=="P-",]$RR*(-1)
  ### Calculate P-value for between-group heterogeneity (asymmetry index)
  test <- rma.mv(RR,Variance_weighted,data=dat1.asymmetry,
                 method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                 mods=~Type,random=list(~1|Title1/Site))
  QM <- test$QM
  QMp <- test$QMp
  
  result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Sensitivity",Treat="All",
                        Micro="All",Type="asymmetry index",QM=QM,QMp=QMp)
  result <- rbind(result,result1) 
  
  ### Different microbe
  for (Micro in unique(dat1$Microbical_type)){
    dat2 <- dat1[(dat1$Microbical_type==Micro),]
    ### Calculate P-value for between-group heterogeneity (across treatment)
    test <- rma.mv(RR,Variance_weighted,data=dat2,
                   method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                   mods=~Type,random=list(~1|Title1/Site))
    QM <- test$QM
    QMp <- test$QMp
    
    result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Sensitivity",Treat="All",
                          Micro=Micro,Type="across treatment",QM=QM,QMp=QMp)
    result <- rbind(result,result1) 
    
    dat2.asymmetry <- dat2
    dat2.asymmetry[dat2.asymmetry$Type=="P-",]$RR <- dat2.asymmetry[dat2.asymmetry$Type=="P-",]$RR*(-1)
    ### Calculate P-value for between-group heterogeneity (asymmetry index)
    test <- rma.mv(RR,Variance_weighted,data=dat2.asymmetry,
                   method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                   mods=~Type,random=list(~1|Title1/Site))
    QM <- test$QM
    QMp <- test$QMp
    
    result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Sensitivity",Treat="All",
                          Micro=Micro,Type="asymmetry index",QM=QM,QMp=QMp)
    result <- rbind(result,result1) 
  }
  
  ### Different treatments
  for (Type in unique(dat1$Type)){
    dat2 <- dat1[(dat1$Type==Type),]
    ### Calculate P-value for between-group heterogeneity (across microbe)
    test <- rma.mv(RR,Variance_weighted,data=dat2,
                   method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                   mods=~Microbical_type,random=list(~1|Title1/Site))
    QM <- test$QM
    QMp <- test$QMp
    
    result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Category="Sensitivity",Treat=Type,
                          Micro="All",Type="across microbe",QM=QM,QMp=QMp)
    result <- rbind(result,result1) 
  }
  result.f <- rbind(result.f,result)
}

write.csv(result.f, file = "2-Total effect-QM.csv")

###### Subgroup analysis: Ecosystem type -------------------------------
result.f <- c()
### For diversity (1:2), homogenization and differentiation (12:13)
for (i in (1:2)){
  dat1 <- data.frame(Site=dat$Site,Title1=dat$Title1,Title2=dat$Title2,Title3=dat$Title3,
                     Category=dat$Biomes,Microbical_type=dat$Microbical_types,Type=dat$Type,
                     Magnitude=dat$Magnitude/100,N_C=dat$N_C,N_T=dat$N_T,
                     Replication_weighted=dat$Replication_weighted,
                     RR=dat[,(i*2+37)],Variance_weighted=dat[,(i*2+38)]) 
  
  ### Response ratio
  dat1 <- na.omit(dat1)
  dat1 <- dat1[!duplicated(dat1[,(ncol(dat1)-1):ncol(dat1)]),]
  
  ### For all microbe
  for (Type in unique(dat1$Type)){
    dat2 <- dat1[(dat1$Type==Type),]
    dat2 <- dat2 %>% group_by(Category) %>% filter(n() >= 3) %>% ungroup()
    dat2 <- dat2 %>% group_by(Category) %>% filter(n_distinct(paste(Title1, Site, sep = "|")) > 1) %>% ungroup()
    
    if (length(unique(dat2$Category)) >= 2){
      result <- c()
      for (class in unique(dat2$Category)){
        dat3 <- dat2[(dat2$Category==class),]
        test <- rma.mv(RR,Variance_weighted,data=dat3,
                       method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                       mods=~1,random=list(~1|Title1/Site))
        
        result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Index="Response ratio",Micro="All",
                              Type=Type,Category=class,Estimate=test$b[1],SE=test$se[1],
                              Z=test$zval,p=test$pval,CL=test$ci.lb,CU=test$ci.ub)
        
        result1$n <- nrow(dat3)
        result1$s <- length(unique(dat3$Title1))
        
        dat.compare <- dat1[dat1$Category==class,]
        if(nrow(dat.compare[(dat.compare$Type!=Type),])>=3){
          test <- rma.mv(RR,Variance_weighted,data=dat.compare,
                         method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                         mods=~Type,random=list(~1|Title1/Site))
          result1$QM.treat <- QM <- test$QM
          result1$QMp.treat <- QMp <- test$QMp  
        }else{result1$QM.treat <- result1$QMp.treat <- NA}
        result <- rbind(result,result1)
      }
      
      ### Calculate P-value for between-group heterogeneity
      test <- rma.mv(RR,Variance_weighted,data=dat2,
                     method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                     mods=~Category,random=list(~1|Title1/Site))
      result$QM.type <- QM <- test$QM
      result$QMp.type <- QMp <- test$QMp    
      result.f <- rbind(result.f,result)
    }else {next}
  }
  
  ### For different microbial taxa
  for (Micro in unique(dat1$Microbical_type)){
    dat2 <- dat1[(dat1$Microbical_type==Micro),]
    for (Type in unique(dat2$Type)){
      dat3 <- dat2[(dat2$Type==Type),]
      dat3 <- dat3 %>% group_by(Category) %>% filter(n() >= 3) %>% ungroup()
      dat3 <- dat3 %>% group_by(Category) %>% filter(n_distinct(paste(Title1, Site, sep = "|")) > 1) %>% ungroup()
      
      if (length(unique(dat3$Category)) >= 2){
        result <- c()
        for (class in unique(dat3$Category)){
          dat4 <- dat3[(dat3$Category==class),]
          test <- rma.mv(RR,Variance_weighted,data=dat4,
                         method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                         mods=~1,random=list(~1|Title1/Site))
          
          result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Index="Response ratio",Micro=Micro,
                                Type=Type,Category=class,Estimate=test$b[1],SE=test$se[1],
                                Z=test$zval,p=test$pval,CL=test$ci.lb,CU=test$ci.ub)
          
          result1$n <- nrow(dat4)
          result1$s <- length(unique(dat4$Title1))
          
          dat.compare <- dat2[dat2$Category==class,]
          if((nrow(dat.compare[(dat.compare$Type!=Type),])>=3)&
             (length(unique(paste0(dat.compare[(dat.compare$Type!=Type),]$Title1,dat.compare[(dat.compare$Type!=Type),]$Site)))>1)){
            test <- rma.mv(RR,Variance_weighted,data=dat.compare,
                           method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                           mods=~Type,random=list(~1|Title1/Site))
            result1$QM.treat <- QM <- test$QM
            result1$QMp.treat <- QMp <- test$QMp  
          }else{result1$QM.treat <- result1$QMp.treat <- NA}
          result <- rbind(result,result1)
        }
        
        ### Calculate P-value for between-group heterogeneity
        test <- rma.mv(RR,Variance_weighted,data=dat3,
                       method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                       mods=~Category,random=list(~1|Title1/Site))
        result$QM.type <- QM <- test$QM
        result$QMp.type <- QMp <- test$QMp   
        result.f <- rbind(result.f,result)
      }else {next}
    }
  }
  
  ### Sensitivity
  dat1$RR <- dat1$RR/abs(dat1$Magnitude)
  dat1$Variance_weighted <- dat1$Variance_weighted/(abs(dat1$Magnitude)*abs(dat1$Magnitude))
  
  ### For all microbe
  for (Type in unique(dat1$Type)){
    dat2 <- dat1[(dat1$Type==Type),]
    dat2 <- dat2 %>% group_by(Category) %>% filter(n() >= 3) %>% ungroup()
    dat2 <- dat2 %>% group_by(Category) %>% filter(n_distinct(paste(Title1, Site, sep = "|")) > 1) %>% ungroup()
    
    if (length(unique(dat2$Category)) >= 2){
      result <- c()
      for (class in unique(dat2$Category)){
        dat3 <- dat2[(dat2$Category==class),]
        test <- rma.mv(RR,Variance_weighted,data=dat3,
                       method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                       mods=~1,random=list(~1|Title1/Site))
        
        result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Index="Sensitivity",Micro="All",
                              Type=Type,Category=class,Estimate=test$b[1],SE=test$se[1],
                              Z=test$zval,p=test$pval,CL=test$ci.lb,CU=test$ci.ub)
        
        result1$n <- nrow(dat3)
        result1$s <- length(unique(dat3$Title1))
        
        dat.compare <- dat1[dat1$Category==class,]
        if(nrow(dat.compare[(dat.compare$Type!=Type),])>=3){
          test <- rma.mv(RR,Variance_weighted,data=dat.compare,
                         method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                         mods=~Type,random=list(~1|Title1/Site))
          result1$QM.treat <- QM <- test$QM
          result1$QMp.treat <- QMp <- test$QMp  
        }else{result1$QM.treat <- result1$QMp.treat <- NA}
        result <- rbind(result,result1)
      }
      
      ### Calculate P-value for between-group heterogeneity
      test <- rma.mv(RR,Variance_weighted,data=dat2,
                     method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                     mods=~Category,random=list(~1|Title1/Site))
      result$QM.type <- QM <- test$QM
      result$QMp.type <- QMp <- test$QMp    
      result.f <- rbind(result.f,result)
    }else {next}
  }
  
  ### For different microbial taxa
  for (Micro in unique(dat1$Microbical_type)){
    dat2 <- dat1[(dat1$Microbical_type==Micro),]
    for (Type in unique(dat2$Type)){
      dat3 <- dat2[(dat2$Type==Type),]
      dat3 <- dat3 %>% group_by(Category) %>% filter(n() >= 3) %>% ungroup()
      dat3 <- dat3 %>% group_by(Category) %>% filter(n_distinct(paste(Title1, Site, sep = "|")) > 1) %>% ungroup()
      
      if (length(unique(dat3$Category)) >= 2){
        result <- c()
        for (class in unique(dat3$Category)){
          dat4 <- dat3[(dat3$Category==class),]
          test <- rma.mv(RR,Variance_weighted,data=dat4,
                         method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                         mods=~1,random=list(~1|Title1/Site))
          
          result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Index="Sensitivity",Micro=Micro,
                                Type=Type,Category=class,Estimate=test$b[1],SE=test$se[1],
                                Z=test$zval,p=test$pval,CL=test$ci.lb,CU=test$ci.ub)
          
          result1$n <- nrow(dat4)
          result1$s <- length(unique(dat4$Title1))
          
          dat.compare <- dat2[dat2$Category==class,]
          if((nrow(dat.compare[(dat.compare$Type!=Type),])>=3)&
             (length(unique(paste0(dat.compare[(dat.compare$Type!=Type),]$Title1,dat.compare[(dat.compare$Type!=Type),]$Site)))>1)){
            test <- rma.mv(RR,Variance_weighted,data=dat.compare,
                           method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                           mods=~Type,random=list(~1|Title1/Site))
            result1$QM.treat <- QM <- test$QM
            result1$QMp.treat <- QMp <- test$QMp  
          }else{result1$QM.treat <- result1$QMp.treat <- NA}
          result <- rbind(result,result1)
        }
        
        ### Calculate P-value for between-group heterogeneity
        test <- rma.mv(RR,Variance_weighted,data=dat3,
                       method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                       mods=~Category,random=list(~1|Title1/Site))
        result$QM.type <- QM <- test$QM
        result$QMp.type <- QMp <- test$QMp   
        result.f <- rbind(result.f,result)
      }else {next}
    }
  }
}

### Percent change
result.f$percent <- (exp(result.f$Estimate) - 1) * 100  
result.f$percent_lower <- (exp(result.f$CL) - 1) * 100 
result.f$percent_upper <- (exp(result.f$CU) - 1) * 100  
write.csv(result.f, file = "3-Subgroup analysis-Ecosystem type.csv")

###### Relative importance of factors -------------------------------
library(rJava)
library(glmulti)
### Depth: 0-30 cm
dat.test <- dat[dat$Depth!="C",]

result.f <- c()
### For diversity (1:2), homogenization and differentiation (12:13)
for (i in (1:2)){
  dat1 <- data.frame(Site=dat.test$Site,Title1=dat.test$Title1,Microbical_type=dat.test$Microbical_types,
                     Type=dat.test$Type,N_C=dat.test$N_C,N_T=dat.test$N_T,
                     Replication_weighted=dat.test$Replication_weighted,
                     RR=dat.test[,(i*2+37)],Variance_weighted=dat.test[,(i*2+38)]) 
  
  ### Factors
  env <- dat.test[,c("MAT","MAP","Clay","Sand","SOC",
                     "TN","CN","pH","Duration","Magnitude")]
  dat2 <- cbind(dat1,env)
  dat2 <- na.omit(dat2)
  dat2 <- dat2[!duplicated(dat2[,8:9]),]
  
  ### All microbe
  for (Type in unique(dat2$Type)){
    dat3 <- dat2[(dat2$Type==Type),]
    model <-  lm(RR ~ MAT+MAP+Clay+Sand+SOC+TN+CN+pH+Duration+Magnitude, data=dat3)
    
    result <- glmulti(model,level=1,crit="aicc")
    result.output <- as.data.frame(coef(result))
    result.output <- data.frame(Estimate=result.output$Estimate, SE=sqrt(result.output$`Uncond. variance`), 
                                Importance=result.output$Importance, row.names=row.names(result.output))
    result.output$z <- result.output$Estimate / result.output$SE
    result.output$p <- 2*pnorm(abs(result.output$z), lower.tail=FALSE)
    names(result.output) <- c("Estimate", "SE", "Importance", "Z", "P")
    result.output$ci.lb <- result.output[[1]] - qnorm(.975) * result.output[[2]]
    result.output$ci.ub <- result.output[[1]] + qnorm(.975) * result.output[[2]]
    result.output <- result.output[order(result.output$Importance, decreasing=TRUE), c(1,2,4:7,3)]
    
    best_model <- result@objects[[1]]
    best_coefficients <- as.data.frame(coef(best_model))
    
    result.output$best_model <- rownames(result.output) %in% rownames(best_coefficients)
    result.output <- data.frame(Factor = colnames(dat.test)[(i*2+37)],Predictor =rownames(result.output),micro = "All",type = Type,result.output)
    rownames(result.output) <- NULL
    result.f <- rbind(result.f,result.output)
  }
  
  ### Bacteria + Fungi
  dat2 <- dat2[dat2$Microbical_type!="Protists",]
  for (Micro in unique(dat2$Microbical_type)){
    dat3 <- dat2[(dat2$Microbical_type==Micro),]
    for (Type in unique(dat3$Type)){
      dat4 <- dat3[(dat3$Type==Type),]
      model <-  lm(RR ~ MAT+MAP+Clay+Sand+SOC+TN+CN+pH+Duration+Magnitude, data=dat4)
      
      result <- glmulti(model,level=1,crit="aicc")
      result.output <- as.data.frame(coef(result))
      result.output <- data.frame(Estimate=result.output$Estimate, SE=sqrt(result.output$`Uncond. variance`), 
                                  Importance=result.output$Importance, row.names=row.names(result.output))
      result.output$z <- result.output$Estimate / result.output$SE
      result.output$p <- 2*pnorm(abs(result.output$z), lower.tail=FALSE)
      names(result.output) <- c("Estimate", "SE", "Importance", "Z", "P")
      result.output$ci.lb <- result.output[[1]] - qnorm(.975) * result.output[[2]]
      result.output$ci.ub <- result.output[[1]] + qnorm(.975) * result.output[[2]]
      result.output <- result.output[order(result.output$Importance, decreasing=TRUE), c(1,2,4:7,3)]
      
      best_model <- result@objects[[1]]
      best_coefficients <- as.data.frame(coef(best_model))
      
      result.output$best_model <- rownames(result.output) %in% rownames(best_coefficients)
      result.output <- data.frame(Factor = colnames(dat.test)[(i*2+37)],Predictor =rownames(result.output),micro = Micro,type = Type,result.output)
      rownames(result.output) <- NULL
      result.f <- rbind(result.f,result.output)
    }
  }
}

write.csv(result.f, file = "4-Factor importance-AIC.csv")

###### Meta-regression -------------------------------
library(MuMIn)
### Depth: 0-30 cm
dat.test <- dat[dat$Depth!="C",]

result.f <- c()
### 1 ### Model selection
### For diversity (1:2), homogenization and differentiation (12:13)
for (i in (1:2)){
  dat1 <- data.frame(Site=dat.test$Site,Title=dat.test$Title1,
                     MAT=dat.test$MAT,Magnitude=dat.test$Magnitude,Duration=dat.test$Duration,
                     Microbical_type=dat.test$Microbical_types,Type=dat.test$Type,
                     N_C=dat.test$N_C,N_T=dat.test$N_T,Replication_weighted=dat.test$Replication_weighted,
                     RR=dat.test[,(i*2+37)],Variance_weighted=dat.test[,(i*2+38)]) 
  dat1$Variance_weighted_1 <- 1/dat1$Variance_weighted
  
  ### Response ratio
  dat1 <- na.omit(dat1)
  dat1 <- dat1[!duplicated(dat1[,(ncol(dat1)-1):ncol(dat1)]),]
  dat1 <- dat1[dat1$Microbical_type!="Protists",]
  
  for (Micro in unique(dat1$Microbical_type)){
    dat2 <- dat1[dat1$Microbical_type==Micro,]
    for (Type in unique(dat2$Type)){
      dat3 <- dat2[dat2$Type==Type,]
      
      test1 <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(Duration)+(1|Title/Site),
                   weights=Variance_weighted_1,dat3)
      test2 <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(Duration)+scale(Magnitude)*scale(Duration)+(1|Title/Site),
                    weights=Variance_weighted_1,dat3)
      test3 <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(Duration)+scale(Magnitude)*scale(Duration)*scale(MAT)+(1|Title/Site),
                    weights=Variance_weighted_1,dat3)
      test4 <- lmer(RR~scale(MAT)+scale(log(Magnitude))+scale(Duration)+(1|Title/Site),
                    weights=Variance_weighted_1,dat3)
      test5 <- lmer(RR~scale(MAT)+scale(log(Magnitude))+scale(Duration)+scale(log(Magnitude))*scale(Duration)+(1|Title/Site),
                    weights=Variance_weighted_1,dat3)
      test6 <- lmer(RR~scale(MAT)+scale(log(Magnitude))+scale(Duration)+scale(log(Magnitude))*scale(Duration)*scale(MAT)+(1|Title/Site),
                    weights=Variance_weighted_1,dat3)
      test7 <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(log(Duration))+(1|Title/Site),
                    weights=Variance_weighted_1,dat3)
      test8 <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(log(Duration))+scale(Magnitude)*scale(log(Duration))+(1|Title/Site),
                    weights=Variance_weighted_1,dat3)
      test9 <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(log(Duration))+scale(Magnitude)*scale(log(Duration))*scale(MAT)+(1|Title/Site),
                    weights=Variance_weighted_1,dat3)
      test10 <- lmer(RR~scale(MAT)+scale(log(Magnitude))+scale(log(Duration))+(1|Title/Site),
                    weights=Variance_weighted_1,dat3)
      test11 <- lmer(RR~scale(MAT)+scale(log(Magnitude))+scale(log(Duration))+scale(log(Magnitude))*scale(log(Duration))+(1|Title/Site),
                    weights=Variance_weighted_1,dat3)
      test12 <- lmer(RR~scale(MAT)+scale(log(Magnitude))+scale(log(Duration))+scale(log(Magnitude))*scale(log(Duration))*scale(MAT)+(1|Title/Site),
                    weights=Variance_weighted_1,dat3)
      
      models <- list(test1 = test1,test2 = test2,test3 = test3,test4 = test4,
                     test5 = test5,test6 = test6,test7 = test7,test8 = test8,
                     test9 = test9,test10 = test10,test11 = test11,test12 = test12)
      
      model_comparison <- data.frame(Model = names(models),
        AIC = sapply(models, AIC),BIC = sapply(models, BIC),logLik = sapply(models, logLik))
      model_comparison$Marginal_R2 <- sapply(models, function(m) r.squaredGLMM(m)[1])
      model_comparison$Conditional_R2 <- sapply(models, function(m) r.squaredGLMM(m)[2])
      model_comparison <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                                     model_comparison)
      rownames(model_comparison) <- NULL
      model_comparison <- model_comparison %>%
        mutate(Best_model = if_else(AIC == min(AIC, na.rm = TRUE), "YES", "NO"))
      result.f <- rbind(result.f,model_comparison)
    }
  }
}  

write.csv(result.f, file = "5-Meta regression-AIC.csv")

### 2 ### Meta-regression - without env
result.f <- c()
### For diversity (1:2), homogenization and differentiation (12:13)
for (i in (1:2)){
  dat1 <- data.frame(Site=dat.test$Site,Title=dat.test$Title1,
                     MAT=dat.test$MAT,Magnitude=dat.test$Magnitude,Duration=dat.test$Duration,
                     Microbical_type=dat.test$Microbical_types,Type=dat.test$Type,
                     N_C=dat.test$N_C,N_T=dat.test$N_T,Replication_weighted=dat.test$Replication_weighted,
                     RR=dat.test[,(i*2+37)],Variance_weighted=dat.test[,(i*2+38)]) 
  dat1$Variance_weighted_1 <- 1/dat1$Variance_weighted
  
  ### Response ratio
  dat1 <- na.omit(dat1)
  dat1 <- dat1[!duplicated(dat1[,(ncol(dat1)-2):ncol(dat1)]),]
  dat1 <- dat1[dat1$Microbical_type!="Protists",]
  
  for (Micro in unique(dat1$Microbical_type)){
    dat2 <- dat1[dat1$Microbical_type==Micro,]
    for (Type in unique(dat2$Type)){
      dat3 <- dat2[dat2$Type==Type,]
      
      if ((colnames(dat.test)[(i*2+37)]=="Richness_RR")){
        if ((Micro=="Fungi")&(Type=="P-")){
         test <- lmer(RR~scale(MAT)+scale(log(Magnitude))+scale(Duration)+(1|Title/Site),
                    weights=Variance_weighted_1,dat3)
         test1 <- anova(test)
         
         result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                              Variables=rownames(coef(summary(test)))[2:4],
                              Estimate=as.numeric(coef(summary(test))[2:4,1]),SE=as.numeric(coef(summary(test))[2:4,2]),
                              F.value=as.numeric(test1$`F value`),p=as.numeric(test1$`Pr(>F)`),
                              n=nrow(dat3),s=length(unique(dat3$Title)))
         result.f <- rbind(result.f,result)
      }else{
        test <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(Duration)+(1|Title/Site),
                     weights=Variance_weighted_1,dat3)
        test1 <- anova(test) 
        
        result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                             Variables=rownames(coef(summary(test)))[2:4],
                             Estimate=as.numeric(coef(summary(test))[2:4,1]),SE=as.numeric(coef(summary(test))[2:4,2]),
                             F.value=as.numeric(test1$`F value`),p=as.numeric(test1$`Pr(>F)`),
                             n=nrow(dat3),s=length(unique(dat3$Title)))
        result.f <- rbind(result.f,result)
      }
      }else if((colnames(dat.test)[(i*2+37)]=="Shannon_RR")){
        if ((Micro=="Bacteria")&(Type=="P-")){
          test <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(log(Duration))+(1|Title/Site),
                       weights=Variance_weighted_1,dat3)
          test1 <- anova(test)
          
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=rownames(coef(summary(test)))[2:4],
                               Estimate=as.numeric(coef(summary(test))[2:4,1]),SE=as.numeric(coef(summary(test))[2:4,2]),
                               F.value=as.numeric(test1$`F value`),p=as.numeric(test1$`Pr(>F)`),
                               n=nrow(dat3),s=length(unique(dat3$Title)))
          result.f <- rbind(result.f,result)
        }else{
          test <- lmer(RR~scale(MAT)+scale(log(Magnitude))+scale(Duration)+(1|Title/Site),
                       weights=Variance_weighted_1,dat3)
          test1 <- anova(test)
          
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=rownames(coef(summary(test)))[2:4],
                               Estimate=as.numeric(coef(summary(test))[2:4,1]),SE=as.numeric(coef(summary(test))[2:4,2]),
                               F.value=as.numeric(test1$`F value`),p=as.numeric(test1$`Pr(>F)`),
                               n=nrow(dat3),s=length(unique(dat3$Title)))
          result.f <- rbind(result.f,result)
        }
      }else if((colnames(dat.test)[(i*2+37)]=="Beta_RR")){
        if ((Micro=="Bacteria")&(Type=="P-")){
          test <- lmer(RR~scale(MAT)+scale(log(Magnitude))+scale(Duration)+(1|Title/Site),
                       weights=Variance_weighted_1,dat3)
          test1 <- anova(test)
          
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=rownames(coef(summary(test)))[2:4],
                               Estimate=as.numeric(coef(summary(test))[2:4,1]),SE=as.numeric(coef(summary(test))[2:4,2]),
                               F.value=as.numeric(test1$`F value`),p=as.numeric(test1$`Pr(>F)`),
                               n=nrow(dat3),s=length(unique(dat3$Title)))
          result.f <- rbind(result.f,result)
        }else if ((Micro=="Fungi")&(Type=="P-")){
          test <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(log(Duration))+(1|Title/Site),
                       weights=Variance_weighted_1,dat3)
          test1 <- anova(test)
          
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=rownames(coef(summary(test)))[2:4],
                               Estimate=as.numeric(coef(summary(test))[2:4,1]),SE=as.numeric(coef(summary(test))[2:4,2]),
                               F.value=as.numeric(test1$`F value`),p=as.numeric(test1$`Pr(>F)`),
                               n=nrow(dat3),s=length(unique(dat3$Title)))
          result.f <- rbind(result.f,result)
        }else{
          test <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(Duration)+(1|Title/Site),
                       weights=Variance_weighted_1,dat3)
          test1 <- anova(test)
          
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=rownames(coef(summary(test)))[2:4],
                               Estimate=as.numeric(coef(summary(test))[2:4,1]),SE=as.numeric(coef(summary(test))[2:4,2]),
                               F.value=as.numeric(test1$`F value`),p=as.numeric(test1$`Pr(>F)`),
                               n=nrow(dat3),s=length(unique(dat3$Title)))
          result.f <- rbind(result.f,result)
        }
      }else if((colnames(dat.test)[(i*2+37)]=="Comp_RR")){
        if ((Micro=="Bacteria")&(Type=="P-")){
          test <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(log(Duration))+scale(Magnitude)*scale(log(Duration))+(1|Title/Site),
                       weights=Variance_weighted_1,dat3)
          test1 <- anova(test)
          
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=rownames(coef(summary(test)))[2:4],
                               Estimate=as.numeric(coef(summary(test))[2:4,1]),SE=as.numeric(coef(summary(test))[2:4,2]),
                               F.value=as.numeric(test1$`F value`[1:3]),p=as.numeric(test1$`Pr(>F)`[1:3]),
                               n=nrow(dat3),s=length(unique(dat3$Title)))
          result.f <- rbind(result.f,result)
        }else if ((Micro=="Bacteria")&(Type=="P+")){
          test <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(log(Duration))+(1|Title/Site),
                       weights=Variance_weighted_1,dat3)
          test1 <- anova(test)
          
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=rownames(coef(summary(test)))[2:4],
                               Estimate=as.numeric(coef(summary(test))[2:4,1]),SE=as.numeric(coef(summary(test))[2:4,2]),
                               F.value=as.numeric(test1$`F value`),p=as.numeric(test1$`Pr(>F)`),
                               n=nrow(dat3),s=length(unique(dat3$Title)))
          result.f <- rbind(result.f,result)
        }else if ((Micro=="Fungi")&(Type=="P-")){
          test <- lmer(RR~scale(MAT)+scale(Magnitude)+scale(Duration)+(1|Title/Site),
                       weights=Variance_weighted_1,dat3)
          test1 <- anova(test)
          
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=rownames(coef(summary(test)))[2:4],
                               Estimate=as.numeric(coef(summary(test))[2:4,1]),SE=as.numeric(coef(summary(test))[2:4,2]),
                               F.value=as.numeric(test1$`F value`),p=as.numeric(test1$`Pr(>F)`),
                               n=nrow(dat3),s=length(unique(dat3$Title)))
          result.f <- rbind(result.f,result)
        }else{
          test <- lmer(RR~scale(MAT)+scale(log(Magnitude))+scale(log(Duration))+(1|Title/Site),
                       weights=Variance_weighted_1,dat3)
          test1 <- anova(test)
          
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=rownames(coef(summary(test)))[2:4],
                               Estimate=as.numeric(coef(summary(test))[2:4,1]),SE=as.numeric(coef(summary(test))[2:4,2]),
                               F.value=as.numeric(test1$`F value`),p=as.numeric(test1$`Pr(>F)`),
                               n=nrow(dat3),s=length(unique(dat3$Title)))
          result.f <- rbind(result.f,result)
        }
      }
    }
  }
}  

write.csv(result.f, file = "5-Meta regression-main.csv")

###### Subgroup analysis: Duration -------------------------------
### Depth: 0-30 cm
dat.test <- dat[dat$Depth!="C",]

result.f <- c()
### For diversity (1:2), homogenization and differentiation (12:13)
for (i in (1:2)){
  dat1 <- data.frame(Site=dat.test$Site,Title1=dat.test$Title1,Title2=dat.test$Title2,Title3=dat.test$Title3,
                     Category=dat.test$Duration_type,Microbical_type=dat.test$Microbical_types,Type=dat.test$Type,
                     Magnitude=dat.test$Magnitude/100,N_C=dat.test$N_C,N_T=dat.test$N_T,Replication_weighted=dat.test$Replication_weighted,
                     RR=dat.test[,(i*2+37)],Variance_weighted=dat.test[,(i*2+38)]) 
 
  ### Response ratio
  dat1 <- na.omit(dat1)
  dat1 <- dat1[!duplicated(dat1[,(ncol(dat1)-1):ncol(dat1)]),]
  
  ### For different microbial taxa
  for (Micro in unique(dat1$Microbical_type)){
    dat2 <- dat1[(dat1$Microbical_type==Micro),]
    for (Type in unique(dat2$Type)){
      dat3 <- dat2[(dat2$Type==Type),]
      dat3 <- dat3 %>% group_by(Category) %>% filter(n() >= 3) %>% ungroup()
      dat3 <- dat3 %>% group_by(Category) %>% filter(n_distinct(paste(Title1, Site, sep = "|")) > 1) %>% ungroup()
      
      if (length(unique(dat3$Category)) >= 2){
        result <- c()
        for (class in unique(dat3$Category)){
          dat4 <- dat3[(dat3$Category==class),]
          test <- rma.mv(RR,Variance_weighted,data=dat4,
                         method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                         mods=~1,random=list(~1|Title1/Site))
          
          result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Index="Response ratio",Micro=Micro,
                                Type=Type,Category=class,Estimate=test$b[1],SE=test$se[1],
                                Z=test$zval,p=test$pval,CL=test$ci.lb,CU=test$ci.ub)
          
          result1$n <- nrow(dat4)
          result1$s <- length(unique(dat4$Title1))
          
          dat.compare <- dat2[dat2$Category==class,]
          if((nrow(dat.compare[(dat.compare$Type!=Type),])>=3)&
             (length(unique(paste0(dat.compare[(dat.compare$Type!=Type),]$Title1,dat.compare[(dat.compare$Type!=Type),]$Site)))>1)){
            test <- rma.mv(RR,Variance_weighted,data=dat.compare,
                           method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                           mods=~Type,random=list(~1|Title1/Site))
            result1$QM.treat <- QM <- test$QM
            result1$QMp.treat <- QMp <- test$QMp  
          }else{result1$QM.treat <- result1$QMp.treat <- NA}
          result <- rbind(result,result1)
        }
        
        ### Calculate P-value for between-group heterogeneity
        test <- rma.mv(RR,Variance_weighted,data=dat3,
                       method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                       mods=~Category,random=list(~1|Title1/Site))
        result$QM.type <- QM <- test$QM
        result$QMp.type <- QMp <- test$QMp   
        result.f <- rbind(result.f,result)
      }else {next}
    }
  }
  
  ### Sensitivity
  dat1$RR <- dat1$RR/abs(dat1$Magnitude)
  dat1$Variance_weighted <- dat1$Variance_weighted/(abs(dat1$Magnitude)*abs(dat1$Magnitude))
  
  ### For different microbial taxa
  for (Micro in unique(dat1$Microbical_type)){
    dat2 <- dat1[(dat1$Microbical_type==Micro),]
    for (Type in unique(dat2$Type)){
      dat3 <- dat2[(dat2$Type==Type),]
      dat3 <- dat3 %>% group_by(Category) %>% filter(n() >= 3) %>% ungroup()
      dat3 <- dat3 %>% group_by(Category) %>% filter(n_distinct(paste(Title1, Site, sep = "|")) > 1) %>% ungroup()
      
      if (length(unique(dat3$Category)) >= 2){
        result <- c()
        for (class in unique(dat3$Category)){
          dat4 <- dat3[(dat3$Category==class),]
          test <- rma.mv(RR,Variance_weighted,data=dat4,
                         method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                         mods=~1,random=list(~1|Title1/Site))
          
          result1 <- data.frame(Factor=colnames(dat)[(i*2+37)],Index="Sensitivity",Micro=Micro,
                                Type=Type,Category=class,Estimate=test$b[1],SE=test$se[1],
                                Z=test$zval,p=test$pval,CL=test$ci.lb,CU=test$ci.ub)
          
          result1$n <- nrow(dat4)
          result1$s <- length(unique(dat4$Title1))
          
          dat.compare <- dat2[dat2$Category==class,]
          if((nrow(dat.compare[(dat.compare$Type!=Type),])>=3)&
             (length(unique(paste0(dat.compare[(dat.compare$Type!=Type),]$Title1,dat.compare[(dat.compare$Type!=Type),]$Site)))>1)){
            test <- rma.mv(RR,Variance_weighted,data=dat.compare,
                           method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                           mods=~Type,random=list(~1|Title1/Site))
            result1$QM.treat <- QM <- test$QM
            result1$QMp.treat <- QMp <- test$QMp  
          }else{result1$QM.treat <- result1$QMp.treat <- NA}
          result <- rbind(result,result1)
        }
        
        ### Calculate P-value for between-group heterogeneity
        test <- rma.mv(RR,Variance_weighted,data=dat3,
                       method="REML",control=list(rel.tol=1e-8,iter.max=1000),
                       mods=~Category,random=list(~1|Title1/Site))
        result$QM.type <- QM <- test$QM
        result$QMp.type <- QMp <- test$QMp   
        result.f <- rbind(result.f,result)
      }else {next}
    }
  }
}

### Percent change
result.f$percent <- (exp(result.f$Estimate) - 1) * 100  
result.f$percent_lower <- (exp(result.f$CL) - 1) * 100 
result.f$percent_upper <- (exp(result.f$CU) - 1) * 100  
write.csv(result.f, file = "3-Subgroup analysis-Duration.csv")

###### Linked with environmental factors -------------------------------
### Depth: 0-30 cm
dat.test <- dat[dat$Depth!="C",]

### 1 ### Meta-regression
result.f <- c()
### For diversity (1:2), homogenization and differentiation (12:13)
for (i in (1:2)){
  dat1 <- data.frame(Site=dat.test$Site,Title=dat.test$Title1,
                     Microbical_type=dat.test$Microbical_types,Type=dat.test$Type,
                     N_C=dat.test$N_C,N_T=dat.test$N_T,Replication_weighted=dat.test$Replication_weighted,
                     RR=dat.test[,(i*2+37)],Variance_weighted=dat.test[,(i*2+38)]) 
  dat1$Variance_weighted_1 <- 1/dat1$Variance_weighted
  
  ### env
  for (j in 3:11){
    dat2 <- cbind(dat1,env=dat.test[,(j*2+37)])
    ### Response ratio
    dat2 <- na.omit(dat2)
    dat2 <- dat2[!duplicated(dat2[,8:9]),]
    dat2 <- dat2[dat2$Microbical_type!="Protists",]
    
    for (Micro in unique(dat2$Microbical_type)){
      dat3 <- dat2[dat2$Microbical_type==Micro,]
      for (Type in unique(dat3$Type)){
        dat4 <- dat3[dat3$Type==Type,]
        
        if(nrow(dat4)>10){
          test <- lmer(env~RR+(1|Title/Site),weights=Variance_weighted_1,dat4)
          test1 <- anova(test)
          
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Method="lmer",Microbe=Micro,Type=Type,
                               Variables=colnames(dat.test)[(j*2+37)],
                               Estimate=as.numeric(coef(summary(test))[2,1]),SE=as.numeric(coef(summary(test))[2,2]),
                               t.F.value=as.numeric(test1["RR",5]),p=as.numeric(test1["RR",6]),
                               n=nrow(dat4),s=length(unique(dat4$Title)))
          result.f <- rbind(result.f,result)
        }else{next}
      }
    }
  }
}  

write.csv(result.f, file = "6-Linked env-Correlation.csv")

### 2 ### piecewiseSEM
library(piecewiseSEM)
library(nlme)
library(dplyr)
library(tidyr)

dat.test <- dat[dat$Depth!="C",]
### Take decreased precipitation (P-) as an example
### Richness
sem_data <- dat.test[, c("Title1", "Site","Biomes","MAT","Duration","Sample_depth","Magnitude",
                         "CN","Type","Microbical_types", "Richness_RR")]
sem_data1 <- sem_data[(sem_data$Type=="P-")&(sem_data$Microbical_types=="Bacteria"),c(1:9,11)]
sem_data1 <- na.omit(sem_data1)
sem_data2 <- sem_data[(sem_data$Type=="P-")&(sem_data$Microbical_types=="Fungi"),c(1:9,11)]
sem_data2 <- na.omit(sem_data2)
sem_data3 <- full_join(sem_data1, sem_data2,
    by = c("Title1", "Site","Biomes","MAT","Duration","Sample_depth","Magnitude","CN","Type"))
colnames(sem_data3)[10:11]<- c('Bac_Richness','Fun_Richness')

### Community differentiation
sem_data <- dat.test[, c("Title1", "Site","Biomes","MAT","Duration","Sample_depth","Magnitude",
                         "CN","Type","Microbical_types", "Comp_RR")]
sem_data1 <- sem_data[(sem_data$Type=="P-")&(sem_data$Microbical_types=="Bacteria"),c(1:9,11)]
sem_data1 <- na.omit(sem_data1)
sem_data2 <- sem_data[(sem_data$Type=="P-")&(sem_data$Microbical_types=="Fungi"),c(1:9,11)]
sem_data2 <- na.omit(sem_data2)
sem_data4 <- full_join(sem_data1, sem_data2,
                       by = c("Title1", "Site","Biomes","MAT","Duration","Sample_depth","Magnitude","CN","Type"))
colnames(sem_data4)[10:11]<- c('Bac_Comp','Fun_Comp')

### MBC
sem_data <- dat.test[, c("Title1", "Site","Biomes","MAT","Duration","Sample_depth","Magnitude",
                         "CN","Type","Microbical_types","MBC_RR")]
sem_data5 <- sem_data[(sem_data$Type=="P-"),c(1:9,11)]
sem_data5 <- na.omit(sem_data5)
sem_data5 <- sem_data5[!duplicated(sem_data5[,c(1:10)]),]

### SOC
sem_data <- dat.test[, c("Title1", "Site","Biomes","MAT","Duration","Sample_depth","Magnitude",
                         "CN","Type","Microbical_types","SOC_RR")]
sem_data6 <- sem_data[(sem_data$Type=="P-"),c(1:9,11)]
sem_data6 <- na.omit(sem_data6)
sem_data6 <- sem_data6[!duplicated(sem_data6[,c(1:10)]),]

### Merged
sem_data7 <- full_join(sem_data3, sem_data4, by = c("Title1", "Site","Biomes","MAT","Duration","Sample_depth","Magnitude",
                                       "CN","Type"))
sem_data7 <- full_join(sem_data7, sem_data5, by = c("Title1", "Site","Biomes","MAT","Duration","Sample_depth","Magnitude",
                                                    "CN","Type"))
sem_data7 <- full_join(sem_data7, sem_data6, by = c("Title1", "Site","Biomes","MAT","Duration","Sample_depth","Magnitude",
                                                    "CN","Type"))
### piecewiseSEM - Pathway result 
### For Bacteria
sem_data8 <- sem_data7[, c("Title1", "Site","Magnitude","Duration","Bac_Richness", "Bac_Comp", "MBC_RR", "SOC_RR")]
names(sem_data8) <- c("Title", "Site","Magnitude","Duration", "Richness", "Comp", "MBC", "SOC")
### For Fungi
sem_data8 <- sem_data7[, c("Title1", "Site","Magnitude","Duration", "Fun_Richness", "Fun_Comp", "MBC_RR", "SOC_RR")]
names(sem_data8) <- c("Title", "Site","Magnitude","Duration", "Richness", "Comp", "MBC", "SOC")

Type = "Decreased"
prefix = "Bacteria"
m1 <- lme(Richness ~ Duration, random = ~1|Title/Site, 
          data = sem_data8, na.action = na.exclude)
m2 <- lme(Comp ~ Duration, random = ~1|Title/Site, 
          data = sem_data8, na.action = na.exclude)
m3 <- lme(SOC ~ Richness + Comp, random = ~1|Title/Site, 
          data = sem_data8, na.action = na.exclude)

sem_model <- psem(m1, m2, m3)
result1 <- fisherC(sem_model)
result2 <- AIC_psem(sem_model)
result3 <- coefs(sem_model, standardize = "scale")
result4  <- summary(sem_model)$R2
result3 <- full_join(result3, result4, by = c("Response"))
result.out <- data.frame(Type=Type,prefix=prefix,Env="SOC",
                        FisherC=result1$Fisher.C,df=result1$df,p=result1$P.Value,AIC=result2$AIC,
                        result3[,c(1:8,13:14)])

###### Upscaling projections -------------------------------
library(raster)
library(terra)

### 1 ### Prediction parameters 
MAT <- rast("MAT_ssp245_2021-2060.tif")
Magnitude <- rast("Magnitude_ssp245_2021-2060.tif")
compareGeom(Magnitude, MAT)
MAT_df <- as.data.frame(MAT, xy = TRUE, cells = FALSE)
Magnitude_df <- as.data.frame(Magnitude, xy = TRUE, cells = FALSE)
Magnitude_df <- Magnitude_df %>%
  left_join(MAT_df, by = c("x", "y"))
Magnitude_df <- na.omit(Magnitude_df)
colnames(Magnitude_df) <- c("lng","lat","Magnitude","MAT")
input <- data.frame(Magnitude_df,Duration="5")
mean(input$Magnitude < 0, na.rm = TRUE) ### Decreased precipitation

### 2 ### Model coefficient
dat.test <- dat[dat$Depth!="C",]
result.f <- c()
### For richness (1), and differentiation (13)
for (i in c(1,13)){
  dat1 <- data.frame(Site=dat.test$Site,Title1=dat.test$Title1,
                     MAT=dat.test$MAT,Magnitude=dat.test$Magnitude,Duration=dat.test$Duration,
                     Microbical_type=dat.test$Microbical_types,Type=dat.test$Type,
                     N_C=dat.test$N_C,N_T=dat.test$N_T,Replication_weighted=dat.test$Replication_weighted,
                     RR=dat.test[,(i*2+37)],Variance_weighted=dat.test[,(i*2+38)]) 
  dat1$Variance_weighted_1 <- 1/dat1$Variance_weighted
  
  ### Response ratio
  dat1 <- na.omit(dat1)
  dat1 <- dat1[!duplicated(dat1[,(ncol(dat1)-2):ncol(dat1)]),]
  dat1 <- dat1[dat1$Microbical_type!="Protists",]
  
  for (Micro in unique(dat1$Microbical_type)){
    dat2 <- dat1[dat1$Microbical_type==Micro,]
    for (Type in unique(dat2$Type)){
      dat3 <- dat2[dat2$Type==Type,]
      
      if ((colnames(dat.test)[(i*2+37)]=="Richness_RR")){
        if ((Micro=="Fungi")&(Type=="P-")){
          test <- lm(RR~MAT+log(Magnitude)+Duration,dat3)
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=c(rownames(coef(summary(test)))[1:4],"sd"),
                               Estimate=c(as.numeric(coef(summary(test))[1:4,1]),summary(test)$sigma))
          result.f <- rbind(result.f,result)
        }else{
          test <- lm(RR~MAT+Magnitude+Duration,dat3)
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=c(rownames(coef(summary(test)))[1:4],"sd"),
                               Estimate=c(as.numeric(coef(summary(test))[1:4,1]),summary(test)$sigma))
          result.f <- rbind(result.f,result)
        }
      }else if((colnames(dat.test)[(i*2+37)]=="Comp_RR")){
        if ((Micro=="Bacteria")&(Type=="P-")){
          test <- lm(RR~MAT+Magnitude+log(Duration)+Magnitude*log(Duration),dat3)
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=c(rownames(coef(summary(test)))[1:5],"sd"),
                               Estimate=c(as.numeric(coef(summary(test))[1:5,1]),summary(test)$sigma))
          result.f <- rbind(result.f,result)
        }else if ((Micro=="Bacteria")&(Type=="P+")){
          test <- lm(RR~MAT+Magnitude+log(Duration),dat3)
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=c(rownames(coef(summary(test)))[1:4],"sd"),
                               Estimate=c(as.numeric(coef(summary(test))[1:4,1]),summary(test)$sigma))
          result.f <- rbind(result.f,result)
        }else if ((Micro=="Fungi")&(Type=="P-")){
          test <- lm(RR~MAT+Magnitude+Duration,dat3)
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=c(rownames(coef(summary(test)))[1:4],"sd"),
                               Estimate=c(as.numeric(coef(summary(test))[1:4,1]),summary(test)$sigma))
          result.f <- rbind(result.f,result)
        }else{
          test <- lm(RR~MAT+log(Magnitude)+log(Duration),dat3)
          result <- data.frame(Factor=colnames(dat.test)[(i*2+37)],Microbe=Micro,Type=Type,
                               Variables=c(rownames(coef(summary(test)))[1:4],"sd"),
                               Estimate=c(as.numeric(coef(summary(test))[1:4,1]),summary(test)$sigma))
          result.f <- rbind(result.f,result)
        }
      }
    }
  }
}  

### 3 ### Output
input$Duration <- Duration <- 5
input$Duration <- as.numeric(input$Duration)
input <- input[abs(input$Magnitude) >= 1 & abs(input$Magnitude) <= 100,]
input <- input[input$MAT >= -10 & input$MAT <= 30,]

output <- c()
for (i in 1:nrow(input)){
  if (input$Magnitude[i]>0){
    Bac.rich_mean <-  result.f$Estimate[1]+result.f$Estimate[2]*input$MAT[i]+result.f$Estimate[3]*abs(input$Magnitude[i])+result.f$Estimate[4]*input$Duration[i]
    #Bac.rich_sim <- rnorm(1000, Bac.rich_mean, result.f$Estimate[5])
    #Bac.rich_mean = (exp(mean(Bac.rich_sim)) - 1) * 100
    #Bac.rich_sd = (((exp(mean(Bac.rich_sim)+sd(Bac.rich_sim)) - 1) * 100)-
    #                 ((exp(mean(Bac.rich_sim)-sd(Bac.rich_sim)) - 1) * 100))/2
    
    Fun.rich_mean <-  result.f$Estimate[11]+result.f$Estimate[12]*input$MAT[i]+result.f$Estimate[13]*abs(input$Magnitude[i])+result.f$Estimate[14]*input$Duration[i]
    #Fun.rich_sim <- rnorm(1000, Fun.rich_mean, result.f$Estimate[15])
    #Fun.rich_mean = (exp(mean(Fun.rich_sim)) - 1) * 100
    #Fun.rich_sd = (((exp(mean(Fun.rich_sim)+sd(Fun.rich_sim)) - 1) * 100)-
    #                 ((exp(mean(Fun.rich_sim)-sd(Fun.rich_sim)) - 1) * 100))/2
    
    Bac.Comp_mean <-  result.f$Estimate[27]+result.f$Estimate[28]*input$MAT[i]+result.f$Estimate[29]*abs(input$Magnitude[i])+result.f$Estimate[30]*log(input$Duration[i])
    #Bac.Comp_sim <- rnorm(1000, Bac.Comp_mean, result.f$Estimate[31])
    #Bac.Comp_mean = (exp(mean(Bac.Comp_sim)) - 1) * 100
    #Bac.Comp_sd = (((exp(mean(Bac.Comp_sim)+sd(Bac.Comp_sim)) - 1) * 100)-
    #                 ((exp(mean(Bac.Comp_sim)-sd(Bac.Comp_sim)) - 1) * 100))/2
    
    Fun.Comp_mean <-  result.f$Estimate[37]+result.f$Estimate[38]*input$MAT[i]+result.f$Estimate[39]*log(abs(input$Magnitude[i]))+result.f$Estimate[40]*log(input$Duration[i])
    #Fun.Comp_sim <- rnorm(1000, Fun.Comp_mean, result.f$Estimate[41])
    #Fun.Comp_mean = (exp(mean(Fun.Comp_sim)) - 1) * 100
    #Fun.Comp_sd = (((exp(mean(Fun.Comp_sim)+sd(Fun.Comp_sim)) - 1) * 100)-
    #                 ((exp(mean(Fun.Comp_sim)-sd(Fun.Comp_sim)) - 1) * 100))/2
  }else{
    Bac.rich_mean <-  result.f$Estimate[6]+result.f$Estimate[7]*input$MAT[i]+result.f$Estimate[8]*abs(input$Magnitude[i])+result.f$Estimate[9]*input$Duration[i]
    #Bac.rich_sim <- rnorm(1000, Bac.rich_mean, result.f$Estimate[10])
    #Bac.rich_mean = (exp(mean(Bac.rich_sim)) - 1) * 100
    #Bac.rich_sd = (((exp(mean(Bac.rich_sim)+sd(Bac.rich_sim)) - 1) * 100)-
    #                 ((exp(mean(Bac.rich_sim)-sd(Bac.rich_sim)) - 1) * 100))/2
    
    Fun.rich_mean <-  result.f$Estimate[16]+result.f$Estimate[17]*input$MAT[i]+result.f$Estimate[18]*log(abs(input$Magnitude[i]))+result.f$Estimate[19]*input$Duration[i]
    #Fun.rich_sim <- rnorm(1000, Fun.rich_mean, result.f$Estimate[20])
    #Fun.rich_mean = (exp(mean(Fun.rich_sim)) - 1) * 100
    #Fun.rich_sd = (((exp(mean(Fun.rich_sim)+sd(Fun.rich_sim)) - 1) * 100)-
    #                 ((exp(mean(Fun.rich_sim)-sd(Fun.rich_sim)) - 1) * 100))/2
    
    Bac.Comp_mean <-  result.f$Estimate[21]+result.f$Estimate[22]*input$MAT[i]+result.f$Estimate[23]*abs(input$Magnitude[i])+result.f$Estimate[24]*log(input$Duration[i])+result.f$Estimate[25]*abs(input$Magnitude[i])*log(input$Duration[i])
    #Bac.Comp_sim <- rnorm(1000, Bac.Comp_mean, result.f$Estimate[26])
    #Bac.Comp_mean = (exp(mean(Bac.Comp_sim)) - 1) * 100
    #Bac.Comp_sd = (((exp(mean(Bac.Comp_sim)+sd(Bac.Comp_sim)) - 1) * 100)-
    #                 ((exp(mean(Bac.Comp_sim)-sd(Bac.Comp_sim)) - 1) * 100))/2
    
    Fun.Comp_mean <-  result.f$Estimate[32]+result.f$Estimate[33]*input$MAT[i]+result.f$Estimate[34]*abs(input$Magnitude[i])+result.f$Estimate[35]*input$Duration[i]
    #Fun.Comp_sim <- rnorm(1000, Fun.Comp_mean, result.f$Estimate[36])
    #Fun.Comp_mean = (exp(mean(Fun.Comp_sim)) - 1) * 100
    #Fun.Comp_sd = (((exp(mean(Fun.Comp_sim)+sd(Fun.Comp_sim)) - 1) * 100)-
    #                 ((exp(mean(Fun.Comp_sim)-sd(Fun.Comp_sim)) - 1) * 100))/2
  }
  
  result1 <- data.frame(input[i,],Bac.rich_mean=Bac.rich_mean,#Bac.rich_sd=Bac.rich_sd,
                        Fun.rich_mean=Fun.rich_mean,#Fun.rich_sd=Fun.rich_sd,
                        Bac.Comp_mean=Bac.Comp_mean,#Bac.Comp_sd=Bac.Comp_sd,
                        Fun.Comp_mean=Fun.Comp_mean)#Fun.Comp_sd=Fun.Comp_sd)
  output <- rbind(output,result1)
}
