Evaluation_TML <- function(algoritmo,env,env.table,occ.treino,occ.teste,N,dir.out,replicas,porcentagem){
  
  options(digits=3)
  
  if(any(algoritmo=="BIO")){
    Env.Bio <- data.frame()
    Thrs.Bio <- data.frame()
  }
  if(any(algoritmo=="MXS")){
    Eval.Mxs <- data.frame()
    Thrs.Mxs <- data.frame()
  }
  if(any(algoritmo=="MXD")){
    Eval.Mxd <- data.frame()
    Thrs.Mxd <- data.frame()
  }
  if(any(algoritmo=="SVR")){
    Eval.Svr <- data.frame()
    Thrs.Svr <- data.frame()
  }
  if(any(algoritmo=="SVC")){
    Eval.Svc <- data.frame()
    Thrs.Svc <- data.frame()
  }
  if(any(algoritmo=="GLM")){
    Eval.Glm <- data.frame()
    Thrs.Glm <- data.frame()
  }
  if(any(algoritmo=="RDF")){
    Eval.Rdf <- data.frame()
    Thrs.Rdf <- data.frame()
  }
  if(any(algoritmo=="MDA")){
    Eval.Mda <- data.frame()
    Thrs.Mda <- data.frame()
  }
  
  for (g in 1:replicas){
    print(paste("Replica...",g,sep=""))
    
    if (porcentagem!=1){
      ocorrencias.teste.rep <- occ.teste[[g]]
    }else{
      ocorrencias.teste.rep <- occ.treino[[g]]
    }
    
    sps <- as.character(unique(ocorrencias.teste.rep[,1]))
    
    for (i in 1:N){
      
      print(as.character(sps[i]))
    
      ocorrencias.teste.sp<-ocorrencias.teste.rep[ocorrencias.teste.rep[,1]==sps[i],]
        
      #Select Pseudo-Absences for Evaluation
      set.seed(N*g)
      pb.eval <- sample(1:nrow(env.table),10000)
      pb.eval<-env.table[pb.eval,c(1,2)]
      colnames(pb.eval) <- c("Long","Lat")
        
      if (any(algoritmo=="BIO")){
        Eval<-dismo::evaluate(ocorrencias.teste.sp[,2:3], pb.eval, Bioclim[[i]], env)
        ROC <- which.max(Eval@TPR + Eval@TNR)
        LPT <- which(round(Eval@t,3)==min(round(Eval@presence,3)))[1]
        N.Pres <- Eval@np
        TSS.ROC <- (Eval@TPR[ROC]+Eval@TNR[ROC])-1
        TSS.LPT <- (Eval@TPR[LPT]+Eval@TNR[LPT])-1
        TPR.ROC <- Eval@TPR[ROC]
        TPR.LPT <- Eval@TPR[LPT]
        TNR.ROC <- Eval@TNR[ROC]
        TNR.LPT <- Eval@TNR[LPT]
        Thr.ROC<- Eval@t[ROC]
        Thr.LPT<- min(Eval@presence)
        AUC <- Eval@auc
        #Create Evaluation table for a species
        Eval.table <- cbind.data.frame(sps[i],N.Pres,TSS.ROC,TPR.ROC,TNR.ROC,Thr.ROC,TSS.LPT,TPR.LPT,TNR.LPT,Thr.LPT,AUC)
        Thrs.table <- cbind.data.frame(sps[i],Thr.ROC,Thr.LPT)
        #Create Evaluation table for a Algorithm
        Eval.Bio <- rbind.data.frame(Eval.Bio,Eval.table)
        Thrs.Bio <- rbind.data.frame(Thrs.Bio,Thrs.table)
      }
      
      if (any(algoritmo=="MXS")){
        Eval<-dismo::evaluate(ocorrencias.teste.sp[,2:3], pb.eval, MaxentS[[i]], env)
        ROC <- which.max(Eval@TPR + Eval@TNR)
        LPT <- which(round(Eval@t,3)==min(round(Eval@presence,3)))[1]
        N.Pres <- Eval@np
        TSS.ROC <- (Eval@TPR[ROC]+Eval@TNR[ROC])-1
        TSS.LPT <- (Eval@TPR[LPT]+Eval@TNR[LPT])-1
        TPR.ROC <- Eval@TPR[ROC]
        TPR.LPT <- Eval@TPR[LPT]
        TNR.ROC <- Eval@TNR[ROC]
        TNR.LPT <- Eval@TNR[LPT]
        Thr.ROC<- Eval@t[ROC]
        Thr.LPT<- min(Eval@presence)
        AUC <- Eval@auc
        #Create Evaluation table for a species
        Eval.table <- cbind.data.frame(sps[i],N.Pres,TSS.ROC,TPR.ROC,TNR.ROC,Thr.ROC,TSS.LPT,TPR.LPT,TNR.LPT,Thr.LPT,AUC)
        Thrs.table <- cbind.data.frame(sps[i],Thr.ROC,Thr.LPT)
        #Create Evaluation table for a Algorithm
        Eval.Mxs <- rbind.data.frame(Eval.Mxs,Eval.table)
        Thrs.Mxs <- rbind.data.frame(Thrs.Mxs,Thrs.table)
      }
      
      if (any(algoritmo=="MXD")){
        Eval<-dismo::evaluate(ocorrencias.teste.sp[,2:3], pb.eval, MaxentD[[i]], env)
        ROC <- which.max(Eval@TPR + Eval@TNR)
        LPT <- which(round(Eval@t,3)==min(round(Eval@presence,3)))[1]
        N.Pres <- Eval@np
        TSS.ROC <- (Eval@TPR[ROC]+Eval@TNR[ROC])-1
        TSS.LPT <- (Eval@TPR[LPT]+Eval@TNR[LPT])-1
        TPR.ROC <- Eval@TPR[ROC]
        TPR.LPT <- Eval@TPR[LPT]
        TNR.ROC <- Eval@TNR[ROC]
        TNR.LPT <- Eval@TNR[LPT]
        Thr.ROC<- Eval@t[ROC]
        Thr.LPT<- min(Eval@presence)
        AUC <- Eval@auc
        #Create Evaluation table for a species
        Eval.table <- cbind.data.frame(sps[i],N.Pres,TSS.ROC,TPR.ROC,TNR.ROC,Thr.ROC,TSS.LPT,TPR.LPT,TNR.LPT,Thr.LPT,AUC)
        Thrs.table <- cbind.data.frame(sps[i],Thr.ROC,Thr.LPT)
        #Create Evaluation table for a Algorithm
        Eval.Mxd <- rbind.data.frame(Eval.Mxd,Eval.table)
        Thrs.Mxd <- rbind.data.frame(Thrs.Mxd,Thrs.table)
      }
      
      if (any(algoritmo=="SVR")){
        Eval<-dismo::evaluate(ocorrencias.teste.sp[,2:3], pb.eval, SVR[[i]], env)
        ROC <- which.max(Eval@TPR + Eval@TNR)
        LPT <- which(round(Eval@t,3)==min(round(Eval@presence,3)))[1]
        N.Pres <- Eval@np
        TSS.ROC <- (Eval@TPR[ROC]+Eval@TNR[ROC])-1
        TSS.LPT <- (Eval@TPR[LPT]+Eval@TNR[LPT])-1
        TPR.ROC <- Eval@TPR[ROC]
        TPR.LPT <- Eval@TPR[LPT]
        TNR.ROC <- Eval@TNR[ROC]
        TNR.LPT <- Eval@TNR[LPT]
        Thr.ROC<- Eval@t[ROC]
        Thr.LPT<- min(Eval@presence)
        AUC <- Eval@auc
        #Create Evaluation table for a species
        Eval.table <- cbind.data.frame(sps[i],N.Pres,TSS.ROC,TPR.ROC,TNR.ROC,Thr.ROC,TSS.LPT,TPR.LPT,TNR.LPT,Thr.LPT,AUC)
        Thrs.table <- cbind.data.frame(sps[i],Thr.ROC,Thr.LPT)
        #Create Evaluation table for a Algorithm
        Eval.Svr <- rbind.data.frame(Eval.Svr,Eval.table)
        Thrs.Svr <- rbind.data.frame(Thrs.Svr,Thrs.table)
      }
      
      if (any(algoritmo=="SVC")){
        Eval<-dismo::evaluate(ocorrencias.teste.sp[,2:3], pb.eval, SVC[[i]], env)
        ROC <- which.max(Eval@TPR + Eval@TNR)
        LPT <- which(round(Eval@t,3)==min(round(Eval@presence,3)))[1]
        N.Pres <- Eval@np
        TSS.ROC <- (Eval@TPR[ROC]+Eval@TNR[ROC])-1
        TSS.LPT <- (Eval@TPR[LPT]+Eval@TNR[LPT])-1
        TPR.ROC <- Eval@TPR[ROC]
        TPR.LPT <- Eval@TPR[LPT]
        TNR.ROC <- Eval@TNR[ROC]
        TNR.LPT <- Eval@TNR[LPT]
        Thr.ROC<- Eval@t[ROC]
        Thr.LPT<- min(Eval@presence)
        AUC <- Eval@auc
        #Create Evaluation table for a species
        Eval.table <- cbind.data.frame(sps[i],N.Pres,TSS.ROC,TPR.ROC,TNR.ROC,Thr.ROC,TSS.LPT,TPR.LPT,TNR.LPT,Thr.LPT,AUC)
        Thrs.table <- cbind.data.frame(sps[i],Thr.ROC,Thr.LPT)
        #Create Evaluation table for a Algorithm
        Eval.Svc <- rbind.data.frame(Eval.Svc,Eval.table)
        Thrs.Svc <- rbind.data.frame(Thrs.Svc,Thrs.table)
      }
      
      if (any(algoritmo=="GLM")){
        Eval<-dismo::evaluate(ocorrencias.teste.sp[,2:3], pb.eval, GLM[[i]], env)
        ROC <- which.max(Eval@TPR + Eval@TNR)
        LPT <- which(round(Eval@t,3)==min(round(Eval@presence,3)))[1]
        N.Pres <- Eval@np
        TSS.ROC <- (Eval@TPR[ROC]+Eval@TNR[ROC])-1
        TSS.LPT <- (Eval@TPR[LPT]+Eval@TNR[LPT])-1
        TPR.ROC <- Eval@TPR[ROC]
        TPR.LPT <- Eval@TPR[LPT]
        TNR.ROC <- Eval@TNR[ROC]
        TNR.LPT <- Eval@TNR[LPT]
        Thr.ROC<- Eval@t[ROC]
        Thr.LPT<- min(Eval@presence)
        AUC <- Eval@auc
        #Create Evaluation table for a species
        Eval.table <- cbind.data.frame(sps[i],N.Pres,TSS.ROC,TPR.ROC,TNR.ROC,Thr.ROC,TSS.LPT,TPR.LPT,TNR.LPT,Thr.LPT,AUC)
        Thrs.table <- cbind.data.frame(sps[i],Thr.ROC,Thr.LPT)
        #Create Evaluation table for a Algorithm
        Eval.Glm <- rbind.data.frame(Eval.Glm,Eval.table)
        Thrs.Glm <- rbind.data.frame(Thrs.Glm,Thrs.table)
      }
      
      if (any(algoritmo=="RDF")){
        Eval<-dismo::evaluate(ocorrencias.teste.sp[,2:3], pb.eval, RDF[[i]], env)
        ROC <- which.max(Eval@TPR + Eval@TNR)
        LPT <- which(round(Eval@t,3)==min(round(Eval@presence,3)))[1]
        N.Pres <- Eval@np
        TSS.ROC <- (Eval@TPR[ROC]+Eval@TNR[ROC])-1
        TSS.LPT <- (Eval@TPR[LPT]+Eval@TNR[LPT])-1
        TPR.ROC <- Eval@TPR[ROC]
        TPR.LPT <- Eval@TPR[LPT]
        TNR.ROC <- Eval@TNR[ROC]
        TNR.LPT <- Eval@TNR[LPT]
        Thr.ROC<- Eval@t[ROC]
        Thr.LPT<- min(Eval@presence)
        AUC <- Eval@auc
        #Create Evaluation table for a species
        Eval.table <- cbind.data.frame(sps[i],N.Pres,TSS.ROC,TPR.ROC,TNR.ROC,Thr.ROC,TSS.LPT,TPR.LPT,TNR.LPT,Thr.LPT,AUC)
        Thrs.table <- cbind.data.frame(sps[i],Thr.ROC,Thr.LPT)
        #Create Evaluation table for a Algorithm
        Eval.Rdf <- rbind.data.frame(Eval.Rdf,Eval.table)
        Thrs.Rdf <- rbind.data.frame(Thrs.Rdf,Thrs.table)
      }
      
      if (any(algoritmo=="MDA")){
        Eval<-dismo::evaluate(ocorrencias.teste.sp[,2:3], pb.eval, MDA[[i]], env)
        ROC <- which.max(Eval@TPR + Eval@TNR)
        LPT <- which(round(Eval@t,3)==min(round(Eval@presence,3)))[1]
        N.Pres <- Eval@np
        TSS.ROC <- (Eval@TPR[ROC]+Eval@TNR[ROC])-1
        TSS.LPT <- (Eval@TPR[LPT]+Eval@TNR[LPT])-1
        TPR.ROC <- Eval@TPR[ROC]
        TPR.LPT <- Eval@TPR[LPT]
        TNR.ROC <- Eval@TNR[ROC]
        TNR.LPT <- Eval@TNR[LPT]
        Thr.ROC<- Eval@t[ROC]
        Thr.LPT<- min(Eval@presence)
        AUC <- Eval@auc
        #Create Evaluation table for a species
        Eval.table <- cbind.data.frame(sps[i],N.Pres,TSS.ROC,TPR.ROC,TNR.ROC,Thr.ROC,TSS.LPT,TPR.LPT,TNR.LPT,Thr.LPT,AUC)
        Thrs.table <- cbind.data.frame(sps[i],Thr.ROC,Thr.LPT)
        #Create Evaluation table for a Algorithm
        Eval.Mda <- rbind.data.frame(Eval.Mda,Eval.table)
        Thrs.Mda <- rbind.data.frame(Thrs.Mda,Thrs.table)
      }
    }
  }
  
  setwd(dir.out)
  
  if(any(algoritmo=="BIO")){
    colnames(Eval.Bio) <- c("Species","N","TSS.ROC","TPR.ROC","TNR.ROC","Thr.ROC","TSS.LPT","TPR.LPT","TNR.LPT","Thr.LPT","AUC")
    Eval.Bio[,1] <- as.character(Eval.Bio[,1])
    Eval.Bio <- Eval.Bio[order(Eval.Bio$Species),]
    colnames(Thrs.Bio) <- c("Species","Thr.ROC","Thr.LPT")
    Thrs.Bio <- Thrs.Bio[order(Thrs.Bio$Species),]
    
    write.xlsx(Eval.Bio,"AvaliacaoBioclim.xls",row.names=F)
    write.xlsx(Thrs.Bio,"ThresholdsBioclim.xls",row.names=F)
  }
  if(any(algoritmo=="MXS")){
    colnames(Eval.Mxs) <- c("Species","N","TSS.ROC","TPR.ROC","TNR.ROC","Thr.ROC","TSS.LPT","TPR.LPT","TNR.LPT","Thr.LPT","AUC")
    Eval.Mxs[,1] <- as.character(Eval.Mxs[,1])
    Eval.Mxs <- Eval.Mxs[order(Eval.Mxs$Species),]
    colnames(Thrs.Mxs) <- c("Species","Thr.ROC","Thr.LPT")
    Thrs.Mxs <- Thrs.Mxs[order(Thrs.Mxs$Species),]
    
    write.xlsx(Eval.Mxs,"AvaliacaoMXS.xls",row.names=F)
    write.xlsx(Thrs.Mxs,"ThresholdsMXS.xls",row.names=F)
  }
  if(any(algoritmo=="MXD")){
    colnames(Eval.Mxd) <- c("Species","N","TSS.ROC","TPR.ROC","TNR.ROC","Thr.ROC","TSS.LPT","TPR.LPT","TNR.LPT","Thr.LPT","AUC")
    Eval.Mxd[,1] <- as.character(Eval.Mxd[,1])
    Eval.Mxd <- Eval.Mxd[order(Eval.Mxd$Species),]
    colnames(Thrs.Mxd) <- c("Species","Thr.ROC","Thr.LPT")
    Thrs.Mxd <- Thrs.Mxd[order(Thrs.Mxd$Species),]
    
    write.xlsx(Eval.Mxd,"AvaliacaoMXD.xls",row.names=F)
    write.xlsx(Thrs.Mxd,"ThresholdsMXD.xls",row.names=F)
  }
  if(any(algoritmo=="SVR")){
    colnames(Eval.Svr) <- c("Species","N","TSS.ROC","TPR.ROC","TNR.ROC","Thr.ROC","TSS.LPT","TPR.LPT","TNR.LPT","Thr.LPT","AUC")
    Eval.Svr[,1] <- as.character(Eval.Svr[,1])
    Eval.Svr <- Eval.Svr[order(Eval.Svr$Species),]
    colnames(Thrs.Svr) <- c("Species","Thr.ROC","Thr.LPT")
    Thrs.Svr <- Thrs.Svr[order(Thrs.Svr$Species),]
    
    write.xlsx(Eval.Svr,"AvaliacaoSVR.xls",row.names=F)
    write.xlsx(Thrs.Svr,"ThresholdsSVR.xls",row.names=F)
  }
  if(any(algoritmo=="SVC")){
    colnames(Eval.Svc) <- c("Species","N","TSS.ROC","TPR.ROC","TNR.ROC","Thr.ROC","TSS.LPT","TPR.LPT","TNR.LPT","Thr.LPT","AUC")
    Eval.Svc[,1] <- as.character(Eval.Svc[,1])
    Eval.Svc <- Eval.Svc[order(Eval.Svc$Species),]
    colnames(Thrs.Svc) <- c("Species","Thr.ROC","Thr.LPT")
    Thrs.Svc <- Thrs.Svc[order(Thrs.Svc$Species),]
    
    write.xlsx(Eval.Svc,"AvaliacaoSVC.xls",row.names=F)
    write.xlsx(Thrs.Svc,"ThresholdsSVC.xls",row.names=F)
  }
  if(any(algoritmo=="GLM")){
    colnames(Eval.Glm) <- c("Species","N","TSS.ROC","TPR.ROC","TNR.ROC","Thr.ROC","TSS.LPT","TPR.LPT","TNR.LPT","Thr.LPT","AUC")
    Eval.Glm[,1] <- as.character(Eval.Glm[,1])
    Eval.Glm <- Eval.Glm[order(Eval.Glm$Species),]
    colnames(Thrs.Glm) <- c("Species","Thr.ROC","Thr.LPT")
    Thrs.Glm <- Thrs.Glm[order(Thrs.Glm$Species),]
    
    write.xlsx(Eval.Glm,"AvaliacaoGLM.xls",row.names=F)
    write.xlsx(Thrs.Glm,"ThresholdsGLM.xls",row.names=F)
  }
  if(any(algoritmo=="RDF")){
    colnames(Eval.Rdf) <- c("Species","N","TSS.ROC","TPR.ROC","TNR.ROC","Thr.ROC","TSS.LPT","TPR.LPT","TNR.LPT","Thr.LPT","AUC")
    Eval.Rdf[,1] <- as.character(Eval.Rdf[,1])
    Eval.Rdf <- Eval.Rdf[order(Eval.Rdf$Species),]
    colnames(Thrs.Rdf) <- c("Species","Thr.ROC","Thr.LPT")
    Thrs.Rdf <- Thrs.Rdf[order(Thrs.Rdf$Species),]
    
    write.xlsx(Eval.Rdf,"AvaliacaoRDF.xls",row.names=F)
    write.xlsx(Thrs.Rdf,"ThresholdsRDF.xls",row.names=F)
  }
  if(any(algoritmo=="MDA")){
    colnames(Eval.Mda) <- c("Species","N","TSS.ROC","TPR.ROC","TNR.ROC","Thr.ROC","TSS.LPT","TPR.LPT","TNR.LPT","Thr.LPT","AUC")
    Eval.Mda[,1] <- as.character(Eval.Mda[,1])
    Eval.Mda <- Eval.Mda[order(Eval.Mda$Species),]
    colnames(Thrs.Mda) <- c("Species","Thr.ROC","Thr.LPT")
    Thrs.Mda <- Thrs.Mda[order(Thrs.Mda$Species),]
    
    write.xlsx(Eval.Mda,"AvaliacaoMDA.xls",row.names=F)
    write.xlsx(Thrs.Mda,"ThresholdsMDA.xls",row.names=F)
  }
  
}