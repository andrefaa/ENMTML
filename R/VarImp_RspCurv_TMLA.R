VarImp_RspCurv <- function(Model,
                           Algorithm,
                           spN,
                           SpDataT,
                           VarColT,
                           Outcome){
  #Function to generate calculations of Variable Importance and Response Curves for Models
  #Written by Andre Andrade
  
  #Parameters:
    #Model: The fitted Model
    #Algorithm: The algorithm used to fit the model
    #spN: Species name
    #SpDataT: Table with environmental and presence-absence informatio
    #VarColT: Vector with variables names
    #Outcome: Outcome predicted by the model, used when varImp is unavailable for the model class
  
  #Create Folder to Save Files
  dir.create(file.path(grep(Algorithm,folders,value=T),"Response Curves & Variable Importance"))
  dir.create(file.path(grep(Algorithm,folders,value=T),"Response Curves & Variable Importance",spN))
  DirV <- file.path(grep(Algorithm,folders,value=T),"Response Curves & Variable Importance",spN)
  
  #Variable Importance & Response Curves
  if(Algorithm %in% c("BIO","DOM","MAH")){
    warning("Variable importance Evaluated using a filter approach")
    V_IMP <- filterVarImp(SpDataTM[,VarColT],Outcome, nonpara = FALSE)
    V_IMP <- V_IMP/sum(V_IMP)
    write.table(V_IMP,file.path(DirV,"VariableImportance.txt"),sep="\t",row.names=T)
    #Response Curves
    for(k in 1:length(VarColT)){
      tiff(file.path(DirV,paste0("Response_Curves_",VarColT[k],".tiff")),width = 3200, height = 3200, units = "px", res = 800)
      response(Model,var=VarColT[k])
      dev.off()
    }
  }
  
  if(Algorithm %in% "ENF"){
    write.table(Model$cor,file.path(DirV,"VariableImportance.txt"),sep="\t",row.names=T)
    #Response Curves
    tiff(file.path(DirV,"ENFA_Axis.tiff"),width = 3200, height = 3200, units = "px", res = 800)
    scatter(Model)
    dev.off()
    for(k in 1:length(VarColT)){
      for(o in 1:2){
        tiff(file.path(DirV,paste0("Response_Curves_",VarColT[k],"_Axis_",o,".tiff")),width = 3200, height = 3200, units = "px", res = 800)
        scatterniche(SpDataT[,VarColT],pr=rep(1,7),xax=k,yax=o)
        dev.off()
      }
    }
  }
  
  if(Algorithm %in% "BRT"){
    V_IMP <- varImp(Model,numTrees = Model$n.trees, scale=T)
    V_IMP <- V_IMP/sum(V_IMP)
    write.table(V_IMP,file.path(DirV,"VariableImportance.txt"),sep="\t",row.names=T)
    #Response Curves
    for(k in 1:length(VarColT)){
      tiff(file.path(DirV,paste0("Response_Curves_",VarColT[k],".tiff")),width = 3200, height = 3200, units = "px", res = 800)
      gbm.plot(Model,variable.no=k,plot.layout=c(1, 1),write.title=F)
      dev.off()
    }
  }
  
  if(Algorithm %in% "RDF"){
    V_IMP <- varImp(Model,type=1,numTrees = Model$ntree, scale=T)
    V_IMP <- round(V_IMP/sum(V_IMP),3)
    write.table(V_IMP,file.path(DirV,"VariableImportance.txt"),sep="\t",row.names=T)
    #Response Curves
    for(k in 1:length(VarColT)){
      tiff(file.path(DirV,paste0("Response_Curves_",VarColT[k],".tiff")),width = 3200, height = 3200, units = "px", res = 800)
      partialPlot(Model,x.var=VarColT[k],pred.data=SpDataT[,VarColT],xlab=VarColT[k],main=NULL)
      dev.off()
    }
  }
  
  if(Algorithm %in% c("MXD","MXS","MLK")){
    warning("Variable importance Evaluated using a filter approach")
    V_IMP <- filterVarImp(SpDataT[,VarColT],Outcome,nonpara = FALSE)
    V_IMP <- V_IMP/sum(V_IMP)
    write.table(V_IMP,file.path(DirV,"VariableImportance.txt"),sep="\t",row.names=T)
    #Response Curves
    for(k in 1:length(VarColT)){
      tiff(file.path(DirV,paste0("Response_Curves_",VarColT[k],".tiff")),width = 3200, height = 3200, units = "px", res = 800)
      plot(Model, vars=VarColT[k],type="logit")
      dev.off()
    }
  }
  
  if(Algorithm %in% "SVM"){
    V_IMP <- varImp(Model,type=1,numTrees = Model$ntree, scale=T)
    V_IMP <- round(V_IMP/sum(V_IMP),3)
    write.table(V_IMP,file.path(DirV,"VariableImportance.txt"),sep="\t",row.names=T)
    #Response Curves
    for(k in 1:length(VarColT)){
      tiff(file.path(DirV,paste0("Response_Curves_",VarColT[k],".tiff")),width = 3200, height = 3200, units = "px", res = 800)
      partialPlot(Model,x.var=VarColT[k],pred.data=SpDataT[,VarColT],xlab=VarColT[k],main=NULL)
      dev.off()
    }
  }
  
  if(Algorithm %in% c("GLM","GAM")){
    V_IMP <- varImp(Model)
    V_IMP <- round(V_IMP/sum(V_IMP),3)
    write.table(V_IMP,file.path(DirV,"VariableImportance.txt"),sep="\t",row.names=T)
    #Response Curves
    for(k in 1:length(VarColT)){
      tiff(file.path(DirV,paste0("Response_Curves_",VarColT[k],".tiff")),width = 3200, height = 3200, units = "px", res = 800)
      response(Model,var=VarColT[k])
      dev.off()
    }
  }
  
  if(Algorithm %in% "GAU"){
    V_IMP <- varImp(Model,type=1,numTrees = Model$ntree, scale=T)
    V_IMP <- round(V_IMP/sum(V_IMP),3)
    write.table(V_IMP,file.path(DirV,"VariableImportance.txt"),sep="\t",row.names=T)
    #Response Curves
    for(k in 1:length(VarColT)){
      tiff(file.path(DirV,paste0("Response_Curves_",VarColT[k],".tiff")),width = 3200, height = 3200, units = "px", res = 800)
      partialPlot(Model,x.var=VarColT[k],pred.data=SpDataT[,VarColT],xlab=VarColT[k],main=NULL)
      dev.off()
    }
  }
}
