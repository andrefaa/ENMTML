VarImp_RspCurv <- function(Model,
                           Algorithm,
                           folders,
                           spN,
                           SpDataT,
                           VarColT,
                           Outcome) {
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
  dir.create(file.path(
    grep(Algorithm, folders, value = T),
    "Response Curves & Variable Importance"
  ))
  DirV <-
    file.path(grep(Algorithm, folders, value = T),
              "Response Curves & Variable Importance")
  
  #Variable Importance & Response Curves
  
  #Bioclim, Domain, Mahalanobis & ENFA----
  if (Algorithm %in% c("BIO", "DOM", "MAH", "ENF")) {
    warning("Variable importance Evaluated using a filter approach")
    V_IMP <- caret::filterVarImp(SpDataT[, VarColT], Outcome, nonpara = FALSE)
    V_IMP <- V_IMP / sum(V_IMP)
    V_IMP <-
      cbind(
        Sp = spN,
        Algorithm = Algorithm,
        Variables = VarColT,
        Importance = V_IMP
      )
    row.names(V_IMP) <- NULL
    if (file.exists(paste(DirV, "/VariableImportance.txt", sep = ""))) {
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    #Response Curves
    if (Algorithm != 'ENF') {
      png(
        file.path(DirV, paste0("ResponseCurves_", spN, ".png")),
        width = 3200,
        height = 3200,
        units = "px",
        res = 800
      )
      dismo::response(Model)
      dev.off()
    } else{
      png(
        file.path(DirV, paste0("ENFAAxis_", spN, ".png")),
        width = 3200,
        height = 3200,
        units = "px",
        res = 800
      )
      ade4::scatter(Model)
      dev.off()
    }
  }
  
  #Maxent----
  if (Algorithm %in% c("MXD", "MXS")) {
    warning("Variable importance Evaluated using a filter approach")
    V_IMP <- caret::filterVarImp(SpDataT[, VarColT], Outcome, nonpara = FALSE)
    V_IMP <- V_IMP / sum(V_IMP)
    V_IMP <-
      cbind(
        Sp = spN,
        Algorithm = Algorithm,
        Variables = VarColT,
        Importance = V_IMP
      )
    row.names(V_IMP) <- NULL
    if (file.exists(paste(DirV, "/VariableImportance.txt", sep = ""))) {
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    #Response Curves
    if (!Algorithm %in% "MLK") {
      png(
        file.path(DirV, paste0("ResponseCurves_", spN, ".png")),
        width = 3200,
        height = 3200,
        units = "px",
        res = 800
      )
      plot(Model, type = "cloglog")
      dev.off()
    } else{
      plot(Model)
    }
  }
  
  
  
  #Boosted Regression Tree----
  if (Algorithm %in% "BRT") {
    V_IMP <- caret::varImp(Model, numTrees = Model$n.trees, scale = T)
    V_IMP <- V_IMP / sum(V_IMP)
    V_IMP <-
      cbind(
        Sp = spN,
        Algorithm = Algorithm,
        Variables = VarColT,
        Importance = V_IMP
      )
    row.names(V_IMP) <- NULL
    if (file.exists(paste(DirV, "/VariableImportance.txt", sep = ""))) {
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    #Response Curves
    png(
      file.path(DirV, paste0("ResponseCurves_", spN, ".png")),
      width = 3200,
      height = 3200,
      units = "px",
      res = 800
    )
    dismo::gbm.plot(Model,
             plot.layout = c(length(Model$var.names) / 3, 3),
             write.title = F)
    dev.off()
  }
  
  #Random Forests----
  if (Algorithm %in% "RDF") {
    V_IMP <- caret::varImp(Model,
                    type = 1,
                    numTrees = Model$ntree,
                    scale = T)
    V_IMP <- round(V_IMP / sum(V_IMP), 3)
    V_IMP <-
      cbind(
        Sp = spN,
        Algorithm = Algorithm,
        Variables = VarColT,
        Importance = V_IMP
      )
    row.names(V_IMP) <- NULL
    if (file.exists(paste(DirV, "/VariableImportance.txt", sep = ""))) {
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    #Response Curves
    png(
      file.path(DirV, paste0("ResponseCurves_", spN, ".png")),
      width = 3200,
      height = 3200,
      units = "px",
      res = 800
    )
    par(mfrow = c(ceiling(nrow(
      Model$importance
    ) / 3), 3))
    for (o in 1:nrow(Model$importance)) {
      randomForest::partialPlot(
        Model,
        x.var = row.names(Model$importance)[o],
        pred.data = SpDataT[, VarColT],
        xlab = row.names(Model$importance)[o],
        main = NULL
      )
    }
    dev.off()
  }
  
  #Support Vector Machine----
  if (Algorithm %in% "SVM") {
    V_IMP <- caret::filterVarImp(SpDataT[, VarColT], Outcome, nonpara = FALSE)
    V_IMP <- V_IMP / sum(V_IMP)
    V_IMP <-
      cbind(
        Sp = spN,
        Algorithm = Algorithm,
        Variables = VarColT,
        Importance = V_IMP
      )
    row.names(V_IMP) <- NULL
    if (file.exists(paste(DirV, "/VariableImportance.txt", sep = ""))) {
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    # #Response Curves
    # png(file.path(DirV,paste0("Response_Curves_",VarColT[k],".png")),width = 3200, height = 3200, units = "px", res = 800)
    # partialPlot(Model,pred.data=SpDataT[,VarColT],main=NULL)
    # dev.off()
  }
  
  #GLM,GAM----
  if (Algorithm %in% c("GLM", "GAM")) {
    V_IMP <- caret::varImp(Model)
    V_IMP <- round(V_IMP / sum(V_IMP), 3)
    V_IMP <-
      data.frame(
        Sp = spN,
        Algorithm = Algorithm,
        Variables = row.names(V_IMP),
        Importance = V_IMP
      )
    row.names(V_IMP) <- NULL
    if (file.exists(paste(DirV, "/VariableImportance.txt", sep = ""))) {
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    #Response Curves
    if (Algorithm %in% "GLM") {
      png(
        file.path(DirV, paste0("ResponseCurves_", spN, ".png")),
        width = 1050 * ceiling(length(Model$coefficients) / 3),
        height = 1050 * ceiling(length(Model$coefficients) / 3),
        units = "px",
        res = 800
      )
      par(mfrow = c(ceiling(length(
        Model$coefficients
      ) / 3), 3))
      gam::plot.Gam(x = Model,
               residuals = F,
               ask = F)
      dev.off()
    } else{
      png(
        file.path(DirV, paste0("ResponseCurves_", spN, ".png")),
        width = 3200,
        height = 3200,
        units = "px",
        res = 800
      )
      par(mfrow = c(ceiling(length(
        Model$var.summary
      ) / 3), 3))
      visreg::visreg(
        Model,
        scale = "response",
        rug = F,
        line = list(lwd = 1),
        plot = T
      )
      dev.off()
    }
  }
  
  #Gaussian----
  if (Algorithm %in% "GAU") {
    V_IMP <- caret::filterVarImp(SpDataT[, VarColT], Outcome, nonpara = FALSE)
    V_IMP <- round(V_IMP / sum(V_IMP), 3)
    V_IMP <-
      cbind(
        Sp = spN,
        Algorithm = Algorithm,
        Variables = row.names(V_IMP),
        Importance = V_IMP
      )
    row.names(V_IMP) <- NULL
    if (file.exists(paste(DirV, "/VariableImportance.txt", sep = ""))) {
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    #Response Curves
    png(
      file.path(DirV, paste0("ResponseCurves_", spN, ".png")),
      width = 3200,
      height = 3200,
      units = "px",
      res = 800
    )
    plot(Model)
    dev.off()
  }
}
