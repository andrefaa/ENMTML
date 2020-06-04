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
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    #Response Curves
  }
  
  #Maxent----
  if (Algorithm %in% c("MXD", "MXS", "MLK")) {
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
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    #Response Curves
  }
  
  
  
  #Boosted Regression Tree----
  if (Algorithm %in% "BRT") {
    V_IMP <- gbm::relative.influence(Model, n.trees = Model$n.trees, scale = T)
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
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    #Response Curves
  }
  
  #Random Forests----
  if (Algorithm %in% "RDF") {
    V_IMP <- randomForest::importance(Model)
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
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    #Response Curves
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
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    # #Response Curves
  }
  
  #GLM,GAM----
  if (any(Algorithm %in% c("GLM", "GAM"))) {
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
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    #Response Curves
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
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t',
        col.names = F
      )
    } else{
      utils::write.table(
        V_IMP,
        paste(DirV, "/VariableImportance.txt", sep = ""),
        append = TRUE,
        quote = FALSE,
        sep = '\t'
      )
    }
    #Response Curves
  }
}
