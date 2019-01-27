# Setup file for survival power analysis


# package loading ---------------------------------------------------------

# Smart installer will check list of packages that are installed, install any
# necessary package that is missing, and load the library:

smartInstallAndLoad <- function(packageVector){
  for(i in 1:length(packageVector)){
    package <- packageVector[i]
    if(!package %in% rownames(installed.packages())){
      install.packages(packageVector[i],repos="http://cran.rstudio.com/",
                       dependencies=TRUE)
    }
  }
  lapply(packageVector, library, character.only = TRUE)
}

# Load and potentially install libraries:

smartInstallAndLoad(c('marked', 'tidyverse', 'stringr', 'lubridate'))

# options -----------------------------------------------------------------

# Options and assigning a library to conflicting functions:

options(stringsAsFactors = FALSE)

select <- 
  dplyr::select

summarize <-
  dplyr::summarize

rename <- 
  dplyr::rename

outPlotsDir <- 'outPlots/'


# functions to simulate encounter histories -------------------------------

# Function to simulate (apparent) survival for a given bird:

simulate_survival <-
  function(phiTrue, nYears) {
    survivalVector = 1
    for(i in 2:nYears){
      survivalVector[i] <-
        ifelse(
          survivalVector[i-1] == 0,
          0,
          rbinom(1, 1,phiTrue))
    }
    survivalVector
  }

# Function to simulate encounter histories for a given bird:

simulate_encounters <-
  function(survivalVector, pTrue) {
    encounterVector = 1
    for(i in 2:length(survivalVector)){
      encounterVector[i] <-
        ifelse(
          survivalVector[i] == 0,
          0,
          rbinom(1, 1,pTrue))
    }
    encounterVector
  }

# Function to generate random survival and encounter histories:

generateRandomCh <- 
  function(
    phiTrue, 
    pTrue, 
    nBirds, 
    nYears) {
    
    encounterList <-
      vector('list', length = nBirds)
    
    
    for(i in 1:nBirds){
      survivalVector <-
        simulate_survival(phiTrue, nYears)
      
      encounterVector <-
        simulate_encounters(survivalVector, pTrue)
      
      encounterList[[i]] <-
        tibble(
          birdId = i,
          survival = paste0(survivalVector, collapse = ''),
          encounter = paste0(encounterVector, collapse = '')
        )
    }
    bind_rows(encounterList)
  }

# Function to clean encounter records, leaving just ch and covariates:

clean_encounters <- 
  function(encounterFrame, simulated = TRUE){
    if(simulated == TRUE){
      encounterFrame <-
        group_ch %>%
        rename(ch = encounter) %>%
        select(-c(birdId, survival))
    } else {
      encounterFrame <-
        group_ch %>%
        rename(ch = survival) %>%
        select(-c(birdId, encounter))
    }
    as.data.frame(encounterFrame)
  }


# functions to generate simulated input values ----------------------------


simulateSurvival_group <- 
  function(pValues, startingPhi, phiDiffs, nBirds, yrs){
  # Make grid of all possible combinations of p and years:
  inputGrid <-
    expand.grid(
      pTrue = pValues,
      phiDiff = phiDiffs/2,
      nBirds = nBirds,
      nYears = yrs
    ) %>%
    mutate(
      phi_group1 = startingPhi + phiDiff,
      phi_group2 = startingPhi - phiDiff
    ) %>%
    select(-phiDiff)
  
  # Simulate for each row of the grid above:
  
  chList <- 
    vector('list', length = nrow(inputGrid))
  
  for(i in 1:nrow(inputGrid)){
    inputs <- inputGrid[i,]
    chList[[i]] <-
      bind_rows(
        generateRandomCh(
          phiTrue = inputs$phi_group1, 
          pTrue = inputs$pTrue, 
          nBirds = inputs$nBirds, 
          nYears = inputs$nYears) %>%
          mutate(group_var = 'group_a'),
        generateRandomCh(
          phiTrue = inputs$phi_group2, 
          pTrue = inputs$pTrue, 
          nBirds = inputs$nBirds, 
          nYears = inputs$nYears) %>%
          mutate(group_var = 'group_b'))
  }
  
}

# model fitting functions -------------------------------------------------


fitModels_group <- function(group_ch, simulated = TRUE) {
  
  # Processed and design data:
  
  encounterData <- 
    clean_encounters(group_ch, simulated) %>%
    as.data.frame %>%
    mutate(group_var = factor(group_var))
  
  procData <-
    encounterData %>%
    process.data
  
  designData <- 
    make.design.data(procData)
  
  # Formulas for phi and p:
  
  Phi.dot <- list(formula = ~1)
  Phi.group_var <- list(formula = ~group_var)
  p.dot <- list(formula = ~1)
  p.group_var <- list(formula = ~group_var)
  
  # Create model list and output:
  
  cml <- 
    create.model.list(c('Phi', 'p'))
  
  crm.wrapper(
    cml,
    data = procData, 
    ddl = designData, 
    hessian = TRUE,
    external = FALSE,
    accumulate = FALSE)
}


