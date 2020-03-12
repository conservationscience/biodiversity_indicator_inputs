##Needs R version 3.3.2 for sds package to work (but .nc versions of functions 
# run on package ncdf4 and should work with later versions)

# https://github.com/conservationscience/model_outputs_to_indicator_inputs.git

# cd C:\\Users\\ssteven\\Dropbox\\Deakin\\Serengeti-analysis\\model_outputs_to_indicator_inputs


# install.packages("packrat")
# install.packages("sds")
# install.packages("dplyr")
# install.packages("raster")
# install.packages("data.table")
# install.packages("readr")
# install.packages("ncdf4")
# install.packages("tidyverse")

# library(packrat)
# #library(sds)
# library(dplyr)
# library(raster)
# library(data.table)
# library(readr)
# library(ncdf4)
# library(tidyverse)
# library(stringr)

#### General Functions #####

#Function to print updates to the console about what the code is doing


Log <- function(...) { 
  if(getOption('MPverbose', TRUE)) cat(...)
}

#List model output files depending on type (cells, grids, massbins, 
#cohort properties, treeheightbins, various text etc)

ListCellOutputFiles <- function(resultsDir){
  
  files <- dir(resultsDir)
  files <- files[grep("BasicOutputs",files)]
  files <- files[grep("Cell",files)]
  
  return(files)
  
}

ListGridOutputFiles <- function(resultsDir){
  
  files <- dir(resultsDir)
  files <- files[grep("GridOutputs",files)]
  
  return(files)
}

ListMassBinsFiles <- function(resultsDir){
  
  files<-dir(resultsDir)
  files<-files[grep("MassBinsOutputs",files)]
  files<-files[grep("Cell",files)]
  
  return(files)
  
}

ListCohortPropertyFiles <- function(resultsDir){
  
  files <- dir(resultsDir)
  files <- files[grep("TrackedCohortProperties",files)]
  
  return(files)
}

ListCohortOutputFiles <- function(resultsDir){
  
  files <- dir(resultsDir)
  files <- files[grep("TrackedCohortsOutputs",files)]
  
  return(files)
}

ListTreeHeightBinsOutputsFiles <- function(resultsDir){
  
  files <- dir(resultsDir)
  files <- files[grep("TreeHeightBinsOutputs",files)]
  
  return(files)
}

ListExtinctionFiles <- function(resultsDir){
  
  files <- dir(resultsDir, pattern = "Extinction")
  return(files)
}

ListGrowthFiles <- function(resultsDir){
  
  files <- dir(resultsDir, pattern = "Growth")
  return(files)
}

ListMaturityFiles <- function(resultsDir){
  
  files <- dir(resultsDir, pattern = "Maturity")
  return(files)
}

ListMetabolismFiles <- function(resultsDir){
  
  files <- dir(resultsDir, pattern = "Metabolism")
  return(files)
}

ListMortalityFiles <- function(resultsDir){
  
  files <- dir(resultsDir, pattern = "Mortality")
  return(files)
}

ListNewCohortFiles <- function(resultsDir){
  
  files <- dir(resultsDir, pattern = "NewCohorts")
  return(files)
}

ListStateFiles <- function(resultsDir){
  
  files <- dir(resultsDir, pattern = "State")
  return(files)
}

ListTrophicFlowsFiles <- function(resultsDir){
  
  files <- dir(resultsDir, pattern = "TrophicFlows")
  return(files)
}

ListModelInputFiles <- function(resultsDir){
  
  files <- dir(resultsDir, pattern = "\\.csv$")
  return(files)
}

#create dimensons for a plot based on number of cells (nCells argument, can be obtained from nlength(files))

.plotDims<-function(nCells){
  widths<-list('1'=8.5/2.54,
               '2'=12.5/2.54,
               '3'=17.5/2.54,
               '4'=12.5/2.54,
               '5'=17.5/2.54,
               '6'=17.5/2.54,
               '7'=17.5/2.54,
               '8'=17.5/2.54,
               '9'=17.5/2.54,
               '10'=17.5/2.54,
               '11'=17.5/2.54,
               '12'=17.5/2.54
  )
  heights<-list('1'=6/2.54,
                '2'=6/2.54,
                '3'=6/2.54,
                '4'=12/2.54,
                '5'=12/2.54,
                '6'=12/2.54,
                '7'=18/2.54,
                '8'=18/2.54,
                '9'=18/2.54,
                '10'=24/2.54,
                '11'=24/2.54,
                '12'=24/2.54
  )
  
  stopifnot(nCells %in% names(widths))
  stopifnot(nCells %in% names(heights))
  
  return(list(width=widths[nCells][[1]],height=heights[nCells][[1]]))
  
}

gridArrange<-function(nCells){
  arrange<-list('1'=c(1,1),
                '2'=c(1,2),
                '3'=c(1,3),
                '4'=c(2,2),
                '5'=c(2,3),
                '6'=c(2,3),
                '7'=c(3,3),
                '8'=c(3,3),
                '9'=c(3,3),
                '10'=c(4,3),
                '11'=c(4,3),
                '12'=c(4,3)
  )
  
  stopifnot(nCells %in% names(arrange))
  
  return(arrange[nCells][[1]])
}


#Plots different functional groups over time (using ncdf4 package)

PlotTemporal.nc<-function(resultsDir,plotName,outDir=NULL,
                       label=NULL,
                       gridSimulation=FALSE,
                       whichCells = NULL,
                       vars=c("autotroph biomass density",
                              "herbivore biomass density",
                              "omnivore biomass density",
                              "carnivore biomass density"),
                       cols=c("#1b9e77",
                              "#66a61e",
                              "#7570b3",
                              "#d95f02"),
                       xlab="Years",ylab=plotName,
                       plotConfidence=TRUE,
                       confidenceInterval=95,
                       returnResults = FALSE){
  
  stopifnot(length(vars)==length(cols)) #checks if vars argument and cols arguments specified in the function arguments are the same length
  
  if ((gridSimulation) & (is.null(whichCells))){ #If you are plotting output from a grid simulation and you haven't specified cells as 'whichCells', print error msg
    stop("Error, if a grid-based simulation, you must specify cells")
  }
  
  if (gridSimulation){
    vars <- gsub(" biomass density","biomass density",vars) #look at the strings in object vars, and replace biomass density with a space in front of it,
    #with just biomass density, no space.
  }
  
  PlantBiomassVariables<-list( # Create list of plant variables with logical values
    "autotroph biomass density" = TRUE,
    "carnivore density" = FALSE,
    "carnivore biomass density" = FALSE,
    "herbivore density" = FALSE,
    "herbivore biomass density" = FALSE,
    "omnivore density" = FALSE,
    "omnivore biomass density" = FALSE,
    "Mean Trophic Level" = FALSE,
    "Max Trophic Index" = FALSE,
    "Generalist herbivoreBiomass" = FALSE,
    "Primary herbivoreBiomass" = FALSE,
    "Secondary herbivoreBiomass" = FALSE,
    "Plantation herbivoreBiomass" = FALSE,
    "Cropland herbivoreBiomass" = FALSE,
    "Pasture herbivoreBiomass" = FALSE,
    "Urban herbivoreBiomass" = FALSE,
    "Generalist omnivoreBiomass" = FALSE,
    "Primary omnivoreBiomass" = FALSE,
    "Secondary omnivoreBiomass" = FALSE,
    "Plantation omnivoreBiomass" = FALSE,
    "Cropland omnivoreBiomass" = FALSE,
    "Pasture omnivoreBiomass" = FALSE,
    "Urban omnivoreBiomass" = FALSE,
    "Generalist carnivoreBiomass" = FALSE,
    "Primary carnivoreBiomass" = FALSE,
    "Secondary carnivoreBiomass" = FALSE,
    "Plantation carnivoreBiomass" = FALSE,
    "Cropland carnivoreBiomass" = FALSE,
    "Pasture carnivoreBiomass" = FALSE,
    "Urban carnivoreBiomass" = FALSE,
    "Generalist herbivoreDensity" = FALSE,
    "Primary herbivoreDensity" = FALSE,
    "Secondary herbivoreDensity" = FALSE,
    "Plantation herbivoreDensity" = FALSE,
    "Cropland herbivoreDensity" = FALSE,
    "Pasture herbivoreDensity" = FALSE,
    "Urban herbivoreDensity" = FALSE,
    "Generalist omnivoreDensity" = FALSE,
    "Primary omnivoreDensity" = FALSE,
    "Secondary omnivoreDensity" = FALSE,
    "Plantation omnivoreDensity" = FALSE,
    "Cropland omnivoreDensity" = FALSE,
    "Pasture omnivoreDensity" = FALSE,
    "Urban omnivoreDensity" = FALSE,
    "Generalist carnivoreDensity" = FALSE,
    "Primary carnivoreDensity" = FALSE,
    "Secondary carnivoreDensity" = FALSE,
    "Plantation carnivoreDensity" = FALSE,
    "Cropland carnivoreDensity" = FALSE,
    "Pasture carnivoreDensity" = FALSE,
    "Urban carnivoreDensity" = FALSE,
    "Primary Landuse Affinity" = FALSE,
    "Secondary Landuse Affinity" = FALSE,
    "Plantation Landuse Affinity" = FALSE,
    "Cropland Landuse Affinity" = FALSE,
    "Pasture Landuse Affinity" = FALSE,
    "Urban Landuse Affinity" = FALSE
    
  )
  
  if (gridSimulation){ # if plotting a grid simulation
    names(PlantBiomassVariables) <- lapply( #rename plant biomass variables
      names(PlantBiomassVariables),function(x) return(
        gsub(" biomass density","biomass density",x))) #remove space from before "biomass"
  }
  
  if (gridSimulation){
    
    LogVariables<-list( #ifplotting a grid simulation, make a list of log variables with logical values
      "Abundance density" = TRUE,
      "Biomass density" = TRUE,
      "autotroph biomass density" = TRUE,
      "carnivore abundance density" = TRUE,
      "carnivore biomass density" = TRUE,
      "herbivore abundance density" = TRUE,
      "herbivore biomass density" = TRUE,
      "omnivore abundance density" = TRUE,
      "omnivore biomass density" = TRUE,
      "Mean Trophic Level" = FALSE,
      "Max Trophic Index" = FALSE
      
    )
    
    GramsBiomassVariables<-list( #ifplotting a grid simulation, make a list of biomass variables with logical values
      "Abundance density" = FALSE,
      "Biomass density" = TRUE,
      "autotroph biomass density" = TRUE,
      "carnivore abundance density" = FALSE,
      "carnivore biomass density" = TRUE,
      "herbivore abundance density" = FALSE,
      "herbivore biomass density" = TRUE,
      "omnivore abundance density" = FALSE,
      "omnivore biomass density" = TRUE,
      "Mean Trophic Level" = FALSE,
      "Max Trophic Index" = FALSE
    )
    
    stopifnot(all(vars %in% names(LogVariables))) #check vars entered as arguments match names in the lists in the functions
    stopifnot(all(vars %in% names(GramsBiomassVariables)))
  }
  
  stopifnot(all(vars %in% names(PlantBiomassVariables)))
  
  if ("SimulationControlParameters.csv" %in% dir(resultsDir)){ # check if results folder holds a simulation parameters file
    initialization <- read.csv(paste(resultsDir,"/SimulationControlParameters.csv",sep="")) #save simulation parameters used to initialize model as object
  } else {
    initialization <- read.csv(paste(resultsDir,"/EcosystemModelInitialisation.csv",sep=""))
  }
  
  Log("Finding Madingley output files\n")
  if (gridSimulation){
    files<-ListGridOutputFiles(resultsDir)
    file.paths <- paste(resultsDir,files,sep = "")
  } else {
    files<-ListCellOutputFiles(resultsDir) 
    file.paths <- paste(resultsDir,files,sep = "")#save the names of cell output files if it is not a grid simulation?
  }
  
  if (!gridSimulation){ #if this isn't a grid simulation, and therefore whichCells is null
    if(!is.null(whichCells)){
      files <- files[sapply(paste("Cell",whichCells-1,sep=""),FUN = function(x) return(grep(x,files)))] # returns the same results  as the previous list cells??
    }
    
    # Find the unique cells in these simulations
    cells.re<-regexpr("Cell[0-9]+",files) #I THINK this searches through the files and selectes ones with the word cells in them
    cells<-as.list(unique(substr(files,cells.re,cells.re+
                                   attr(cells.re,"match.length")-1))) # makes a list of cell names from 1 to length of cells simulated (ie removes duplicates created from multiple simulations)
    
  }
  
  if (gridSimulation) {
    # Find the simulation numbers
    sims.re<-regexpr("_[0-9]+",files) # searches for files with word sims in them?
    sims<-as.list(unique(substr(files,sims.re,sims.re+
                                  attr(sims.re,"match.length")-1))) #returns _0_ ??? Not sure if you'd get something different if there were multiple sims?
    
  } else {
    # Find the simulation numbers
    sims.re<-regexpr("_[0-9]+_",files) # Might need to update this if you have >10 sims?
    sims<-as.list(unique(substr(files,sims.re,sims.re+
                                  attr(sims.re,"match.length")-1))) # as above
  }
  
  if(is.null(label)){
    label<-unique(substr(files,1,sims.re-1)) #makes a label with file and simulation number?
    stopifnot(length(label)==1) #checks the label is only one word long
    label<-label[1] #subsets label to only include first word
  } else {
    label <- paste("BasicOutputs_",label,sep="")
  }
  
  if (!gridSimulation) Log(paste("Found results for ",length(cells)," cells\n",sep="")) #prints number of cells to the console
  Log(paste("Found results for ",length(sims)," simulations\n",sep="")) #prints number of simulations to the console
  
  Log("Getting basic information about simulations\n")
  # if (gridSimulation){
  #   sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sims[1],
  #                   ".nc",sep="")
  # } else {
  #   sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sims[1],cells[1], # creates path to sds output
  #                   ".nc",sep="")
  # }
  # 
  data <- nc_open(file.paths)
  
  times <- ncvar_get(data,"Time step") #gets monthly timesteps
  
  y <- times/12 #converts to number of years
  
  y <- ceiling(y) #idk what this part is for
  
  years <- unique(y) #maybe to simplify into one time series if you have multiple cells and sims?
  
  Log("Initializing plot\n")
  if (gridSimulation){
    dims <- .plotDims(length(whichCells))
  } else {
    dims<-.plotDims(length(cells))
  }
  
  if(!is.null(outDir)){
    pdf(paste(outDir,plotName,".pdf",sep=""),
        width = dims$width,height = dims$height) #this creates a pdf file in the output directory but can't seem to open it - perhaps this is a setup file?
  }
  
  if (gridSimulation) cells <- as.list(1:length(whichCells))
  
  par(mfrow=gridArrange(length(cells))) #setting up parameters for plotting
  par(mar=c(3,3.5,0.5,8.5))
  par(tck=-0.01)
  par(las=1)
  
  Log("Plotting\n")
  ret <- list()
  
  lapply(cells,FUN=function(cell){ 
    
    # Create matrices to hold the results for each specified variable
    allResults<-list()
    for (i in 1:length(vars)){
      allResults[[i]]<-matrix(data = NA,nrow = length(sims),
                              ncol = length(times))
    }
    names(allResults)<-vars
  
    # Loop over simulations in the ensemble
    s<-1
    for (sim in sims){
      if (gridSimulation){
        sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sim,".nc",sep="")
      } else {
        sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sim,cell,
                        ".nc",sep="")
      }
      
    # Populate the results matrices
      for (var in vars){
        if (gridSimulation){
          allResults[var][[1]][s,]<-ncvar_get(data,var)[,whichCells[[cell]][1],whichCells[[cell]][2]]
        } else {
          allResults[var][[1]][s,]<-ncvar_get(data,var)
          
        }
      }
      s<-s+1
    }
    
    for (var in vars){
      if (gridSimulation){
        if (LogVariables[var][[1]]){
          allResults[[var]] <- exp(allResults[[var]])-1
        }
        if (GramsBiomassVariables[var][[1]]){
          allResults[[var]] <- allResults[[var]]/1000.0
        }
        
        return(allResults)
        
      }
      # if (PlantBiomassVariables[var][[1]]){
      #  lp.ratio <- .GetLPRatios(resultsDir)[as.integer(gsub("Cell","",cell))+1]
      #  allResults[[var]] <- allResults[[var]]/lp.ratio
      # }
    }
    
    yearAvgs<-lapply(allResults,FUN = function(resultsMatrix){
      return(
        apply(resultsMatrix,1,FUN = function(simResults){
          tapply(simResults,y,mean)
        })
      )
    })
    
    # For each variable get the median and confidence limits across simulations
    summaryValues<-lapply(yearAvgs,FUN = function(x,ret){
      medianValues<-apply(x,1,median)
      upperValues<-apply(x,1,quantile,
                         probs=c(1-(1-(confidenceInterval/100))/2),
                         na.rm=TRUE)
      lowerValues<-apply(x,1,quantile,
                         probs=c((1-(confidenceInterval/100))/2),
                         na.rm=TRUE)
      
      medianValues[medianValues==0]<-NA
      upperValues[upperValues==0]<-NA
      upperValues[is.na(medianValues)]<-NA
      lowerValues[lowerValues==0]<-NA
      lowerValues[is.na(medianValues)]<-NA
      
      return(list(median=medianValues,
                  lower=lowerValues,
                  upper=upperValues))
    })
    
    maxVal<-max(unlist(lapply(summaryValues,function(x){
      return(max(x$upper,na.rm=TRUE))}
    )))
    minVal<-min(unlist(lapply(summaryValues,function(x){
      return(min(x$lower,na.rm=TRUE))
    })))
    
    par(mgp=c(1.6,0.3,0))
    plot(999999999999999,999999999999999,xlim=range(years),
         ylim=c(minVal,maxVal),log="y",
         xlab=xlab,yaxt="n",ylab=NA,bty="l")
    par(mgp=c(2.6,0.3,0))
    axis(2)
    title(ylab=ylab)
    
    if (plotConfidence){
      mapply(FUN=function(summaryValues,cols){
        X.Vec<-c(years,max(years),rev(years),min(years))
        Y.Vec<-c(summaryValues$lower,tail(summaryValues$upper,1),
                 rev(summaryValues$upper),summaryValues$lower[1])
        X.Vec<-X.Vec[!is.na(Y.Vec)]
        Y.Vec<-Y.Vec[!is.na(Y.Vec)]
        polygon(X.Vec,Y.Vec,col=paste(cols,"33",sep=""),border=NA)
      },summaryValues,cols)
    }
    
    mapply(FUN = function(summaryValues,cols){
      points(years,summaryValues$median,type="l",col=cols)
    },summaryValues,as.list(cols))
    
    if (returnResults){
      ret[[cell]] <<- summaryValues
    }
    
    legend(x = range(years)[2]+diff(range(years))*0.07,
           y = 10^(log10(minVal)+0.7*(log10(maxVal)-log10(minVal))),
           legend = gsub(" biomass density","",vars),
           lty=1,col = cols,xpd=TRUE,bty = "n")
    
}
    ) 
  
  if (!is.null(outDir)) invisible(dev.off())
  
  if (returnResults){
    return(ret)
  }
  
}


#Plots different functional groups over time using sds package

PlotTemporal.sds <-function(resultsDir,plotName,outDir=NULL,
                       label=NULL,
                       gridSimulation=FALSE,
                       whichCells = NULL,
                       vars=c("autotroph biomass density",
                              "herbivore biomass density",
                              "omnivore biomass density",
                              "carnivore biomass density"),
                       cols=c("#1b9e77",
                              "#66a61e",
                              "#7570b3",
                              "#d95f02"),
                       xlab="Years",ylab=plotName,
                       plotConfidence=TRUE,
                       confidenceInterval=95,
                       returnResults = FALSE){
  
  stopifnot(length(vars)==length(cols)) #checks if vars argument and cols arguments specified in the function arguments are the same length
  
  if ((gridSimulation) & (is.null(whichCells))){ #If you are plotting output from a grid simulation and you haven't specified cells as 'whichCells', print error msg
    stop("Error, if a grid-based simulation, you must specify cells")
  }
  
  if (gridSimulation){
    vars <- gsub(" biomass density","biomass density",vars) #look at the strings in object vars, and replace biomass density with a space in front of it,
    #with just biomass density, no space.
  }
  
  PlantBiomassVariables<-list( # Create list of plant variables with logical values
    "autotroph biomass density" = TRUE,
    "carnivore density" = FALSE,
    "carnivore biomass density" = FALSE,
    "herbivore density" = FALSE,
    "herbivore biomass density" = FALSE,
    "omnivore density" = FALSE,
    "omnivore biomass density" = FALSE,
    "Mean Trophic Level" = FALSE,
    "Max Trophic Index" = FALSE,
    "Generalist herbivoreBiomass" = FALSE,
    "Primary herbivoreBiomass" = FALSE,
    "Secondary herbivoreBiomass" = FALSE,
    "Plantation herbivoreBiomass" = FALSE,
    "Cropland herbivoreBiomass" = FALSE,
    "Pasture herbivoreBiomass" = FALSE,
    "Urban herbivoreBiomass" = FALSE,
    "Generalist omnivoreBiomass" = FALSE,
    "Primary omnivoreBiomass" = FALSE,
    "Secondary omnivoreBiomass" = FALSE,
    "Plantation omnivoreBiomass" = FALSE,
    "Cropland omnivoreBiomass" = FALSE,
    "Pasture omnivoreBiomass" = FALSE,
    "Urban omnivoreBiomass" = FALSE,
    "Generalist carnivoreBiomass" = FALSE,
    "Primary carnivoreBiomass" = FALSE,
    "Secondary carnivoreBiomass" = FALSE,
    "Plantation carnivoreBiomass" = FALSE,
    "Cropland carnivoreBiomass" = FALSE,
    "Pasture carnivoreBiomass" = FALSE,
    "Urban carnivoreBiomass" = FALSE,
    "Generalist herbivoreDensity" = FALSE,
    "Primary herbivoreDensity" = FALSE,
    "Secondary herbivoreDensity" = FALSE,
    "Plantation herbivoreDensity" = FALSE,
    "Cropland herbivoreDensity" = FALSE,
    "Pasture herbivoreDensity" = FALSE,
    "Urban herbivoreDensity" = FALSE,
    "Generalist omnivoreDensity" = FALSE,
    "Primary omnivoreDensity" = FALSE,
    "Secondary omnivoreDensity" = FALSE,
    "Plantation omnivoreDensity" = FALSE,
    "Cropland omnivoreDensity" = FALSE,
    "Pasture omnivoreDensity" = FALSE,
    "Urban omnivoreDensity" = FALSE,
    "Generalist carnivoreDensity" = FALSE,
    "Primary carnivoreDensity" = FALSE,
    "Secondary carnivoreDensity" = FALSE,
    "Plantation carnivoreDensity" = FALSE,
    "Cropland carnivoreDensity" = FALSE,
    "Pasture carnivoreDensity" = FALSE,
    "Urban carnivoreDensity" = FALSE,
    "Primary Landuse Affinity" = FALSE,
    "Secondary Landuse Affinity" = FALSE,
    "Plantation Landuse Affinity" = FALSE,
    "Cropland Landuse Affinity" = FALSE,
    "Pasture Landuse Affinity" = FALSE,
    "Urban Landuse Affinity" = FALSE
    
  )
  
  if (gridSimulation){ # if plotting a grid simulation
    names(PlantBiomassVariables) <- lapply( #rename plant biomass variables
      names(PlantBiomassVariables),function(x) return(
        gsub(" biomass density","biomass density",x))) #remove space from before "biomass"
  }
  
  if (gridSimulation){
    
    LogVariables<-list( #ifplotting a grid simulation, make a list of log variables with logical values
      "Abundance density" = TRUE,
      "Biomass density" = TRUE,
      "autotroph biomass density" = TRUE,
      "carnivore abundance density" = TRUE,
      "carnivore biomass density" = TRUE,
      "herbivore abundance density" = TRUE,
      "herbivore biomass density" = TRUE,
      "omnivore abundance density" = TRUE,
      "omnivore biomass density" = TRUE,
      "Mean Trophic Level" = FALSE,
      "Max Trophic Index" = FALSE
      
    )
    
    GramsBiomassVariables<-list( #ifplotting a grid simulation, make a list of biomass variables with logical values
      "Abundance density" = FALSE,
      "Biomass density" = TRUE,
      "autotroph biomass density" = TRUE,
      "carnivore abundance density" = FALSE,
      "carnivore biomass density" = TRUE,
      "herbivore abundance density" = FALSE,
      "herbivore biomass density" = TRUE,
      "omnivore abundance density" = FALSE,
      "omnivore biomass density" = TRUE,
      "Mean Trophic Level" = FALSE,
      "Max Trophic Index" = FALSE
    )
    
    stopifnot(all(vars %in% names(LogVariables))) #check vars entered as arguments match names in the lists in the functions
    stopifnot(all(vars %in% names(GramsBiomassVariables)))
  }
  
  stopifnot(all(vars %in% names(PlantBiomassVariables)))
  
  if ("SimulationControlParameters.csv" %in% dir(resultsDir)){ # check if results folder holds a simulation parameters file
    initialization <- read.csv(paste(resultsDir,"/SimulationControlParameters.csv",sep="")) #save simulation parameters used to initialize model as object
  } else {
    initialization <- read.csv(paste(resultsDir,"/EcosystemModelInitialisation.csv",sep=""))
  }
  
  Log("Finding Madingley output files\n")
  if (gridSimulation){
    files<-ListGridOutputFiles(resultsDir)
  } else {
    files<-ListCellOutputFiles(resultsDir) #save the names of cell output files if it is not a grid simulation?
  }
  
  if (!gridSimulation){ #if this isn't a grid simulation, and therefore whichCells is null
    if(!is.null(whichCells)){
      files <- files[sapply(paste("Cell",whichCells-1,sep=""),FUN = function(x) return(grep(x,files)))] # returns the same results  as the previous list cells??
    }
    
    # Find the unique cells in these simulations
    cells.re<-regexpr("Cell[0-9]+",files) #I THINK this searches through the files and selectes ones with the word cells in them
    cells<-as.list(unique(substr(files,cells.re,cells.re+
                                   attr(cells.re,"match.length")-1))) # makes a list of cell names from 1 to length of cells simulated (ie removes duplicates created from multiple simulations)
    
  }
  
  if (gridSimulation) {
    # Find the simulation numbers
    sims.re<-regexpr("_[0-9]+",files) # searches for files with word sims in them?
    sims<-as.list(unique(substr(files,sims.re,sims.re+
                                  attr(sims.re,"match.length")-1))) #returns _0_ ??? Not sure if you'd get something different if there were multiple sims?
    
  } else {
    # Find the simulation numbers
    sims.re<-regexpr("_[0-9]+_",files) 
    sims<-as.list(unique(substr(files,sims.re,sims.re+
                                  attr(sims.re,"match.length")-1))) # as above
  }
  
  if(is.null(label)){
    label<-unique(substr(files,1,sims.re-1)) #makes a label with file and simulation number?
    stopifnot(length(label)==1) #checks the label is only one word long
    label<-label[1] #subsets label to only include first word
  } else {
    label <- paste("BasicOutputs_",label,sep="")
  }
  
  if (!gridSimulation) Log(paste("Found results for ",length(cells)," cells\n",sep="")) #prints number of cells to the console
  Log(paste("Found results for ",length(sims)," simulations\n",sep="")) #prints number of simulations to the console
  
  Log("Getting basic information about simulations\n")
  if (gridSimulation){
    sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sims[1],
                    ".nc",sep="")
  } else {
    sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sims[1],cells[1], # creates path to sds output
                    ".nc",sep="")
  }
  
  data <- open.sds(sds.path) #returns the sds key 
  
  times <- get.sds(data,"Time step") #gets monthly timesteps
  
  y <- times/12 #converts to number of years
  
  y <- ceiling(y) #idk what this part is for
  
  years <- unique(y) #maybe to simplify into one time series if you have multiple cells and sims?
  
  Log("Initializing plot\n")
  if (gridSimulation){
    dims <- .plotDims(length(whichCells))
  } else {
    dims<-.plotDims(length(cells))
  }
  
  if(!is.null(outDir)){
    pdf(paste(outDir,plotName,".pdf",sep=""),
        width = dims$width,height = dims$height) #this creates a pdf file in the output directory but can't seem to open it - perhaps this is a setup file?
  }
  
  if (gridSimulation) cells <- as.list(1:length(whichCells))
  
  par(mfrow=gridArrange(length(cells))) #setting up parameters for plotting
  par(mar=c(3,3.5,0.5,8.5))
  par(tck=-0.01)
  par(las=1)
  
  Log("Plotting\n")
  ret <- list()
  
  lapply(cells,FUN=function(cell){
    
    # Create matrices to hold the results for each specified variable
    allResults<-list()
    for (i in 1:length(vars)){
      allResults[[i]]<-matrix(data = NA,nrow = length(sims),
                              ncol = length(times))
    }
    names(allResults)<-vars
    
    # Loop over simulations in the ensemble
    s<-1
    for (sim in sims){
      if (gridSimulation){
        sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sim,".nc",sep="")
      } else {
        sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sim,cell,
                        ".nc",sep="")
      }
      
      data<-open.sds(sds.path)
      
      # Populate the results matrices
      for (var in vars){
        if (gridSimulation){
          allResults[var][[1]][s,]<-get.sds(data,var)[,whichCells[[cell]][1],whichCells[[cell]][2]]
        } else {
          allResults[var][[1]][s,]<-get.sds(data,var)
          
        }
      }
      s<-s+1
    }
    
    for (var in vars){
      if (gridSimulation){
        if (LogVariables[var][[1]]){
          allResults[[var]] <- exp(allResults[[var]])-1
        }
        if (GramsBiomassVariables[var][[1]]){
          allResults[[var]] <- allResults[[var]]/1000.0
        }
        
        return(allResults)
        
      }
      # if (PlantBiomassVariables[var][[1]]){
      #  lp.ratio <- .GetLPRatios(resultsDir)[as.integer(gsub("Cell","",cell))+1]
      #  allResults[[var]] <- allResults[[var]]/lp.ratio
      # }
    }
    
    yearAvgs<-lapply(allResults,FUN = function(resultsMatrix){
      return(
        apply(resultsMatrix,1,FUN = function(simResults){
          tapply(simResults,y,mean)
        })
      )
    })#, #ADDED COMMA
    
    # For each variable get the median and confidence limits across simulations
    summaryValues<-lapply(yearAvgs,FUN = function(x,ret){
      medianValues<-apply(x,1,median)
      upperValues<-apply(x,1,quantile,
                         probs=c(1-(1-(confidenceInterval/100))/2),
                         na.rm=TRUE)
      lowerValues<-apply(x,1,quantile,
                         probs=c((1-(confidenceInterval/100))/2),
                         na.rm=TRUE)
      
      medianValues[medianValues==0]<-NA
      upperValues[upperValues==0]<-NA
      upperValues[is.na(medianValues)]<-NA
      lowerValues[lowerValues==0]<-NA
      lowerValues[is.na(medianValues)]<-NA
      
      return(list(median=medianValues,
                  lower=lowerValues,
                  upper=upperValues))
    })#, #ADDED COMMA
    
    maxVal<-max(unlist(lapply(summaryValues,function(x){
      return(max(x$upper,na.rm=TRUE))}
    )))#, #ADDED COMMA
    minVal<-min(unlist(lapply(summaryValues,function(x){
      return(min(x$lower,na.rm=TRUE))
    })))#, #ADDED COMMA
    
    par(mgp=c(1.6,0.3,0))#, #ADDED COMMA
    plot(999999999999999,999999999999999,xlim=range(years),
         ylim=c(minVal,maxVal),log="y",
         xlab=xlab,yaxt="n",ylab=NA,bty="l")#, #ADDED COMMA
    par(mgp=c(2.6,0.3,0))#, #ADDED COMMA
    axis(2)#, #ADDED COMMA
    title(ylab=ylab)#, #ADDED COMMA
    
    if (plotConfidence){
      mapply(FUN=function(summaryValues,cols){
        X.Vec<-c(years,max(years),rev(years),min(years))
        Y.Vec<-c(summaryValues$lower,tail(summaryValues$upper,1),
                 rev(summaryValues$upper),summaryValues$lower[1])
        X.Vec<-X.Vec[!is.na(Y.Vec)]
        Y.Vec<-Y.Vec[!is.na(Y.Vec)]
        polygon(X.Vec,Y.Vec,col=paste(cols,"33",sep=""),border=NA)
      },summaryValues,cols)
    }
    
    mapply(FUN = function(summaryValues,cols){
      points(years,summaryValues$median,type="l",col=cols)
    },summaryValues,as.list(cols))#, #ADDED COMMA
    
    if (returnResults){
      ret[[cell]] <<- summaryValues
    }
    
    legend(x = range(years)[2]+diff(range(years))*0.07,
           y = 10^(log10(minVal)+0.7*(log10(maxVal)-log10(minVal))),
           legend = gsub(" biomass density","",vars),
           lty=1,col = cols,xpd=TRUE,bty = "n")
    
  }
  ) 
  
  if (!is.null(outDir)) invisible(dev.off())
  
  if (returnResults){
    return(ret)
  }
  
}

#Function to extract the data from the nc outputs for each cell simulation

extract.data <- function(file.paths, files){
  
  #Use the file paths to extract the sds keys you need to get the data
  data <- lapply(file.paths, open.sds)
  
  #Get a list of data content
  schema <- lapply(data, schema.sds)
  
  #Names list so we know which cell and simulation the output belongs to
  names(schema) <- file.paths
  
  #Get a list of variables these outputs contain
  schema1 <- schema[[1]]
  vars <- names(schema1$variables)
  
  #Create empty parent.list for each cell+simulation, each containing a variable.list
  #containing the output data for 29 variables
  parent.list <- list()
  variable.list <- list()
  
  for (j in seq_along(data)) {#Loop over each cell for each simulation
    cell.sim <- data[[j]]
    names(cell.sim) <- files[[j]]####
    for (i in seq_along(vars)) {
      variable.list[[i]] <- get.sds(cell.sim,vars[[i]]) #Loop over each variable
    }
    names(variable.list) <- vars #make sure variables are names
    parent.list[[j]] <- variable.list
  }
  
  names(parent.list) <- files#### #make sure each cell simulation is named
  return(parent.list)
}

#Version 2 that uses ncdf4 package not sds

extract.nc.data <- function(file.paths, files){
  
  #Use the file paths to extract the sds keys you need to get the data
  data <- lapply(file.paths, nc_open)
  
  #Get a list of data content
  #schema <- lapply(data, schema.sds)
  
  #Names list so we know which cell and simulation the output belongs to
  names(data) <- files
  
  #Get a list of variables these outputs contain
  schema1 <- data[[1]]
  vars <- names(schema1$var)
  
  #Create empty parent.list for each cell+simulation, each containing a variable.list
  #containing the output data for 29 variables
  parent.list <- list()
  variable.list <- list()
  
  for (j in seq_along(data)) {#Loop over each cell for each simulation
    cell.sim <- data[[j]]
    #names(cell.sim) <- files[[1]]####
    for (i in seq_along(vars)) {
      variable.list[[i]] <- ncvar_get(cell.sim,vars[[i]]) #Loop over each variable
    }
    names(variable.list) <- vars #make sure variables are names
    parent.list[[j]] <- variable.list
  }
  
  names(parent.list) <- files#### #make sure each cell simulation is named
  return(parent.list)
  
  nc_close(data) #Check this doesn't stop the function from working
}



# Function to pull out relevant elements from raw madingley outputs

filter.by.pattern <- function(pattern, your.list){
  
    names(your.list) %>% 
    str_detect(pattern) %>%
    keep(your.list, .)
  
}

remove.burn.in <- function(data, burnin) {
  
  data[,(burnin + 1):ncol(data)]
  
}

# Function to pull out specific variables and remove burn-in period from Madingley MassBinsOutput

subset.massbins <- function(data, pattern, burnin) {
  
  Temp1 <- list()
  Temp2 <- list()
  Temp3 <- list()
  Output <- list()
  
  for (i in seq_along(data)){
    
    replicate <- data[[i]]
    
    for (j in seq_along(pattern)){
      
      Temp1[[j]] <- filter.by.pattern(pattern[j], replicate)
      
    }
    
    Temp2[[i]] <- Temp1
  }
  
  # Loop over replicates, and groups in replicates, to remove burn in timesteps
  
  for (i in seq_along(Temp2)){
    
    replicate <- Temp2[[i]]
    
    for (j in seq_along(replicate)){
      
      Temp3[[j]] <- remove.burn.in(replicate[[j]][[1]], burnin)
    }
    
    names(Temp3) <- pattern
    
    Output[[i]] <- Temp3
    
  }
  Output
  
}


## Output function - this function is dependent on the extract.data function.
#It searches the results directory, searches for output type you want (i.e. massbins,
#extinction etc., processes it and returns a data frame).  You can also use
#it with lapply() to process all files in the directory at once

output.function <- function(output.type, resultsDir) {
  
  if (output.type == "BasicOutputs") {
    
    files <- ListCellOutputFiles(resultsDir)
    file.paths <- paste(resultsDir,files,sep = "")
    results <- extract.nc.data(file.paths, files)
    names(results) <- files
  }
  
  else if (output.type == "GridOutputs") {
    files <- ListGridOutputFiles(resultsDir)
    results <- files
    names(results) <- files
    
  }
  
  else if (output.type == "MassBinsOutputs") {
    files <- ListMassBinsFiles(resultsDir)
    file.paths <- paste(resultsDir,files,sep = "")
    results <- extract.nc.data(file.paths, files)
    names(results) <- files
    
  }
  
  else if (output.type == "TrackedCohortOutputs") {
    files <- ListCohortOutputFiles(resultsDir)
    file.paths <- paste(resultsDir,files,sep = "")
    results <- extract.nc.data(file.paths, files)
    names(results) <- files
    
  }
  
  else if (output.type == "TreeHeightBinsOutputs") {
    files <- ListTreeHeightBinsOutputsFiles(resultsDir)
    file.paths <- paste(resultsDir,files,sep = "")
    results <- extract.nc.data(file.paths, files)
    names(results) <- files
    
  }
  
  else if (output.type == "Extinction") {
    
    files <- ListExtinctionFiles(resultsDir)
    
    if (length(files) == 0) {
      
      message(paste("no files of type", output.type, "exist in this results directory")) 
      
    } else {
      
      file.paths <- paste(resultsDir,files,sep = "")
      results <- lapply(file.paths, read_tsv)
      names(results) <- files
      
    }
  }
  
  else if (output.type == "Growth") {
    
    files <- ListGrowthFiles(resultsDir)
    
    if (length(files) == 0) {
      
      message(paste("no files of type", output.type, "exist in this results directory")) 
      
    } else {
      
      file.paths <- paste(resultsDir,files,sep = "")
      results <- lapply(file.paths, read_tsv)
      names(results) <- files
      
    }
  }
  
  else if (output.type == "Maturity") {
    
    files <- ListMaturityFiles(resultsDir)
    
    if (length(files) == 0) {
      
      message(paste("no files of type", output.type, "exist in this results directory")) 
      
    } else {
      
      file.paths <- paste(resultsDir,files,sep = "")
      results <- lapply(file.paths, read_tsv)
      names(results) <- files
      
    }
  }
  
  else if (output.type == "Metabolism") {
    files <- ListMetabolismFiles(resultsDir)
    
    if (length(files) == 0) {
      
      message(paste("no files of type", output.type, "exist in this results directory")) 
      
    } else {
      
      file.paths <- paste(resultsDir,files,sep = "")
      results <- lapply(file.paths, read_tsv)
      names(results) <- files
      
    }
  }
  
  else if (output.type == "Mortality") {
    files <- ListMortalityFiles(resultsDir)
    
    if (length(files) == 0) {
      
      message(paste("no files of type", output.type, "exist in this results directory")) 
      
    } else {
      
      file.paths <- paste(resultsDir,files,sep = "")
      results <- lapply(file.paths, read_tsv)
      names(results) <- files
     
    }
  }
  
  else if (output.type == "Newcohorts") {
    
    files <- ListNewCohortFiles(resultsDir)
    
    if (length(files) == 0) {
      
      message(paste("no files of type", output.type, "exist in this results directory")) 
      
    } else {
      
      file.paths <- paste(resultsDir,files,sep = "")
      results <- lapply(file.paths, read_tsv)
      names(results) <- files
      
    }
  }
  
  else if (output.type == "State") {
    
    files <- ListStateFiles(resultsDir)
    
    if (length(files) == 0) {
      
      message(paste("no files of type", output.type, "exist in this results directory")) 
      
    } else {
      
      file.paths <- paste(resultsDir,files,sep = "")
      results <- lapply(file.paths, read_tsv)
      names(results) <- files
      
    }
  }
  
  else if (output.type == "TrackedCohortProperties") {
    
    files <- ListCohortPropertyFiles(resultsDir)
    
    if (length(files) == 0) {
      
      message(paste("no files of type", output.type, "exist in this results directory")) 
      
    } else {
      
      file.paths <- paste(resultsDir,files,sep = "")
      results <- lapply(file.paths, read_tsv)
      names(results) <- files
      
    }
  }
  
  else if (output.type == "TrophicFlows") {
    
    files <- ListTrophicFlowsFiles(resultsDir)
    
    if (length(files) == 0) {
      
      message(paste("no files of type", output.type, "exist in this results directory")) 
      
    } else {
      
      file.paths <- paste(resultsDir,files,sep = "")
      results <- lapply(file.paths, read_tsv)
      names(results) <- files
      
    }
  }
  
  else if (output.type %in% c("CohortFunctionalGroupDefinitions",
                              "EcologicalParameters","EnvironmentalDataLayers",
                              "FileLocationParameters","OutputControlParameters",
                              "PlantCohortFunctionalGroupDefinitions","Scenarios",
                              "SimulationControlParameters","SpecificLocations",
                              "StockFunctionalGroupDefinitions")) {
    
    files <- ListModelInputFiles(resultsDir)
    results <- files
    names(results) <- files
  }
  
  else {
    
    message("You have not entered a valid output type for this results directory") 
    
  }
  
  #return(files) #print the file names of the output type you specified (change this to dataframes eventually)
  
  if (length(files > 0)) {
    
    cells.re<-regexpr("Cell[0-9]+",files) #I THINK this searches through the files and selectes ones with the word cells in them
    cells<-as.list(unique(substr(files,cells.re,cells.re+
                                   attr(cells.re,"match.length")-1)))
    
    sims.re<-regexpr("_[0-9]+_",files)
    sims<-as.list(unique(substr(files,sims.re,sims.re+
                                  attr(sims.re,"match.length")-1)))
    
    
    Log(paste("Found ",output.type, " results for ",length(cells)," cells\n",sep="")) #prints number of cells to the console
    Log(paste("Found ", output.type, " results for ",length(sims)," simulations\n",sep=""))
    
  } 
  #names(files) <- output.type
  
  if (length(files > 0)) {
    
    return(results)
    
  } else {
    
    rm(files) #This part doesn't work, still adds an empty 
    #character element
  }
  
}

#Save processed output

#Add system time to object name in case you process more than one set of outputs per day

save.and.date <- function(your.object, filename){
  
  time <- as.character(Sys.time())
  time <- gsub(":", "-", time)
  objectname <- paste(time, filename,".rds",sep = "")
  save(your.object, file = paste(output_file_path,objectname, sep = "/" ))
  
}

PlotCells <- function(inputs,outDir=NULL,simulation_number,map){
  
  require(rgdal)
  
  if (map == "Africa") {
    
    outline<-readOGR("N:/Quantitative-Ecology/Indicators-Project/Madingley_setup_data/MadingleyPlots/data","Africa",verbose = FALSE)
    
  } else {
    
    outline<-readOGR("N:/Quantitative-Ecology/Indicators-Project/Madingley_setup_data/MadingleyPlots/data","World",verbose = FALSE)
  }
  
  locations <- read.csv(paste(inputs,"/SpecificLocations.csv",sep=""))
  
  if (!is.null(outDir)){
    
    pdf(paste(outDir, "/", simulation_number, "_CellsMap.pdf",sep=""),
        width = 12.5/2.54,height=6.25/2.54)
  }
  
  
  par(mar=c(0,0,0,0))
  
  plot(outline)
  
  points(locations$Longitude,locations$Latitude,pch=16,col="red")
  
  text(locations$Longitude,locations$Latitude,1:dim(locations)[1],pos=2)
  
  if (!is.null(outDir)) invisible(dev.off())
  
}

save_plot <- function(inputs, outputs, simulation_number, map) {
  
  require(rgdal)
  
  existing_plots <- grepl("biomass.tiff", list.files(outputs))
  plotting_complete <- any(existing_plots = TRUE)
  
  if (!plotting_complete) {
  
  files <- list.files(inputs)
  scenario <- sub("BasicOutputs_", "", files[1])
  scenario <- sub("_0_Cell0.nc", "", scenario)
  
  # Save plot of coordinates of modelled grid cell
  PlotCells(inputs, outputs, simulation_number, map)
 
   #Set up tiff file save for plot
  
  plotName <- paste(simulation_number, "_", scenario, "_biomass",".tiff",sep="")
  tiff(file = (paste(outputs,plotName, sep = "/")), units ="in", width=10, height=5, res=200)
  
  # Plot biomass over time
  biomassPlot <- PlotTemporal.nc(inputs, label = scenario, "Temporal Biomass Plot")
  
  print(paste("plot of ",simulation_number, "_", scenario, " output biomass saved to output folder", sep = ""))

  #turn off devices
 
  dev.off()
  
  } else {
    
    folder_name <- basename(outputs)
    
    print(paste(folder_name, "outputs have already been plotted and saved", sep = " "))
    
  }
  
}


