  #########################################################################################################################
  #########################################################################################################################
  ###### PROJECT:        DBNs
  ###### NAME:           ValidateTM.R
  ###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
  ###### AFFILIATION:    Florida International University
  ###### DESCRIPTION:    Converts the graphml network into an adjacency matrix
  ######                 Also outputs an adjacency matrix with just the validated interactions
  ######                 Calculates the hypergeometric distribution values and other stats
  ######                 Does a binomial test
  #########################################################################################################################
  #########################################################################################################################
  

  rm(list=ls())
  library(scales)
  library(stringr)
  library("stats")
  library(poisbinom)
  library(plotrix)
  
  options("scipen"=100, "digits"=4)
  
  TMDirectionCounts = T
  TMSignCounts = F
  
  clinical = tolower(c("Week sample obtained"))
  folder = "" #"Demo" #"finalFiguresBootscore" #"finalFiguresBootscore"
  setwd(folder)

  nameOfNetwork = paste(folder,"Networks/HEGTM/TMGHCombinations_GeneReference/","human_ibd_microbiota_dbn_sample_noalignment_sr14d_top100x100_genes_reference_dbnIntraBoot.graphml",sep="")
  nameOfNetwork = paste(folder,"Networks/HEGTM/TMGHCombinations_GeneReference/","human_ibd_microbiota_dbn_sample_alignment_sr14d_top100x100_genes_reference_dbnIntraBoot.graphml",sep="")
  nameOfNetwork = paste(folder,"Networks/HEGTM/TMGHCombinations_GeneReference/","human_ibd_microbiota_genes_metabolites_dbn_sample_alignment_sr14d_top100x100_host_genes_reference_dbnIntraBoot.graphml",sep="")
  nameOfNetwork = paste(folder,"Networks/HEGTM/TMGHCombinations_GeneReference/","human_ibd_microbiota_genes_metabolites_dbn_sample_noalignment_sr14d_top100x100_host_genes_reference_dbnIntraBoot.graphml",sep="")
  
    network = tryCatch({
    readChar(nameOfNetwork, file.info(nameOfNetwork)$size)
  }, warning = function(e) {
    ""
  }, error = function(e) {
    ""
  })
  if (network == "")
    next
  


  
  #####################################################################
  ############# LOAD INFORMATION FROM  GRAPHML NETWORK ############# 
  #####################################################################
  
  #READ NETWORK
  network = gsub("t0", "ti", gsub("tn", "ti+1", tolower(unlist(strsplit(x =network, "\r?\n")))))
  network = gsub("\n<data key=\"key_bootscore\">[0-9]*\\.[0-9]*</data>\n","",network)
  newNetwork = network[grep(pattern = paste(".*target=\".*_ti\\+1\">","",sep=""), x = network, ignore.case = T)]
  newNetwork = gsub("s__", "", gsub("\\(", "\\\\(", gsub("\\)", "\\\\)", gsub("\\[", "\\\\]", gsub("\\[", "\\\\]", gsub("_ti", "", gsub("\\+1", "", gsub("\\[", "[", gsub("\\]", "]",  gsub("\">", "", gsub("m__", "", gsub("<edge id=\"e\\d*\" source=\"", "", newNetwork))))))))))))
  
  #Remove everything that has to do with genes
  newNetwork = newNetwork[setdiff(1:length(newNetwork),grep(".*g__.*",newNetwork))]
  newNetwork = newNetwork[setdiff(1:length(newNetwork),grep(".*hg.__.*",newNetwork))]
  
  #get the names of taxa
  namesOfTaxa = network[grep("<node id=\\\"s__.*",network)]
  namesOfTaxa = unique(gsub("_ti", "", gsub("\\+1", "", gsub("\\\"/>", "", gsub("<node id=\\\"", "", gsub("s__", "", namesOfTaxa))))))
  taxaNewNames = gsub("_"," ",namesOfTaxa)
  
  #get the names of metabolites
  namesOfMetabolites = network[grep("<node id=\\\"m__.*",network)]
  namesOfMetabolites = unique(gsub("_ti", "", gsub("\\+1", "", gsub("\\\"/>", "", gsub("<node id=\\\"", "", gsub("m__", "", namesOfMetabolites))))))
  namesOfMetabolites = gsub("_"," ",namesOfMetabolites)
  
  #get the weight of the edges
  edgeWeight = network[grep("<data key=\\\"key_weight.*",network)]
  edgeWeight = as.numeric(gsub("\"", "", gsub("<data key=\\\"key_weight\\\">", "", gsub("</data>", "", edgeWeight)))) 
  
  #get the bootscore of the edges
  edgeBootscore = network[grep("<data key=\\\"key_bootscore.*",network)]
  edgeBootscore = as.numeric(gsub("\"", "", gsub("<data key=\\\"key_bootscore\\\">", "", gsub("</data>", "", edgeBootscore)))) 
  
  
  #Extract the information from the network
  pairs = c()
  for (i in 1:length(newNetwork)){
    edge =  c(trimws(strsplit(newNetwork[i],"\" target=\"")[[1]], which ="both"),edgeWeight[i],edgeBootscore[i])
    pairs = rbind(pairs,edge)
  }
  uniqueNames = gsub("\\\\", "", gsub("_unclassified", "", gsub("_noname", "", unique(c(pairs[,1],pairs[,2])))))
  
  #Replace the new taxa names into the final all names
  for (i in 1:length(namesOfTaxa)){
    for (j in 1:length(uniqueNames)){
      if (namesOfTaxa[i] == uniqueNames[j]){
        uniqueNames[j] = taxaNewNames[i]
      }
    }
  }

  #Replace the new taxa genus names into the map pairs fil
  for (i in 1:length(namesOfTaxa)){
    for (j in 1:nrow(pairs)){
      if (namesOfTaxa[i] == as.character(pairs[j,1])){
        pairs[j,1] = taxaNewNames[i]
      }
      if (namesOfTaxa[i] == as.character(pairs[j,2])){
        pairs[j,2] = taxaNewNames[i]
      }
    }
  }

  
  #####################################################################
  ############# CREATE ADJACENCY MATRIX FROM GRAPHML ############# 
  #####################################################################
  
  adjacencyMatrix <- function(uniqueNames, pairs){
    output = matrix(0, nrow = length(uniqueNames), ncol = length(uniqueNames)) 
    colnames(output) = rownames(output) = uniqueNames
    
    #Fill in the adjacency matrix
    for (i in 1:nrow(pairs)){
      # output[grep(paste("^",pairs[i,1],"$",sep=""),uniqueNames),grep(paste("^",pairs[i,2],"$",sep=""),uniqueNames)] = 1 + output[grep(paste("^",pairs[i,1],"$",sep=""),uniqueNames),grep(paste("^",pairs[i,2],"$",sep=""),uniqueNames)]
      output[grep(paste("^",pairs[i,1],sep=""),uniqueNames),grep(paste("^",pairs[i,2],sep=""),uniqueNames)] = pairs[i,3]
    }
    #output graphml as csv
    # write.table(output,paste(nameOfNetwork,".csv",sep=""),sep=",",quote=T,col.names=NA)
    return (output)
  }
  
  
  #####################################################################
  ############# CREATE VALIDATED  ############# 
  #####################################################################
  
    ######### FILTER OUT EVERYTHING BUT TAXA - METABOLITES
  getJustTMInteractions = function(TMDirectionCounts,pairs, taxaNewNames, namesOfMetabolites ){
    pairsTM = c() #Only focus on taxa-metabolite
    if (TMDirectionCounts){
      for (i in nrow(pairs):1){
        if ((pairs[i,1] %in% taxaNewNames && pairs[i,2] %in% namesOfMetabolites && pairs[i,3] > 0)){
          pairsTM = rbind(pairsTM,pairs[i,])
        }
      }
    }else{
      for (i in nrow(pairs):1){
        if ((pairs[i,1] %in% taxaNewNames && pairs[i,2] %in% namesOfMetabolites) || pairs[i,1] %in% namesOfMetabolites && pairs[i,2] %in% taxaNewNames){
          pairsTM = rbind(pairsTM,pairs[i,])
        }
      }
    }
    pairsTM = unique(pairsTM)
  return(pairsTM)
  }
  
    
    
  # Calculate the metabolites that have interactions with our taxa
  retrieveRelevantMet = function(validTM, namesOfMetabolites ){
    relevantMet = c()
    for (metList in validTM$ValidMetabolitesHMDB){
      relevantMet = paste(relevantMet,metList,sep="|")
    }
    relevantMet = unique(strsplit(relevantMet,("\\|"))[[1]])
    relevantMet = relevantMet[-which(relevantMet %in% "")]
    #Intersect with the metabolites in the network
    relevantMet = intersect(relevantMet,namesOfMetabolites)
    return (relevantMet)
  }
  
  
    #VALIDATE INTERACTIONS USING KEGG DATABASE FILE
  computeValidatedInteractions = function(output,uniqueNames, clinical, TMDirectionCounts, validTM, relevantMet){
    found = c()
    numberOfTMValidated = 0
    numberOfTMNotValidated = 0
    interactionsValidated = c()
    interactionsNOTValidated = c()
    taxaNotFound = c()
    metabolitesNotInteresting = c()
    outputValidated = output
    totalTMInteractions = 0
    for (i in 1:nrow(output)){
      for (j in 1:ncol(output)){
        if (output[i,j] != 0){
          # if (uniqueNames[j] %in% "hmdb03357"){
          #   print(uniqueNames[i])
          #   print(uniqueNames[j])
          #   print("here")
          # }
  
          
          #If any one is clinical, no interested
          if (uniqueNames[i] %in% clinical || uniqueNames[j] %in% clinical){
            # print(paste(uniqueNames[i], "->" , uniqueNames[i] , "NOT Validated"))
            next
          }
          #If it is taxa-taxa, continue
          if (!grepl("hmdb", uniqueNames[i]) && !grepl("hmdb", uniqueNames[j])){
            # print(paste(uniqueNames[i], "->" , uniqueNames[i] , "NOT Validated"))
            next
          }
          #If it is metabolite-metabolite, continue
          if (grepl("hmdb", uniqueNames[i]) && grepl("hmdb", uniqueNames[j])){
            # print(paste(uniqueNames[i], "->" , uniqueNames[i] , "NOT Validated"))
            next
          }
          #Assign the taxa and the metabolite to the correct variable
          if(grepl("hmdb", uniqueNames[i])){
            metabolite = uniqueNames[i]
            taxa = uniqueNames[j]
            if(TMDirectionCounts ){
              print("Wrong direction")
              print(paste(uniqueNames[i], "->" , uniqueNames[j] , "NOT Validated"))
              next
            }
          }else{
            metabolite = uniqueNames[j]
            taxa = uniqueNames[i]
          }
    
          coefficient = as.numeric(output[i,j])
          if (TMDirectionCounts && TMSignCounts && coefficient <0){
            print("Wrong sign")
            print(paste(uniqueNames[i], "->" , uniqueNames[j] , "NOT Validated"))
            next
          }
          
          #Check if the taxa is found
          indexTaxa = grep(tolower(taxa), tolower(validTM$SpeciesName)) #CAREFUL WITH FRAGMENTFOUNDINKEGG
          if (length(indexTaxa) == 0){
            # print(paste("TAXA NOT FOUND IN KEGG DATABASE:",taxa))
            taxaNotFound = c(taxaNotFound,taxa)
            next
          }
          
          #Check if the metabolite is relevant
          if (!metabolite %in% relevantMet){
            # print(paste("Metabolite not relevant: ",metabolite))
            metabolitesNotInteresting = c(metabolitesNotInteresting,metabolite)
            next
          }
  
          #see if the metabolite is related to the taxa in the database
          if (length(grep(metabolite, validTM$ValidMetabolitesHMDB[indexTaxa]))>0){
            print(paste("VALIDATED: ",uniqueNames[i],uniqueNames[j]))
            interactionsValidated = rbind(interactionsValidated,c(taxa,metabolite, output[i,j]))
            found= cbind(found,paste(taxa,metabolite,sep=" - "))
            outputValidated[i,j] = output[i,j]
            numberOfTMValidated = numberOfTMValidated+1
            }else{
              print(paste("Not validated: ",uniqueNames[i],uniqueNames[j]))
              numberOfTMNotValidated = numberOfTMNotValidated+1
              outputValidated[i,j] = 0
              interactionsNOTValidated = rbind(interactionsNOTValidated,c(taxa,metabolite, output[i,j]))
          }
        }
      }
    }
    taxaNotFound = unique(taxaNotFound)
    metabolitesNotInteresting = unique(metabolitesNotInteresting)
    print(paste("There are:", numberOfTMNotValidated,"NOT Validated interactions by the DB"))
    print(paste("There are:", numberOfTMValidated,"Validated interactions by the DB"))
  
    return(list(numberOfTMValidated = numberOfTMValidated, numberOfTMNotValidated = numberOfTMNotValidated, interactionsValidated = interactionsValidated, interactionsNOTValidated = interactionsNOTValidated, found = found))
  }

  
  #####################################################################
  ############# VALIDATED INTERACTIONS FOR THIS SET OF TAXA ############# 
  #####################################################################
  
  #compute the total number of possible interactions║
  computeNumberOfPossibleInteractions  = function(validTM, taxaNewNames, found){
    
    validTMAux = validTM
    indexTaxaValid = c()
    metabolites = c()
    for (i in 1:length(taxaNewNames)){
      index =  grep(taxaNewNames[i], validTMAux$SpeciesName)
      #get the taxa that are relevant to us
      indexTaxaValid = c(indexTaxaValid,index)
      #get the metabolites that are relevant to us
      metabolites = unlist(strsplit(validTMAux[index,]$ValidMetabolitesHMDB,("\\|")))
      if (length(intersect(metabolites,namesOfMetabolites))!= length(metabolites)){
        # print(index)
      }
      validTMAux$ValidMetabolitesHMDB[index] = paste(intersect(metabolites,namesOfMetabolites), collapse = '|')
    }
    validMetabolites = validTMAux[indexTaxaValid,]$ValidMetabolitesHMDB
    strValidMetabolites = paste(validMetabolites, collapse = '|')
    # write.table(validTMAux[indexTaxaValid,],paste("validTMAux.csv",sep=""),sep=",",quote=F,col.names=NA)
    
    #Now compute the total number of validated interactions in the database, EXCEPT for taxa and metabolites not in the networks
    validMetabolites = strsplit(paste(validMetabolites, collapse = '|'),"\\|")[[1]]
    blank = which(validMetabolites %in% "")
    if (length(blank)>0)
      validMetabolites = validMetabolites[-blank]
    nTotalValidatedInteractionsOfJustOurData = length(validMetabolites)
    nTotalValidatedInteractions = length(validMetabolites)
    
    # #Now compute the total number of validated interactions in the database, even for taxa and metabolites not in the networks
    validMetabolitesInteraction = validTMAux$ValidMetabolitesHMDB
    validMetabolitesInteraction = strsplit(paste(validMetabolitesInteraction, collapse = '|'),"\\|")[[1]]
    blank = which(validMetabolitesInteraction %in% "")
    if (length(blank)>0)
      validMetabolitesInteraction = validMetabolitesInteraction[-blank]
    # nTotalValidatedInteractions = length(validMetabolitesInteraction)
    
    
    #For every validated interaction in the network, count the number of taxa with which that metabolite interacts with
    count = c()
    for (interaction in found ){
      metabolite = strsplit(interaction[[1]]," - ")[[1]][2]
      count = c(count,sum(str_count(validMetabolitesInteraction, metabolite)))
    }
    if (length(count) > 0) hist(count)
    
    return(count)
  }
  
  
  ###############################################
  ############################################### BINOMIAL
  ###############################################
  binomialTest <- function(relevantMet, validTM, taxaNewNames, namesOfTaxa, namesOfMetabolites, numberOfTMValidated, numberOfTMNotValidated, interactionsValidated, interactionsNOTValidated){
    print("DOING BINOMIAL")
   
     # Calculate the probability of a metabolite being associated with a taxa, based on the validation file for
    # our metabolites and our taxa
    #Iterate over the mets
    probabilityTable = c()
    for (met in relevantMet){
      possibleTaxaInteractions = 0
      actualTaxaInteractions = 0
      #For every valid interaction
      for (i in 1:nrow(validTM)){
        #If it is one of our taxa
        if (validTM$SpeciesName[i] %in% taxaNewNames){
          possibleTaxaInteractions = possibleTaxaInteractions+1
          if (met %in% strsplit(validTM$ValidMetabolitesHMDB[i],"\\|")[[1]]){
            actualTaxaInteractions = actualTaxaInteractions+1
          }
        }
      }
      # possibleTaxaInteractions = length(taxaNewNames)
      # probabilityTable = rbind(probabilityTable,c(met,actualTaxaInteractions,possibleTaxaInteractions,actualTaxaInteractions/possibleTaxaInteractions))
      probabilityTable = rbind(probabilityTable,c(met,actualTaxaInteractions,possibleTaxaInteractions,actualTaxaInteractions/length(taxaNewNames)))
      }
    colnames(probabilityTable) = c("Metabolite","actualTaxaInteractions","possibleTaxaInteractions","probability")
   
    print("probabilityTable") 
    print(probabilityTable) 
    print("interactionsValidated")
    print(interactionsValidated)
    print("interactionsNOTValidated")
    print(interactionsNOTValidated)
    
    if (length(interactionsValidated) == 0){
      probBin = 1
    }else{
      #calculate probabilities for the trials
      allInteractions = rbind(interactionsValidated,interactionsNOTValidated)
      probabilityVector = c()
      for (i in 1:nrow(allInteractions)){
        met = allInteractions[i,2]
        for (j in 1:nrow(probabilityTable)){
          if (probabilityTable[j,1] %in% met){
            probabilityVector = c(probabilityVector,probabilityTable[j,4])
          }
        }
      }
      probabilityVector = as.numeric(probabilityVector)
      # probBin = dpoisbinom(numberOfTMValidated,probabilityVector)
      probBin = tryCatch({
        ppoisbinom(numberOfTMValidated,probabilityVector,lower_tail=F)
      }, error = function(e) {
        1
      })

    }
    
    FN = 0    # fn = not in the network, but in the database
    for (row in  validTM$ValidMetabolitesHMDB){
      FN = FN + str_count(row,"hmdb")
    }
    precision = numberOfTMValidated / (numberOfTMValidated+numberOfTMNotValidated) # (tp)/(tp+fp)
    recall = numberOfTMValidated / (numberOfTMValidated+FN)# (tp)/(tp+fn)

    
    writeBinomial = c(length(namesOfTaxa),possibleTaxaInteractions,length(namesOfMetabolites),length(relevantMet),numberOfTMValidated, numberOfTMNotValidated,precision, recall, round(probBin,6))
    names(writeBinomial) = c("Taxa","TaxaFoundInValidationDB","Met","MetFoundInValidationDB", "TP","FP","Precision","Recall","Pr(validated>=TP)")
    print(writeBinomial)
    return(writeBinomial)
  }
  

  
  ##################### FUNCTION CALLS
  
  # folder = "C:/Users/danir/Google Drive/RESEARCH/DBNs/Heterogeneous/Validation/VMH/" #"Demo" #"finalFiguresBootscore" #"finalFiguresBootscore"
  # setwd(folder)
  
  
  # validTM = read.csv(paste(folder,"validtTaxaToMetabolitesFileHashSignifPos.csv",sep=""), sep = "$")
  validTM = read.csv(paste(folder,"MIMOSA/DBNMimosa/","validtTaxaToMetabolitesFileHashAll.csv",sep=""), sep = "$")
  # validTM = read.csv(paste(folder,"validtTaxaToMetabolitesFileHashVMH.csv",sep=""), sep = "$")
  # validTM = read.csv(paste(folder,"validtMetabolitesToTaxaFileHashVMH.csv",sep=""), sep = "$")
  
  validTM$SpeciesName = tolower(as.character(validTM$Taxa))
  validTM$ValidMetabolitesKEGG = tolower(as.character(validTM$KEGGCompounds))
  validTM$ValidMetabolitesHMDB = tolower(as.character(validTM$HMDBMetabolites))
  
  relevantMet = retrieveRelevantMet(validTM, namesOfMetabolites)
  
  # pairsTM = getJustTMInteractions(TMDirectionCounts,pairs, taxaNewNames, namesOfMetabolites)
  output = adjacencyMatrix(uniqueNames, pairs)
  
  # computeNumberOfPossibleInteractions (validTM, taxaNewNames, found)
  
  validationResult = computeValidatedInteractions(output,uniqueNames, clinical, TMDirectionCounts, validTM, relevantMet)
    numberOfTMValidated = validationResult$numberOfTMValidated
    numberOfTMNotValidated = validationResult$numberOfTMNotValidated
    interactionsValidated = validationResult$interactionsValidated
    interactionsNOTValidated = validationResult$interactionsNOTValidated
    found = validationResult$found
  
  binomialTest(relevantMet, validTM, taxaNewNames, namesOfTaxa, namesOfMetabolites, numberOfTMValidated, numberOfTMNotValidated, interactionsValidated, interactionsNOTValidated)
    
    

  ######
  # Uncomment this to do the precission-recall curve
  binomialTotal = c()
  for (t in seq(0.05,1,0.05)){
    print(t)

    outputThresholdPairsTM = adjacencyMatrix(uniqueNames, pairs[as.numeric(pairs[,4]) >= t,])
    validationResult = computeValidatedInteractions(outputThresholdPairsTM,uniqueNames, clinical, TMDirectionCounts, validTM, relevantMet)
      numberOfTMValidated = validationResult$numberOfTMValidated
      numberOfTMNotValidated = validationResult$numberOfTMNotValidated
      interactionsValidated = validationResult$interactionsValidated
      interactionsNOTValidated = validationResult$interactionsNOTValidated
      found = validationResult$found

    binomial = binomialTest(relevantMet, validTM, taxaNewNames, namesOfTaxa, namesOfMetabolites, numberOfTMValidated, numberOfTMNotValidated, interactionsValidated, interactionsNOTValidated)
    binomialTotal = rbind(binomialTotal,c(binomial,t))

  }
  colnames(binomialTotal)[10] = "t"
  
  outputFile = "outBinomial.csv"
  write.table(binomialTotal, outputFile, sep=",", append=TRUE, row.names = F,quote=F)
  
  #add a row to the output file
  # if(!file.exists(outputFile)){
  #   write.csv(binomialTotal, outputFile,row.names = F,quote=F) 
  # }else{
  #   write.table(binomialTotal, outputFile, sep=",", append=TRUE, col.names=FALSE, row.names = F,quote=F)
  # }
  
  binomialTotal = as.data.frame(binomialTotal)
  write.table(binomialTotal, file = "",append = TRUE, sep = "\t",  quote=FALSE)
  
  plot(binomialTotal$Precision, binomialTotal$Recall, type="b",main = "Precision-Recall", xlab = "Precision", ylab="Recall")
  
  # gap.plot( binomialTotal$t,binomialTotal$`Pr(validated>=TP)`,ytics = c(0,0.05,0.1,0.2,1), gap=c(0.25,0.98), type="b",main = "Cumulative distribution function", xlab = "t", ylab="Pr(validated>=TP)")
  plot( binomialTotal$t,binomialTotal$`Pr(validated>=TP)`, type="b",main = "Cumulative distribution function", xlab = "t", ylab="Pr(validated>=TP)")
  
  TMDirectionCounts
  TMSignCounts 
  nameOfNetwork  
  
  
  
  
  
  
  