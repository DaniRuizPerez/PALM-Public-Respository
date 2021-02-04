  #########################################################################################################################
  #########################################################################################################################
  ###### PROJECT:        DBNs
  ###### NAME:           ExportNetworksToList.R
  ###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
  ###### AFFILIATION:    Florida International University
  ###### DESCRIPTION:    Exports the .graphml files of one folder to a list of edges sorted by importance
  #########################################################################################################################
  #########################################################################################################################
  

  rm(list=ls())
  library(scales)
  library(stringr)
  library("stats")
  library(poisbinom)
  library(plotrix)
  
  options("scipen"=100, "digits"=4)
  

  folder = ""
  setwd(folder)
  files = list.files(path = folder, pattern = ".graphml")
  
  for (f in files){

      nameOfNetwork = paste(folder,f,sep="/")
  
      network = tryCatch({
      readChar(nameOfNetwork, file.info(nameOfNetwork)$size)
    }, warning = function(e) {
      ""
    }, error = function(e) {
      ""
    })
    if (network == "")
      next
    
    network = str_replace(network,"0.0.","0.")
  
    
    #####################################################################
    ############# LOAD INFORMATION FROM  GRAPHML NETWORK ############# 
    #####################################################################
    
    #READ NETWORK
    network =  gsub("\\(", "", gsub("\\)", "",gsub("\\[", "", gsub("\\]","",gsub("t0", "ti", gsub("tn", "ti+1", tolower(unlist(strsplit(x =network, "\r?\n")))))))))
    network = gsub("\n<data key=\"key_bootscore\">[0-9]*\\.[0-9]*</data>\n","",network)
    newNetwork = network[grep(pattern = paste(".*target=\".*_ti\\+1\">","",sep=""), x = network, ignore.case = T)]
    newNetwork = gsub("s__", "", gsub("\\(", "", gsub("\\)", "", gsub("\\[", "", gsub("\\]", "", gsub("_ti", "", gsub("\\+1", "",  gsub("\">", "", gsub("g__", "", gsub("<edge id=\"e\\d*\" source=\"", "", newNetwork))))))))))
    
  
    #get the weight of the edges
    edgeWeight = network[grep("<data key=\\\"key_weight.*",network)]
    edgeWeight = as.numeric(gsub("\"", "", gsub("<data key=\\\"key_weight\\\">", "", gsub("</data>", "", edgeWeight)))) 
    
    #get the bootscore of the edges
    edgeBootscore = network[grep("<data key=\\\"key_bootscore.*",network)]
    edgeBootscore = as.numeric(gsub("\"", "", gsub("<data key=\\\"key_bootscore\\\">", "", gsub("</data>", "", edgeBootscore)))) 
    
    
    #Extract the information from the network
    pairs = c()
    for (i in 1:length(newNetwork)){
      edge =  c(trimws(strsplit(newNetwork[i],"\" target=\"")[[1]], which ="both"),as.numeric(edgeWeight[i]),as.numeric(edgeBootscore[i]))
      pairs = rbind(pairs,edge)
    }
  
    
    colnames(pairs) = c("source","target", "weigh", "bootscore")
    pairs = as.data.frame(pairs)
  
    pairs = pairs[order(abs(as.numeric(as.character(pairs[,3]))*as.numeric(as.character(pairs[,4]))),decreasing=T ),]
    
    write.table(pairs, paste(strsplit(f,".graphml")[[1]],"_list.txt"), sep=",", append=F, row.names = F,quote=F)
    
  }
 
  
  
  