#########################################################################################################################
#########################################################################################################################
###### PROJECT:        DBNs
###### NAME:           NormalizeEdges.R
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
###### AFFILIATION:    Florida International University
###### 
###### DESCRIPTION:    This file normalizes the weight of the files in a folder by updating them
#########################################################################################################################
#########################################################################################################################




library(scales)
library(stringr)
library(igraph)
options("scipen"=100, "digits"=4)

folder = "" #"Demo" #"finalFiguresBootscore" #"finalFiguresBootscore"
setwd(folder)

files = list.files(path = folder, pattern = ".graphml")

for (k in 1:1){
  for (f in files){
    
    nameOfNetwork = paste(folder,f,sep="/")
    abundances = unlist(read.csv('../../meanAbundanceiHMPAlignmentHost.txt',sep="\t"))
    
    networkO = readChar(nameOfNetwork, file.info(nameOfNetwork)$size)
    # We make sure that the weights are in the correct format
    networkO = gsub("(\\d\\.)\\d+\\.","\\1",networkO)
  
    file = file(f, open = "wb")
    writeChar(object = networkO, con = file,  eos = NULL)
    close(file)
  
    
    network= unlist(strsplit(x =networkO, "\r\n"))
    network= unlist(strsplit(x =network, "\n"))
    
  
    ################## EDGE WIDTH
    weights = str_extract(network[grep("key=\"key_weight",network)], "\\-*\\d*\\.*\\d+")

    ########### WE NORMALIZE THE WEIGHTS BASED ON w_i \average{x_i} / \sum_j  abs(w_j) \average{x_j}
    
    library(igraph)
    edgeList = igraph::get.edgelist(igraph::read.graph(nameOfNetwork, format = "graphml"))
    edgeAttr = igraph::get.edge.attribute(igraph::read.graph(nameOfNetwork, format = "graphml"))
    nodeAttr = igraph::get.vertex.attribute(igraph::read.graph(nameOfNetwork, format = "graphml"))
    
    
    for (i in 1:length(edgeAttr$weight)){
      if ( edgeAttr$weight[i] > 0){
        edgeAttr$weight[i] = min(edgeAttr$weight[i],100)
      }else {
        edgeAttr$weight[i] = max(edgeAttr$weight[i],-100)
      }
    }

    
    
    #This function returns the abundance of a node based on the index of the incoming node
    abundanceOfTaxa <- function(nameOfIncomingNode,abundances){
      # we look for the abundance of that parent node
      abundance = abundances[which(names(abundances) %in% substr(nameOfIncomingNode,1,nchar(nameOfIncomingNode)-3))] 
      if (length(abundance)==0){
        #the node doesn't have abundance, it is clinical. We make it to be an average of the other abundances
        abundance = mean(abundances)
      }
      return (abundance)
    }
    
    normalizedWeights =c()
    #Iterate over all posible children. to normalize incoming edges to each of them
    for (i in 1:length(nodeAttr[[1]])){
      indexOfIncomingEdges = which(edgeList[,2]==i)
      indexOfOutgoingEdges = edgeList[indexOfIncomingEdges,1]
      
      #We need to calculate the sum of the multiplication of each weight times the mean abundance for all incoming edges of this children
      sum = 0
      for (j in 1:length(indexOfIncomingEdges)){
        nameOfIncomingNode = nodeAttr$id[indexOfOutgoingEdges[j]]
        abundance = abundanceOfTaxa(nameOfIncomingNode,abundances)
        sum = sum + abs(edgeAttr$weight[indexOfIncomingEdges[j]]*abundance)
      }
      
      #Iterate over all incoming edges for this children
      for (j in 1:length(indexOfIncomingEdges)){
        nameOfIncomingNode = nodeAttr$id[indexOfOutgoingEdges[j]]
        abundance = abundanceOfTaxa(nameOfIncomingNode,abundances)
        normalizedWeights = as.vector(c(normalizedWeights,edgeAttr$weight[indexOfIncomingEdges[j]]*abundance/sum))
      }
    }
    
  
    
    normalizedWeights = as.character(round(normalizedWeights,10))
    
    #we update the weights in the file
    for(i in 1:length(weights)){
      networkO = gsub(weights[i], normalizedWeights[i], networkO)
    }
    
    networkO = gsub("--","-",networkO)
    
    #save it and take care of the null character that appears
    writeChar(networkO,nameOfNetwork,nchars = nchar(networkO))
    r = readBin(nameOfNetwork, raw(), file.info(nameOfNetwork)$size)
    r[r==as.raw(0)] = as.raw(0x20) ## replace with 0x20 = <space>
    writeBin(r, nameOfNetwork)
  
  
  }
}