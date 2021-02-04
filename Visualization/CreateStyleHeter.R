#########################################################################################################################
#########################################################################################################################
###### PROJECT:        DBNs
###### NAME:           CreateStyle.R
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
###### AFFILIATION:    Florida International University
###### 
###### DESCRIPTION:    This file creates the style.xml file needed to visualize it in Cytoscape. It takes as input a network
######                , a file with the mean abundance and a base style file to update.
#########################################################################################################################
#########################################################################################################################


library(scales)
library(stringr)
options("scipen"=100, "digits"=4)


#nameOfNetwork = "DemoFigure.graphml"#list.files(pattern =".*\\.graphml")
nameOfNetwork = "/Networks/human_ibd_microbiota_genes_metabolites_dbn_sample_alignment_sr14d_top100x100_filtered_dbnIntraBoot100boots.graphml"#list.files(pattern =".*\\.graphml")
nameOfNetwork = "/Networks/human_ibd_microbiota_genes_metabolites_dbn_sample_alignment_sr14d_top100x100_host_genes_reference_dbnIntraBoot.graphml"#list.files(pattern =".*\\.graphml")
nameOfNetwork = paste(folder,nameOfNetwork,sep="/")
abundances = unlist(read.csv('meanAbundanceiHMPAlignmentHost.txt',sep="\t"))

# READ STYLE BASE
f = readChar("styleBase.xml", file.info("styleBase.xml")$size)

split = strsplit(f, "<dependency name=\"nodeSizeLocked\" value=\"true\"/>", fixed = T)
before = unlist(split)[1]
after = unlist(split)[2]

#READ NETWORK
network = readChar(nameOfNetwork, file.info(nameOfNetwork)$size)
network= unlist(strsplit(x =network, "\r\n"))
network= unlist(strsplit(x =network, "\n"))

#CLASSIFY ATTRIBUTES BASED ON TYPE
namesAll = str_extract(network[grep("<node id=",network)], "((s|g|m|hg.)__)?.+\\_ti(\\+1)?")
namesAll = gsub("<node id=\"", "", namesAll)
namesAll = gsub("_ti", "", namesAll)
namesAll = gsub("\\+1", "", namesAll)
# namesAll = gsub("\\[", "", namesAll)
# namesAll = gsub("\\]", "", namesAll)
namesAll = unique(namesAll)
prefix = substr(namesAll,1,3)
#namesAll = substr(namesAll,4,1000)



taxa = namesAll[which(prefix %in% "s__")]
host = c(namesAll[which(prefix %in% "hgi")],namesAll[which(prefix %in% "hgr")],namesAll[which(prefix %in% "hgs")])
genes = namesAll[which(prefix %in% "g__")]
metabolites = namesAll[which(prefix %in% "m__")]
clinical =  c("Week sample obtained_ti","Week sample obtained_ti+1") #namesAll[! (prefix %in% "g__"| prefix %in% "s__" | prefix %in% "m__")]

all = c(host,taxa, genes, metabolites, clinical)

################## NODE COLOR -> timepoint
aux ="\n            <visualProperty name=\"NODE_FILL_COLOR\" default=\"#ff9232\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">"
for (i in 1:length(all)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"#4d93a8\" attributeValue=\"", paste(all[i],"_ti",sep=""),"\"/>",sep="")
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
f = paste(before,"<dependency name=\"nodeSizeLocked\" value=\"true\"/>",aux,after,sep= "")




############### SIZE -> incoming edges
# maxAcum = 1
# for (i in 1:length(all))
#   maxAcum = max(maxAcum, length(grep(pattern = paste(".*target=\"",all[i],"_ti+1","\".*",sep=""), x = network, ignore.case = T)))
# 
# aux =paste(aux,"\n            <visualProperty name=\"NODE_SIZE\" default=\"40\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
# for (i in 1:length(all)){
#   aux = paste(aux,"\n                    <discreteMappingEntry value=\"",(20*length(grep(pattern = paste(".*target=\"",all[i],"_ti","\".*",sep=""), x = network, ignore.case = T))/maxAcum+40),"\" attributeValue=\"",namesAll[i],"_ti","\"/>",sep="")
#   aux = paste(aux,"\n                    <discreteMappingEntry value=\"",(20*length(grep(pattern = paste(".*target=\"",all[i],"_ti\\+1","\".*",sep=""), x = network, ignore.case = T))/maxAcum+40),"\" attributeValue=\"",namesAll[i],"_ti+1","\"/>",sep="")
# }
# aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
# f = paste(before,"<dependency name=\"nodeSizeLocked\" value=\"true\"/>",aux,after,sep= "")



############### TRANSPARENCY -> abundance
abundancesScaled = c(rescale(abundances[1:length(host)],c(90,255)),
rescale(abundances[(1+length(host)):(length(taxa)+length(host))],c(90,255)),
rescale(abundances[(1+length(host)+length(taxa)):(length(genes)+length(taxa)+length(host))],c(90,255)),
rescale(abundances[(1+length(host)+length(taxa)+length(genes)):(length(genes)+length(taxa)+length(host)+length(metabolites))],c(90,255)))

abundancesScaled = c(200,abundancesScaled)
aux =paste(aux,"\n            <visualProperty name=\"NODE_TRANSPARENCY\" default=\"255\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
for (i in 1:length(namesAll)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"",round(abundancesScaled[i],1),"\" attributeValue=\"",all[i],"_ti","\"/>",sep="")
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"",round(abundancesScaled[i],1),"\" attributeValue=\"",all[i],"_ti+1","\"/>",sep="")
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
f = paste(before,"<dependency name=\"nodeSizeLocked\" value=\"true\"/>",aux,after,sep= "")



############### SHAPE -> type of data

aux =paste(aux,"\n            <visualProperty name=\"NODE_SHAPE\" default=\"ELLIPSE\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
for (i in 1:length(genes)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"DIAMOND\" attributeValue=\"",genes[i],"_ti","\"/>",sep="")
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"DIAMOND\" attributeValue=\"",genes[i],"_ti+1","\"/>",sep="")
}
for (i in 1:length(host)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"HEXAGON\" attributeValue=\"",host[i],"_ti","\"/>",sep="")
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"HEXAGON\" attributeValue=\"",host[i],"_ti+1","\"/>",sep="")
}
for (i in 1:length(clinical)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"TRIANGLE\" attributeValue=\"",clinical[i],"","\"/>",sep="")
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"TRIANGLE\" attributeValue=\"",clinical[i],"","\"/>",sep="")
}
for (i in 1:length(metabolites)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"RECTANGLE\" attributeValue=\"",metabolites[i],"_ti","\"/>",sep="")
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"RECTANGLE\" attributeValue=\"",metabolites[i],"_ti+1","\"/>",sep="")
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
f = paste(before,"<dependency name=\"nodeSizeLocked\" value=\"true\"/>",aux,after,sep= "")


split = strsplit(f, "<dependency name=\"arrowColorMatchesEdge\" value=\"false\"/>", fixed = T)
before = unlist(split)[1]
after = unlist(split)[2]



################## EDGE LINE TYPE -> intra-inter
aux = "\n            <visualProperty name=\"EDGE_LINE_TYPE\" default=\"SOLID\">\n                <discreteMapping attributeType=\"string\" attributeName=\"shared name\">"
for (i in 1:length(all)){
   for (j in i:length(all)){
      aux = paste(aux,"\n                    <discreteMappingEntry value=\"EQUAL_DASH\" attributeValue=\"",all[i],"_ti+1 (-) ", all[j],"_ti+1","\"/>",sep="")
      aux = paste(aux,"\n                    <discreteMappingEntry value=\"EQUAL_DASH\" attributeValue=\"",all[i],"_ti (-) ", all[j],"_ti","\"/>",sep="")
  }
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")


################## EDGE PAINT
# aux =paste(aux,"\n            <visualProperty name=\"EDGE_UNSELECTED_PAINT\" default=\"#CC0033\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
# for (i in 1:length(names)){
#     aux = paste(aux,"\n                    <discreteMappingEntry value=\"#FF9999\" attributeValue=\"",names[i],"_ti (-) ", names[i],"_ti+1","\"/>",sep="")
#   }
# 
# aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")



################  EDGE TRANSPARENCY based on bootscore. Make the self loops more transparent
# aux =paste(aux,"\n            <visualProperty name=\"EDGE_TRANSPARENCY\" default=\"255\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
# for (i in 1:length(all)){
#   aux = paste(aux,"\n                    <discreteMappingEntry value=\"70\" attributeValue=\"",all[i],"_ti (-) ", all[i],"_ti+1","\"/>",sep="")
# }
# aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")




# ################## EDGE TRANSPARENCY
weights = str_extract(network[grep("key=\"key_bootScore",network)], "\\-*\\d+\\.+\\d+")
weights = weights[!is.na(weights)]
maxAbsWeight =  max(abs(as.numeric(weights)))
minAbsWeight =  min(abs(as.numeric(weights)))

# ######## For edge coefficient
aux =paste(aux,"\n<visualProperty name=\"EDGE_TRANSPARENCY\" default=\"2.0\">\n")
aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"bootScore\">\n")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"50.0\" greaterValue=\"50.0\" equalValue=\"50.0\" attrValue=\"",minAbsWeight,"\"/>\n",sep="")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"50.0\" greaterValue=\"255.0\" equalValue=\"255.0\" attrValue=\"",maxAbsWeight,"\"/>\n",sep="")
aux =paste(aux,"  </continuousMapping>\n")
aux =paste(aux,"  </visualProperty>\n")




# ######## make self loops Invisible:
aux =paste(aux,"\n<visualProperty default=\"true\" name=\"EDGE_VISIBLE\">")
aux =paste(aux,"\n<discreteMapping attributeName=\"name\" attributeType=\"string\">")
for (i in 1:length(all)){
  aux = paste(aux,"\n <discreteMappingEntry attributeValue=\"",all[i],"_ti (-) ", all[i],"_ti+1\""," value=\"false\"/>",sep="")
}
aux =paste(aux,"  \n</discreteMapping>\n")
aux =paste(aux,"  </visualProperty>\n")




# ######## For edge coefficient
# aux =paste(aux,"\n<visualProperty name=\"EDGE_TRANSPARENCY\" default=\"2.0\">\n")
# aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"bootScore\">\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"130.0\" greaterValue=\"130.0\" equalValue=\"130.0\" attrValue=\"",minAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"130.0\" greaterValue=\"255.0\" equalValue=\"255.0\" attrValue=\"",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  </continuousMapping>\n")
# aux =paste(aux,"  </visualProperty>\n")



################## EDGE WIDTH
weights = str_extract(network[grep("key=\"key_weight",network)], "\\-*\\d+\\.+\\d+")
weights = weights[!is.na(weights)]

maxAbsWeight =  max(abs(as.numeric(weights)))
medianAbsWeight =  median(abs(as.numeric(weights)))
# maxAbsWeight = 150 #when we hide intra edges we need this



######## For edge coefficient normalized
aux =paste(aux,"\n<visualProperty name=\"EDGE_WIDTH\" default=\"2.0\">\n")
aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"weight\">\n")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"15.0\" greaterValue=\"15.0\" equalValue=\"15.0\" attrValue=\"-",1,"\"/>\n",sep="")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"1.0\" equalValue=\"1.0\" attrValue=\"0.0\"/>\n")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"15.0\" greaterValue=\"15.0\" equalValue=\"15.0\" attrValue=\"",1,"\"/>\n",sep="")
aux =paste(aux,"  </continuousMapping>\n")
aux =paste(aux,"  </visualProperty>\n")



# ######## For edge coefficient
# aux =paste(aux,"\n<visualProperty name=\"EDGE_WIDTH\" default=\"2.0\">\n")
# aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"weight\">\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"10.0\" greaterValue=\"15.0\" equalValue=\"15.0\" attrValue=\"-",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"10.0\" equalValue=\"10.0\" attrValue=\"-",medianAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"1.0\" equalValue=\"1.0\" attrValue=\"0.0\"/>\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"10.0\" equalValue=\"10.0\" attrValue=\"",medianAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"10.0\" greaterValue=\"15.0\" equalValue=\"15.0\" attrValue=\"",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  </continuousMapping>\n")
# aux =paste(aux,"  </visualProperty>\n")
  
########### For edge confidence
# aux =paste(aux,"<visualProperty name=\"EDGE_WIDTH\" default=\"2.0\">\n")
# aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"weight\">\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"16.0\" greaterValue=\"20.0\" equalValue=\"20.0\" attrValue=\"-",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"8.0\" greaterValue=\"16.0\" equalValue=\"16.0\" attrValue=\"-",80.0,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"8.0\" greaterValue=\"8.0\" equalValue=\"8.0\" attrValue=\"-",20.0,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"1.0\" equalValue=\"1.0\" attrValue=\"0.0\"/>\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"8.0\" greaterValue=\"8.0\" equalValue=\"8.0\" attrValue=\"",20.0,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"8.0\" greaterValue=\"16.0\" equalValue=\"16.0\" attrValue=\"",80.0,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"16.0\" greaterValue=\"20.0\" equalValue=\"20.0\" attrValue=\"",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  </continuousMapping>\n")
# aux =paste(aux,"  </visualProperty>\n")


fileConn<-file("style.xml")
writeLines(paste(before,"<dependency name=\"arrowColorMatchesEdge\" value=\"false\"/>",aux,after,sep=""), fileConn)
close(fileConn)
