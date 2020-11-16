#Plots a network of patient linkages with edges weighted and shaded by number of SNPs

args <- commandArgs(TRUE)

library(igraph)
library(GGally)
library(network)
library(intergraph)
library(sna)
library(ggplot2)

#Import the patient linkages
patientLinkages <- read.table(args[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#Output file prefix
outFile <- args[2]

#Extract the linkages below the required SNP threshold
patientLinkagesCutoff <- patientLinkages[which(patientLinkages$Minimum.SNP.distance <= 38),]

#Extract the unique patients
uniquePatients <- unique(c(patientLinkagesCutoff$Patient1, patientLinkagesCutoff$Patient2))

#Will be filled with the patient information to label the network
patientInformation <- data.frame(matrix("", ncol = 4, nrow = length(uniquePatients)), stringsAsFactors = FALSE)
names(patientInformation) <- c("Patient", "Country", "Location", "CF_status")

#Iterate through the patients and add their information to patientInformation
for (patient in 1:length(uniquePatients)) {
  patientSamples <- patientLinkagesCutoff[which(patientLinkagesCutoff$Patient1 == uniquePatients[patient] | patientLinkagesCutoff$Patient2 == uniquePatients[patient]),]
  patientInformation$Patient[patient] <- uniquePatients[patient]
  if (patientSamples$Patient1[1] == uniquePatients[patient]) {
    patientInformation$Country[patient] <- patientSamples$Patient1.country[1]
    patientInformation$Location[patient] <- patientSamples$Patient1.location1[1]
    patientInformation$CF_status[patient] <- patientSamples$Patient1.CF.status[1]} else {
    if (patientSamples$Patient2[1] == uniquePatients[patient]) {
      patientInformation$Country[patient] <- patientSamples$Patient2.country[1]
      patientInformation$Location[patient] <- patientSamples$Patient2.location1[1]
      patientInformation$CF_status[patient] <- patientSamples$Patient2.CF.status[1]}}
  if (patientInformation$Location[patient] == "") {
    patientInformation$Location[patient] <- patientInformation$Country[patient]}}

#Create a network from the patients, assumes that the first 2 columns in the data frame are the patient identifiers
clusters <- graph_from_data_frame(patientLinkagesCutoff, directed = FALSE)
clustersWithVertices <- graph_from_data_frame(patientLinkagesCutoff, directed = FALSE, vertices = patientInformation)

#Add a CF status attribute that will be filled with CF status
clustersWithCFStatus <- set_vertex_attr(clusters, "CF_status", value = "grey")

#Iterate through the nodes and add their CF status to the attribute
for (node in V(clustersWithCFStatus)) {
  if (patientInformation$CF_status[which(patientInformation$Patient == names(V(clusters))[node])] == "CF") {
    V(clustersWithCFStatus)$CF_status[node] <- "dodgerblue"} else {
    if (patientInformation$CF_status[which(patientInformation$Patient == names(V(clusters))[node])] == "Non-CF") {
      V(clustersWithCFStatus)$CF_status[node] <- "orange"}}}

#Set a SNP distance scale to size the edges
E(clustersWithCFStatus)$SNP.distance <- (39 - E(clustersWithCFStatus)$Minimum.SNP.distance)/20
E(clustersWithVertices)$SNP.distance <- (39 - E(clustersWithVertices)$Minimum.SNP.distance)/20

#Set a SNP distance scale to set alpha for the edges
E(clustersWithCFStatus)$SNP.distance.alpha <- (39 - E(clustersWithCFStatus)$Minimum.SNP.distance)/39
E(clustersWithVertices)$SNP.distance.alpha <- (39 - E(clustersWithVertices)$Minimum.SNP.distance)/39

pdf(paste(outFile, ".cf.status.pdf", sep = ""))
clusterPlot <- ggnet2(clusters, size = 2, color = V(clustersWithCFStatus)$CF_status)
print(clusterPlot)
dev.off()

pdf(paste(outFile, ".cf.status.centre.pdf", sep = ""))
clusterPlot <- ggnet2(clusters, size = 2, color = V(clustersWithVertices)$Location) + scale_color_discrete() + guides(color = FALSE)
print(clusterPlot)
dev.off()

pdf(paste(outFile, ".cf.status.country.pdf", sep = ""))
clusterPlot <- ggnet2(clusters, size = 2, color = V(clustersWithVertices)$Country) + scale_color_discrete() + scale_color_manual(values = c("blue", "magenta", "steelblue", "red", "purple", "black", "goldenrod", "grey", "orange", "tan", "coral", "dodgerblue", "red"))
print(clusterPlot)
dev.off()