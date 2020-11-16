#Bootstraps a mutation spectrum to generate a distribution of proportions that can be compared with an observed mutation spectrum
#to determine statistical significance

arg <- commandArgs(TRUE)

bootstrapSpectrum <- read.table(arg[1],sep="\t",header=TRUE) #Import the distribution that will be bootstrapped
mutationSpectrum <- read.table(arg[2],sep="\t",header=TRUE) #Import the distribution that will be compared with the bootstraps
numberBootstraps <- as.numeric(arg[3])
outFile <- arg[4]

numberSamples <- sum(mutationSpectrum[,2]) #Will be the number of mutations that are sampled

bootstrapMutations <- data.frame(matrix(0,ncol=(numberBootstraps+1),nrow=length(bootstrapSpectrum[,1]))) #Will be filled with the proportion of each type of mutation in each bootstrap
names(bootstrapMutations) <- c("Mutation",paste("Bootstrap",1:numberBootstraps,sep=""))
bootstrapMutations[,1] <- bootstrapSpectrum[,1]

mutations <- c() #Will be filled with the mutations that can be sampled
for (sample in 1:length(bootstrapSpectrum[,1])) { #Iterate through the mutations
  mutations <- c(mutations,rep(toString(bootstrapSpectrum[sample,1]),bootstrapSpectrum[sample,2]))}

for (bootstrap in 1:numberBootstraps) {
  bootstrapProportion <- data.frame(matrix(0,ncol=3,nrow=length(bootstrapSpectrum[,1])))
  names(bootstrapProportion) <- c("Mutation","Number.of.mutations","Proportion.of.mutations")
  bootstrapProportion[,1] <- bootstrapSpectrum[,1]
  bootstrapSample <- sample(mutations,numberSamples,replace=FALSE)
  for (sequence in 1:length(bootstrapProportion[,1])) { #Iterate through the mutations
    bootstrapProportion[sequence,2] <- length(which(bootstrapSample == toString(bootstrapProportion[sequence,1])))}
  bootstrapProportion[,3] <- bootstrapProportion[,2]/numberSamples
  bootstrapMutations[,(bootstrap+1)] <- bootstrapProportion[,3]}

write.table(bootstrapMutations,paste(outFile,".bootstrap.proportions.txt",sep=""),quote=FALSE,row.names=FALSE,sep="\t")

totalMutations <- data.frame(matrix(0,ncol=5,nrow=length(bootstrapSpectrum[,1]))) #Will be filled with the distribution of the bootstraps of each mutation
names(totalMutations) <- c("Mutation","Proportion","Mean.bootstrap.proportion","Bootstrap.sd","Significant")
totalMutations[,1] <- bootstrapSpectrum[,1]

for (mutationType in 1:length(bootstrapMutations[,1])) { #Iterate through the mutations
  totalMutations[mutationType,2] <- mutationSpectrum[mutationType,3]
  totalMutations[mutationType,3] <- mean(as.matrix(bootstrapMutations[mutationType,2:length(bootstrapMutations[1,])]))
  totalMutations[mutationType,4] <- sd(as.matrix(bootstrapMutations[mutationType,2:length(bootstrapMutations[1,])]))
  if (mutationSpectrum[mutationType,3] > (totalMutations[mutationType,3] + (2*totalMutations[mutationType,4]))) { #Check if the proportion is greater than expected
    totalMutations[mutationType,5] <- "Greater"} else {
    if (mutationSpectrum[mutationType,3] < (totalMutations[mutationType,3] - (2*totalMutations[mutationType,4]))) { #Check if the proportion is lower than expected
      totalMutations[mutationType,5] <- "Lesser"} else {
      totalMutations[mutationType,5] <- "No"}}}

write.table(totalMutations,paste(outFile,"mutation.comparison.txt"),quote=FALSE,row.names=FALSE,sep="\t")