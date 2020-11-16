#Combines the mutation spectra from a series of clusters
#To run:

arg <- commandArgs(TRUE)

mutationSpectra <- arg #Will be a list of mutation spectra

mutationSpectrum1 <- read.table(arg[1],sep="\t",header=TRUE)
mutationSpectrum <- data.frame(matrix(0,ncol=3,nrow=length(mutationSpectrum1[,1]))) #Will be filled with the combined mutational spectra
names(mutationSpectrum) <- c("Mutation","Number.of.mutations","Proportion.of.mutations")
mutationSpectrum[,1] <- mutationSpectrum1[,1]

for (sample in 1:length(mutationSpectra)) {
  mutation <- read.table(mutationSpectra[sample],sep="\t",header=TRUE) #Import the mutational spectrum
  for (mutationSample in 1:length(mutation[,1])) { #Iterate through the mutations
    mutationSpectrum[mutationSample,2] <- mutationSpectrum[mutationSample,2] + mutation[mutationSample,2]}}

totalMutations <- sum(mutationSpectrum[,2]) #Assign to the total number of mutations in the spectra
mutationSpectrum[,3] <- (mutationSpectrum[,2]/totalMutations)

write.table(mutationSpectrum,"combined.mutational.spectrum.txt",quote=FALSE,row.names=FALSE,sep="\t")