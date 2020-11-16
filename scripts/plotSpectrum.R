#Plots the contextual mutation spectrum, run on the output from extractMutationSpectrumContext.py

arg <- commandArgs(TRUE)

library(ggplot2)

mutationSpectrum <- read.table(arg[1],sep="\t",header=TRUE) #Import the mutation spectrum
outFile <- arg[2]

mutationSpectrum$Mutation.factor <- factor(mutationSpectrum$Mutation,levels=mutationSpectrum$Mutation)

totalMutations <- sum(mutationSpectrum[,2]) #Assign to the total number of mutations
mutationSpectrum$Proportion.of.mutations = mutationSpectrum$Number.of.mutations/totalMutations #Assign to the proportion of mutations of each type
mutationSpectrum$Type <- c(rep("cyan",16),rep("black",16),rep("red",16),rep("grey",16),rep("green",16),rep("magenta",16))

pdf(outFile)
mutationPlot <- ggplot(mutationSpectrum,aes(x=Mutation.factor,y=Proportion.of.mutations)) + theme_classic() + geom_bar(stat="identity",fill=mutationSpectrum$Type)
print(mutationPlot)
dev.off()