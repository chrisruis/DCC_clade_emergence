#Plots a phylogenetic tree, colouring branches by a given label, colour conversion file contains 2 columns with no header, column 1 is the
#label and column 2 is the colour
#Plots the tree in time units given a most recent sequence date in format Year:Month:Day, e.g. 2010:01:01
#To run: RScript plotPhylogenyColour.R -t $BEASTTree -l $Label -c $ColourConversion -d $MostRecentSampleDate -o $OutputTree

arg <- commandArgs(TRUE)

library(treeio)
library(ggtree)
library(ggpubr)
library(stringr)

hh <- paste(unlist(arg), collapse=' ') #Combine the arguments into a vector

listOfOptions <- unlist(strsplit(hh,'-'))[-1] #Separate the arguments by name

optionArguments <- list()

for (i in 1:length(listOfOptions)) { #Assign the arguments to their options
  optionArguments[[i]] <- strsplit(listOfOptions[i],' ')[[1]][-1]
  names(optionArguments)[i] <- strsplit(listOfOptions[i],' ')[[1]][1]}

phylogeny <- read.beast(optionArguments$t) #Import the BEAST tree to be plotted
date <- str_replace_all(optionArguments$d,":","-") #The date of the most recent sample, as Year:Month:Day
label <- optionArguments$l #The label that will be used to colour the tree
colourConversion <- read.table(optionArguments$c,sep="\t",header=FALSE) #Import the conversion between the label and colours
outFile <- optionArguments$o

colours <- list() #Will be filled with the colours of each label
for (sample in 1:length(colourConversion[,1])) {
  colours[toString(colourConversion[sample,1])] = toString(colourConversion[sample,2])}

pdf(outFile)
phylogenyPlot <- ggtree(phylogeny,mrsd=date,aes(color=location)) + theme_tree2() + scale_x_continuous(limits=c(1900,2020),breaks=c(1900,1920,1940,1960,1980,2000,2020)) + scale_color_manual(values = colours)
print(phylogenyPlot)
dev.off()