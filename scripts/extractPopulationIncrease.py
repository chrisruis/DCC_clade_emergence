#Calculates the distribution of dates at which multiple population increases have occurred from a Bayesian skyline plot
#To run:

from cStringIO import StringIO
from Bio import Phylo as p
from operator import itemgetter
import argparse

def getDistribution(distributionFile): #This function takes a log file and returns the header region
	distribution = distributionFile.readlines()
	positions = []
	for line in distribution:
		if not "tree STATE" in line:
			positions.append(line)
	positions.pop(-1)
	return positions

def getGroupSizes(logFile): #This function takes a log file and returns the columns corresponding to the GroupSizes
	logInformation = logFile.readlines()
	positions = []
	for i, identifier in enumerate(logInformation[0].strip().split("\t")):
		if "GroupSizes" in identifier:
			positions.append(i)
	return positions

def getPopulationSizes(logFile): #This function takes a log file and returns the columns corresponding to the PopulationSizes
	logInformation = logFile.readlines()
	positions = []
	for i,identifier in enumerate(logInformation[0].strip().split("\t")):
		if "PopSizes" in identifier:
			positions.append(i)
	return positions

def getNodeHeights(tree,date): #This function takes a tree and returns node heights
	nodeHeights = []
	for clade in tree.get_nonterminals():
		nodeHeights.append(float(date)-(float(max(tree.depths().iteritems(),key=itemgetter(1))[1]) - float(tree.depths()[clade])))
	sortedNodeHeights = sorted(nodeHeights)
	return sortedNodeHeights

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("posterior_trees", help = "Nexus file containing trees")
    parser.add_argument("log_file", help = "Log file with header region removed")
    #parser.add_argument("nexus", help = "Posterior distribution of trees from BEAST without any changes")
    parser.add_argument("-p", help = "Minimum percentage change for a significant increase, default 10", default = "10")
    parser.add_argument("-l", help = "Date of the latest sample")
    parser.add_argument("-n", help = "Print update every nth tree. Default is 1000", default = "1000")
    parser.add_argument("-d", help = "Optional to set date before which the increase in population size needs to have occurred in the sample", default = "10000")
    parser.add_argument("-o", help = "Prefix of output file")
    args = parser.parse_args()

    outFile = open(args.o,"w")

    with open(args.log_file) as fileobject:
		GroupSizesPositions = getGroupSizes(fileobject)

    with open(args.log_file) as fileobject:
		PopulationSizesPositions = getPopulationSizes(fileobject)

    with open(args.log_file) as fileobject:
		logfile = fileobject.readlines()[1:]

    j = 0

    populationIncreases = [] #Will be filled with the dates of the population increase in each sample
    
    with open(args.posterior_trees) as fileobject:
        for line,logline in zip(fileobject,logfile):
            if j %int(args.n) == 0:
                print "Tree", j
            j += 1
            phylogeny = p.read(StringIO(line),"newick")
            nodeHeight = getNodeHeights(phylogeny,args.l)
            groupSizes = [logline.strip().split("\t")[i] for i in GroupSizesPositions][::-1] #Assign to the column of each of the GroupSizes
            populationSizes = [logline.strip().split("\t")[i] for i in PopulationSizesPositions][::-1] #Assign to the column of each of the PopulationSizes
            samplePopulationIncreases = [] #Will be filled with the dates at which the population increases in the sample
            for populationSize in range(len(populationSizes)-1): #Iterate through the population sizes in the MCMC step
                if float(populationSizes[populationSize+1]) > (float(populationSizes[populationSize])+(float(populationSizes[populationSize])*(float(args.p)/100.0))): #Check whether the population size has increased by more than the required amount
                    if populationSize == 0: #Check if the population size increases after the first group
                        nodeNumber = groupSizes[0]
                        samplePopulationIncreases.append(nodeHeight[int(nodeNumber)-1])
                    else: #Sum the nodes until the population increase
                        nodeNumber = sum([int(k) for k in groupSizes[0:(populationSize+1)]])
                        samplePopulationIncreases.append(nodeHeight[int(nodeNumber)-1])
            outFile.write("\t".join(str(l) for l in samplePopulationIncreases) + "\n")
    
    outFile.close()