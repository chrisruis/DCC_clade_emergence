#Extracts the mutational spectrum from a PAML rst file and divides the mutations into within patient, branch leading to patient, internal branch and terminal branch
#Gets the context from the vcf file from gubbins and the reference file used for initial mapping
#Labels the tree to identify within patient clusters using a sample information file
#Differs from mutationalSpectrum.py in the contextual mutations combined
#This script combines mutations in the reverse complement base triplet so the impact of the mutation on either strand is indistinguishable
#To run:

from cStringIO import StringIO
from Bio import AlignIO
from Bio import Phylo as p
from collections import OrderedDict
import argparse

def convertVCFPositions(vcf): #This function takes a VCF file and returns the conversion between alignment position and genome position as a dictionary with the alignment positions as the keys
    alignmentPosition = 1 #Will be incremented with each position in the vcf file
    positionDict = {} #Will be filled with each of the alignment positions and the corresponding genome position
    for position in vcf: #Iterate through the lines the vcf file
        if position.strip()[0] != "#": #Check if the line is a comment
            positionDict[alignmentPosition] = position.strip().split("\t")[1]
            alignmentPosition += 1
    return positionDict

def complement(base): #This function takes a base and returns the complement base
    if base == "A":
        return "T"
    elif base == "C":
        return "G"
    elif base == "G":
        return "C"
    elif base == "T":
        return "A"

def getTreeRst(rst): #This function takes an rst file and returns the tree with node labels
    for i,line in enumerate(rst):
        if line == "tree with node labels for Rod Page's TreeView\n":
            return rst[i+1].strip().replace(" ","")

#Adds a label to each tip in a given tree
def labelRstTree(tree, sampleInformation, column):
    #Identify the column of the label of interest
    for i,name in enumerate(sampleInformation[0].strip().split("\t")):
        if name == column:
            labelLocation = i 
    
    #Sample names as keys, labels as values
    sampleDict = {}

    #Iterate through the samples and add to sampleDict
    for sample in sampleInformation[1:]:
        sampleDict[sample.strip().split("\t")[0]] = sample.strip().split("\t")[labelLocation]
    
    #Iterate through the tips and append the patient ID
    for tip in tree.get_terminals():
        #Check if the sample name contains _, if so need to split differently to obtain the tip name after the PAML ID
        if len(str(tip).split("_")) > 2:
            tip.name = tip.name + "____" + sampleDict["_".join(str(tip).split("_")[1:])]
        else:
            tip.name = tip.name + "____" + sampleDict[str(tip).split("_")[1]]
    
    return tree

def getBranchMutations(rst): #This function takes an rst file and returns a dictionary with branches as keys and the changes along the branch as values
    branchDict = {}
    for i, line in enumerate(rst): #Iterate through the rst file
        if line[0:6] == "Branch": #Check if the line starts a mutation section
            branchDict[line.strip().split()[2]] = [] #Will be filled with each of the mutations along the branch
            blankLines = 0 #Will be incremented with each blank line, break with second blank line
            for rstLine in rst[i+1:]: #Iterate through the remaining rst lines until the second blank line which ends the mutation section
                if rstLine == "\n":
                    if blankLines == 1:
                        break
                    blankLines += 1
                else:
                    position = int(rstLine.strip().split(" ")[0]) #Assign to the alignment position of the mutation
                    base = rstLine.strip().split(" ")[1] #Assign to the base at the start of the branch
                    mutant = rstLine.strip().split(" ")[4] #Assign to the base at the end of the branch
                    if base != "-" and mutant != "-":
                        branchDict[line.strip().split()[2]].append([base,position,mutant])
    return branchDict

#Checks if all samples in a clade have the same label after ____ in their names
def getMonophyletic(clade):
    #Labels in the clade
    labels = []

    #Iterate through tips in the clade and add to labels
    for tip in clade.get_terminals():
        labels.append(str(tip).split("____")[-1])
    
    #Check if the clade contains a single label
    if len(set(labels)) == 1:
        return True
    else:
        return False

def getParentNode(tree,clade): #This function takes a tree and a clade within the tree and returns the identifier of the upstream clade
    if len(tree.get_path(clade)) >= 2: #Check if the clade is one downstream of the root
        node = tree.get_path(clade)[-2].confidence
    else:
        node = tree.root.confidence
    return node

def getMutationSpectrum():
    mutation = OrderedDict()
    mutation["ACAA"] = 0
    mutation["ACAC"] = 0
    mutation["ACAG"] = 0
    mutation["ACAT"] = 0
    mutation["CCAA"] = 0
    mutation["CCAC"] = 0
    mutation["CCAG"] = 0
    mutation["CCAT"] = 0
    mutation["GCAA"] = 0
    mutation["GCAC"] = 0
    mutation["GCAG"] = 0
    mutation["GCAT"] = 0
    mutation["TCAA"] = 0
    mutation["TCAC"] = 0
    mutation["TCAG"] = 0
    mutation["TCAT"] = 0
    mutation["ACGA"] = 0
    mutation["ACGC"] = 0
    mutation["ACGG"] = 0
    mutation["ACGT"] = 0
    mutation["CCGA"] = 0
    mutation["CCGC"] = 0
    mutation["CCGG"] = 0
    mutation["CCGT"] = 0
    mutation["GCGA"] = 0
    mutation["GCGC"] = 0
    mutation["GCGG"] = 0
    mutation["GCGT"] = 0
    mutation["TCGA"] = 0
    mutation["TCGC"] = 0
    mutation["TCGG"] = 0
    mutation["TCGT"] = 0
    mutation["ACTA"] = 0
    mutation["ACTC"] = 0
    mutation["ACTG"] = 0
    mutation["ACTT"] = 0
    mutation["CCTA"] = 0
    mutation["CCTC"] = 0
    mutation["CCTG"] = 0
    mutation["CCTT"] = 0
    mutation["GCTA"] = 0
    mutation["GCTC"] = 0
    mutation["GCTG"] = 0
    mutation["GCTT"] = 0
    mutation["TCTA"] = 0
    mutation["TCTC"] = 0
    mutation["TCTG"] = 0
    mutation["TCTT"] = 0
    mutation["ATAA"] = 0
    mutation["ATAC"] = 0
    mutation["ATAG"] = 0
    mutation["ATAT"] = 0
    mutation["CTAA"] = 0
    mutation["CTAC"] = 0
    mutation["CTAG"] = 0
    mutation["CTAT"] = 0
    mutation["GTAA"] = 0
    mutation["GTAC"] = 0
    mutation["GTAG"] = 0
    mutation["GTAT"] = 0
    mutation["TTAA"] = 0
    mutation["TTAC"] = 0
    mutation["TTAG"] = 0
    mutation["TTAT"] = 0
    mutation["ATCA"] = 0
    mutation["ATCC"] = 0
    mutation["ATCG"] = 0
    mutation["ATCT"] = 0
    mutation["CTCA"] = 0
    mutation["CTCC"] = 0
    mutation["CTCG"] = 0
    mutation["CTCT"] = 0
    mutation["GTCA"] = 0
    mutation["GTCC"] = 0
    mutation["GTCG"] = 0
    mutation["GTCT"] = 0
    mutation["TTCA"] = 0
    mutation["TTCC"] = 0
    mutation["TTCG"] = 0
    mutation["TTCT"] = 0
    mutation["ATGA"] = 0
    mutation["ATGC"] = 0
    mutation["ATGG"] = 0
    mutation["ATGT"] = 0
    mutation["CTGA"] = 0
    mutation["CTGC"] = 0
    mutation["CTGG"] = 0
    mutation["CTGT"] = 0
    mutation["GTGA"] = 0
    mutation["GTGC"] = 0
    mutation["GTGG"] = 0
    mutation["GTGT"] = 0
    mutation["TTGA"] = 0
    mutation["TTGC"] = 0
    mutation["TTGG"] = 0
    mutation["TTGT"] = 0

    return(mutation)

#Creates a dictionary with patients as keys and CF status as values
def getPatientStatus(cfStatus):
    patientStatus = {}

    #Iterate through the patients and add to patientStatus
    for patient in cfStatus[1:]:
        patientStatus[patient.strip().split(",")[0]] = patient.strip().split(",")[3]
    
    return(patientStatus)

#Identifies the CF status of a within patient clade
def getLabel(clade, patientStatus):
    labels = [] #Will be filled with the labels of the tips in the clade
    for tip in clade.get_terminals():
        labels.append(patientStatus[str(tip).split("____")[-1]])
    
    #Check if the clade contains a single label
    if len(set(labels)) == 1:
        return(labels[0])
    
    else:
        print("Clade does not contain a single label")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", help = "PAML rst file")
    parser.add_argument("-v", help = "VCF file from gubbins which contains the conversion from alignment position to genome position")
    parser.add_argument("-s", help = "Reference sequence that was mapped against")
    parser.add_argument("-i", help = "Sample information file that contains information to label the tree, names must match those in the tree")
    parser.add_argument("-c", help = "Name of the column containing the patient label, the names in this column should not contain _, default Patient", default = "Patient")
    parser.add_argument("-p", help = "patients.with.multiple.samples.summary.csv containing patient CF status")
    parser.add_argument("-o", help = "Prefix of output file")
    args = parser.parse_args()

    #Import the PAML rst file
    rst = open(args.r).readlines()
    #Import the VCF file
    vcf = open(args.v).readlines()
    #Import the reference sequence
    reference = AlignIO.read(args.s,"fasta")
    #Import the sample information
    sampleInformation = open(args.i).readlines()
    #Import the patient CF statuses
    cfStatus = open(args.p).readlines()

    #Open the output spectrum files
    outFileCF = open(args.o + ".cf.within.patient.mutation.spectrum.txt","w")
    outFileCF.write("Mutation\tNumber.of.mutations\n")
    outFileNonCF = open(args.o + ".non.cf.within.patient.mutation.spectrum.txt","w")
    outFileNonCF.write("Mutation\tNumber.of.mutations\n")
    outFileUnknown = open(args.o + ".unknown.cf.status.within.patient.mutation.spectrum.txt","w")
    outFileUnknown.write("Mutation\tNumber.of.mutations\n")
    outFileLeading = open(args.o + ".leading.to.patient.mutation.spectrum.txt","w")
    outFileLeading.write("Mutation\tNumber.of.mutations\n")
    outFileInternal = open(args.o + ".internal.branch.mutation.spectrum.txt","w")
    outFileInternal.write("Mutation\tNumber.of.mutations\n")
    outFileTerminal = open(args.o + ".terminal.branch.mutation.spectrum.txt","w")
    outFileTerminal.write("Mutation\tNumber.of.mutations\n")

    #Create the mutational spectra that will be filled with mutations in the respective class
    withinCFPatient = getMutationSpectrum()
    withinNonCFPatient = getMutationSpectrum()
    withinUnknownPatient = getMutationSpectrum()
    internalMutation = getMutationSpectrum()
    terminalMutation = getMutationSpectrum()
    leadingPatientMutation = getMutationSpectrum()

    #Alignment positions as keys and the genome positions as values
    vcfPositions = convertVCFPositions(vcf)
    #Extract the tree with node labels from the rst file
    rstTree = p.read(StringIO(getTreeRst(rst)),"newick")
    #Label each of the tips in the tree with the corresponding patient identifier
    rstTreeLabelled = labelRstTree(rstTree,sampleInformation,args.c)
    #Assign to the mutations that occurred along each branch in the tree, branches as keys, mutations as values
    branchMutations = getBranchMutations(rst)
    #Patients as keys, CF status as values
    patientStatus = getPatientStatus(cfStatus)

    #Add a label to determine if the clade has already been analysed
    for clade in rstTreeLabelled.find_clades():
        clade.analysed = "No"

    #Iterate through the clades, check which category they are in and add their mutations to the respective spectrum
    for clade in rstTreeLabelled.find_clades():
        #Do not analyse the root
        if len(rstTreeLabelled.get_path(clade)) != 0:
            #Check if the clade contains multiple samples from a single patient, if so its internal branches are within patient
            if getMonophyletic(clade) and (clade.is_terminal() == False):
                #Do not analyse previously analysed clades
                if clade.analysed == "No":
                    #The CF status of the patient
                    patientLabel = getLabel(clade, patientStatus)
                    #Check if the clade contains nonterminal branches, a clade with 2 samples from the same patient will not
                    if len(clade.get_nonterminals()) > 1:
                        #Iterate through the internal branches in the clade and add the mutations to the respective spectrum
                        for subclade in clade.get_nonterminals()[1:]:
                            #The branch name
                            branch = str(getParentNode(rstTreeLabelled,subclade)) + ".." + str(subclade.confidence)
                            #Iterate through the mutations, get their contexts and add to the spectrum
                            for location in branchMutations[branch]:
                                #The alignment position of the mutation
                                position = location[1]
                                #The genome position of the mutation
                                genomePosition = int(vcfPositions[position])
                                #The base at the start of the branch
                                base = location[0]
                                #The base at the end of the branch
                                mutant = location[2]
                                #Upstream base in the genome
                                upstreamNucleotide = reference[0].seq[genomePosition-2]
                                #Downstream base in the genome
                                downstreamNucleotide = reference[0].seq[genomePosition]
                                #Check if all bases are nucleotides
                                if base != "-" and mutant != "-" and upstreamNucleotide != "N" and upstreamNucleotide != "-" and downstreamNucleotide != "N" and downstreamNucleotide != "-":
                                    #Check the CF status of the patient
                                    if patientLabel == "CF":
                                        #Check if the mutation is included
                                        if (upstreamNucleotide + base + mutant + downstreamNucleotide) in withinCFPatient:
                                            withinCFPatient[upstreamNucleotide + base + mutant + downstreamNucleotide] += 1
                                        else:
                                            withinCFPatient[complement(downstreamNucleotide) + complement(base) + complement(mutant) + complement(upstreamNucleotide)] += 1
                                    elif patientLabel == "Non-CF":
                                        #Check if the mutation is included
                                        if (upstreamNucleotide + base + mutant + downstreamNucleotide) in withinNonCFPatient:
                                            withinNonCFPatient[upstreamNucleotide + base + mutant + downstreamNucleotide] += 1
                                        else:
                                            withinNonCFPatient[complement(downstreamNucleotide) + complement(base) + complement(mutant) + complement(upstreamNucleotide)] += 1
                                    else:
                                        #Check if the mutation is included
                                        if (upstreamNucleotide + base + mutant + downstreamNucleotide) in withinUnknownPatient:
                                            withinUnknownPatient[upstreamNucleotide + base + mutant + downstreamNucleotide] += 1
                                        else:
                                            withinUnknownPatient[complement(downstreamNucleotide) + complement(base) + complement(mutant) + complement(upstreamNucleotide)] += 1
                        #Change the label of each subclade in the clade so they are not analysed again
                        for subcladeLabel in clade.get_nonterminals():
                            subcladeLabel.analysed = "Yes"
                    #Iterate through the terminal branches in the within patient clade
                    for terminalBranch in clade.get_terminals():
                        #The branch name
                        branch = str(getParentNode(rstTreeLabelled,terminalBranch)) + ".." + str(terminalBranch).split("_")[0]
                        #Iterate through the mutations, get their contexts and add to the spectrum
                        for location in branchMutations[branch]:
                            #The alignment position of the mutation
                            position = location[1]
                            #The genome position of the mutation
                            genomePosition = int(vcfPositions[position])
                            #The base at the start of the branch
                            base = location[0]
                            #The base at the end of the branch
                            mutant = location[2]
                            #Upstream base in the genome
                            upstreamNucleotide = reference[0].seq[genomePosition-2]
                            #Downstream base in the genome
                            downstreamNucleotide = reference[0].seq[genomePosition]
                            #Check if all bases are nucleotides
                            if base != "-" and mutant != "-" and upstreamNucleotide != "N" and upstreamNucleotide != "-" and downstreamNucleotide != "N" and downstreamNucleotide != "-":
                                #Check the CF status of the patient
                                if patientLabel == "CF":
                                    #Check if the mutation is included
                                    if (upstreamNucleotide + base + mutant + downstreamNucleotide) in withinCFPatient:
                                        withinCFPatient[upstreamNucleotide + base + mutant + downstreamNucleotide] += 1
                                    else:
                                        withinCFPatient[complement(downstreamNucleotide) + complement(base) + complement(mutant) + complement(upstreamNucleotide)] += 1
                                elif patientLabel == "Non-CF":
                                    #Check if the mutation is included
                                    if (upstreamNucleotide + base + mutant + downstreamNucleotide) in withinNonCFPatient:
                                        withinNonCFPatient[upstreamNucleotide + base + mutant + downstreamNucleotide] += 1
                                    else:
                                        withinNonCFPatient[complement(downstreamNucleotide) + complement(base) + complement(mutant) + complement(upstreamNucleotide)] += 1
                                else:
                                    #Check if the mutation is included
                                    if (upstreamNucleotide + base + mutant + downstreamNucleotide) in withinUnknownPatient:
                                        withinUnknownPatient[upstreamNucleotide + base + mutant + downstreamNucleotide] += 1
                                    else:
                                        withinUnknownPatient[complement(downstreamNucleotide) + complement(base) + complement(mutant) + complement(upstreamNucleotide)] += 1
                    #Change the label of each tip in the clade so they are not analysed again
                    for terminalBranchLabel in clade.get_terminals():
                        terminalBranchLabel.analysed = "Yes"
                    #The name of the branch leading to the within patient cluster to identify the leading to patient mutations
                    cladeBranch = str(getParentNode(rstTreeLabelled,clade)) + ".." + str(clade.confidence)
                    #Iterate through the mutations on the branch leading to the within patient cluster and add to the spectrum
                    for location in branchMutations[cladeBranch]:
                        #The alignment position of the mutation
                        position = location[1]
                        #The genome position of the mutation
                        genomePosition = int(vcfPositions[position])
                        #The base at the start of the branch
                        base = location[0]
                        #The base at the end of the branch
                        mutant = location[2]
                        #Upstream base in the genome
                        upstreamNucleotide = reference[0].seq[genomePosition-2]
                        #Downstream base in the genome
                        downstreamNucleotide = reference[0].seq[genomePosition]
                        #Check if all bases are nucleotides
                        if base != "-" and mutant != "-" and upstreamNucleotide != "N" and upstreamNucleotide != "-" and downstreamNucleotide != "N" and downstreamNucleotide != "-": #Check if the mutation involves two bases
                            #Check if the mutation is included
                            if (upstreamNucleotide + base + mutant + downstreamNucleotide) in leadingPatientMutation:
                                leadingPatientMutation[upstreamNucleotide + base + mutant + downstreamNucleotide] += 1
                            else:
                                leadingPatientMutation[complement(downstreamNucleotide) + complement(base) + complement(mutant) + complement(upstreamNucleotide)] += 1
                    #Change the label of the clade to analysed
                    clade.analysed = "Yes"

            #Check if the clade is internal and is not the ancestor of a monophyletic within patient cluster
            elif clade.is_terminal() == False and clade.analysed == "No":
                #The branch name
                branch = str(getParentNode(rstTreeLabelled,clade)) + ".." + str(clade.confidence)
                #Iterate through the mutations on the branch and add to the spectrum
                for location in branchMutations[branch]:
                    #The alignment position of the mutation
                    position = location[1]
                    #The genome position of the mutation
                    genomePosition = int(vcfPositions[position])
                    #The base at the start of the branch
                    base = location[0]
                    #The base at the end of the branch
                    mutant = location[2]
                    #Upstream base in the genome
                    upstreamNucleotide = reference[0].seq[genomePosition-2]
                    #Downstream base in the genome
                    downstreamNucleotide = reference[0].seq[genomePosition]
                    #Check if all bases are nucleotides
                    if base != "-" and mutant != "-" and upstreamNucleotide != "N" and upstreamNucleotide != "-" and downstreamNucleotide != "N" and downstreamNucleotide != "-":
                        #Check if the mutation is included
                        if (upstreamNucleotide + base + mutant + downstreamNucleotide) in internalMutation:
                            internalMutation[upstreamNucleotide + base + mutant + downstreamNucleotide] += 1
                        else:
                            internalMutation[complement(downstreamNucleotide) + complement(base) + complement(mutant) + complement(upstreamNucleotide)] += 1
            
            #Check if the clade is a tip branch that hasn't been analysed
            elif clade.is_terminal() and clade.analysed == "No":
                #The branch name
                branch = str(getParentNode(rstTreeLabelled,clade)) + ".." + str(clade).split("_")[0]
                #Iterate through the mutations on the branch and add to the spectrum
                for location in branchMutations[branch]:
                    #The alignment position of the mutation
                    position = location[1]
                    #The genome position of the mutation
                    genomePosition = int(vcfPositions[position])
                    #The base at the start of the branch
                    base = location[0]
                    #The base at the end of the branch
                    mutant = location[2]
                    #Upstream base in the genome
                    upstreamNucleotide = reference[0].seq[genomePosition-2]
                    #Downstream base in the genome
                    downstreamNucleotide = reference[0].seq[genomePosition]
                    #Check if all bases are nucleotides
                    if base != "-" and mutant != "-" and upstreamNucleotide != "N" and upstreamNucleotide != "-" and downstreamNucleotide != "N" and downstreamNucleotide != "-":
                        #Check if the mutation is included
                        if (upstreamNucleotide + base + mutant + downstreamNucleotide) in terminalMutation:
                            terminalMutation[upstreamNucleotide + base + mutant + downstreamNucleotide] += 1
                        else:
                            terminalMutation[complement(downstreamNucleotide) + complement(base) + complement(mutant) + complement(upstreamNucleotide)] += 1
    
    #Write the spectra
    for nucleotideCF in withinCFPatient:
        outFileCF.write(nucleotideCF + "\t" + str(withinCFPatient[nucleotideCF]) + "\n")
    for nucleotideNonCF in withinNonCFPatient:
        outFileNonCF.write(nucleotideNonCF + "\t" + str(withinNonCFPatient[nucleotideNonCF]) + "\n")
    for nucleotideUnknown in withinUnknownPatient:
        outFileUnknown.write(nucleotideUnknown + "\t" + str(withinUnknownPatient[nucleotideUnknown]) + "\n")
    for leadingNucleotide in leadingPatientMutation:
        outFileLeading.write(leadingNucleotide + "\t" + str(leadingPatientMutation[leadingNucleotide]) + "\n")
    for internalNucleotide in internalMutation:
        outFileInternal.write(internalNucleotide + "\t" + str(internalMutation[internalNucleotide]) + "\n")
    for terminalNucleotide in terminalMutation:
        outFileTerminal.write(terminalNucleotide + "\t" + str(terminalMutation[terminalNucleotide]) + "\n")
    
    outFileCF.close()
    outFileNonCF.close()
    outFileUnknown.close()
    outFileLeading.close()
    outFileInternal.close()
    outFileTerminal.close()