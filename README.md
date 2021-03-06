# DCC_clade_emergence
Data and scripts used in DCC clade emergence manuscript

Repository DOI: [![DOI](https://zenodo.org/badge/313015200.svg)](https://zenodo.org/badge/latestdoi/313015200)

The data within the data folder are:
alignments - contains alignments for each FastBAPS cluster. Named as subspecies_clusterNumber so abscessus_2.fasta is the alignment for abscessus subspecies FastBAPS cluster 2

trees - phylogenetic trees for each FastBAPS cluster

temporal_trees - BEAST MCC trees for each DCC

BEAST_XMLs - XML files used in BEAST analyses. These files are as follows:
DCC1.xml - informed substitution rate prior from TempEst analysis, no phylogeography. Used to reconstruct temporal and population history
DCC1_uniform.xml - uniform substitution rate prior between 10^-5 and 10^-9, no phylogeography
DCC1_date_randomisation1.xml - dates randomised, uniform substitution rate prior. Used for date randomisation test
DCC1_phylogeography.xml - asymmetric phylogeography
DCC1_phylogeography_sample1.xml - asymmetric phylogeography with over-represented continents downsampled

mutational_spectra - counts of context-specific mutations used in mutational spectrum analyses. non_DCC_internal_spectrum.csv is the mutation counts for the non-DCC internal branches. DCC_internal_spectrum.csv is the mutation counts for the DCC internal branches. within_CF_patient_spectrum.csv is the mutation counts for within CF patient mutations

sequence_metadata.xlsx contains metadata for all sequences included in the study

The scripts in the scripts folder are:
mutationalSpectrum.v1.py - used to reconstruct the mutational spectrum from a PAML rst file

combineMutationalSpectra.R - combines multiple mutational spectra into a single spectrum

bootstrapMutationSpectrum.R - calculates the significance of differences in mutational spectra using a permutation test

plotSpectrum.R - plots a mutational spectrum

plotMCCTrees.R - plots a MCC tree coloured by inferred location

mapDCCMigrationsMap.R - plots supported DCC migration rates on a world map

plotTransmissionNetwork.R - plots a transmission network at a given SNP cutoff

extractPopulationIncrease.py - used to extract the dates of increases in relative genetic diversity from a posterior distribution of trees
