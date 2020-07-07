# Hybridization_and_stress
For the data analysis of this paper we provide:
	A) 9 Tables in .csv format
	B) 2 scripts (For main and supplementary figures)


A) The tables contain the data compiled in « .csv » (raw data can be requested from the author). Each table represents an experiment and is associated with the figures. Below we summarize the headings of the tables (some are common among tables). The tables are:
-Fig_1.csv
-Fig_2.csv
-Fig_3.csv
-Fig_S2.csv
-Fig_S3A.csv
-Fig_S3B.csv
-Fig_S4.csv
-Fig_S6.csv
-Fig_S7.csv


Headeds of the tables:
	-well: well disposition in the plate
	-day: day of the experimental evolution (1-21)
	-plate: plate number 
	-tecan: instrument used
	-time: at which measurements was taken (in seconds)
	-od: optical density at 595 nm
	-strain: genotype
	-NQO_c: concentration of the chemical
	-rep/repi:  replicate number
	-hour: time in hours
	-NQO: whether the well contains 4-NQO or not
	-gen: number of generations elapsed (200 corresponds to 100)
	-evol: condition of evolution
	-corr: correlation between the growth rate (r) and the maximum population size (K) parameters
	-note_fit: fit of the growthcurver package for each growth curve
	-rval: the maximum growth rate (r)
	-max_size: the maximal population size (K)
	-check: check the ones are under 0.2 OD
	-Final_perc (“Division” and “Final columns” were used to get this column): percentage of growth rate compared to same strain in normal conditions
	-Final_perc_check: removing from the percentage of growth rate compared to same strain in normal conditions in the ones which are under 0.2 OD

The headers that are not mentioned are not used during the scripts presented.

B) 2 scripts:
	-Main figures (« Review_Main_figures_paper_CBR.R » )
	-Supplementary figures (« Review_Supp_figures_paper_CBR.R » )

INPUT: These scripts require that the user adds the path to the directory containing the data files in setwd(""). Afterwards, the scripts can be executed without additional arguments because they use the relative paths to each of the data sets.

OUTPUT: The user will receive as output all the figures in "pdf" and "jpg" format in the directory containing the data tables.
