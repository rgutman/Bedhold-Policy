# Bedhold-Policy
Includes code for the paper A Bayesian Procedure for Estimating the Causal Effects of Nursing Home Bed-hold Policy


This project includes several files:

1. simulatedDataSetForPaper.csv - This is a file that include simulated data. The variables are as follows:
    * acute_qtr - Number of acute cases observed at quarter
    * died_qtr - Number of deaths observed in a quarter
    * nhmod - Number of residents' months in a quarter
    * facid - facility ID
    * bedhold- an indicator representing if a bedhold policy was enacted in that quarter for that facility
    * tobed - quintile of facility in regards to total number of beds 5-highest 1-lowest
    * BedHoldChange - an indicator whether the facility was in a state that enacted bedhold policy
    * state - a state ID for the facility
    * year - year of the ovbservation
    * quarter- quarter of the observation within year.

2) simulatedDataCode.R - R code to run the main analysis. It generates results table and three graphs that are reported in the supplementary material of the paper. It reads file simulatedDataSetForPaper.csv so the file should be available.

3) 
