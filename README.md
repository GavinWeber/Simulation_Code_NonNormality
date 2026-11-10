# Repository Guide
All files included in this repository are simulation code or results for the paper "Operating Characteristics of Welch’s t-test, the Brunner-Munzel test and the C2 test Under Severe Non-Normality". The files are organized as follows: 

### Simulation Code
The file Simulation_R_Code.R is the main file for running the Monte Carlo simulation described in the paper. It was structured to run on a Linux-based server that has R installed. This file's purpose includes: 
- Generating custom functions that simulate the real-world distributions described in the paper
- Defining parameters and controls for the simulation procedure
- Running the Monte Carlo simulation
- Compiling the results
- Generating the images and tables to analyze the results
Thus, anyone wishing to reproduce the complete results should use this code.

The paper notes that this simulation was ran several times with different seeds. These seeds were 123, 13, and 6098. 


### Simulation Results
The full simulation takes several hours to run on a server, therefore we have included the results from several runs of the simulation that can be used in conjunction with the code to reproduce the figures shown in the paper. These files are located in the Results folder in this repository. 
The following files can be found in the Results folder: 
- Simulation_Results.rds (the results of the simulation using seed 123, can be split up into individual datasets for each distribution using list2env() in R)
- Appendix_Codebook.xlsx (an Excel file containing an explanation of each column in the other Appendix Excel files)
- Appendix_AllResults_1.xlsx (an Excel file containing the results of each iteration of Welch's t-test, the Brunner-Munzel test, and the C2 test under both null and alternate hypotheses for all distributions generated using seed 123)
- Appendix_AllResults_2.xlsx (same as Appendix_AllResults_1.xlsx, but distributions were generated using seed 13)
- Appendix_AllResults_3.xlsx (same as Appendix_AllResults_1.xlsx, but distributions were generated using seed 6098)
