# REE_pH
This is the repository for the REE-pH proxy in the paper Lin and Catling (2023).

Speciation.py produces the plots for speciation.
(PHREEQC SCRIPT) is the script for PHREEQC. Please refer to https://www.usgs.gov/software/phreeqc-version-3 to install and run the software. The sciprt can be directly used in the input. The output should be similar to 

REE_pH_proxy.py is the script for the REE-pH proxy model calibration mentioned in Figure 3c of Lin and Catling (2023). It generates the linear regression model of seawter based on the average values of REE concentration and pH from seawater measurements.

# How to use: 
To estimate seawater pH based on REE conctration in rocks, please import the REE concentrations in .csv format. Please refer to the test files for the spreadsheet format. 







