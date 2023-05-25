# REE_pH
This is the repository for the REE-pH proxy in the paper Lin and Catling (in review).
All scripts are written in python, except the REE_speciation.pqi. REE_speciation.pqi is the input file for PHREEQC, a USGS software avaliable at https://www.usgs.gov/software/phreeqc-version-3. PHREEQC is only required to run the PHREEQC script with the modeling parameters modified from (Schijf et al. 2015).

## How to use: 
To estimate seawater pH based on REE conctration in rocks, please import the REE concentrations in .csv format. Please refer to the test files in the Test folder for the spreadsheet format. To read in csv

## Speciation.py
Speciation.py visualizes the change of REE concentrations as the pH changes in seawater. This script recreates Supplementary Figure 2, the Heavy/ Light REE ratio chance from pH 5.7 to 8.7. The input file is ree_speciation_new.csv is from the PHREEQC model using REE_speciation.pqi. 

## REE_speciation.pqi 
REE_speciation.pqi is the input file for PHREEQC. After installing PHREEQC, you can directly use this file. The first line of the input file refers to the location of the database on the user end. You can change it to the loca file where you've installed the PHREEQC software. After finishing the run, PHREEQC will produce multiple output files, and the selected_output_all.sel. The REE speciation data from ree_speciation_new.csv (the input csv for Speciation.py) is a subset from selected_output_all.sel.

## REE_pH_proxy.py
REE_pH_proxy.py is the script for the REE-pH proxy model calibration mentioned in Figure 3c of the manuscript. It generates the linear regression model of seawater based on the average values of REE concentration and pH from seawater measurements.

For the pH data, you will need GLODAPv2.2022_Merged_Master_File from the GLODAPv2 website: https://www.glodap.info/index.php/merged-and-adjusted-data-product-v2-2022/ (Lauvset et al. 2022).

The values in  DREE.csv, molar mass_mgmol.csv, atomic_number.csv and are interchangable. We use Post-Archean Australian Shale (PAAS; Pourmand et al. 2012) as the REE concentration to normalize the REE conctrations. You can change the DREE.csv file to other shale normalizing standards by changing the values on the spreadsheet.


## REE_pH_estimates.py
In the code, we use the limestone REE concentration from Toyama and Terakada (2019) as the example (REE_limestone_toyama_noFJ1YK1_ugg.csv). You can change the csv files to other formations to estiamte pH values from corresponding literature. 


## Contact & Publications
If you have questions about the code please email the corresponding author: pcl1225@uw.edu





