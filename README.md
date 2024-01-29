# REE-pH

**This is the repository for REE-pH, the REE-pH proxy in the paper by Lin and Catling (in review). This work has not been published in a peer-reviewed journal yet. Please email me before you use the code to make sure that the pH values are generated correctly (pcl1225@uw.edu)**

All scripts are written in Python, except REE_speciation.pqi. The Python scripts can be run independently to produce results without PHREEQC. If you would like to use the PHREEQC model, please refer to the PHREEQC section below.

## How to use: 
To calibrate the pH proxy using modern seawater data, you will need GLODAPv2.2022_Merged_Master_File.csv from the GLODAPv2 website: https://www.glodap.info/index.php/merged-and-adjusted-data-product-v2-2022/ (Lauvset et al. 2022). Please download the CSV and place it in the pH_data folder.

To estimate seawater pH based on REE concentration in rocks, please import the REE concentrations of marine carbonates in CSV format. Please refer to the CSV files in the REE_data folder for the spreadsheet format. In the code, we use the limestone samples from Toyama and Terakado (2019) as the example (REE_limestone_toyama_noFJ1YK1_ugg.csv). You can change the CSV files to other formations to estimate different pH values. 

The CSV files in the coeff folder can be replaced by other coefficient values. You can change them to values from other literature. For example, we use Post-Archean Australian Shale (PAAS; Pourmand et al. 2012) to normalize the REE concentrations. The REE concentrations for PAAS are in paas.csv. You can change the paas.csv file to other shale normalizing standards by importing other spreadsheets as input.

## PHREEQC:

REE_speciation.pqi is the input file for PHREEQC, a USGS software available at https://www.usgs.gov/software/phreeqc-version-3. The pqi file is only for PHREEQC with the REE scavenging model modified from (Schijf et al. 2015). The dat file is the database file for the equations (Pitzer).

After installing PHREEQC, you can directly use REE_speciation.pqi. The first line of the file refers to the path of the database on the user end. You can change it to the local path where you've installed the PHREEQC software. 
Directly put the file in PHREEQC, and you should be able to load the pqi file. After running the model, PHREEQC will produce multiple output files, including selected_output_all.sel. The REE speciation data from ree_speciation.csv (the input data for speciation.py) is a subset of trivalent ions (REE3+), monocarbonato complex (REECO3+), and dicarbonato complex (REE(CO3)2-) from selected_output_all.sel.


## Contact & Publications
If you have questions about the code, please email the corresponding author: pcl1225@uw.edu





