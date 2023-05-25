# REE_pH

**This is the repository for the REE-pH proxy in the paper Lin and Catling (in review). This work has not been published in a peer-reviewed journal yet. Please do not use this model in a paper that will be published. Please mmail me if you have questions (pcl1225@uw.edu)**

All scripts are written in python, except REE_speciation.pqi. The python scipts can be run indepdently to produce results without PHREEQC. If you would like to use the PHREEQC model, please refer to the PHREEQC section below.

## How to use: 
To calibrate the pH proxy using modern seawater data, you will need GLODAPv2.2022_Merged_Master_File.csv from the GLODAPv2 website: https://www.glodap.info/index.php/merged-and-adjusted-data-product-v2-2022/ (Lauvset et al. 2022). Please download the CSV and place it in the pH_data folder.

To estimate seawater pH based on REE conctration in rocks, please import the REE concentrations of marine carbonates in CSV format. Please refer to the CSV files in the REE_data folder for the spreadsheet format. In the code, we use the limestone samples from Toyama and Terakada (2019) as the example (REE_limestone_toyama_noFJ1YK1_ugg.csv). You can change the CSV files to other formations to estiamte pH values from corresponding literature. 

The CSV files in the coeff folder can be replaced by other coefficient values. You can change them to values from other literature. For example, we use Post-Archean Australian Shale (PAAS; Pourmand et al. 2012) as the REE concentration to normalize the REE conctrations. The REE contrations for PAAS are in paas.csv. You can change the paas.csv file to other shale normalizing standards by importing another spreadsheet as an input.

## PHREEQC:

REE_speciation.pqi is the input file for PHREEQC, a USGS software avaliable at https://www.usgs.gov/software/phreeqc-version-3. PHREEQC is only required to run the PHREEQC script with the REE scavenging model modified from (Schijf et al. 2015).

REE_speciation.pqi is the input file for PHREEQC. After installing PHREEQC, you can directly use this file. The first line of the input file refers to the location of the database on the user end. You can change it to the local path where you've installed the PHREEQC software. After running the model, PHREEQC will produce multiple output files, including selected_output_all.sel. The REE speciation data from ree_speciation.csv (the input csv for speciation.py) is a subset from selected_output_all.sel.


## Contact & Publications
If you have questions about the code please email the corresponding author: pcl1225@uw.edu





