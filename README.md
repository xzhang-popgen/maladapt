# maladapt

#######MALADAPT PIPELINE#######
The development of MaLAdapt requires the following steps:

#1. Create simulations of archaic introgression as training data

SLiM templates for three types of DFEs (additive, recessive, partial recessive) can be found in the slim/ directory

Under the python/ directory, the following scripts should be used under the numeric order:
	1slim_simulations.py #create the simulations, extract features in 50kb overlapping windows
	2compile_simulations.py #combine all simulations from different DFE simulations into a large dataframe
	3adding_exon+r_info_simulations.py #add the separately computed exon density and recombination rate information to the aforementioned combined dataframe
	4downsample_class-ratio_simulations.py #downsample the non-AI class windows (2:1 ratio)

For running the slim simulation from bash script (UCLA hoffman2 cluster), a sample job script can be found in "sample_bash_command.sh"

For the non-AI selective sweep simulations, the slim and python scripts can be found under the nonAIsweep/ directory in the slim/ directory

#2. Train MaLAdapt model
The script used for training MaLAdapt can be found in the python/ directory ("5trainMaLAdapt.py")


#3. Assessment of MaLAdapt performance under parameter-misspecified scenarios

slim templates and python simulation scripts can be found in misspec_sims/ directory


#4. Apply MaLAdapt to 1000 genomes empirical data

The empirical/ directory includes the python scripts used to compute features from 1KG populations, add the exon density/recombination rate features, and the application of MaLAdapt prediction

All scripts should be used under the numeric order


#5. Additional files on Google Drive (https://drive.google.com/drive/folders/10r8e5WbhcgAIjC0DVmIe4saVYODRgFCO?usp=share_link)

	1. A pretrained model of MaLAdapt (published version in Zhang et al. 2023 MBE)
	2. 1000 genomic segments used in MaLAdapt training simulations


#######MALADAPT PIPELINE REQUIRED MODULES#######
Required python3 modules:

	1. msprime==0.7.0
	2. tskit==0.1.5
	3. pyslim==0.401

Required slim version: 3.2.0

