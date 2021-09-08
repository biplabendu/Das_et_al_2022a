#!/bin/bash
#SBATCH -J data_structure 
#SBATCH -t 4:0:0
#SBATCH --mem=5G

### Data Structure
# This script creates the directories to organize the data. It is important to use this created directory structure to execute the other scripts. 
#
#
## Here we defined the absolute path to where the directory needs to be created and create a directory for the experiment
# example, in this case Time-course experiment 6 (TC6)
path = /home/uu_bio_fg/rbrouns/data/ant_fungus # don’t add a / at the end of the path
experiment_name = TC6
mkdir $path/TC6/
#
# In the experiment directory, create a directory for storing the data, scripts and the results
mkdir  $path/TC6/data/
mkdir  $path/TC6/scripts
mkdir  $path/TC6/results/
#
# In both the data and result directory, make a directory for the species we are analyzing. 
# in this case O. camp-florani.
species_name = ophio
mkdir  $path/TC6/data/ophio/
mkdir  $path/TC6/results/ophio/
#
### END