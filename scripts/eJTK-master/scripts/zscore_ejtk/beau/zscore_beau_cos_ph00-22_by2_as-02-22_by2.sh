# This script runs rhythmicity analysis (eJTK)
# 
#################################################
## Do the following before you run this script ##
#################################################
##
# 01. Go to:
	# "~/Documents/GitHub/Das_et_al_2022a/scripts/eJTK-master/"
# 02. Run:
	# chmod 755 eJTK-CalcP.py
	# cd bin
	# python setup.py build_ext --inplace
	# cd ..
## 
## status: DONE
#################################################
## Find more information in the protocol for eJTK
## File name: Running_eJTK_on_Mac_BD.pdf
#################################################

## 01. HOUSEKEEPING
##
# path to the directory with normalized gene expression
species=ophio_cflo
billu_path=Documents/GitHub
roos_path=Dropbox/Ant-fungus/02_git/Git_Das_folder2
data_path=~/${roos_path}/Das_et_al_2022a/results/normalized_gene_exp/zscore/${species}/
# specify the data file 
#################################################
	# species: ophio_cflo
	# dataset: ld (z-score)
data=${species}_zscores_noNAs.txt
#################################################
# Define the period you want to analyse
period_number=24
# specify the identifier for the output files 
##			 # (Billu's file naming system for period of 8h is cos08. However, the ref_file for this period is named period8. So the if statement is build to get the exact same naming system Billu has for the files)
if [[ ${period_number} = 8 ]]; then
	output_name=cos0${period_number}_ph0022by2_as0222by2_zscore
else
	output_name=cos${period_number}_ph0022by2_as0222by2_zscore
fi
# specify the name of the output directory (will create one, if absent)
output_dir=ejtk_output/${species}/zscore/


## 02. Setup rhythmicity-test-conditions
#
# Specify the arguments for the eJTK-CalcP.py function

args=(
	# data file with normalized gene expression (fpkm/rpkm or z-score)
	-f ${data_path}${data} \
	# test waveform (default: cosine)
	-w ref_files/waveform_cosine.txt \
	# periodicity to test
	-p ref_files/period${period_number}.txt \
	# find phases at specified ZTs
	-s ref_files/phases_00-22_by2.txt \
	# allowed asymmetries for rhythmic waveform
	-a ref_files/asymmetries_02-22_by2.txt \
	# identifier for the output files
	-x ${output_name}
)

## 03. Activate conda environment 
##
# change directory to the folder containing eJTK-Calc.py
cd ~/${roos_path}/Das_et_al_2022a/scripts/eJTK-master/
#
# Activate the conda environment with python-version2
. /Users/roos_brouns/opt/anaconda3/bin/activate && conda activate /Users/roos_brouns/opt/anaconda3/envs/py2
#
##
## 04. Run the rhythmicity analysis
##
./eJTK-CalcP.py "${args[@]}"
##
## 05. Move output to a separate folder
##
cd $data_path
mv *${output_name}* ../../../${output_dir}
# else
# 	mkdir $output_dir
# 	mv *${output_name}* $output_dir
# fi
#
#
### END