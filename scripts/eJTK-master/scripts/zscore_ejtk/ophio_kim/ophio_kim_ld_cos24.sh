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
##
## 01. HOUSEKEEPING
##
# path to the directory with normalized gene expression
data_path=~/Documents/GitHub/Das_et_al_2022a/results/normalized_gene_exp/zscore/ophio_kim/
# specify the data file 
	# species: ophio_kim
	# dataset: ld
data=ophio_kim_LD_zscores_noNAs.txt
# specify the identifier for the output files
output_name=cos24_ph0020by4_as0420by4
# specify the name of the output directory (will create one, if absent)
output_dir=ejtk_output
#
##
## 02. Setup rhythmicity-test-conditions
##
# Specify the arguments for the eJTK-CalcP.py function
#
args=(
	# data file with normalized gene expression (fpkm/rpkm or z-score)
	-f ${data_path}${data} \
	# test waveform (default: cosine)
	-w ref_files/waveform_cosine.txt \
	# periodicity to test
	-p ref_files/period24.txt \
	# find phases at specified ZTs
	-s ref_files/phases_00-20_by4.txt \
	# allowed asymmetries for rhythmic waveform
	-a ref_files/asymmetries_04-20_by4.txt \
	# identifier for the output files
	-x ${output_name}
)
#
##
## 03. Activate conda environment 
##
# change directory to the folder containing eJTK-Calc.py
cd ~/Documents/GitHub/Das_et_al_2022a/scripts/eJTK-master/
#
# Activate the conda environment with python-version2
. /opt/anaconda3/bin/activate && conda activate /opt/anaconda3/envs/ejtk_py2
#
##
## 04. Run the rhythmicity analysis
##
./eJTK-CalcP.py "${args[@]}"
##
## 05. Move output to a separate folder
##
cd $data_path
if [ -d "$output_dir" ]; then
    mv *${output_name}* $output_dir
else
	mkdir $output_dir
	mv *${output_name}* $output_dir
fi



	