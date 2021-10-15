### Housekeeping
# 
# specify path to folder
path = '/Users/roos_brouns/Dropbox/Ant-fungus/02_scripts/Git_Das_folder2/Das_et_al_2022a'
#
# specify species
species = 'ophio_cflo'
#
# load modules
import sqlite3 # to connect to database
import pandas as pd # data anaylis handeling
import numpy as np  # following three are for plotting
import matplotlib.pyplot as plt
import seaborn as sns
import colorsys
import matplotlib
#
# Add defenition for color palette making for plotting. 
# source: source: https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
def scale_lightness(rgb, scale_l):
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

### Load in Data
#
# Connect to the database
conn = sqlite3.connect(f'{path}/data/databases/new_TC6_fungal_data.db')
#
# read data from DB into dataframe  
exp_val = pd.read_sql_query(f"SELECT * from {species}_fpkm", conn)
#
### Clean data
#
# drop 'start' and 'end' columns
exp_val.drop(['start','end'], axis=1, inplace=True)

### Part 01. Plot expressioin value of modules
#
# A. Get the data you want to plot
#
# Load in the data with the gene ID's and their assigned module
all_modules = pd.read_csv(f'{path}/results/networks/{species}_gene_IDs_and_module_identity.csv')
# Clean the data by dropping exsessive columns
all_modules.drop(['start', 'end'], axis=1, inplace=True)
#
#
# Input the colors of the modules you want to plot
# source : https://pynative.com/python-accept-list-input-from-user/
print("Attention: If the code does not work and a error message is given, this means that the colors you gave are not defined as modules")
input_string = input('Enter elements of a list separated by space ')
print("\n")
module_list = input_string.split()
# print list
print('giving input colors', module_list)

if len(module_list) == 1:
    # Define the module you want to plot
    module_name = module_list[0]
    # select all the data of that module
    module_data = all_modules.loc[all_modules['module_identity'] == module_name]
    # select the gene ID's in of the module
    module_IDs = pd.DataFrame(module_data['gene_ID_robin'])
    # get the expression values of the selected gene ID's in the module
    exp_val_module = module_IDs.merge(exp_val, on='gene_ID_robin', how='left')
    #

    # B. Plot the data
    # 
    # Transform the data so we can plot
    t_exp_val_module = exp_val_module.T
    t_exp_val_module.drop(['gene_ID_robin','gene_ID_ncbi'], axis=0, inplace=True)
    #
    # Make the color palette for plotting
    color = matplotlib.colors.ColorConverter.to_rgb(module_name)
    rgbs = [scale_lightness(color, scale) for scale in [0.3, .6, 1, 1.1, 1.25]]
    # Show the palette color
    # sns.palplot(rgbs) 
    #
    # Set background and the color palette
    sns.set_theme()
    sns.set_palette(rgbs)
    # Plot the gene expression values against the time, without legend and with titles
    ax = t_exp_val_module.plot(legend=False)
    ax.set_title(species[0].upper()+f'{species[1:len(species)]}-{module_name} gene expression', fontsize=15, fontweight='bold')
    ax.set_xlabel('Time point')
    ax.set_ylabel('Expression value (FPKM)')
    #
    ### Done.

elif len(module_list) > 1:
    
    for color in module_list:

        # Define the module you want to plot
        module_name = color
        # select all the data of that module
        module_data = all_modules.loc[all_modules['module_identity'] == module_name]
        # select the gene ID's in of the module
        module_IDs = pd.DataFrame(module_data['gene_ID_robin'])
        # get the expression values of the selected gene ID's in the module
        exp_val_module = module_IDs.merge(exp_val, on='gene_ID_robin', how='left')
        #

        # B. Plot the data
        # 
        # Transform the data so we can plot
        t_exp_val_module = exp_val_module.T
        t_exp_val_module.drop(['gene_ID_robin','gene_ID_ncbi'], axis=0, inplace=True)
        #
        # Make the color palette for plotting
        color = matplotlib.colors.ColorConverter.to_rgb(module_name)
        rgbs = [scale_lightness(color, scale) for scale in [0.3, .6, 1, 1.1, 1.25]]
        # Show the palette color
        # sns.palplot(rgbs)
        #
        # set background and the color palette
        sns.set_theme()
        sns.set_palette(rgbs)
        #
        # Plot the gene expression values against the time, without legend and with titles
        ax = t_exp_val_module.plot(legend=False)
        ax.set_title(species[0].upper()+f'{species[1:len(species)]}-{module_name} gene expression', fontsize=15, fontweight='bold')
        ax.set_xlabel('Time point')
        ax.set_ylabel('Expression value (FPKM)')
        #
        ### Done.

else:
    print('No input is given')