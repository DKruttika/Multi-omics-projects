import sys
import os
sys.executable

import scipy.io
base_dir = '/path/box/folder'
filename = str(sys.argv[1])
filename = str(filename.split('.')[0])

# load the .mat file
mat = scipy.io.loadmat(base_dir + filename + '_norm_intensity.mat')

# extract key 'header'
values = [item for sublist in mat['header'][0] for item in sublist]
header = values
print(header)

# convert the 'X' to a dataframe
import pandas as pd
df = pd.DataFrame(mat['X'])

# how to access every elements of array of arrays and assign as column header of data frame df
df.columns = header

df.head()





# reading in csv file
import pandas as pd
df_csv = pd.read_csv(base_dir + filename + '.dat')
df_csv.head()
# add last column from df_csv to df as a new column
df['cell_type'] = df_csv.iloc[:,-1]
df['cell_number'] = df_csv.iloc[:,0]
df['ROI_ID'] = df_csv.iloc[:,5]

# subset Object_ID if value in column cell_subtyping_v1_d2 is 'Tumor cells'
df_filtered = df[df['ROI_ID'].str.contains('Tumor')]
df_filtered = df_filtered[df_filtered['cell_type'] == 'Tumor cells']

# remove the last three column from df
df_filtered = df_filtered.iloc[:,:-3]
df_filtered.head()

# calculate column means except for the last three columns and store in dataframe
df_means = df_filtered.mean(axis=0)
df_means = pd.DataFrame(df_means)
df_means.columns = ['mean_intensity']

# give name to index column
df_means.index.name = 'cell_type'


# parse file name to get the patient ID and then paste together
sample_name = '_'.join(filename.split('_')[6:9])

df_means['filename'] = sample_name

# write out the data frame to a text file
df_means.to_csv(filename + '_norm_intensity_celltype.csv', index=True)
