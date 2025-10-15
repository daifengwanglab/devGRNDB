import pandas as pd
import numpy as np
import ast
import re

type = 'Telencephalon'
linType = 'InN'
species = 'human'

file = species+'/'+type+'/'+type+'_'+linType


# Load ctx file
ctx = pd.read_csv(file+'_ctx.csv')  # Should contain a column like "TF" and "enrichment_target_genes"
ctx = ctx.drop([0,1])
ctx = ctx[['Unnamed: 0','Enrichment.6']]
# Adjust the column names if needed
ctx.columns = ['TF', 'TargetGenes']
df = ctx

final_df = pd.DataFrame()
for idx,row in df.iterrows():
    tf = row["TF"]
    x = row["TargetGenes"]
    if x[-1]!="]":
        x = x+"]"

    # Replace np.float64(...) with the float value
    s_clean = re.sub(r'np\.float64\((.*?)\)', r'\1', x)

    # Convert string to list of tuples
    tuple_list = ast.literal_eval(s_clean)


    temp = pd.DataFrame(tuple_list, columns=['TG', 'Value'])
    temp['TF'] = tf
    final_df = pd.concat([final_df, temp])

print(final_df.shape)
final_df = final_df[["TF","TG","Value"]]
df_unique = final_df.drop_duplicates()
print(df_unique.shape)
df_unique.to_csv(file+'_ctx_cleaned.csv')

