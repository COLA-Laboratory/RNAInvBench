import pandas as pd

df = pd.read_csv("antaRNA_GCHotknot.out", header=None)
pd.set_option("display.max_columns", None)
print(df.head())
df['structure'] = df[0].str.split(':').str[0].str.strip()
df['structure'] = df['structure'].str.replace('{','').str.replace("'","")

df2 = pd.read_csv("full_sorted_pk.csv", delimiter=",", header=None)
print(df2)
real_structs = df2[2].to_list()

structures = df['structure'].to_list()
# Calculate the amount of correct structures pls?
counter = 0
for p_struct, r_struct in zip(structures, real_structs):
    if p_struct == r_struct:
        counter += 1

print("correct: ", counter)

print("total: ", len(real_structs))
print(len(structures))