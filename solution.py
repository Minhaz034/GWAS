import pandas as pd
import numpy as np
import csv
from scipy.stats import fisher_exact

data_path = './1002154424.csv'
df = pd.read_csv(data_path)


results = []

# Iterate over each row (SNP)
for index, row in df.iterrows():
    contingency_table = [
        [row['Case_Num_C_Allele'], row['Case_Num_T_Allele']],
        [row['Control_Num_C_Allele'], row['Control_Num_T_Allele']]
    ]
    # Perform Fisher's exact test
    odds_ratio, p_value = fisher_exact(contingency_table, alternative='two-sided')

    # Check significance (You may choose alternative='greater' based on your research questions)
    is_significant = p_value < 5e-8

    # Append to results
    results.append([row['SNP'], p_value, is_significant])

# Write results to a CSV file
with open('results.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['SNP', 'p_value', 'is_significant'])
    csvwriter.writerows(results)


print(len(contingency_table))
