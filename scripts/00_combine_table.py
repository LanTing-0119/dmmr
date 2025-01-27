import pandas as pd

# Load the two Excel files into pandas DataFrames
file1 = '/Users/lanting/Project/dmmr/data/申请数据病理明细 new.xlsx'  # Replace with your actual file path
file2 = '/Users/lanting/Project/dmmr/data/msih_dmmr_crc240807final.xlsx'  # Replace with your actual file path

# Read the Excel files
df1 = pd.read_excel(file1)  # Read the first file
df2 = pd.read_excel(file2)  # Read the second file

# Merge the DataFrames based on the common columns '姓名' in df1 and 'Name' in df2
merged_df = pd.merge(df2, df1[['姓名', '分子编号']], left_on='Name', right_on='姓名', how='right')

# Drop the redundant '姓名' column if you no longer need it after the merge
merged_df = merged_df.drop(columns=['姓名'])


# Reorder the columns to put '分子编号' as the first column
merged_df = merged_df[['分子编号'] + [col for col in merged_df.columns if col != '分子编号']]

# Remove rows with NaN values in '分子编号' (or other columns as needed)
merged_df = merged_df.dropna(subset=['分子编号'])

# Reset the index if necessary
merged_df = merged_df.reset_index(drop=True)

# Save the merged DataFrame to a new Excel file (optional)
merged_df.to_excel('/Users/lanting/Project/dmmr/data/merged_output.xlsx', index=False)

# Display the first few rows to verify the merge
print(merged_df.head())
