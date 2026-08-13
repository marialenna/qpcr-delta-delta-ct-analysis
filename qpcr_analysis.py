import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Load raw data (handling European/Excel semicolon separators)
raw_file = "hospital.csv"
df = pd.read_csv(raw_file, sep=';')

# Clean column names (strip spaces, handle potential naming variations)
df.columns = df.columns.str.strip()

# Print Raw Data
print("--- Raw Hospital Data ---")
print(df)
print("\n")

# Detect correct column names dynamically
target_col = [c for c in df.columns if 'Target' in c][0]
gapdh_col = [c for c in df.columns if 'GAPDH' in c or '18' in str(df[c].iloc[0]) or 'Ct' in c][-1]

# Convert strings to numbers in case Excel saved them with commas (e.g. 22,4)
df[target_col] = pd.to_numeric(df[target_col].astype(str).str.replace(',', '.'), errors='coerce')
df[gapdh_col] = pd.to_numeric(df[gapdh_col].astype(str).str.replace(',', '.'), errors='coerce')

# 2. Clean data (remove missing values and Ct > 35)
df_clean = df.dropna(subset=[target_col, gapdh_col]).copy()
df_clean = df_clean[(df_clean[target_col] <= 35.0) & (df_clean[gapdh_col] <= 35.0)]

print("--- Cleaned Data ---")
print(df_clean)
print("\n")

# 3. Calculate 2^(-Delta Delta Ct)
df_clean['Delta_Ct'] = df_clean[target_col] - df_clean[gapdh_col]
control_mean = df_clean[df_clean['Group'] == 'Control']['Delta_Ct'].mean()
df_clean['Delta_Delta_Ct'] = df_clean['Delta_Ct'] - control_mean
df_clean['Fold_Change'] = 2 ** (-df_clean['Delta_Delta_Ct'])

# 4. Plot
summary = df_clean.groupby('Group')['Fold_Change'].agg(['mean', 'std']).reset_index()
plt.figure(figsize=(6, 5))
plt.bar(summary['Group'], summary['mean'], yerr=summary['std'], capsize=5, color=['#a8a8a8', '#2b5c8f'], edgecolor='black')
plt.title('Relative Gene Expression (Hospital Data)', fontsize=12, fontweight='bold')
plt.ylabel('Fold Change (Normalized to GAPDH)', fontsize=10)
plt.grid(axis='y', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig("qpcr_fold_change_plot.png", dpi=300)
print("SUCCESS: Plot saved as 'qpcr_fold_change_plot.png'")