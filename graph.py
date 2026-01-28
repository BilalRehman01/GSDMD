import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# Load the dataset
df = pd.read_csv('Qsar_affinity.csv')

# Rename and split the combined first column
df.rename(columns={'SMILES,cats_euclidean_dist': 'temp_col'}, inplace=True)
df[['SMILES', 'cats_euclidean_dist']] = df['temp_col'].str.rsplit(',', n=1, expand=True)

# Convert columns to numeric types, coercing errors
df['cats_euclidean_dist'] = pd.to_numeric(df['cats_euclidean_dist'], errors='coerce')
df['Affinity'] = pd.to_numeric(df['Affinity'], errors='coerce')

# Drop rows with missing values that may have resulted from coercion
df.dropna(subset=['cats_euclidean_dist', 'Affinity'], inplace=True)

# --- Plotting ---
# CORRECTED: Use a more standard and compatible style name
plt.style.use('seaborn-whitegrid')
fig, ax = plt.subplots(figsize=(8, 6))

# Create a regression plot
sns.regplot(
    x='cats_euclidean_dist',
    y='Affinity',
    data=df,
    ax=ax,
    scatter_kws={'alpha': 0.6, 's': 50, 'edgecolor': 'w', 'linewidth': 0.7},
    line_kws={'color': 'red', 'linewidth': 2}
)

# --- Calculate R-squared ---
# Extracting the data for linear regression
x = df['cats_euclidean_dist']
y = df['Affinity']

# Perform linear regression to get the R-squared value
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
r_squared = r_value**2

# --- Final Touches ---
# Add R-squared value to the plot
ax.text(0.05, 0.95, f'$R^2 = {r_squared:.2f}$', transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

# Set titles and labels
ax.set_title('Binding Affinity vs. QSAR Score', fontsize=16, fontweight='bold')
ax.set_xlabel('QSAR Score (cats_euclidean_dist)', fontsize=12)
ax.set_ylabel('Binding Affinity', fontsize=12)

# Display the plot
plt.tight_layout()
plt.savefig('corrected_qsar_plot.png')
