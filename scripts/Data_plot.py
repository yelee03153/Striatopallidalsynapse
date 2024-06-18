import pandas as pd 
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import numpy as np

# Function to add statistical annotations
def add_stat_annotation(ax, x_positions, y_max, p_value):
    x1, x2 = x_positions
    y, h, col = y_max, 0.1 * y_max, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    ax.text((x1 + x2) * 0.5, y + h, f"{'*' if p_value < 0.05 else 'ns'}", ha='center', va='bottom', color=col)

# Create an empty dictionary to hold the dataframes
sheets_data = {}
excel_data=pd.ExcelFile("data/6-OHDA PPR trace_ Additional Analysis.xlsx")

# Define the measurements based on the provided image
measurements = [
    "PPR", "Peak amplitude", "Time of peak", "Area(pA*ms)", "Half-width (ms)",
    " Time of maximum rise slope (ms) 500ms=Stimulation apply time", 
    " Time of maximum decay slope (ms) 500ms=Stimulation apply time",
    "Rise time", "Rise slope", "Decay time", "Decay slope"
]

# Create an empty dictionary to hold the dataframes
sheets_data = {}
sheet_names=["DL","VL","DM","VM"]
# Process each sheet in the Excel file
for sheet_name in sheet_names:
    # Read the sheet into a dataframe
    df = excel_data.parse(sheet_name)
    
    # Filter rows where the second column indicates "Base", "QP", and "After"
    df_base = df[df.iloc[:, 1] == "Base"]
    df_during = df[df.iloc[:, 1] == "QP"]
    df_washout = df[df.iloc[:, 1] == "After"]
    
    # Create a multi-index for the dataframe
    tuples = [(sheet_name, measurement) for measurement in measurements]
    index = pd.MultiIndex.from_tuples(tuples, names=['Sheet', 'Measurement'])
    
    # Create the multi-index dataframe with the desired columns
    columns = pd.MultiIndex.from_product([measurements, ['Baseline', 'During', 'Wash-out']], names=['Measurement', 'Phase'])
    multi_index_df = pd.DataFrame(index=range(len(df_base)), columns=columns)
    
    # Fill the multi-index dataframe with the data
    for i, measurement in enumerate(measurements):
        # Extracting the data for each measurement
        base_data = df_base[measurement].values
        during_data = df_during[measurement].values
        washout_data = df_washout[measurement].values

        # Check if the data has the correct shape
            # Fill the dataframe with the extracted data
        multi_index_df[(measurement, 'Baseline')] = base_data
        multi_index_df[(measurement, 'During')] =during_data
        multi_index_df[(measurement, 'Wash-out')] = washout_data

    sheets_data[sheet_name] = multi_index_df

# Concatenate all the dataframes into a single dataframe
final_multi_index_df = pd.concat(sheets_data.values(), keys=sheets_data.keys())


# Perform the statistical tests and store results
stat_results = {}

for measurement in measurements:
    data = final_multi_index_df[measurement]
    
    # Perform statistical tests
    stat_results[measurement] = {}
    for condition1, condition2 in [('Baseline', 'During'), ('During', 'Wash-out'), ('Baseline', 'Wash-out')]:
        group1 = data[condition1].dropna()
        group2 = data[condition2].dropna()
        
        if len(group1) > 0 and len(group2) > 0:
            stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
            stat_results[measurement][(condition1, condition2)] = p_value

# Now plot the data with statistical annotations
num_plots = len(measurements)
ncols = 2
nrows = (num_plots + ncols - 1) // ncols  # Calculate rows needed for given columns

# Now plot the data with statistical annotations
fig, axes = plt.subplots(nrows, ncols, figsize=(15, 5 * nrows), sharex=True)

# Flatten axes array for easy iteration
axes = axes.flatten()


for ax, measurement in zip(axes, measurements):
    data = final_multi_index_df[measurement].dropna()
    
    # Plotting each paired sample for Baseline, During, and Wash-out with grey lines
    for idx in data.index:
        ax.plot(['Baseline', 'During', 'Wash-out'], data.loc[idx], color='lightgrey', marker='o')
    
    # Plotting dots for each phase
    ax.scatter(['Baseline'] * len(data), data['Baseline'], color='blue', label='Baseline')
    ax.scatter(['During'] * len(data), data['During'], color='orange', label='During')
    ax.scatter(['Wash-out'] * len(data), data['Wash-out'], color='red', label='Wash-out')
    
    ax.set_title(measurement)
    ax.set_ylabel('Value')
    ax.set_xlabel('Phase')
    
    # Calculate y_max for annotation placement
    y_max = max(np.nanmax(pd.to_numeric(data['Baseline'], errors='coerce')),
                np.nanmax(pd.to_numeric(data['During'], errors='coerce')),
                np.nanmax(pd.to_numeric(data['Wash-out'], errors='coerce')))
    
    # Add statistical annotations if significant
    for (condition1, condition2), p_value in stat_results[measurement].items():
        if p_value < 0.05:
            x_positions = [0, 1] if condition1 == 'Baseline' and condition2 == 'During' else \
                          [1, 2] if condition1 == 'During' and condition2 == 'Wash-out' else [0, 2]
            add_stat_annotation(ax, x_positions, y_max, p_value)

plt.tight_layout()
plt.savefig("results/plots/6OHDA_quantif.pdf")
