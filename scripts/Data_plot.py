import pandas as pd 
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import numpy as np
import seaborn as sns
# Function to add statistical annotations
def add_stat_annotation(ax, x_positions, y_max, p_value):
    x1, x2 = x_positions
    y, h, col = y_max, 0.1 * y_max, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    ax.text((x1 + x2) * 0.5, y + h, f"{'*' if p_value < 0.05 else 'ns'}", ha='center', va='bottom', color=col)

def add_stat_annotation_corss(ax, x_positions, y_max, p_value, comparison_label=""):
    x1, x2 = x_positions
    y, h, col = y_max, 0.1 * y_max, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    label = f"{comparison_label}\n{'*' if p_value < 0.05 else 'ns'}"
    ax.text((x1 + x2) * 0.5, y + h, label, ha='center', va='bottom', color=col)


# Create an empty dictionary to hold the dataframes
sheets_data = {}
excel_data=pd.ExcelFile("data/10uM QP PPR trace_ Additional Analysis.xlsx")
file_paths = {
    "10uM QP": "data/10uM QP PPR trace_ Additional Analysis.xlsx",
    "6-OHDA": "data/6-OHDA PPR trace_ Additional Analysis.xlsx",
    "10uM CdCl2": "data/10uM CdCl2 PPR trace_ Additional Analysis.xlsx"
}

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
plt.savefig("results/plots/QUIN_quantif.pdf")

combined_data = {}

# Loop through each file and read the data
for file_label, file_path in file_paths.items():
    excel_data = pd.ExcelFile(file_path)
    sheets_data = {}
    print(file_path)
    # Process each sheet in the Excel file
    for sheet_name in sheet_names:
        # Read the sheet into a dataframe
        df = excel_data.parse(sheet_name)
        
        # Filter rows where the second column indicates "Base", "QP", and "After"
        df_base = df[df.iloc[:, 1] == "Base"]
        df_during = df[df.iloc[:, 1].isin(["QP", "During"])] if any(df.iloc[:, 1].isin(["QP", "During"])) else pd.DataFrame()
        df_washout = df[df.iloc[:, 1] == "After"] if "After" in df.iloc[:, 1].values else pd.DataFrame()
        
        # Create a multi-index dataframe with the desired columns
        columns = pd.MultiIndex.from_product([measurements, ['Baseline', 'During', 'Wash-out']], names=['Measurement', 'Phase'])
        multi_index_df = pd.DataFrame(index=range(len(df_base)), columns=columns)
        
        # Fill the multi-index dataframe with the data
        for measurement in measurements:
            # Extracting the data for each measurement
            base_data = df_base[measurement].values
            during_data = df_during[measurement].values if measurement in df_during.columns else np.nan
            washout_data = df_washout[measurement].values if measurement in df_washout.columns else np.nan

            # Fill the dataframe with the extracted data
            multi_index_df[(measurement, 'Baseline')] = base_data
            if during_data is not np.nan:
                multi_index_df[(measurement, 'During')] =during_data
            if washout_data is not np.nan:
                multi_index_df[(measurement, 'Wash-out')] = washout_data

        sheets_data[sheet_name] = multi_index_df

    # Concatenate all the sheets into one dataframe
    file_multi_index_df = pd.concat(sheets_data.values(), keys=sheets_data.keys(), names=['Sheet', 'Index'])
    
    # Store the dataframe with a label for the file
    combined_data[file_label] = file_multi_index_df

# Concatenate all dataframes from each file into a single dataframe
final_combined_df = pd.concat(combined_data.values(), keys=combined_data.keys(), names=['File', 'Sheet', 'Index'])

# Perform statistical tests and store results
stat_results = {}

# Loop through each measurement and compare across conditions and files
for measurement in measurements:
    data = final_combined_df[measurement]
    stat_results[measurement] = {}

    # Perform intra-file statistical tests for conditions
    for file_label in file_paths.keys():
        file_data = data.loc[file_label]
        for condition1, condition2 in [('Baseline', 'During'), ('During', 'Wash-out'), ('Baseline', 'Wash-out')]:
            group1 = file_data[condition1].dropna()
            group2 = file_data[condition2].dropna()

            if len(group1) > 0 and len(group2) > 0:
                stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
                stat_results[measurement][(file_label, condition1, condition2)] = p_value

    # Perform inter-file statistical tests for each condition
    for condition in ['Baseline', 'During', 'Wash-out']:
        groups = [data.loc[file_label, condition].dropna() for file_label in file_paths.keys()]
        for i, group1 in enumerate(groups):
            for j, group2 in enumerate(groups):
                if i < j and len(group1) > 0 and len(group2) > 0:
                    file_label1 = list(file_paths.keys())[i]
                    file_label2 = list(file_paths.keys())[j]
                    stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
                    stat_results[measurement][(condition, file_label1, file_label2)] = p_value

# Now plot the data with statistical annotations
num_plots = len(measurements)
ncols = 2
nrows = (num_plots + ncols - 1) // ncols  # Calculate rows needed for given columns

fig, axes = plt.subplots(nrows, ncols, figsize=(15, 5 * nrows), sharex=True)
axes = axes.flatten()

# Create colors for each file
colors = {'10uM QP': 'blue', '6-OHDA': 'green', '10uM CdCl2': 'red'}



for measurement in measurements:
    # Create a new figure with 12 subplots arranged in a 3x4 grid
    fig, axes = plt.subplots(4, 3, figsize=(20, 12))  # 3 rows and 4 columns grid
    fig.suptitle(f'Comparison of {measurement} Across Conditions and Regions', fontsize=20)
    global_min = float('inf')
    global_max = float('-inf')

    for file_label in file_paths.keys():
        for region in sheet_names:
            for phase in ['Baseline', 'During', 'Wash-out']:
                    # Check if the measurement and phase exist in the DataFrame
                    if measurement in final_combined_df.columns.get_level_values('Measurement') and phase in final_combined_df.columns.get_level_values('Phase'):
                        # Access data using .loc to specify levels of interest in MultiIndex
                        phase_values =  final_combined_df[measurement].loc[file_label].loc[region][phase].dropna()
                        if not phase_values.empty:
                            global_min = min(global_min, phase_values.min())
                            global_max = max(global_max, phase_values.max())
    # Flatten axes for easier iteration
    axes = axes.flatten()

    # Define the phases and corresponding index positions in the grid
    phases = ['Baseline', 'During', 'Wash-out']

    # Create a dictionary to map each region and phase to a subplot index
    subplot_index = {}
    for row, region in enumerate(sheet_names):
        for col, phase in enumerate(phases):
            subplot_index[(region, phase)] = row * len(phases) + col

    # Iterate through each region and phase to populate subplots
    for region in sheet_names:
        for phase in phases:
            ax_index = subplot_index[(region, phase)]
            ax = axes[ax_index]

            # Prepare data for swarm plot for each region and phase
            phase_data = []
            for file_label in file_paths.keys():
                try:
                    # Check if the measurement and phase exist in the DataFrame
                    if measurement in final_combined_df.columns.get_level_values('Measurement') and phase in final_combined_df.columns.get_level_values('Phase'):
                        # Access data using .loc to specify levels of interest in MultiIndex
                        phase_values =  final_combined_df[measurement].loc[file_label].loc[region][phase].dropna()
                        if not phase_values.empty:
                            phase_data.append(pd.DataFrame({
                                'Value': phase_values,
                                'File': file_label,
                            }))
                except KeyError:
                    print(f"No data available for file: {file_label}, region: {region}, phase: {phase}")

            # If no data exists for this phase, turn off the subplot
            if len(phase_data) == 0 or all([df.empty for df in phase_data]):
                ax.axis('off')  # Turn off this subplot if no data
                continue

            # Combine data from all files for this phase and region
            phase_data_df = pd.concat(phase_data)
            sns.violinplot(x='File', y='Value', data=phase_data_df, ax=ax,  palette='pastel')

            # Create swarm plot for the given phase and region
            sns.swarmplot(x='File', y='Value', data=phase_data_df, ax=ax, color='k', alpha=0.7)

            ax.set_ylim(global_min - 0.1 * abs(global_max - global_min), global_max + 0.1 * abs(global_max - global_min))  # Add padding to min and max

            # Set subplot titles and labels
            ax.set_title(f'{region} - {phase}')
            ax.set_ylabel(measurement)

    #plt.show()

    # Adjust layout and save the figure
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust for suptitle
    plt.savefig(f"results/plots/{measurement}_Region_Phase_Comparison.pdf")
    plt.close()
    final_combined_df.to_pickle("results/features_all_data.pkl")