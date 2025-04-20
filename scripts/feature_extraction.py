# -----------------------------------------------------------------------------
# Feature Extraction Script
# -----------------------------------------------------------------------------
# Description:
# Extracts electrophysiological features from Excel data across experimental
# conditions and striatal regions, organizing them into a MultiIndex DataFrame and
# saving to a pickle file for downstream modeling and analysis.

import pandas as pd
import numpy as np

# Define input Excel files for each experimental condition
file_paths = {
    "10uM QP": "data/10uM QP PPR trace_ Additional Analysis.xlsx",
    "6-OHDA": "data/6-OHDA PPR trace_ Additional Analysis.xlsx",
    "10uM CdCl2": "data/10uM CdCl2 PPR trace_ Additional Analysis.xlsx"
}

# Measurements to extract
measurements = [
    "PPR", "Peak amplitude", "Time of peak", "Area(pA*ms)", "Half-width (ms)",
    " Time of maximum rise slope (ms) 500ms=Stimulation apply time", 
    " Time of maximum decay slope (ms) 500ms=Stimulation apply time",
    "Rise time", "Rise slope", "Decay time", "Decay slope"
]

# Brain regions (sheet names)
sheet_names = ["DL", "VL", "DM", "VM"]

# Where the extracted features will be stored
combined_data = {}

# Loop through all Excel files and extract features from each
for file_label, file_path in file_paths.items():
    excel_data = pd.ExcelFile(file_path)
    sheets_data = {}

    for sheet_name in sheet_names:
        df = excel_data.parse(sheet_name)

        # Subsets based on experimental phase labels
        df_base = df[df.iloc[:, 1] == "Base"]
        df_during = df[df.iloc[:, 1].isin(["QP", "During"])]
        df_washout = df[df.iloc[:, 1] == "After"]

        # Prepare the multi-index DataFrame for features
        columns = pd.MultiIndex.from_product([measurements, ['Baseline', 'During', 'Wash-out']], names=['Measurement', 'Phase'])
        multi_index_df = pd.DataFrame(index=range(len(df_base)), columns=columns)

        for measurement in measurements:
            base_data = df_base[measurement].values if measurement in df_base.columns else np.full(len(df_base), np.nan)

            if not df_during.empty and measurement in df_during.columns:
                during_data = df_during[measurement].values
            else:
                during_data = np.full(len(df_base), np.nan)

            if not df_washout.empty and measurement in df_washout.columns:
                washout_data = df_washout[measurement].values
            else:
                washout_data = np.full(len(df_base), np.nan)

            multi_index_df[(measurement, 'Baseline')] = base_data
            multi_index_df[(measurement, 'During')] = during_data
            multi_index_df[(measurement, 'Wash-out')] = washout_data

        sheets_data[sheet_name] = multi_index_df

    # Combine all sheets into one DataFrame per file
    file_multi_index_df = pd.concat(sheets_data.values(), keys=sheets_data.keys(), names=['Sheet', 'Index'])
    combined_data[file_label] = file_multi_index_df

# Combine data from all files into a final MultiIndex DataFrame
final_combined_df = pd.concat(combined_data.values(), keys=combined_data.keys(), names=['File', 'Sheet', 'Index'])

# Save as pickle file for reuse
final_combined_df.to_pickle("results/features/features_all_data.pkl")
print("✅ Feature data saved to: results/features/features_all_data.pkl")

