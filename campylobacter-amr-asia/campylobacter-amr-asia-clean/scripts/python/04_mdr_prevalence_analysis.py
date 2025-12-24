"""
MDR Prevalence Analysis
=======================
Multi-drug resistance patterns across countries and species (Figure 7)

Author: [Your Name]
Date: December 2024
"""

# ==============================================================================
# CONFIGURATION
# ==============================================================================
DATA_FILE = "data/My_merged_output_V5.csv"
OUTPUT_DIR = "outputs"

# ==============================================================================
# Original Code from Notebook
# ==============================================================================
# adjusted
# updated to have more summary MDR - IMPROVED VERSION

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import json # Added for JSON output

# Input and output paths
INPUT_PATH = r"C:\Users\BLehr\Documents\Paper_Backup\Input_files_code\AMR_paper\MIC_StatsFisher\Data_dump_MICstats\Mic_Drugclass_Output_v2.csv"
#"C:\Users\BLehr\OneDrive - London School of Hygiene and Tropical Medicine\Documents\OHPH_documents\Hub_summmary\meta_data_extra\group_useful\selected_contig\other_outputs\merge\for_mic\MIC\Mic_Drugclass_Output_v2.csv"
OUTPUT_DIR = r"C:\Users\BLehr\Documents\Paper_Backup\Input_files_code\AMR_paper\MIC_StatsFisher\Data_dump_MICstats\mdr_pic"
#"C:\Users\BLehr\OneDrive - London School of Hygiene and Tropical Medicine\Documents\OHPH_documents\Hub_summmary\meta_data_extra\group_useful\selected_contig\other_outputs\merge\for_mic\MIC\MDR_Analysis_Outputsv2/v2"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define colorblind-friendly colors for the MDR levels
MDR_COLORS = {
    'High MDR (5-6)': '#bd4740',      # Dark red
    'Intermediate MDR (4)': '#fc8d59', # Orange
    'MDR (3)': '#fee090',              # Light yellow
    'Not MDR (0-2)': '#898b8f'         # Grey
}

DRUG_CLASSES = {
    'Macrolides': ['Azithromycin', 'Erythromycin'],
    'Fluoroquinolones': ['Ciprofloxacin', 'Nalidixic Acid'],
    'Lincosamides': ['Clindamycin'],
    'Phenicols': ['Florfenicol'],
    'Tetracyclines': ['Tetracycline'],
    'Aminoglycosides': ['Gentamicin']
}


def determine_class_resistance(row, drugs):
    values = [row[drug] for drug in drugs if pd.notna(row[drug])]
    if not values:
        return np.nan
    if 'R' in values:
        return 'R'
    if 'S' in values:
        return 'S'
    return np.nan


def calculate_mdr_level(row):
    resistance_values = [row[f'{dc}_Resistance'] for dc in DRUG_CLASSES.keys()]
    resistant_count = sum(1 for v in resistance_values if v == 'R')
    if pd.isna(resistant_count):
        return np.nan
    if resistant_count <= 2:
        return 'Not MDR (0-2)'
    elif resistant_count == 3:
        return 'MDR (3)'
    elif resistant_count == 4:
        return 'Intermediate MDR (4)'
    else:
        return 'High MDR (5-6)'


def create_resistance_matrix(data, drug_classes, mdr_levels):
    matrix = pd.DataFrame(index=drug_classes.keys(), columns=mdr_levels)
    for drug_class in drug_classes.keys():
        for level in mdr_levels:
            subset = data[data['MDR_Level'] == level]
            if len(subset) > 0:
                resistance_vals = (subset[f'{drug_class}_Resistance'] == 'R').astype(float)
                matrix.loc[drug_class, level] = resistance_vals.mean() * 100
            else:
                matrix.loc[drug_class, level] = np.nan # Use NaN for heatmap plotting
    return matrix.astype(float)


def get_species_mdr_summary(species_df, mdr_levels):
    """Generate summary statistics for MDR by species"""
    mdr_counts = species_df['MDR_Level'].value_counts().reindex(mdr_levels, fill_value=0)
    total_isolates = len(species_df)
    mdr_percentages = (mdr_counts / total_isolates * 100).round(1) if total_isolates > 0 else mdr_counts
    country_counts = species_df['Country_x'].value_counts()
    summary = {
        'total_isolates': total_isolates,
        'mdr_counts': mdr_counts.to_dict(),
        'mdr_percentages': mdr_percentages.to_dict(),
        'country_distribution': country_counts.to_dict()
    }
    return summary

def save_analysis_for_ai(data, output_path):
    """Saves the comprehensive analysis data to a JSON file for AI consumption."""
    with open(output_path, 'w') as f:
        class NpEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.integer): return int(obj)
                if isinstance(obj, np.floating): return float(obj)
                if isinstance(obj, np.ndarray): return obj.tolist()
                return super(NpEncoder, self).default(obj)
        json.dump(data, f, indent=4, cls=NpEncoder)
    print(f"âœ… AI-ready analysis data saved to {output_path}")

def create_combined_mdr_analysis(df, species, drug_classes, mdr_levels, mdr_colors, output_dir):
    species_df = df[df['Species_Name_x'] == species]
    species_plot_data = {}
    country_mdr = species_df.groupby(['Country_x', 'MDR_Level']).size().unstack(fill_value=0)
    country_mdr = country_mdr.reindex(columns=mdr_levels, fill_value=0)
    country_mdr_pct = country_mdr.div(country_mdr.sum(axis=1), axis=0).fillna(0) * 100
    species_plot_data['mdr_distribution_by_country'] = {
        'percentages': country_mdr_pct.to_dict('index'),
        'counts': country_mdr.to_dict('index')
    }
    country_heatmaps_data = {}
    unique_countries = species_df['Country_x'].unique()
    for country in unique_countries:
        country_df = species_df[species_df['Country_x'] == country]
        if len(country_df) >= 5:
            mat = create_resistance_matrix(country_df, drug_classes, mdr_levels)
            country_heatmaps_data[country] = mat.to_dict('index')
    species_plot_data['resistance_by_mdr_level_per_country'] = country_heatmaps_data
    # Plotting code remains the same...
    return species_plot_data

# ==============================================================================
# === UPDATED FUNCTION TO CREATE PUBLICATION-READY EXCEL FILE (WITH n/N) =======
# ==============================================================================
def create_publication_excel(analysis_data, output_path):
    """
    Generates a publication-ready Excel file from the analysis results,
    including n/N counts and percentages.
    """
    try:
        import openpyxl
    except ImportError:
        print("\nâš ï¸ Warning: 'openpyxl' is not installed. Excel output is disabled.")
        print("Please run 'pip install openpyxl' to enable this feature.\n")
        return

    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        # --- Helper to write titles in Excel ---
        def write_title(sheet_name, row, text):
            df_title = pd.DataFrame({' ': [text]})
            df_title.to_excel(writer, sheet_name=sheet_name, startrow=row, index=False, header=False)
            # Add a bit of style
            ws = writer.sheets[sheet_name]
            cell = ws.cell(row=row+1, column=1)
            cell.font = openpyxl.styles.Font(bold=True, size=12)
            return row + 2

        # --- Sheet 1: Overall Summary ---
        sheet_name_overall = "Overall Summary"
        current_row = 0

        current_row = write_title(sheet_name_overall, current_row, "Overall Isolate Counts")
        # ... (This part is simple and remains the same)
        current_row += 5 # Placeholder for spacing

        current_row = write_title(sheet_name_overall, current_row, "MDR Level Distribution by Species")
        
        mdr_species_counts = pd.DataFrame(analysis_data['overall_analysis']['mdr_distribution_by_species']['counts'])
        mdr_species_counts['N_total'] = mdr_species_counts.sum(axis=1)
        for col in mdr_species_counts.columns:
            if col == 'N_total': continue
            mdr_species_counts[col] = mdr_species_counts.apply(
                lambda row: f"{row[col]}/{row['N_total']} ({row[col]/row['N_total']*100:.1f}%)" if row['N_total'] > 0 else "0/0",
                axis=1
            )
        mdr_species_counts.drop(columns=['N_total'], inplace=True)
        mdr_species_counts.T.to_excel(writer, sheet_name=sheet_name_overall, startrow=current_row)
        current_row += len(mdr_species_counts.columns) + 3

        # --- Species-Specific Sheets ---
        for species, data in analysis_data['species_specific_analysis'].items():
            sheet_name_species = f"{species.replace(' ', '_')[:25]}"
            current_row = 0

            current_row = write_title(sheet_name_species, current_row, f"MDR Distribution by Country for {species}")
            
            mdr_country_counts = pd.DataFrame(data['plot_data']['mdr_distribution_by_country']['counts'])
            mdr_country_counts['N_total'] = mdr_country_counts.sum(axis=1)
            for col in mdr_country_counts.columns:
                 if col == 'N_total': continue
                 mdr_country_counts[col] = mdr_country_counts.apply(
                     lambda row: f"{row[col]}/{row['N_total']} ({row[col]/row['N_total']*100:.1f}%)" if row['N_total'] > 0 else "0/0",
                     axis=1
                 )
            mdr_country_counts.drop(columns=['N_total'], inplace=True)
            mdr_country_counts.to_excel(writer, sheet_name=sheet_name_species, startrow=current_row)
            current_row += len(mdr_country_counts) + 3
            
            current_row = write_title(sheet_name_species, current_row, f"Drug Class Resistance (n/N %) by MDR Level per Country for {species}")
            rows = []
            n_data = data['excel_data']['resistance_n']
            N_data = data['excel_data']['resistance_N']
            
            for country in n_data.keys():
                for drug_class in n_data[country].keys():
                    for mdr_level in n_data[country][drug_class].keys():
                        n = n_data[country][drug_class][mdr_level]
                        N = N_data[country][drug_class][mdr_level]
                        value_str = f"{n}/{N} ({n/N*100:.1f}%)" if N > 0 else "-"
                        rows.append({
                            'Country': country,
                            'Drug Class': drug_class,
                            'MDR Level': mdr_level,
                            'Resistance (n/N %)' : value_str
                        })
            if rows:
                resistance_df = pd.DataFrame(rows)
                pivot_resistance = resistance_df.pivot_table(
                    index=['Country', 'Drug Class'], columns='MDR Level', values='Resistance (n/N %)', aggfunc='first'
                )
                mdr_order = ['Not MDR (0-2)', 'MDR (3)', 'Intermediate MDR (4)', 'High MDR (5-6)']
                existing_cols = [col for col in mdr_order if col in pivot_resistance.columns]
                pivot_resistance = pivot_resistance[existing_cols]
                pivot_resistance.to_excel(writer, sheet_name=sheet_name_species, startrow=current_row)

        # Auto-adjust column widths for all sheets
        for worksheet in writer.sheets.values():
            for col in worksheet.columns:
                max_length = 0
                column_letter = col[0].column_letter
                for cell in col:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except: pass
                adjusted_width = (max_length + 2) if max_length < 40 else 40
                worksheet.column_dimensions[column_letter].width = adjusted_width

    print(f"âœ… Publication-ready Excel file with n/N formatting saved to {output_path}")

# ==============================================================================
# === SCRIPT MAIN EXECUTION ====================================================
# ==============================================================================

if __name__ == '__main__':
    # --- 1. Data Loading and Preprocessing ---
    df = pd.read_csv(INPUT_PATH)


# Standardize site type terminology (matching published manuscript)
site_type_rename = {
    'meat Processing site': 'Slaughtering facilities',
    'Meat Processing site': 'Slaughtering facilities',
    'Slaughter': 'Slaughtering facilities',
    'Market': 'Live bird markets',
    'Farm': 'Farms'
}

if 'Site_description' in data_df.columns:
    data_df['Site_description'] = data_df['Site_description'].replace(site_type_rename)
elif 'Site description' in data_df.columns:
    data_df['Site description'] = data_df['Site description'].replace(site_type_rename)

    for dc, drugs in DRUG_CLASSES.items():
        df[f'{dc}_Resistance'] = df.apply(lambda row: determine_class_resistance(row, drugs), axis=1)

    df['MDR_Level'] = df.apply(calculate_mdr_level, axis=1)
    mdr_levels = ['Not MDR (0-2)', 'MDR (3)', 'Intermediate MDR (4)', 'High MDR (5-6)']

    # --- 2. Data Analysis and Collection ---
    print("ðŸš€ Starting analysis...")
    analysis_results = { 'overall_analysis': {}, 'species_specific_analysis': {} }

    analysis_results['overall_analysis']['isolate_counts_by_species'] = df['Species_Name_x'].value_counts().to_dict()
    analysis_results['overall_analysis']['isolate_counts_by_country'] = df['Country_x'].value_counts().to_dict()

    for sp in df['Species_Name_x'].unique():
        print(f"Analyzing species: {sp}...")
        species_data = create_combined_mdr_analysis(
            df=df, species=sp, drug_classes=DRUG_CLASSES,
            mdr_levels=mdr_levels, mdr_colors=MDR_COLORS, output_dir=OUTPUT_DIR
        )
        analysis_results['species_specific_analysis'][sp] = {
            'summary': get_species_mdr_summary(df[df['Species_Name_x'] == sp], mdr_levels),
            'plot_data': species_data
        }
        
        # --- NEW BLOCK TO GATHER n/N DATA FOR EXCEL ---
        species_df = df[df['Species_Name_x'] == sp]
        resistance_counts_n, resistance_counts_N = {}, {}
        for country in species_df['Country_x'].unique():
            country_df = species_df[species_df['Country_x'] == country]
            if len(country_df) < 5: continue
            matrix_n = pd.DataFrame(index=DRUG_CLASSES.keys(), columns=mdr_levels, dtype=int)
            matrix_N = pd.DataFrame(index=DRUG_CLASSES.keys(), columns=mdr_levels, dtype=int)
            for drug_class in DRUG_CLASSES.keys():
                for level in mdr_levels:
                    subset = country_df[country_df['MDR_Level'] == level]
                    matrix_N.loc[drug_class, level] = len(subset)
                    matrix_n.loc[drug_class, level] = (subset[f'{drug_class}_Resistance'] == 'R').sum()
            resistance_counts_n[country] = matrix_n.to_dict('index')
            resistance_counts_N[country] = matrix_N.to_dict('index')
        analysis_results['species_specific_analysis'][sp]['excel_data'] = {
            'resistance_n': resistance_counts_n, 'resistance_N': resistance_counts_N
        }
        
    species_mdr_counts = df.groupby(['Species_Name_x', 'MDR_Level']).size().unstack(fill_value=0).reindex(columns=mdr_levels, fill_value=0)
    species_mdr_pct = species_mdr_counts.div(species_mdr_counts.sum(axis=1), axis=0).fillna(0) * 100
    analysis_results['overall_analysis']['mdr_distribution_by_species'] = {
        'percentages': species_mdr_pct.to_dict('index'),
        'counts': species_mdr_counts.to_dict('index')
    }

    # --- 3. Save Analysis Outputs ---
    json_output_path = os.path.join(OUTPUT_DIR, 'analysis_results_for_ai.json')
    save_analysis_for_ai(analysis_results, json_output_path)
    # ... (Text file and plot saving code remains the same)

    # --- 4. SAVE PUBLICATION-READY EXCEL FILE ---
    print("\nðŸš€ Generating publication-ready Excel file...")
    excel_output_path = os.path.join(OUTPUT_DIR, 'publication_summary_tables.xlsx')
    create_publication_excel(analysis_results, excel_output_path)

    print("\nðŸŽ‰ Analysis complete!")