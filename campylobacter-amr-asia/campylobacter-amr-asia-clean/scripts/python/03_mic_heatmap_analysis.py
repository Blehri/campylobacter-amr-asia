"""
MIC Heatmap Analysis
====================
Interactive heatmaps of phenotypic resistance with ECOFF comparisons (Figure 6)

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
#MIC HEATMAP WITH REFORMATED EXCEL OUTPUT
# =============================================================================
# 0) Imports
# =============================================================================
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint          # NEW (CI)
import itertools
from pathlib import Path

# =============================================================================
# 1) Define Antibiotic Classes
# =============================================================================
ANTIBIOTIC_CLASSES = {
    'Aminoglycosides': ['Gentamicin'],
    'Fluoroquinolones': ['Ciprofloxacin', 'Nalidixic Acid'],
    'Lincosamides': ['Clindamycin'],
    'Macrolides': ['Azithromycin', 'Erythromycin'],
    'Phenicols'     : ['Florfenicol'],
    'Tetracyclines' : ['Tetracycline']
}

# =============================================================================
# 2) Calculate Resistance Rates, CIs, and Perform Statistical Analysis
# =============================================================================
def process_and_analyze_resistance(df, grouping_vars=['Country', 'Species_Name']):
    """
    Summarises antibiotic class resistance, adds 95 % Wilson score CIs,
    performs pairwise Fisher tests with BH correction, and
    returns DataFrames for the heat-map and a detailed statistical report.
    """
    df = df.copy()

    # Fill missing grouping variables with "Unknown"
    for col in grouping_vars:
        if col not in df.columns:
            print(f"Warning: Column '{col}' missing in dataset.")
            df[col] = 'Unknown'
        else:
            df[col] = df[col].fillna('Unknown')

    # Create a unique group identifier
    df['group'] = df[grouping_vars].apply(
        lambda row: ' - '.join(row.values.astype(str)),
        axis=1
    )

    # --- Part A: Calculate Resistance Rates (+ 95 % CI) for Each Group ---
    summary_rows = []
    for group_key, group_data in df.groupby(grouping_vars):
        row_dict = dict(zip(grouping_vars,
                            group_key if isinstance(group_key, tuple) else [group_key]))
        row_dict['Total_Samples'] = len(group_data)

        for class_name, antibiotics_list in ANTIBIOTIC_CLASSES.items():
            available_cols = [ab for ab in antibiotics_list if ab in group_data.columns]

            if not available_cols:
                row_dict[f'{class_name}_Rate']          = np.nan
                row_dict[f'{class_name}_N']             = 0
                row_dict[f'{class_name}_Resistant_Count'] = 0
                ### NEW (CI): add blank CI columns
                row_dict[f'{class_name}_CI_Lower']      = np.nan
                row_dict[f'{class_name}_CI_Upper']      = np.nan
                continue

            # Any isolate tested for *any* drug in the class counts as 'tested'
            valid_isolates = group_data[available_cols].notna().any(axis=1)
            n_tested = valid_isolates.sum()

            if n_tested == 0:
                row_dict[f'{class_name}_Rate']          = np.nan
                row_dict[f'{class_name}_N']             = 0
                row_dict[f'{class_name}_Resistant_Count'] = 0
                row_dict[f'{class_name}_CI_Lower']      = np.nan
                row_dict[f'{class_name}_CI_Upper']      = np.nan
            else:
                resistant_mask = (
                    group_data[available_cols].isin(['R']).any(axis=1)
                    & valid_isolates
                )
                n_resistant = resistant_mask.sum()
                rate = (n_resistant / n_tested) * 100.0

                ### NEW (CI): compute Wilson CI and convert to %
                ci_low, ci_high = proportion_confint(
                    count=n_resistant,
                    nobs=n_tested,
                    alpha=0.05,
                    method='wilson'
                )
                ci_low  *= 100
                ci_high *= 100

                row_dict[f'{class_name}_Rate']          = rate
                row_dict[f'{class_name}_N']             = int(n_tested)
                row_dict[f'{class_name}_Resistant_Count'] = int(n_resistant)
                row_dict[f'{class_name}_CI_Lower']      = ci_low
                row_dict[f'{class_name}_CI_Upper']      = ci_high

        summary_rows.append(row_dict)

    resistance_summary_df = pd.DataFrame(summary_rows)

    # Friendly combined label
    if 'Country' in resistance_summary_df.columns and 'Species_Name' in resistance_summary_df.columns:
        resistance_summary_df['Country_Species'] = (
            resistance_summary_df['Country'] + " - " + resistance_summary_df['Species_Name']
        )

    # --- Part B: Pairwise Fisher Tests with BH Correction (S10-style format) ---
    stats_results = []
    temp_df = resistance_summary_df.copy()
    
    # Get unique countries and species
    countries = temp_df['Country'].unique()
    species_list = temp_df['Species_Name'].unique()

    for species in species_list:
        species_data = temp_df[temp_df['Species_Name'] == species]
        
        for class_name in ANTIBIOTIC_CLASSES.keys():
            p_values = []
            comparison_data = []
            
            # Compare all country pairs for this species
            for c1, c2 in itertools.combinations(countries, 2):
                d1_df = species_data[species_data['Country'] == c1]
                d2_df = species_data[species_data['Country'] == c2]
                
                if d1_df.empty or d2_df.empty:
                    continue
                
                d1 = d1_df.iloc[0]
                d2 = d2_df.iloc[0]
                
                resistant1 = int(d1[f'{class_name}_Resistant_Count'])
                total1 = int(d1[f'{class_name}_N'])
                rate1 = d1[f'{class_name}_Rate']
                ci_lower1 = d1[f'{class_name}_CI_Lower']
                ci_upper1 = d1[f'{class_name}_CI_Upper']
                
                resistant2 = int(d2[f'{class_name}_Resistant_Count'])
                total2 = int(d2[f'{class_name}_N'])
                rate2 = d2[f'{class_name}_Rate']
                ci_lower2 = d2[f'{class_name}_CI_Lower']
                ci_upper2 = d2[f'{class_name}_CI_Upper']
                
                if total1 == 0 or total2 == 0:
                    continue
                
                cont_table = [[resistant1, total1 - resistant1],
                              [resistant2, total2 - resistant2]]
                _, p_val = fisher_exact(cont_table)
                
                p_values.append(p_val)
                comparison_data.append({
                    'Comparison_Type': 'Country',
                    'Species': species,
                    'Country_1': c1,
                    'Country_2': c2,
                    'Antibiotic_Class': class_name,
                    'Country_1_n/N': f"{resistant1}/{total1}",
                    'Country_1_Rate': round(rate1, 2),
                    'Country_1_95%CI': f"({round(ci_lower1, 2)}-{round(ci_upper1, 2)})",
                    'Country_2_n/N': f"{resistant2}/{total2}",
                    'Country_2_Rate': round(rate2, 2),
                    'Country_2_95%CI': f"({round(ci_lower2, 2)}-{round(ci_upper2, 2)})",
                    'Rate_Difference': round(rate1 - rate2, 2),
                    'P_Value': p_val
                })
            
            # Apply BH correction for this antibiotic class and species
            if p_values:
                _, pvals_bh, _, _ = multipletests(p_values, alpha=0.05,
                                                  method='fdr_bh')
                for i, comp in enumerate(comparison_data):
                    comp['P_Value_BH_Corrected'] = pvals_bh[i]
                    comp['Significant_BH'] = pvals_bh[i] < 0.05
                    stats_results.append(comp)
    
    # Also add species comparisons within countries
    for country in countries:
        country_data = temp_df[temp_df['Country'] == country]
        
        for class_name in ANTIBIOTIC_CLASSES.keys():
            p_values = []
            comparison_data = []
            
            for s1, s2 in itertools.combinations(species_list, 2):
                d1_df = country_data[country_data['Species_Name'] == s1]
                d2_df = country_data[country_data['Species_Name'] == s2]
                
                if d1_df.empty or d2_df.empty:
                    continue
                
                d1 = d1_df.iloc[0]
                d2 = d2_df.iloc[0]
                
                resistant1 = int(d1[f'{class_name}_Resistant_Count'])
                total1 = int(d1[f'{class_name}_N'])
                rate1 = d1[f'{class_name}_Rate']
                ci_lower1 = d1[f'{class_name}_CI_Lower']
                ci_upper1 = d1[f'{class_name}_CI_Upper']
                
                resistant2 = int(d2[f'{class_name}_Resistant_Count'])
                total2 = int(d2[f'{class_name}_N'])
                rate2 = d2[f'{class_name}_Rate']
                ci_lower2 = d2[f'{class_name}_CI_Lower']
                ci_upper2 = d2[f'{class_name}_CI_Upper']
                
                if total1 == 0 or total2 == 0:
                    continue
                
                cont_table = [[resistant1, total1 - resistant1],
                              [resistant2, total2 - resistant2]]
                _, p_val = fisher_exact(cont_table)
                
                p_values.append(p_val)
                comparison_data.append({
                    'Comparison_Type': 'Species',
                    'Country': country,
                    'Species_1': s1,
                    'Species_2': s2,
                    'Antibiotic_Class': class_name,
                    'Species_1_n/N': f"{resistant1}/{total1}",
                    'Species_1_Rate': round(rate1, 2),
                    'Species_1_95%CI': f"({round(ci_lower1, 2)}-{round(ci_upper1, 2)})",
                    'Species_2_n/N': f"{resistant2}/{total2}",
                    'Species_2_Rate': round(rate2, 2),
                    'Species_2_95%CI': f"({round(ci_lower2, 2)}-{round(ci_upper2, 2)})",
                    'Rate_Difference': round(rate1 - rate2, 2),
                    'P_Value': p_val
                })
            
            # Apply BH correction
            if p_values:
                _, pvals_bh, _, _ = multipletests(p_values, alpha=0.05,
                                                  method='fdr_bh')
                for i, comp in enumerate(comparison_data):
                    comp['P_Value_BH_Corrected'] = pvals_bh[i]
                    comp['Significant_BH'] = pvals_bh[i] < 0.05
                    stats_results.append(comp)

    stats_summary_df = pd.DataFrame(stats_results)
    
    # Sort by p-value within each comparison type and antibiotic class
    if not stats_summary_df.empty:
        stats_summary_df = stats_summary_df.sort_values(
            by=['Comparison_Type', 'Antibiotic_Class', 'P_Value']
        )

    return resistance_summary_df, stats_summary_df

# =============================================================================
# 3) Heat-map (unchanged apart from doc-string)
# =============================================================================
def create_resistance_heatmap(
    resistance_df,
    row_label_col='Country_Species',
    output_file=None,
    image_width=1200,
    image_height=900,
    image_scale=2
):
    """
    Creates a formatted heat-map. Cell annotations remain the rounded %
    (CIs are exported in Excel, not shown on the plot).
    """
    df = resistance_df.copy()

    # --- Pivot Data ---
    rate_cols = [f'{c}_Rate' for c in ANTIBIOTIC_CLASSES if f'{c}_Rate' in df.columns]
    pivot = df.pivot_table(index=row_label_col, values=rate_cols)
    pivot.columns = [c.replace('_Rate', '') for c in pivot.columns]

    # --- Sort Rows ---
    sort_df = df[[row_label_col, 'Country', 'Species_Name']].drop_duplicates()
    sort_df['Species_Order'] = sort_df['Species_Name'].apply(
        lambda s: 0 if 'jejuni' in str(s).lower() else (1 if 'coli' in str(s).lower() else 2)
    )
    sort_df = sort_df.sort_values(by=['Country', 'Species_Order'])
    pivot = pivot.reindex(index=sort_df[row_label_col].tolist())

    # --- Annotations ---
    ann = []
    for i, row in enumerate(pivot.index):
        for j, col in enumerate(pivot.columns):
            val = pivot.iloc[i, j]
            txt = f'{val:.0f}%' if pd.notna(val) else ''
            txt_col = 'white' if val > 60 else 'black'
            ann.append(dict(x=col, y=row, text=txt, showarrow=False,
                            font=dict(color=txt_col, size=12)))

    # --- Heat-map ---
    fig = go.Figure(data=go.Heatmap(
        z=pivot.values,
        x=pivot.columns,
        y=pivot.index,
        colorscale='Blues',
        zmin=0, zmax=100,
        colorbar=dict(
            title=dict(text="<b>Resistance Rate (%)</b>", font=dict(color='black')),
            tickfont=dict(color='black')
        )
    ))

    fig.update_layout(
        title=dict(text="<b>Resistance Rates Heat-map by Country and Species</b>",
                   font=dict(color='black', size=20), x=0.5),
        xaxis_title="<b>Antibiotic Class</b>",
        yaxis_title="<b>Country - Species</b>",
        font=dict(color="black"),
        yaxis=dict(autorange='reversed',
                   tickfont=dict(color='black', size=12)),
        xaxis=dict(tickfont=dict(color='black', size=12)),
        annotations=ann,
        height=image_height,
        width=image_width,
        margin=dict(l=280, r=50, t=100, b=100)
    )

    if output_file:
        pio.write_image(fig, output_file, scale=image_scale)
        print(f"Saved heat-map to '{output_file}'")

    return fig

# =============================================================================
# 4) Format Excel Output with S10-style Statistical Analysis
# =============================================================================
def save_formatted_excel(res_summary, stats_summary, excel_out):
    """
    Save Excel with formatted statistical analysis sheet similar to S10 style
    """
    with pd.ExcelWriter(excel_out, engine='xlsxwriter') as writer:
        # Save resistance rates summary
        res_summary.to_excel(writer, sheet_name='Resistance_Rates_Summary', index=False)
        
        # Format statistical analysis for country comparisons
        if not stats_summary.empty:
            # Country comparisons
            country_comps = stats_summary[stats_summary['Comparison_Type'] == 'Country'].copy()
            if not country_comps.empty:
                # Reorder columns for better readability
                country_cols = [
                    'Comparison_Type', 'Species', 'Country_1', 'Country_2', 
                    'Antibiotic_Class', 'Country_1_n/N', 'Country_1_Rate',
                    'Country_1_95%CI', 'Country_2_n/N', 'Country_2_Rate', 
                    'Country_2_95%CI', 'Rate_Difference', 'P_Value', 
                    'P_Value_BH_Corrected', 'Significant_BH'
                ]
                country_comps = country_comps[country_cols]
                country_comps.to_excel(writer, sheet_name='Country_Comparisons', index=False)
            
            # Species comparisons
            species_comps = stats_summary[stats_summary['Comparison_Type'] == 'Species'].copy()
            if not species_comps.empty:
                species_cols = [
                    'Comparison_Type', 'Country', 'Species_1', 'Species_2',
                    'Antibiotic_Class', 'Species_1_n/N', 'Species_1_Rate',
                    'Species_1_95%CI', 'Species_2_n/N', 'Species_2_Rate', 
                    'Species_2_95%CI', 'Rate_Difference', 'P_Value',
                    'P_Value_BH_Corrected', 'Significant_BH'
                ]
                species_comps = species_comps[species_cols]
                species_comps.to_excel(writer, sheet_name='Species_Comparisons', index=False)
        
        # Get the workbook and add formatting
        workbook = writer.book
        
        # Define formats
        header_format = workbook.add_format({
            'bold': True,
            'text_wrap': True,
            'valign': 'top',
            'fg_color': '#D7E4BD',
            'border': 1
        })
        
        number_format = workbook.add_format({'num_format': '0.00'})
        percent_format = workbook.add_format({'num_format': '0.00'})
        pval_format = workbook.add_format({'num_format': '0.0000E+00'})
        
        # Format each sheet
        for sheet_name in writer.sheets:
            worksheet = writer.sheets[sheet_name]
            worksheet.set_row(0, 30, header_format)  # Header row
            
            # Auto-adjust column widths
            if sheet_name == 'Country_Comparisons':
                worksheet.set_column('A:A', 15)  # Comparison_Type
                worksheet.set_column('B:B', 25)  # Species
                worksheet.set_column('C:D', 15)  # Countries
                worksheet.set_column('E:E', 18)  # Antibiotic_Class
                worksheet.set_column('F:F', 10)  # Country_1_n/N
                worksheet.set_column('G:G', 12, percent_format)  # Country_1_Rate
                worksheet.set_column('H:H', 15)  # Country_1_95%CI
                worksheet.set_column('I:I', 10)  # Country_2_n/N
                worksheet.set_column('J:J', 12, percent_format)  # Country_2_Rate
                worksheet.set_column('K:K', 15)  # Country_2_95%CI
                worksheet.set_column('L:L', 15, percent_format)  # Rate_Difference
                worksheet.set_column('M:N', 12, pval_format)  # P-values
                worksheet.set_column('O:O', 12)  # Significant
            elif sheet_name == 'Species_Comparisons':
                worksheet.set_column('A:A', 15)  # Comparison_Type
                worksheet.set_column('B:B', 15)  # Country
                worksheet.set_column('C:D', 25)  # Species
                worksheet.set_column('E:E', 18)  # Antibiotic_Class
                worksheet.set_column('F:F', 10)  # Species_1_n/N
                worksheet.set_column('G:G', 12, percent_format)  # Species_1_Rate
                worksheet.set_column('H:H', 15)  # Species_1_95%CI
                worksheet.set_column('I:I', 10)  # Species_2_n/N
                worksheet.set_column('J:J', 12, percent_format)  # Species_2_Rate
                worksheet.set_column('K:K', 15)  # Species_2_95%CI
                worksheet.set_column('L:L', 15, percent_format)  # Rate_Difference
                worksheet.set_column('M:N', 12, pval_format)  # P-values
                worksheet.set_column('O:O', 12)  # Significant

# =============================================================================
# 5) Main
# =============================================================================
def main():
    # --- Paths ---
    input_file = (
        r"C:\Users\BLehr\Documents\Paper_Backup\OHPH_documents\Hub_summary\meta_data_extra\group_useful\selected_contig\other_outputs\merge\for_mic\MIC\updated_file_mic_v2.csv"
    )
    output_dir = Path(
        r"C:\Users\BLehr\Documents\Paper_Backup\OHPH_documents\Hub_summary\meta_data_extra\group_useful\selected_contig\other_outputs\merge\for_mic\MIC\heatmap"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    heatmap_png = output_dir / "resistance_heatmap_formatted.png"
    excel_out   = output_dir / "resistance_stats_summary_S10_format.xlsx"

    # --- Load and process ---
    print(f"Reading data from: {input_file}")
    df = pd.read_csv(input_file)


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

    print("Calculating resistance rates, 95 % CIs and statistical analysisâ€¦")
    res_summary, stats_summary = process_and_analyze_resistance(
        df, grouping_vars=['Country', 'Species_Name']
    )

    # Abbreviate species names for heat-map row labels
    res_summary['Country_Species'] = (
        res_summary['Country_Species']
        .str.replace('Campylobacter jejuni', 'C. jejuni', regex=False)
        .str.replace('Campylobacter coli',   'C. coli',   regex=False)
    )

    # --- Save Excel with S10-style formatting ---
    print(f"Saving Excel report with S10-style formatting to {excel_out} â€¦")
    save_formatted_excel(res_summary, stats_summary, excel_out)
    print("Excel report saved âœ”")

    # --- Save heat-map ---
    print(f"Building heat-map â†’ {heatmap_png} â€¦")
    fig = create_resistance_heatmap(
        res_summary,
        output_file=heatmap_png,
        image_width=1300,
        image_height=800,
        image_scale=3
    )
    print("Heat-map saved âœ”")

    # Display in interactive window (optional)
    fig.show()
    print("\nâœ… All tasks complete!")

# --- Entry Point ---
if __name__ == "__main__":
    main()