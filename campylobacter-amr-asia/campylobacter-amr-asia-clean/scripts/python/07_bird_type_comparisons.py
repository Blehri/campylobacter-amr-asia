"""
Bird Type Statistical Comparisons
=================================
Statistical tests comparing resistance between bird types

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
#bird type
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests
import numpy as np
import os

# ------------------------------------------------------------------------------
# 1. Load and Prepare the Data
# ------------------------------------------------------------------------------
# Define file paths at the top
input_file_path = r"DATA_FILE  # Original path: C:\Users\BLehr\...\My_merged_output_V5.csv"
output_directory = r"C:\Users\BLehr\OneDrive - London School of Hygiene and Tropical Medicine\Documents\OHPH_documents\Hub_summmary\meta_data_extra\group_useful\selected_contig\other_outputs\merge\outputs\Birdtype_ARG\updatev2"

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

data_df = pd.read_csv(input_file_path, encoding='ISO-8859-1')


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

# Clean up and rename columns for clarity and consistency
data_df.columns = data_df.columns.str.strip()
data_df.rename(columns={
    'update_poultry_types_y': 'Bird_Type',
    'Species_Name': 'Species' # Using a consistent name
}, inplace=True)


# ------------------------------------------------------------------------------
# 2. Define Analysis Parameters
# ------------------------------------------------------------------------------
countries_of_interest = ['Bangladesh', 'India', 'Vietnam', 'Sri Lanka']
bird_types = ['Desi', 'Sonali', 'Exotic broiler', 'Hybrid']
species_of_interest = ['Campylobacter jejuni', 'Campylobacter coli']

antibiotic_categories = {
    'ARGs': [
        "AAC(6')-Ie-APH(2'')-Ia bifunctional protein", "APH(2'')-If", "APH(2'')-Ig", "APH(3')-IIIa",
        "APH(3')-VIIa", "ANT(4')-Ia", "aad(6)", "ErmA", "ErmB", "OXA-184", "OXA-185", "OXA-193", "OXA-447",
        "OXA-449", "OXA-465", "OXA-576", "OXA-605", "OXA-617", "OXA-622", "OXA-623", "OXA-625",
        "OXA-631", "OXA-632", "OXA-635", "OXA-638", "cmeA", "cmeB", "cmeC", "cmeR", "tet(L)", "tet(O)", "tet(O/M/O)", "tet(W)", "tet(M)",
        "Campylobacter coli chloramphenicol acetyltransferase", "Limosilactobacillus reuteri cat-TC", "fexA", "lnuC", "lnuP", "SAT-4", "optrA",  
        "vatE", "vatH", "mel", "T86I", "A2075G"
    ],
    'Antibiotic Resistance Gene Classes': [
        'aminoglycoside bifunctional resistance protein', 
        "APH(2'')", "APH(3')", "ANT(4')", 'ANT(6)', 'Erm 23S ribosomal RNA methyltransferase', 'OXA beta-lactamase',
        'resistance-nodulation-cell division (RND) antibiotic efflux pump', 
        'major facilitator superfamily (MFS) antibiotic efflux pump', 
        'tetracycline-resistant ribosomal protection protein', 
        'chloramphenicol acetyltransferase (CAT)', 
        'lincosamide nucleotidyltransferase (LNU)',
        'streptothricin acetyltransferase (SAT)',
        'streptogramin vat acetyltransferase',
        'msr-type ABC-F protein',
        'Miscellaneous ABC-F subfamily ATP-binding cassette ribosomal protection proteins', 
        'fluoroquinolone resistant gyrA', 
        '23S rRNA with mutation conferring resistance to macrolide antibiotics',
    ],
    'Drug Classes': [
        'AminoAB', 'StreptoAB MacroAB LincosAB', 'Ceph Penam Carba', 'Ceph FusiAB MacroAB FluoroAB', 'TetraAB', 'PhenAB', 
        'LincosAB', 'NucleoAB', 'StreptoAB', 'StreptoAB MacroAB', 'PhenAB OxazoAB', 'FluoroAB', 'MacroAB', 
    ]
}


# ------------------------------------------------------------------------------
# 3. Statistical Comparison Functions
# ------------------------------------------------------------------------------
def perform_statistical_test(table, category_name):
    """Helper function to perform Chi-squared or Fisher's exact test."""
    p_value = np.nan
    test_used = 'N/A'
    
    use_fisher = category_name == 'ARGs'
    
    if not use_fisher:
        try:
            _, p_value_chi, _, expected = chi2_contingency(table)
            if (expected < 5).any():
                use_fisher = True
                test_used = 'Fisher (Fallback)'
            else:
                p_value = p_value_chi
                test_used = 'Chi-squared'
        except (ValueError, ZeroDivisionError):
            use_fisher = True

    if use_fisher:
        try:
            _, p_value = fisher_exact(table, alternative='two-sided')
            # Only set test_used if it hasn't been set by fallback logic
            if test_used == 'N/A':
                 test_used = 'Fisher'
        except ValueError:
            p_value = np.nan

    return p_value, test_used

def compare_resistance_across_countries(df, countries, bird_types_list, species_list, categories):
    """Compares two countries for the same bird type."""
    all_results = []
    for category_name, antibiotics in categories.items():
        for antibiotic in antibiotics:
            for bird_type in bird_types_list:
                for specie in species_list:
                    for i, country1 in enumerate(countries):
                        for country2 in countries[i + 1:]:
                            data1 = df[(df['Country'] == country1) & (df['Bird_Type'] == bird_type) & (df['Species'] == specie)]
                            data2 = df[(df['Country'] == country2) & (df['Bird_Type'] == bird_type) & (df['Species'] == specie)]

                            if data1.empty or data2.empty or antibiotic not in df.columns:
                                continue

                            table = np.array([
                                [data1[antibiotic].sum(), max(0, data1.shape[0] - data1[antibiotic].sum())],
                                [data2[antibiotic].sum(), max(0, data2.shape[0] - data2[antibiotic].sum())]
                            ])

                            if np.all(table == 0): continue
                            
                            p_value, test_used = perform_statistical_test(table, category_name)
                            n1, N1 = int(table[0, 0]), int(table[0].sum())
                            n2, N2 = int(table[1, 0]), int(table[1].sum())

                            all_results.append({
                                'Category': category_name, 'Bird Type': bird_type, 'Species': specie, 'Antibiotic': antibiotic,
                                'Country1': country1, 'Count1 (n/N)': f"{n1}/{N1}", 'Percentage1': (n1 / N1 * 100) if N1 > 0 else 0,
                                'Country2': country2, 'Count2 (n/N)': f"{n2}/{N2}", 'Percentage2': (n2 / N2 * 100) if N2 > 0 else 0,
                                'P-Value': p_value, 'Test Used': test_used
                            })

    results_df = pd.DataFrame(all_results)
    if not results_df.empty and 'P-Value' in results_df.columns:
        p_values = results_df['P-Value'].dropna()
        if not p_values.empty:
            _, q_values, _, _ = multipletests(p_values, method='fdr_bh')
            results_df['Corrected P-Value (q-value)'] = np.nan
            results_df.loc[p_values.index, 'Corrected P-Value (q-value)'] = q_values

    col_order = ['Category', 'Bird Type', 'Species', 'Antibiotic', 'Country1', 'Count1 (n/N)', 'Percentage1',
                 'Country2', 'Count2 (n/N)', 'Percentage2', 'P-Value', 'Corrected P-Value (q-value)', 'Test Used']
    return results_df[[col for col in col_order if col in results_df.columns]]


def compare_resistance_within_countries(df, countries, bird_types_list, species_list, categories):
    """Compares two bird types within the same country."""
    all_results = []
    for category_name, antibiotics in categories.items():
        for antibiotic in antibiotics:
            for country in countries:
                for specie in species_list:
                    for i, bird_type1 in enumerate(bird_types_list):
                        for bird_type2 in bird_types_list[i + 1:]:
                            data1 = df[(df['Country'] == country) & (df['Bird_Type'] == bird_type1) & (df['Species'] == specie)]
                            data2 = df[(df['Country'] == country) & (df['Bird_Type'] == bird_type2) & (df['Species'] == specie)]

                            if data1.empty or data2.empty or antibiotic not in df.columns:
                                continue

                            table = np.array([
                                [data1[antibiotic].sum(), max(0, data1.shape[0] - data1[antibiotic].sum())],
                                [data2[antibiotic].sum(), max(0, data2.shape[0] - data2[antibiotic].sum())]
                            ])

                            if np.all(table == 0): continue
                            
                            p_value, test_used = perform_statistical_test(table, category_name)
                            n1, N1 = int(table[0, 0]), int(table[0].sum())
                            n2, N2 = int(table[1, 0]), int(table[1].sum())

                            all_results.append({
                                'Category': category_name, 'Country': country, 'Species': specie, 'Antibiotic': antibiotic,
                                'Bird Type1': bird_type1, 'Count1 (n/N)': f"{n1}/{N1}", 'Percentage1': (n1 / N1 * 100) if N1 > 0 else 0,
                                'Bird Type2': bird_type2, 'Count2 (n/N)': f"{n2}/{N2}", 'Percentage2': (n2 / N2 * 100) if N2 > 0 else 0,
                                'P-Value': p_value, 'Test Used': test_used
                            })
                            
    results_df = pd.DataFrame(all_results)
    if not results_df.empty and 'P-Value' in results_df.columns:
        p_values = results_df['P-Value'].dropna()
        if not p_values.empty:
            _, q_values, _, _ = multipletests(p_values, method='fdr_bh')
            results_df['Corrected P-Value (q-value)'] = np.nan
            results_df.loc[p_values.index, 'Corrected P-Value (q-value)'] = q_values
            
    col_order = ['Category', 'Country', 'Species', 'Antibiotic', 'Bird Type1', 'Count1 (n/N)', 'Percentage1',
                 'Bird Type2', 'Count2 (n/N)', 'Percentage2', 'P-Value', 'Corrected P-Value (q-value)', 'Test Used']
    return results_df[[col for col in col_order if col in results_df.columns]]


# ------------------------------------------------------------------------------
# 4. Run Comparisons and Save Results to Separate Excel Files
# ------------------------------------------------------------------------------
print("Running across-country comparisons (by bird type)...")
results_across = compare_resistance_across_countries(
    data_df, countries_of_interest, bird_types, species_of_interest, antibiotic_categories
)

print("Running within-country comparisons (by bird type)...")
results_within = compare_resistance_within_countries(
    data_df, countries_of_interest, bird_types, species_of_interest, antibiotic_categories
)

# Define two separate file paths for the bird type analysis
output_path_across = os.path.join(output_directory, 'across_country_comparisons_birdType.xlsx')
output_path_within = os.path.join(output_directory, 'within_country_comparisons_birdType.xlsx')

# Save each result set to its own Excel file
print(f"\nSaving across-country results to: {output_path_across}")
results_across.to_excel(output_path_across, engine='xlsxwriter', index=False)

print(f"Saving within-country results to: {output_path_within}")
results_within.to_excel(output_path_within, engine='xlsxwriter', index=False)

print("\nAnalysis complete. All results have been saved to separate files. âœ…")