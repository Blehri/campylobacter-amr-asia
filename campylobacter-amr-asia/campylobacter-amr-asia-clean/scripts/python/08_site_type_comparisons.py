"""
Site Type Statistical Comparisons
=================================
Statistical tests comparing resistance between site types (Table S7)

Author: [Burhan Lehri]
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
#site type
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests
import numpy as np
import os

# ------------------------------------------------------------------------------
# 1. Load the data
# ------------------------------------------------------------------------------
# Define file paths at the top
input_file_path = r"DATA_FILE  # Original path: C:\Users\BLehr\...\My_merged_output_V5.csv"
output_directory = r"C:\Users\BLehr\OneDrive - London School of Hygiene and Tropical Medicine\Documents\OHPH_documents\Hub_summmary\meta_data_extra\group_useful\selected_contig\other_outputs\merge\outputs\Sitetype_ARG\adj\merged"

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

# Clean up column names
data_df.columns = data_df.columns.str.strip()
data_df.rename(columns={'Site description': 'Site_description', 'Species Name': 'Species_Name'}, inplace=True)


# ------------------------------------------------------------------------------
# 2. Define countries, site types, species, and antibiotic categories
# ------------------------------------------------------------------------------
countries_of_interest = ['Bangladesh', 'India', 'Vietnam', 'Sri Lanka']
site_types = ['Farm', 'Market', 'Meat processing facility']
species_of_interest = ['Campylobacter jejuni', 'Campylobacter coli']

antibiotic_categories = {
    'ARGs': [
        "AAC(6')-Ie-APH(2'')-Ia bifunctional protein", "APH(2'')-If", "APH(2'')-Ig", "APH(3')-IIIa",
        "APH(3')-VIIa", "ANT(4')-Ia", "aad(6)", "ErmA", "ErmB", "OXA-184", "OXA-185", "OXA-193", "OXA-447",
        "OXA-449", "OXA-465", "OXA-576", "OXA-605", "OXA-617", "OXA-622", "OXA-623", "OXA-625",
        "OXA-631", "OXA-632", "OXA-635", "OXA-638", "cmeA", "cmeB", "cmeC", "cmeR", "tet(L)", "tet(O)",
        "tet(O/M/O)", "tet(W)", "tet(M)", "Campylobacter coli chloramphenicol acetyltransferase",
        "Limosilactobacillus reuteri cat-TC", "fexA", "lnuC", "lnuP", "SAT-4", "optrA", "vatE",
        "vatH", "mel", "T86I", "A2075G"
    ],
    'Antibiotic Resistance Gene Classes': [
        'aminoglycoside bifunctional resistance protein',
        "APH(2'')", "APH(3')", "ANT(4')", 'ANT(6)',
        'Erm 23S ribosomal RNA methyltransferase', 'OXA beta-lactamase',
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
        'AminoAB', 'StreptoAB MacroAB LincosAB', 'Ceph Penam Carba',
        'Ceph FusiAB MacroAB FluoroAB', 'TetraAB', 'PhenAB',
        'LincosAB', 'NucleoAB', 'StreptoAB', 'StreptoAB MacroAB',
        'PhenAB OxazoAB', 'FluoroAB', 'MacroAB'
    ]
}

# ------------------------------------------------------------------------------
# 3. Statistical comparison functions
# ------------------------------------------------------------------------------
def perform_statistical_test(table, category_name):
    """Helper function to perform Chi-squared or Fisher's exact test."""
    p_value = np.nan
    test_used = 'N/A'
    
    # Use Fisher's for ARGs or if any expected count is < 5
    use_fisher = category_name == 'ARGs'
    
    if not use_fisher:
        try:
            _, p_value_chi, _, expected = chi2_contingency(table)
            if (expected < 5).any():
                use_fisher = True  # Fallback to Fisher
            else:
                p_value = p_value_chi
                test_used = 'Chi-squared'
        except (ValueError, ZeroDivisionError):
            use_fisher = True # Fallback if Chi-squared fails

    if use_fisher:
        try:
            _, p_value = fisher_exact(table, alternative='two-sided')
            test_used = 'Fisher'
        except ValueError:
            p_value = np.nan # Test fails (e.g., all zeros in a row/column)

    return p_value, test_used

def compare_resistance_across_countries(df, countries, site_types, species_list, categories):
    """Compare resistance across countries, returning a detailed DataFrame."""
    all_results = []
    for category_name, antibiotics in categories.items():
        for antibiotic in antibiotics:
            for i, country1 in enumerate(countries):
                for country2 in countries[i + 1:]:
                    for site_type in site_types:
                        for specie in species_list:
                            data1 = df[(df['Country'] == country1) & (df['Site_description'] == site_type) & (df['Species_Name'] == specie)]
                            data2 = df[(df['Country'] == country2) & (df['Site_description'] == site_type) & (df['Species_Name'] == specie)]

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
                                'Category': category_name, 'Site Type': site_type, 'Species': specie, 'Antibiotic': antibiotic,
                                'Country1': country1, 'Count1 (n/N)': f"{n1}/{N1}", 'Percentage1': (n1 / N1 * 100) if N1 > 0 else 0,
                                'Country2': country2, 'Count2 (n/N)': f"{n2}/{N2}", 'Percentage2': (n2 / N2 * 100) if N2 > 0 else 0,
                                'P-Value': p_value, 'Test Used': test_used
                            })

    results_df = pd.DataFrame(all_results)
    if not results_df.empty:
        p_values = results_df['P-Value'].dropna()
        if not p_values.empty:
            _, q_values, _, _ = multipletests(p_values, method='fdr_bh')
            results_df['Corrected P-Value (q-value)'] = np.nan
            results_df.loc[p_values.index, 'Corrected P-Value (q-value)'] = q_values

    # Reorder columns for clarity
    col_order = ['Category', 'Site Type', 'Species', 'Antibiotic', 'Country1', 'Count1 (n/N)', 'Percentage1',
                 'Country2', 'Count2 (n/N)', 'Percentage2', 'P-Value', 'Corrected P-Value (q-value)', 'Test Used']
    return results_df[[col for col in col_order if col in results_df.columns]]


def compare_resistance_within_countries(df, countries, site_types_list, species_list, categories):
    """Compare resistance within countries across site types, returning a detailed DataFrame."""
    all_results = []
    for category_name, antibiotics in categories.items():
        for antibiotic in antibiotics:
            for country in countries:
                for specie in species_list:
                    for i, site_type1 in enumerate(site_types_list):
                        for site_type2 in site_types_list[i + 1:]:
                            data1 = df[(df['Country'] == country) & (df['Site_description'] == site_type1) & (df['Species_Name'] == specie)]
                            data2 = df[(df['Country'] == country) & (df['Site_description'] == site_type2) & (df['Species_Name'] == specie)]

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
                                'Site Type1': site_type1, 'Count1 (n/N)': f"{n1}/{N1}", 'Percentage1': (n1 / N1 * 100) if N1 > 0 else 0,
                                'Site Type2': site_type2, 'Count2 (n/N)': f"{n2}/{N2}", 'Percentage2': (n2 / N2 * 100) if N2 > 0 else 0,
                                'P-Value': p_value, 'Test Used': test_used
                            })
                            
    results_df = pd.DataFrame(all_results)
    if not results_df.empty:
        p_values = results_df['P-Value'].dropna()
        if not p_values.empty:
            _, q_values, _, _ = multipletests(p_values, method='fdr_bh')
            results_df['Corrected P-Value (q-value)'] = np.nan
            results_df.loc[p_values.index, 'Corrected P-Value (q-value)'] = q_values
            
    # Reorder columns for clarity
    col_order = ['Category', 'Country', 'Species', 'Antibiotic', 'Site Type1', 'Count1 (n/N)', 'Percentage1',
                 'Site Type2', 'Count2 (n/N)', 'Percentage2', 'P-Value', 'Corrected P-Value (q-value)', 'Test Used']
    return results_df[[col for col in col_order if col in results_df.columns]]


# ------------------------------------------------------------------------------
# 4. Run comparisons and save results to SEPARATE Excel files
# ------------------------------------------------------------------------------
print("Running across-country comparisons...")
results_across = compare_resistance_across_countries(
    data_df, countries_of_interest, site_types, species_of_interest, antibiotic_categories
)

print("Running within-country comparisons...")
results_within = compare_resistance_within_countries(
    data_df, countries_of_interest, site_types, species_of_interest, antibiotic_categories
)

# Define two separate file paths
output_path_across = os.path.join(output_directory, 'across_country_comparisons_siteType.xlsx')
output_path_within = os.path.join(output_directory, 'within_country_comparisons_siteType.xlsx')

# Save each result set to its own file
print(f"\nSaving across-country results to: {output_path_across}")
results_across.to_excel(output_path_across, engine='xlsxwriter', index=False)

print(f"Saving within-country results to: {output_path_within}")
results_within.to_excel(output_path_within, engine='xlsxwriter', index=False)

print("\nAnalysis complete. All results have been saved to separate files. âœ…")
