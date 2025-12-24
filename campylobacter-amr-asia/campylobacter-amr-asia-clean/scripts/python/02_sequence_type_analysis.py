"""
Sequence Type Analysis
======================
ST-specific resistance distributions and statistical comparisons (Figure 5)

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

# v502 species and country ST
#enhanced stats raw data
# #ST code: Requires minimum 5 samples per group
#adj style THISONE
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import to_hex
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests
import os

################################################################################
#                           FILE I/O AND DATA READING                          #
################################################################################

# Ensure the output directory exists
output_dir = r"OUTPUT_DIR  # Original: C:\Users\BLehr\Desktop\analysis_output"
#"C:\Users\BLehr\OneDrive - London School of Hygiene and Tropical Medicine\Documents\OHPH_documents\Hub_summmary\meta_data_extra\group_useful\selected_contig\other_outputs\merge\outputs\ST\ST_country\adj\n\N_V2"
os.makedirs(output_dir, exist_ok=True)

# Load the CSV file
file_path = r"DATA_FILE  # Original path: C:\Users\BLehr\...\My_merged_output_V6.csv"
data = pd.read_csv(file_path)


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

# Filter the data for the countries of interest
countries_of_interest = ['Bangladesh', 'India', 'Vietnam', 'Sri Lanka']
data_filtered = data[data['Country'].isin(countries_of_interest)]

################################################################################
#                      BASIC CLEAN-UP & COLUMN ADJUSTMENTS                     #
################################################################################

# Replace '-' in 'ST' with 'Uncharacterized'
data_filtered['ST'] = data_filtered['ST'].replace('-', 'Uncharacterized')

# Create a new column to specify the species type
data_filtered['Species_Type'] = data_filtered['Species_Name'].apply(
    lambda x: 'jejuni' if 'jejuni' in str(x).lower() else 'coli' if 'coli' in str(x).lower() else 'Unknown'
)

# Ensure uncharacterized STs are differentiated by their species
data_filtered.loc[data_filtered['ST'] == 'Uncharacterized', 'ST'] = (
    data_filtered['Species_Type'] + '_Uncharacterized'
)

# Replace spaces in column names with underscores for easier access
data_filtered.columns = data_filtered.columns.str.replace(' ', '_')

def replace_spaces_in_list(lst):
    return [item.replace(' ', '_') for item in lst]

################################################################################
#                          DEFINE ANTIBIOTIC CATEGORIES                        #
################################################################################

antibiotic_categories = {
    'ARGs': [
        "AAC(6')-Ie-APH(2'')-Ia_bifunctional_protein", "APH(2'')-If", "APH(2'')-Ig", "APH(3')-IIIa",
        "APH(3')-VIIa", "ANT(4')-Ia", "aad(6)", "ErmA", "ErmB","OXA-184", "OXA-185", "OXA-193", "OXA-447",
        "OXA-449", "OXA-465", "OXA-576", "OXA-605", "OXA-617", "OXA-622", "OXA-623", "OXA-625",
        "OXA-631", "OXA-632", "OXA-635", "OXA-638", "cmeA", "cmeB", "cmeC", "cmeR", "tet(L)", "tet(O)", 
        "tet(O/M/O)", "tet(W)", "tet(M)",
        "Campylobacter_coli_chloramphenicol_acetyltransferase", 
        "Limosilactobacillus_reuteri_cat-TC", 
        "fexA", "lnuC", "lnuP", "SAT-4", "optrA",  
        "vatE", "vatH", "mel", "vgaD",  # added "vgaD" so it doesn't cause KeyError
        "T86I", "A2075G"
    ],
    'Antibiotic_Resistance_Gene_Classes': replace_spaces_in_list([
        'aminoglycoside bifunctional resistance protein', 
        "APH(2'')", "APH(3')", "ANT(4')", 'ANT(6)', 'Erm 23S ribosomal RNA methyltransferase', 
        'OXA beta-lactamase',
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
    ]),
    'Drug_Classes': replace_spaces_in_list([
        'AminoAB', 'StreptoAB MacroAB LincosAB', 'Ceph Penam Carba', 'Ceph FusiAB MacroAB FluoroAB', 
        'TetraAB', 'PhenAB', 'LincosAB', 'NucleoAB', 'StreptoAB', 
        'StreptoAB MacroAB', 'PhenAB OxazoAB','FluoroAB', 'MacroAB'
    ])
}

################################################################################
#                             CUSTOM COLOR PALETTES                             #
################################################################################

custom_colors = {
    'Drug_Classes': {
        'AminoAB': '#FFA07A',
        'Ceph_Penam_Carba': '#DCDCDC',
        'TetraAB': '#FF0000',
        'NucleoAB': '#FFC0CB',
        'PhenAB_OxazoAB': '#800080',
        'FluoroAB': '#1E90FF',
        'MacroAB': '#6e2c00',
        'Ceph_FusiAB_MacroAB_FluoroAB': '#00CED1',
        'LincosAB': '#FFFF00',
        'PhenAB': '#008000',
        'StreptoAB_MacroAB_LincosAB': '#FFFFE0',
        'StreptoAB': '#8A2BE2',  
        'StreptoAB_MacroAB': '#9370DB'
    },
    'Antibiotic_Resistance_Gene_Classes': {
        "aminoglycoside_bifunctional_resistance_protein": '#FF8C00',
        "APH(2'')": '#FFA07A',
        "APH(3')": '#FF7F50',
        "ANT(4')": '#FF69B4',
        "ANT(6)": '#FF4500',
        "resistance-nodulation-cell_division_(RND)_antibiotic_efflux_pump": '#00CED1',
        "chloramphenicol_acetyltransferase_(CAT)": '#008000',
        "lincosamide_nucleotidyltransferase_(LNU)": '#FFFF00',
        "Erm_23S_ribosomal_RNA_methyltransferase": '#FFFFE0',
        "OXA_beta-lactamase": '#DCDCDC',
        "tetracycline-resistant_ribosomal_protection_protein": '#DC143C',
        "major_facilitator_superfamily_(MFS)_antibiotic_efflux_pump": '#00FF7F',
        "streptothricin_acetyltransferase_(SAT)": '#FFC0CB',
        "Miscellaneous_ABC-F_subfamily_ATP-binding_cassette_ribosomal_protection_proteins": '#800080',
        "fluoroquinolone_resistant_gyrA": '#1E90FF',
        "23S_rRNA_with_mutation_conferring_resistance_to_macrolide_antibiotics": '#6e2c00',
        "streptogramin_vat_acetyltransferase": '#8A2BE2',
        "msr-type_ABC-F_protein": '#9370DB'
    },
    'ARGs': {
        "APH(2'')-If": '#FFA07A',
        "AAC(6')-Ie-APH(2'')-Ia_bifunctional_protein": '#FF8C00',
        "APH(2'')-Ig": '#FFA500',
        "APH(3')-IIIa": '#FF7F50',
        "APH(3')-VIIa": '#FF6347',
        "aad(6)": '#FF4500',
        "ANT(4')-Ia": '#edc37e',
        "cmeA": '#00CED1',
        "cmeB": '#20B2AA',
        "cmeC": '#40E0D0',
        "cmeR": '#48D1CC',
        "OXA-184": '#DCDCDC',
        "OXA-185": '#D3D3D3',
        "OXA-193": '#C0C0C0',
        "OXA-447": '#A9A9A9',
        "OXA-449": '#808080',
        "OXA-465": '#696969',
        "OXA-576": '#778899',
        "OXA-605": '#708090',
        "OXA-617": '#696969',
        "OXA-622": '#2F4F4F',
        "OXA-623": '#DCDCDC',
        "OXA-625": '#D3D3D3',
        "OXA-631": '#C0C0C0',
        "OXA-632": '#A9A9A9',
        "OXA-635": '#808080',
        "OXA-638": '#696969',
        "T86I": '#1E90FF',
        "lnuP": '#FFFF00',
        "lnuC": '#FFD700',
        "A2075G": '#6e2c00',
        "SAT-4": '#FFC0CB',
        "Campylobacter_coli_chloramphenicol_acetyltransferase": '#008000',
        "Limosilactobacillus_reuteri_cat-TC": '#00FA9A',
        "fexA": '#00FF7F',
        "optrA": '#800080',
        "ErmA": '#FFFFE0',
        "ErmB": '#FFFACD',
        "tet(L)": '#FF0000',
        "tet(O)": '#DC143C',
        "tet(O/M/O)": '#B22222',
        "tet(M)": '#cd4c50',
        "tet(W)": '#8B0000',
        "vatE": '#BA55D3',
        "vatH": '#9370DB',
        "vgaD": '#8A2BE2',
        "mel": '#9370DB'
    }
}

def replace_spaces_in_dict_keys(d):
    return {k.replace(' ', '_'): v for k, v in d.items()}

for cat in custom_colors:
    custom_colors[cat] = replace_spaces_in_dict_keys(custom_colors[cat])

################################################################################
#               STATISTICAL COMPARISONS (ACROSS and WITHIN STs)                #
################################################################################

def compare_resistance_across_sequence_types(df, countries, species_list, categories):
    """Perform statistical comparison across STs within same country & species."""
    all_results = []
    for category_name, antibiotics in categories.items():
        for antibiotic in antibiotics:
            for country in countries:
                for specie in species_list:
                    st_list = df[(df['Country'] == country) & (df['Species_Type'] == specie)]['ST'].unique()
                    for i, st1 in enumerate(st_list):
                        for st2 in st_list[i + 1:]:
                            data1 = df[(df['Country'] == country) & (df['Species_Type'] == specie) & (df['ST'] == st1)][antibiotic]
                            data2 = df[(df['Country'] == country) & (df['Species_Type'] == specie) & (df['ST'] == st2)][antibiotic]
                            
                            if data1.empty or data2.empty:
                                continue
                            if data1.count() < 5 or data2.count() < 5:
                                continue  # skip if < 5 samples per group

                            # Binary data => use Fisher's
                            if category_name == 'ARGs':
                                table = [
                                    [data1.sum(), data1.count() - data1.sum()],
                                    [data2.sum(), data2.count() - data2.sum()]
                                ]
                                # sanity check
                                if any(x < 0 for row in table for x in row):
                                    continue
                                if any(sum(row) == 0 for row in table):
                                    continue
                                _, p_value = fisher_exact(table, alternative='two-sided')
                                test_used = 'Fisher'
                            else:
                                # Chi-squared or fallback
                                table = np.array([
                                    [data1.sum(), data1.count() - data1.sum()],
                                    [data2.sum(), data2.count() - data2.sum()]
                                ])
                                if any(x < 0 for row in table for x in row):
                                    continue
                                if any(sum(row) == 0 for row in table):
                                    continue
                                try:
                                    _, p_value, _, _ = chi2_contingency(table)
                                    test_used = 'Chi-squared'
                                except ValueError:
                                    _, p_value = fisher_exact(table, alternative='two-sided')
                                    test_used = 'Fisher (Fallback)'

                            # Add the actual counts to the results
                            all_results.append({
                                'Category': category_name,
                                'Country': country,
                                'Species': specie,
                                'Antibiotic': antibiotic,
                                'Sequence Type 1': st1,
                                'Sequence Type 2': st2,
                                'ST1_Count': data1.count(),
                                'ST2_Count': data2.count(),
                                'ST1_Positive': data1.sum(),
                                'ST2_Positive': data2.sum(),
                                'P-Value': p_value,
                                'Test Used': test_used
                            })

    results_df = pd.DataFrame(all_results)
    if not results_df.empty:
        _, corrected_p_values, _, _ = multipletests(results_df['P-Value'], method='fdr_bh')
        results_df['Corrected P-Value'] = corrected_p_values
    else:
        results_df['Corrected P-Value'] = []
    return results_df

def compare_resistance_within_sequence_type(df, countries, species_list, categories):
    """Perform statistical comparison across countries for the SAME ST & species."""
    all_results = []
    for category_name, antibiotics in categories.items():
        for antibiotic in antibiotics:
            for specie in species_list:
                st_list = df[df['Species_Type'] == specie]['ST'].unique()
                for st in st_list:
                    # only compare if ST is present in at least 2 countries
                    countries_with_st = df[(df['Species_Type'] == specie) & (df['ST'] == st)]['Country'].unique()
                    if len(countries_with_st) < 2:
                        continue
                    for i, country1 in enumerate(countries):
                        for country2 in countries[i + 1:]:
                            if country1 not in countries_with_st or country2 not in countries_with_st:
                                continue
                            data1 = df[(df['Country'] == country1) & (df['Species_Type'] == specie) & (df['ST'] == st)][antibiotic]
                            data2 = df[(df['Country'] == country2) & (df['Species_Type'] == specie) & (df['ST'] == st)][antibiotic]
                            
                            if data1.empty or data2.empty:
                                continue
                            if data1.count() < 5 or data2.count() < 5:
                                continue

                            if category_name == 'ARGs':
                                table = [
                                    [data1.sum(), data1.count() - data1.sum()],
                                    [data2.sum(), data2.count() - data2.sum()]
                                ]
                                if any(x < 0 for row in table for x in row):
                                    continue
                                if any(sum(row) == 0 for row in table):
                                    continue
                                _, p_value = fisher_exact(table, alternative='two-sided')
                                test_used = 'Fisher'
                            else:
                                table = np.array([
                                    [data1.sum(), data1.count() - data1.sum()],
                                    [data2.sum(), data2.count() - data2.sum()]
                                ])
                                if any(x < 0 for row in table for x in row):
                                    continue
                                if any(sum(row) == 0 for row in table):
                                    continue
                                try:
                                    _, p_value, _, _ = chi2_contingency(table)
                                    test_used = 'Chi-squared'
                                except ValueError:
                                    _, p_value = fisher_exact(table, alternative='two-sided')
                                    test_used = 'Fisher (Fallback)'

                            all_results.append({
                                'Category': category_name,
                                'Species': specie,
                                'Sequence Type': st,
                                'Antibiotic': antibiotic,
                                'Country 1': country1,
                                'Country 2': country2,
                                'Country1_Count': data1.count(),
                                'Country2_Count': data2.count(),
                                'Country1_Positive': data1.sum(),
                                'Country2_Positive': data2.sum(),
                                'P-Value': p_value,
                                'Test Used': test_used
                            })

    results_df = pd.DataFrame(all_results)
    if not results_df.empty:
        _, corrected_p_values, _, _ = multipletests(results_df['P-Value'], method='fdr_bh')
        results_df['Corrected P-Value'] = corrected_p_values
    else:
        results_df['Corrected P-Value'] = []
    return results_df

################################################################################
#                CALCULATING PERCENTAGES FOR ARG/CLASS PRESENCE                #
################################################################################

def calculate_percentages(df, countries, species_list, categories):
    """Calculate % presence of each antibiotic or gene-class for each ST."""
    percentages_dict = {}
    # New dictionary to store raw counts
    raw_counts_dict = {}
    
    for country in countries:
        for specie in species_list:
            for category_name, antibiotics in categories.items():
                st_list = df[(df['Country'] == country) & (df['Species_Type'] == specie)]['ST'].unique()
                percentages = pd.DataFrame(0, index=st_list, columns=antibiotics)
                
                # Create DataFrames for numerators and denominators
                numerators = pd.DataFrame(0, index=st_list, columns=antibiotics)
                denominators = pd.DataFrame(0, index=st_list, columns=antibiotics)
                
                for st in st_list:
                    st_data = df[
                        (df['Country'] == country) &
                        (df['Species_Type'] == specie) &
                        (df['ST'] == st)
                    ]
                    if st_data.empty:
                        continue
                    for antibiotic in antibiotics:
                        if antibiotic not in st_data.columns:
                            continue
                        if category_name == 'ARGs':
                            presence = (st_data[antibiotic] > 0).sum()
                            total = st_data[antibiotic].count()
                            numerators.loc[st, antibiotic] = presence
                            denominators.loc[st, antibiotic] = total
                            percentages.loc[st, antibiotic] = (presence / total) * 100 if total > 0 else np.nan
                        else:
                            # For non-ARG categories, we're storing total usage counts
                            usage = st_data[antibiotic].sum()
                            numerators.loc[st, antibiotic] = usage
                            # We'll calculate the denominator (total for all antibiotics) later
                            percentages.loc[st, antibiotic] = usage
                
                if category_name != 'ARGs':
                    # For non-ARG categories, calculate percentages based on total usage
                    row_sums = percentages.sum(axis=1)
                    for st in st_list:
                        for antibiotic in antibiotics:
                            # Store the total as denominator for each antibiotic
                            denominators.loc[st, antibiotic] = row_sums[st] if st in row_sums else 0
                    
                    # Now calculate the percentages
                    percentages = percentages.div(row_sums, axis=0).fillna(0) * 100
                
                percentages_dict[(country, specie, category_name)] = percentages
                raw_counts_dict[(country, specie, category_name)] = {'numerators': numerators, 'denominators': denominators}
    
    return percentages_dict, raw_counts_dict

def add_percentages(significant_df, percentages_dict, raw_counts_dict):
    """Add Percentage1 and Percentage2 for across ST comparisons."""
    percentage1_list = []
    percentage2_list = []
    numerator1_list = []
    numerator2_list = []
    denominator1_list = []
    denominator2_list = []
    n1_list = []
    n2_list = []

    for _, row in significant_df.iterrows():
        country = row['Country']
        species = row['Species']
        category_name = row['Category']
        antibiotic = row['Antibiotic']
        st1 = row['Sequence Type 1']
        st2 = row['Sequence Type 2']

        percentages = percentages_dict.get((country, species, category_name), pd.DataFrame())
        raw_counts = raw_counts_dict.get((country, species, category_name), {'numerators': pd.DataFrame(), 'denominators': pd.DataFrame()})
        
        # Extract values
        percentage1 = percentages.at[st1, antibiotic] if not percentages.empty and st1 in percentages.index else np.nan
        percentage2 = percentages.at[st2, antibiotic] if not percentages.empty and st2 in percentages.index else np.nan
        
        # Get raw counts
        numerator1 = raw_counts['numerators'].at[st1, antibiotic] if not raw_counts['numerators'].empty and st1 in raw_counts['numerators'].index else np.nan
        numerator2 = raw_counts['numerators'].at[st2, antibiotic] if not raw_counts['numerators'].empty and st2 in raw_counts['numerators'].index else np.nan
        denominator1 = raw_counts['denominators'].at[st1, antibiotic] if not raw_counts['denominators'].empty and st1 in raw_counts['denominators'].index else np.nan
        denominator2 = raw_counts['denominators'].at[st2, antibiotic] if not raw_counts['denominators'].empty and st2 in raw_counts['denominators'].index else np.nan
        
        # Format n= strings
        n1 = f"n={int(numerator1)}/{int(denominator1)}" if not np.isnan(numerator1) and not np.isnan(denominator1) else "n=NA"
        n2 = f"n={int(numerator2)}/{int(denominator2)}" if not np.isnan(numerator2) and not np.isnan(denominator2) else "n=NA"

        percentage1_list.append(percentage1)
        percentage2_list.append(percentage2)
        numerator1_list.append(numerator1)
        numerator2_list.append(numerator2)
        denominator1_list.append(denominator1)
        denominator2_list.append(denominator2)
        n1_list.append(n1)
        n2_list.append(n2)

    significant_df['Percentage1'] = percentage1_list
    significant_df['Percentage2'] = percentage2_list
    significant_df['Numerator1'] = numerator1_list
    significant_df['Numerator2'] = numerator2_list
    significant_df['Denominator1'] = denominator1_list
    significant_df['Denominator2'] = denominator2_list
    significant_df['n1'] = n1_list
    significant_df['n2'] = n2_list
    
    return significant_df

def add_percentages_within(significant_df, percentages_dict, raw_counts_dict):
    """Add Percentage1 and Percentage2 for within ST (across countries)."""
    percentage1_list = []
    percentage2_list = []
    numerator1_list = []
    numerator2_list = []
    denominator1_list = []
    denominator2_list = []
    n1_list = []
    n2_list = []

    for _, row in significant_df.iterrows():
        country1 = row['Country 1']
        country2 = row['Country 2']
        species = row['Species']
        category_name = row['Category']
        antibiotic = row['Antibiotic']
        st = row['Sequence Type']

        percentages1 = percentages_dict.get((country1, species, category_name), pd.DataFrame())
        percentages2 = percentages_dict.get((country2, species, category_name), pd.DataFrame())
        
        raw_counts1 = raw_counts_dict.get((country1, species, category_name), {'numerators': pd.DataFrame(), 'denominators': pd.DataFrame()})
        raw_counts2 = raw_counts_dict.get((country2, species, category_name), {'numerators': pd.DataFrame(), 'denominators': pd.DataFrame()})

        percentage1 = percentages1.at[st, antibiotic] if not percentages1.empty and st in percentages1.index else np.nan
        percentage2 = percentages2.at[st, antibiotic] if not percentages2.empty and st in percentages2.index else np.nan
        
        # Get raw counts
        numerator1 = raw_counts1['numerators'].at[st, antibiotic] if not raw_counts1['numerators'].empty and st in raw_counts1['numerators'].index else np.nan
        numerator2 = raw_counts2['numerators'].at[st, antibiotic] if not raw_counts2['numerators'].empty and st in raw_counts2['numerators'].index else np.nan
        denominator1 = raw_counts1['denominators'].at[st, antibiotic] if not raw_counts1['denominators'].empty and st in raw_counts1['denominators'].index else np.nan
        denominator2 = raw_counts2['denominators'].at[st, antibiotic] if not raw_counts2['denominators'].empty and st in raw_counts2['denominators'].index else np.nan
        
        # Format n= strings
        n1 = f"n={int(numerator1)}/{int(denominator1)}" if not np.isnan(numerator1) and not np.isnan(denominator1) else "n=NA"
        n2 = f"n={int(numerator2)}/{int(denominator2)}" if not np.isnan(numerator2) and not np.isnan(denominator2) else "n=NA"

        percentage1_list.append(percentage1)
        percentage2_list.append(percentage2)
        numerator1_list.append(numerator1)
        numerator2_list.append(numerator2)
        denominator1_list.append(denominator1)
        denominator2_list.append(denominator2)
        n1_list.append(n1)
        n2_list.append(n2)

    significant_df['Percentage1'] = percentage1_list
    significant_df['Percentage2'] = percentage2_list
    significant_df['Numerator1'] = numerator1_list
    significant_df['Numerator2'] = numerator2_list
    significant_df['Denominator1'] = denominator1_list
    significant_df['Denominator2'] = denominator2_list
    significant_df['n1'] = n1_list
    significant_df['n2'] = n2_list
    
    return significant_df

################################################################################
#                           MAIN STATS PIPELINE RUN                            #
################################################################################

species_list = ['jejuni', 'coli']

# 1) Compare across STs (within each country & species)
results_across_sequence_types = compare_resistance_across_sequence_types(
    data_filtered, countries_of_interest, species_list, antibiotic_categories
)

# 2) Compare within ST (across countries)
results_within_sequence_type = compare_resistance_within_sequence_type(
    data_filtered, countries_of_interest, species_list, antibiotic_categories
)

# 3) Calculate percentages and raw counts
percentages_dict, raw_counts_dict = calculate_percentages(
    data_filtered, countries_of_interest, species_list, antibiotic_categories
)

# 4) Extract significant results (q < 0.05) for across ST comparisons
significant_results_across_sequence_types = results_across_sequence_types[
    results_across_sequence_types['Corrected P-Value'] < 0.05
]
significant_results_across_sequence_types = add_percentages(
    significant_results_across_sequence_types, percentages_dict, raw_counts_dict
)

# 5) Extract significant results (q < 0.05) for within ST comparisons
significant_results_within_sequence_type = results_within_sequence_type[
    results_within_sequence_type['Corrected P-Value'] < 0.05
]
significant_results_within_sequence_type = add_percentages_within(
    significant_results_within_sequence_type, percentages_dict, raw_counts_dict
)

# 6) Save all results to Excel
excel_file_path_all = os.path.join(output_dir, 'all_p_values_across_sequence_types.xlsx')
excel_file_path_significant = os.path.join(output_dir, 'significant_p_values_across_sequence_types.xlsx')
results_across_sequence_types.to_excel(excel_file_path_all, index=False, sheet_name='AllPValues')
significant_results_across_sequence_types.to_excel(excel_file_path_significant, index=False, sheet_name='SignificantPValues')

excel_file_path_all_within = os.path.join(output_dir, 'all_p_values_within_sequence_types.xlsx')
excel_file_path_significant_within = os.path.join(output_dir, 'significant_p_values_within_sequence_types.xlsx')
results_within_sequence_type.to_excel(excel_file_path_all_within, index=False, sheet_name='AllPValues')
significant_results_within_sequence_type.to_excel(excel_file_path_significant_within, index=False, sheet_name='SignificantPValues')

print(f'All p-values (across ST) -> {excel_file_path_all}')
print(f'Significant p-values (across ST) -> {excel_file_path_significant}')
print(f'All p-values (within ST) -> {excel_file_path_all_within}')
print(f'Significant p-values (within ST) -> {excel_file_path_significant_within}')

# 7) Save CSVs of significant results by species
for specie in species_list:
    significant_species = significant_results_across_sequence_types[significant_results_across_sequence_types['Species'] == specie]
    output_path = os.path.join(output_dir, f'significant_results_{specie}.csv')
    significant_species.to_csv(output_path, index=False)
    print(f'Significant results for {specie} saved to {output_path}')

    # Also save the corresponding percentages with raw counts
    percentages_species = {k: v for k, v in percentages_dict.items() if k[1] == specie}
    raw_counts_species = {k: v for k, v in raw_counts_dict.items() if k[1] == specie}
    
    for (country, _, category_name), percentages_df in percentages_species.items():
        raw_counts = raw_counts_species.get((country, specie, category_name))
        
        # Create a combined DataFrame with percentages and raw counts
        combined_df = percentages_df.copy()
        combined_df = combined_df.reset_index().rename(columns={'index': 'ST'})
        
        # Add n= columns
        for antibiotic in combined_df.columns:
            if antibiotic == 'ST':
                continue
                
            # Add numerator and denominator columns
            numerator_col = f"{antibiotic}_numerator"
            denominator_col = f"{antibiotic}_denominator"
            n_col = f"{antibiotic}_n"
            
            # Extract raw counts
            if raw_counts:
                combined_df[numerator_col] = raw_counts['numerators'][antibiotic].reindex(combined_df['ST']).values
                combined_df[denominator_col] = raw_counts['denominators'][antibiotic].reindex(combined_df['ST']).values
                
                # Create n= format
                combined_df[n_col] = combined_df.apply(
                    lambda row: f"n={int(row[numerator_col])}/{int(row[denominator_col])}" 
                    if not np.isnan(row[numerator_col]) and not np.isnan(row[denominator_col]) 
                    else "n=NA", 
                    axis=1
                )
        
        percentage_output_path = os.path.join(output_dir, f'percentages_{country}_{specie}_{category_name}.csv')
        combined_df.to_csv(percentage_output_path, index=False)
        print(f'Percentages with raw counts for {country}, {specie}, {category_name} saved to {percentage_output_path}')

################################################################################
#                DATA SUBSET FOR PLOTTING (EXCLUDE UNCHARACTERIZED)            #
################################################################################

data_filtered_plots = data_filtered[~data_filtered['ST'].str.contains('Uncharacterized')]

################################################################################
#                GLOBAL DICTIONARIES & LISTS FOR MANHATTAN PLOTS               #
################################################################################

# (A) Renaming dictionary for display of ARG gene names
rename_dict = {
    "AAC(6')-Ie-APH(2'')-Ia_bifunctional_protein": "aac(6')-Ie-aph(2'')-Ia",
    "APH(2'')-If": "aph(2'')-If",
    "APH(2'')-Ig": "aph(2'')-Ig",
    "APH(3')-IIIa": "aph(3')-IIIa",
    "APH(3')-VIIa": "aph(3')-VIIa",
    "aad(6)": "aad(6)",
    "ANT(4')-Ia": "ant(4')-Ia",
    "cmeA": "cmeA", 
    "cmeB": "cmeB",
    "cmeC": "cmeC",
    "cmeR": "cmeR",
    "OXA-184": "oxa-184",
    "OXA-185": "oxa-185",
    "OXA-193": "oxa-193",
    "OXA-447": "oxa-447",
    "OXA-449": "oxa-449",
    "OXA-465": "oxa-465",
    "OXA-576": "oxa-576",
    "OXA-605": "oxa-605",
    "OXA-617": "oxa-617",
    "OXA-622": "oxa-622",
    "OXA-623": "oxa-623",
    "OXA-625": "oxa-625",
    "OXA-631": "oxa-631",
    "OXA-632": "oxa-632",
    "OXA-635": "oxa-635",
    "OXA-638": "oxa-638",
    "T86I": "T86I",
    "lnuP": "lnuP",
    "lnuC": "lnuC",
    "A2075G": "A2075G",
    "SAT-4": "sat-4",
    "Campylobacter_coli_chloramphenicol_acetyltransferase": "Camp cat",
    "Limosilactobacillus_reuteri_cat-TC": "L. reuteri cat-TC",
    "fexA": "fexA",
    "optrA": "optrA",
    "ErmA": "ErmA",
    "ErmB": "ErmB",
    "tet(L)": "tet(L)",
    "tet(O)": "tet(O)",
    "tet(O/M/O)": "tet(O/M/O)",
    "tet(M)": "tet(M)",
    "tet(W)": "tet(W)",
    "vatE": "vatE",
    "vatH": "vatH",
    "vgaD": "vgaD",
    "mel": "mel"
}

# (B) Mapping ARGs -> Drug classes
arg_to_drug_class = {
    "AAC(6')-Ie-APH(2'')-Ia_bifunctional_protein": 'AminoAB',
    "APH(2'')-If": 'AminoAB',
    "APH(2'')-Ig": 'AminoAB',
    "APH(3')-IIIa": 'AminoAB',
    "APH(3')-VIIa": 'AminoAB',
    "ANT(4')-Ia": 'AminoAB',
    "aad(6)": 'AminoAB',

    "ErmA": 'StreptoAB_MacroAB_LincosAB',
    "ErmB": 'StreptoAB_MacroAB_LincosAB',

    "OXA-184": 'Ceph_Penam_Carba',
    "OXA-185": 'Ceph_Penam_Carba',
    "OXA-193": 'Ceph_Penam_Carba',
    "OXA-447": 'Ceph_Penam_Carba',
    "OXA-449": 'Ceph_Penam_Carba',
    "OXA-465": 'Ceph_Penam_Carba',
    "OXA-576": 'Ceph_Penam_Carba',
    "OXA-605": 'Ceph_Penam_Carba',
    "OXA-617": 'Ceph_Penam_Carba',
    "OXA-622": 'Ceph_Penam_Carba',
    "OXA-623": 'Ceph_Penam_Carba',
    "OXA-625": 'Ceph_Penam_Carba',
    "OXA-631": 'Ceph_Penam_Carba',
    "OXA-632": 'Ceph_Penam_Carba',
    "OXA-635": 'Ceph_Penam_Carba',
    "OXA-638": 'Ceph_Penam_Carba',

    "cmeA": 'Ceph_FusiAB_MacroAB_FluoroAB',
    "cmeB": 'Ceph_FusiAB_MacroAB_FluoroAB',
    "cmeC": 'Ceph_FusiAB_MacroAB_FluoroAB',
    "cmeR": 'Ceph_FusiAB_MacroAB_FluoroAB',

    "tet(L)": 'TetraAB',
    "tet(O)": 'TetraAB',
    "tet(O/M/O)": 'TetraAB',
    "tet(W)": 'TetraAB',
    "tet(M)": 'TetraAB',

    "Campylobacter_coli_chloramphenicol_acetyltransferase": 'PhenAB',
    "Limosilactobacillus_reuteri_cat-TC": 'PhenAB',
    "fexA": 'PhenAB',

    "lnuC": 'LincosAB',
    "lnuP": 'LincosAB',

    "SAT-4": 'NucleoAB',

    "optrA": 'PhenAB_OxazoAB',

    "vatE": 'StreptoAB',
    "vatH": 'StreptoAB',
    "vgaD": 'StreptoAB', 
    "mel": 'StreptoAB_MacroAB',

    "T86I": 'FluoroAB',
    "A2075G": 'MacroAB'
}

# (C) List of drug classes in a preferred order
drug_classes_order = [
    'AminoAB',
    'Ceph_Penam_Carba',
    'Ceph_FusiAB_MacroAB_FluoroAB',
    'TetraAB',
    'LincosAB',
    'NucleoAB',
    'StreptoAB',
    'StreptoAB_MacroAB',
    'StreptoAB_MacroAB_LincosAB',
    'PhenAB',
    'PhenAB_OxazoAB',
    'FluoroAB',
    'MacroAB',
]

# Build a dictionary grouping ARGs by drug class:
drug_class_to_args = {}
for arg in antibiotic_categories['ARGs']:
    dc = arg_to_drug_class.get(arg, 'Unknown')
    drug_class_to_args.setdefault(dc, []).append(arg)

# Create an ordered list of ARGs grouped by these classes
args_grouped = []
drug_class_positions = []
current_pos = 0
for dc in drug_classes_order:
    args_in_class = drug_class_to_args.get(dc, [])
    args_grouped.extend(args_in_class)
    if len(args_in_class) > 0:
        # store the start/end indices for that drug class
        drug_class_positions.append((dc, current_pos, current_pos + len(args_in_class)))
        current_pos += len(args_in_class)

n_genes = len(args_grouped)

################################################################################
#                          PLOTTING FUNCTIONS                                  #
################################################################################

import matplotlib.cm as cm
from matplotlib.colors import to_hex
from matplotlib.patches import Patch

# Define pastel colors for drug-class backgrounds
drug_class_pastels = {
    'AminoAB': '#FFE5D9',
    'Ceph_Penam_Carba': '#F0F0F0',
    'Ceph_FusiAB_MacroAB_FluoroAB': '#E0F7FA',
    'TetraAB': '#FFE5E5',
    'LincosAB': '#FFFDE7',
    'NucleoAB': '#FCE4EC',
    'StreptoAB': '#F3E5F5',
    'StreptoAB_MacroAB': '#EDE7F6',
    'StreptoAB_MacroAB_LincosAB': '#FFFEF0',
    'PhenAB': '#E8F5E9',
    'PhenAB_OxazoAB': '#F5E5F5',
    'FluoroAB': '#E3F2FD',
    'MacroAB': '#F5E5DC',
}

def create_legend_figure(labels, handles, title, output_path):
    fig_legend = plt.figure(figsize=(6, 4))
    plt.legend(handles, labels, loc='center', frameon=False)
    plt.title(title)
    plt.axis('off')
    plt.savefig(output_path, dpi=800, bbox_inches='tight')
    plt.close(fig_legend)

def create_stacked_bar_charts(
    df, country, species_list, categories, custom_colors, output_dir,
    figsize=(14, 10),
    bar_font_size=12,
    tick_label_size=20,
    legend_font_size=20,
    y_axis_font_size=18,
    subplot_label_size=24,
    min_pct_to_show=5.0,
    min_samples=5,
    show_percentages=True,
    round_percentages=True
):
    """
    Enhanced horizontal stacked bar charts with publication-quality styling
    """
    dfp = df[~df['ST'].str.contains('Uncharacterized')]

    for sp in species_list:
        sts_all = dfp[(dfp['Country']==country)&(dfp['Species_Type']==sp)]['ST']
        sts = [st for st, cnt in sts_all.value_counts().items() if cnt >= min_samples]
        
        if not sts:
            print(f"No STs with â‰¥{min_samples} samples for {country}-{sp}. Skipping bar chart.")
            continue
            
        fig, ax = plt.subplots(figsize=figsize)
        
        dcats = categories['Drug_Classes']
        pct = pd.DataFrame(0, index=sts, columns=dcats)

        sizes = {}
        for st in sts:
            sub = dfp[(dfp['Country']==country)&(dfp['Species_Type']==sp)&(dfp['ST']==st)]
            sizes[st] = len(sub)
            for dc in dcats:
                pct.loc[st, dc] = sub[dc].sum() if dc in sub.columns else 0

        totals = pct.sum(axis=1)
        pct = pct.div(totals, axis=0).fillna(0)*100

        bottoms = np.zeros(len(pct))
        for dc in dcats:
            vals = pct[dc].values
            col = custom_colors['Drug_Classes'].get(dc, '#CCCCCC')
            bars = ax.barh(np.arange(len(pct)), vals, left=bottoms, color=col, label=dc)
            
            if show_percentages:
                for i,(v,b) in enumerate(zip(vals,bottoms)):
                    if v >= min_pct_to_show:
                        if round_percentages:
                            display_val = f"{round(v)}%"
                        else:
                            display_val = f"{v:.1f}%"
                        ax.text(b+v/2, i, display_val, 
                                ha='center', va='center', 
                                fontsize=bar_font_size, 
                                fontweight='bold')
            bottoms += vals

        ax.set_yticks(np.arange(len(pct)))
        ax.set_yticklabels([f"{st}\n(n={sizes[st]})" for st in pct.index], 
                           fontsize=y_axis_font_size, fontweight='bold')
        ax.set_xlabel('Percentage', fontsize=y_axis_font_size, fontweight='bold')
        ax.set_title(f"{country} â€” {sp} Drug Class Distribution (â‰¥{min_samples} samples)", 
                     fontsize=subplot_label_size, fontweight='bold', pad=20)
        ax.tick_params(axis='x', labelsize=tick_label_size)
        ax.tick_params(axis='y', labelsize=y_axis_font_size, pad=15)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels,
                  loc='center left',
                  bbox_to_anchor=(1.02, 0.5),
                  fontsize=legend_font_size,
                  frameon=True, fancybox=True, shadow=True)

        plt.tight_layout()
        out = os.path.join(output_dir, f"barchart_{country}_{sp}.png")
        plt.savefig(out, dpi=1200, bbox_inches='tight', facecolor='white')
        plt.close(fig)

def create_manhattan_plot(
    df, results_df, country, species, categories, custom_colors, output_dir,
    figure_size=(18, 12),
    marker_size=400,
    legend_markerscale=2.5,
    x_label_font_size=26,
    y_label_font_size=20,
    tick_label_size=20,
    legend_font_size=24,
    dashed_line_color='grey',
    min_samples=5,
    show_x_title=True,
    subplot_label_size=32
):
    """
    Enhanced ST-level Manhattan plot with publication-quality styling
    """
    dfp = df[~df['ST'].str.contains('Uncharacterized')]

    # FILTER FOR MINIMUM SAMPLES
    sts_all = dfp[(dfp['Country']==country)&(dfp['Species_Type']==species)]['ST']
    sts = [st for st, cnt in sts_all.value_counts().items() if cnt >= min_samples]
    
    if not sts:
        print(f"No STs with â‰¥{min_samples} samples for {country}-{species}. Skipping Manhattan plot.")
        return
    
    cmap = cm.get_cmap('tab20', len(sts))
    st_colors = {st: to_hex(cmap(i)) for i,st in enumerate(sts)}

    fig, ax = plt.subplots(figsize=figure_size)
    ax.set_facecolor('whitesmoke')

    # Add pastel drug-class background bands
    for dc, start, end in drug_class_positions:
        col = drug_class_pastels.get(dc, '#F5F5F5')
        ax.axvspan(start-0.5, end-0.5, color=col, alpha=0.6, linewidth=0, zorder=0)

    # Scatter plot for each ST
    for st in sts:
        sub = dfp[(dfp['Country']==country)&(dfp['Species_Type']==species)&(dfp['ST']==st)]
        presence = [
            (sub[gene].sum()/sub[gene].count()*100) if gene in sub.columns and sub[gene].count()>0 else 0
            for gene in args_grouped
        ]
        ax.scatter(
            np.arange(n_genes),
            presence,
            color=st_colors[st],
            s=marker_size,
            label=f"{st} (n={len(sub)})",
            edgecolors='black',
            linewidth=1.0,
            zorder=5
        )

    # Add significance stars
    for i,g in enumerate(args_grouped):
        rows = results_df[
            (results_df['Country']==country)&
            (results_df['Species']==species)&
            (results_df['Antibiotic']==g)
        ]
        if any((rows['P-Value']<0.05)&(rows['Corrected P-Value']<0.05)):
            ax.scatter(i, 106, marker='*', s=marker_size*0.8, color='black', zorder=10)

    # Styling
    ax.set_xticks(np.arange(n_genes))
    ax.set_xticklabels([rename_dict.get(g,g) for g in args_grouped],
                       rotation=90, fontsize=x_label_font_size, fontweight='bold')
    
    # Add vertical lines between drug classes
    for dc, start, end in drug_class_positions:
        if start > 0:
            ax.axvline(x=start-0.5, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    
    ax.set_yticks(np.arange(0,101,10))
    ax.tick_params(axis='y', labelsize=tick_label_size)
    ax.set_ylim(0,110)
    ax.axhline(100, linestyle='--', color=dashed_line_color, linewidth=0.5)

    if show_x_title:
        ax.set_xlabel('AMR Genes grouped by Drug Class', fontsize=x_label_font_size, fontweight='bold')
    ax.set_ylabel('Percentage Presence', fontsize=y_label_font_size, fontweight='bold')

    # Title
    ax.set_title(f'{country} - {species} (â‰¥{min_samples} samples)', 
                 fontsize=subplot_label_size, fontweight='bold', pad=20)

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 0:
        ax.legend(handles, labels,
                  loc='center left',
                  bbox_to_anchor=(1.02, 0.5),
                  fontsize=legend_font_size,
                  markerscale=legend_markerscale,
                  frameon=True, fancybox=True, shadow=True)

    plt.tight_layout()
    out = os.path.join(output_dir, f"manhattan_{country}_{species}.png")
    plt.savefig(out, dpi=1200, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def create_combined_manhattan_plots(
    df, results_df, species, categories, custom_colors, output_dir,
    figsize=(50, 32),
    x_label_font_size=26,
    y_label_font_size=20,
    tick_label_size=20,
    legend_font_size=24,
    subplot_label_size=32,
    marker_size=400,
    legend_markerscale=2.5,
    min_samples=5,
    show_x_title=False
):
    """
    Enhanced 2Ã—2 grid of Manhattan plots with publication-quality styling
    """
    fig = plt.figure(figsize=figsize)
    grid = plt.GridSpec(2, 2, figure=fig, hspace=0.5, wspace=0.4)

    positions = np.arange(n_genes)
    country_positions = {
        'Bangladesh': (0, 0, 'A'),
        'India':      (0, 1, 'B'),
        'Vietnam':    (1, 0, 'C'),
        'Sri Lanka':  (1, 1, 'D')
    }

    for country, (r, c, label) in country_positions.items():
        ax = fig.add_subplot(grid[r, c])
        ax.set_facecolor('whitesmoke')

        # FILTER FOR MINIMUM SAMPLES
        sts_all = df[(df['Country']==country)&(df['Species_Type']==species)]['ST']
        st_list = [st for st, cnt in sts_all.value_counts().items() if cnt >= min_samples]

        if not st_list:
            ax.text(0.5, 0.5, f'No STs with â‰¥{min_samples} samples', 
                    transform=ax.transAxes, ha='center', va='center', fontsize=20,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"))
            ax.set_xticks([])
            ax.set_yticks([])
            ax.text(0.5, 1.12, f'{label}) {country}', transform=ax.transAxes,
                    ha='center', fontsize=subplot_label_size, fontweight='bold')
            continue

        # Add pastel drug-class backgrounds
        for dc, start, end in drug_class_positions:
            col = drug_class_pastels.get(dc, '#F5F5F5')
            ax.axvspan(start-0.5, end-0.5, color=col, alpha=0.6, linewidth=0, zorder=0)

        # Build % table
        pct = pd.DataFrame(0, index=args_grouped, columns=st_list)
        for st in st_list:
            sub = df[
                (df['Country']==country) &
                (df['Species_Type']==species) &
                (df['ST']==st)
            ]
            for gene in args_grouped:
                if gene in sub.columns and sub[gene].count()>0:
                    pct.at[gene, st] = sub[gene].sum()/sub[gene].count()*100

        # Plot each ST
        cmap = cm.get_cmap('tab20', len(st_list))
        st_colors = {st: to_hex(cmap(i)) for i, st in enumerate(st_list)}
        
        for st in st_list:
            sub = df[
                (df['Country']==country) &
                (df['Species_Type']==species) &
                (df['ST']==st)
            ]
            vals = pct[st].values
            if not np.any(vals > 0):
                continue
            ax.scatter(
                positions,
                vals,
                color=st_colors[st],
                s=marker_size,
                label=f"{st} (n={len(sub)})",
                edgecolors='black',
                linewidth=1.0,
                zorder=5
            )

        # Significance stars
        for i, gene in enumerate(args_grouped):
            rows = results_df[
                (results_df['Country']==country) &
                (results_df['Species']==species) &
                (results_df['Antibiotic']==gene)
            ]
            if any((rows['P-Value']<0.05)&(rows['Corrected P-Value']<0.05)):
                ax.scatter(i, 106, color='black', marker='*', s=marker_size*0.8, zorder=10)

        # Styling
        ax.set_xticks(positions)
        ax.set_xticklabels([rename_dict.get(g,g) for g in args_grouped],
                           rotation=90, fontsize=x_label_font_size, fontweight='bold')
        ax.set_yticks(np.arange(0,101,10))
        ax.tick_params(axis='y', labelsize=tick_label_size)
        ax.set_ylim(0,110)
        ax.axhline(100, linestyle='--', color='grey', linewidth=0.5)
        
        # Vertical lines between drug classes
        for dc, start, end in drug_class_positions:
            if start > 0:
                ax.axvline(x=start-0.5, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)

        if r == 1 and show_x_title:
            ax.set_xlabel('AMR Genes grouped by Drug Class', fontsize=x_label_font_size, fontweight='bold')
        if c == 0:
            ax.set_ylabel('Percentage Presence', fontsize=y_label_font_size, fontweight='bold')
        
        ax.text(0.5, 1.12, f'{label}) {country}', transform=ax.transAxes,
                ha='center', fontsize=subplot_label_size, fontweight='bold')

        # Individual legend for each subplot - positioned to the right
        handles, labels = ax.get_legend_handles_labels()
        if len(handles) > 0:
            ax.legend(handles, labels,
                      loc='center left', bbox_to_anchor=(1.02, 0.5),
                      fontsize=legend_font_size,
                      markerscale=legend_markerscale,
                      frameon=True, fancybox=True, shadow=True)

    plt.suptitle(f"Manhattan Plots â€” {species} (â‰¥{min_samples} samples)", 
                 fontsize=subplot_label_size+4, y=0.95, fontweight='bold')
    
    plt.subplots_adjust(top=0.90, bottom=0.05, left=0.06, right=0.98, hspace=0.35, wspace=0.4)
    out = os.path.join(output_dir, f"combined_manhattan_{species}.png")
    plt.savefig(out, dpi=1200, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def create_combined_stacked_bar_charts(
    df, species, categories, custom_colors, output_dir,
    figsize=(28, 28),
    tick_label_size=18,
    legend_font_size=20,
    y_axis_font_size=18,
    subplot_label_size=28,
    bar_font_size=12,
    min_samples=5,
    show_percentages=True,
    min_pct_to_show=5.0,
    round_percentages=True
):
    """
    Enhanced 2Ã—2 grid of stacked bar charts with publication-quality styling
    """
    fig = plt.figure(figsize=figsize)
    grid = plt.GridSpec(2, 2, figure=fig, hspace=0.5, wspace=0.5)

    country_positions = {
        'Bangladesh': (0, 0, 'E'),
        'India':      (0, 1, 'F'),
        'Vietnam':    (1, 0, 'G'),
        'Sri Lanka':  (1, 1, 'H')
    }

    for country, (r, c, label) in country_positions.items():
        ax = fig.add_subplot(grid[r, c])
        
        sts_all = df[(df['Country']==country)&(df['Species_Type']==species)]['ST']
        st_list = [st for st, cnt in sts_all.value_counts().items() if cnt >= min_samples]

        if not st_list:
            ax.text(0.5, 0.5, f'No STs with â‰¥{min_samples} samples', 
                    transform=ax.transAxes, ha='center', va='center', fontsize=18,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"))
            ax.set_xticks([])
            ax.set_yticks([])
            ax.text(-0.1, 1.15, f'{label}) {country}', transform=ax.transAxes,
                    fontsize=subplot_label_size, fontweight='bold')
            continue

        dcats = categories['Drug_Classes']
        pct = pd.DataFrame(0, index=st_list, columns=dcats)

        for st in st_list:
            sub = df[
                (df['Country']==country)&
                (df['Species_Type']==species)&
                (df['ST']==st)
            ]
            for dc in dcats:
                pct.loc[st, dc] = sub[dc].sum() if dc in sub.columns else 0

        pct = pct.div(pct.sum(axis=1), axis=0).fillna(0)*100

        bottom = np.zeros(len(pct))
        for dc in dcats:
            vals = pct[dc].values
            col = custom_colors['Drug_Classes'].get(dc, '#CCCCCC')
            bars = ax.barh(np.arange(len(pct)), vals, left=bottom, color=col, label=dc)
            
            if show_percentages:
                for i,(v,b) in enumerate(zip(vals,bottom)):
                    if v >= min_pct_to_show:
                        if round_percentages:
                            display_val = f"{round(v)}%"
                        else:
                            display_val = f"{v:.1f}%"
                        ax.text(b+v/2, i, display_val,
                                ha='center', va='center',
                                fontsize=bar_font_size,
                                fontweight='bold')
            bottom += vals

        ax.set_yticks(np.arange(len(pct)))
        ax.set_yticklabels([f"{st}\n(n={int(df[(df['Country']==country)&(df['Species_Type']==species)&(df['ST']==st)].shape[0])})"
                             for st in pct.index],
                            fontsize=y_axis_font_size, fontweight='bold')
        ax.set_xlabel('Percentage', fontsize=tick_label_size, fontweight='bold')
        ax.tick_params(axis='x', labelsize=tick_label_size)
        ax.tick_params(axis='y', labelsize=y_axis_font_size, pad=15)
        
        ax.text(-0.1, 1.15, f'{label}) {country}', transform=ax.transAxes,
                fontsize=subplot_label_size, fontweight='bold')

    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 0:
        ncols = min(6, len(handles))
        fig.legend(handles, labels,
                   loc='center',
                   bbox_to_anchor=(0.5, 0.02),
                   ncol=ncols,
                   fontsize=legend_font_size,
                   frameon=True,
                   fancybox=True,
                   shadow=True)

    plt.subplots_adjust(top=0.95, bottom=0.18, left=0.12, right=0.95)

    out = os.path.join(output_dir, f"combined_barcharts_{species}.png")
    plt.savefig(out, dpi=1200, bbox_inches='tight', facecolor='white')
    plt.close(fig)

################################################################################
#                          CREATE & SAVE THE PLOTS                             #
################################################################################

# 1) Stacked bar charts per country & species (WITH FILTERING)
for country in countries_of_interest:
    create_stacked_bar_charts(
        data_filtered_plots, country, species_list,
        antibiotic_categories, custom_colors, output_dir,
        min_samples=5, show_percentages=True, round_percentages=True
    )

# 2) Single Manhattan per country & species (WITH FILTERING)
for country in countries_of_interest:
    for sp in species_list:
        create_manhattan_plot(
            data_filtered_plots, results_across_sequence_types,
            country, sp, antibiotic_categories, custom_colors, output_dir,
            marker_size=400, legend_markerscale=2.5,
            min_samples=5
        )

# 3) Combined 2Ã—2 Manhattan (WITH FILTERING)
for sp in species_list:
    create_combined_manhattan_plots(
        data_filtered_plots, results_across_sequence_types,
        sp, antibiotic_categories, custom_colors, output_dir,
        min_samples=5, show_x_title=False
    )

# 4) Combined 2Ã—2 Bar charts (WITH FILTERING)
for sp in species_list:
    create_combined_stacked_bar_charts(
        data_filtered_plots, sp, antibiotic_categories, custom_colors, output_dir,
        min_samples=5, show_percentages=True, min_pct_to_show=5.0
    )

print("All plots generated and saved to:", output_dir)
print(f"All plots now consistently filter for â‰¥5 samples, matching the statistical analysis.")

#st code last