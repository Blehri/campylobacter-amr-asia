"""
ARG Resistance Visualizations
=============================
Stacked bar charts showing ARG prevalence by country and species (Figures 3-4)

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
# updated v3
#updated country 08-07-25
# Version 502 - Complete Adjusted Code country and species SIZE ADJUSTED
import os
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# 1. Define input CSV path
# ------------------------------------------------------------------------------
file_path = r"DATA_FILE  # Original path: C:\Users\BLehr\...\My_merged_output_V5.csv"

# ------------------------------------------------------------------------------
# 2. Load data and clean columns
# ------------------------------------------------------------------------------
data_df = pd.read_csv(file_path, encoding='ISO-8859-1')
data_df.columns = data_df.columns.str.replace(r'\s+', ' ', regex=True)

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


# ------------------------------------------------------------------------------
# 3. Define interest groups & rename mappings
# ------------------------------------------------------------------------------
countries_of_interest = ['Bangladesh', 'India', 'Vietnam', 'Sri Lanka']
species_of_interest = ['Campylobacter jejuni', 'Campylobacter coli']

rename_species = {
    'Campylobacter jejuni': 'C. jejuni',
    'Campylobacter coli': 'C. coli'
}

# Antibiotic categories
antibiotic_categories = {
    'ARGs': [
        "AAC(6')-Ie-APH(2'')-Ia bifunctional protein", "APH(2'')-If", "APH(2'')-Ig", "APH(3')-IIIa",
        "APH(3')-VIIa", "ANT(4')-Ia", "aad(6)", "ErmA", "ErmB", "OXA-184", "OXA-185", "OXA-193", "OXA-447",
        "OXA-449", "OXA-465", "OXA-576", "OXA-605", "OXA-617", "OXA-622", "OXA-623", "OXA-625",
        "OXA-631", "OXA-632", "OXA-635", "OXA-638", "cmeA", "cmeB", "cmeC", "cmeR", "tet(L)", "tet(O)", 
        "tet(O/M/O)", "tet(W)", "tet(M)",
        "Campylobacter coli chloramphenicol acetyltransferase", "Limosilactobacillus reuteri cat-TC", 
        "fexA", "lnuC", "lnuP", "SAT-4", "optrA", "vatE", "vatH", "mel", "T86I", "A2075G"
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
        'AminoAB', 'StreptoAB MacroAB LincosAB', 'Ceph Penam Carba', 'Ceph FusiAB MacroAB FluoroAB', 
        'TetraAB', 'PhenAB', 'LincosAB', 'NucleoAB', 'StreptoAB', 'StreptoAB MacroAB', 
        'PhenAB OxazoAB', 'FluoroAB', 'MacroAB', 
    ]
}

# Custom colors for each category (as in your code)
custom_colors = {
    'Drug Classes': {
        'AminoAB': '#FFA07A',
        'Ceph Penam Carba': '#DCDCDC',
        'TetraAB': '#FF0000',
        'NucleoAB': '#FFC0CB',
        'PhenAB OxazoAB': '#800080',
        'FluoroAB': '#1E90FF',
        'MacroAB': '#6e2c00',
        'Ceph FusiAB MacroAB FluoroAB': '#00CED1',
        'LincosAB': '#FFFF00',
        'PhenAB': '#008000',
        'StreptoAB MacroAB LincosAB': '#FFFFE0',
        'StreptoAB': '#8A2BE2',
        'StreptoAB MacroAB': '#9370DB'
    },
    'Antibiotic Resistance Gene Classes': {
        "aminoglycoside bifunctional resistance protein": '#FF8C00',
        "APH(2'')": '#FFA07A',
        "APH(3')": '#FF7F50',
        "ANT(4')": '#FF69B4',
        "ANT(6)": '#FF4500',
        "resistance-nodulation-cell division (RND) antibiotic efflux pump": '#00CED1',
        "chloramphenicol acetyltransferase (CAT)": '#008000',
        "lincosamide nucleotidyltransferase (LNU)": '#FFFF00',
        "Erm 23S ribosomal RNA methyltransferase": '#FFFFE0',
        "OXA beta-lactamase": '#DCDCDC',
        "tetracycline-resistant ribosomal protection protein": '#DC143C',
        "major facilitator superfamily (MFS) antibiotic efflux pump": '#00FF7F',
        "streptothricin acetyltransferase (SAT)": '#FFC0CB',
        "Miscellaneous ABC-F subfamily ATP-binding cassette ribosomal protection proteins": '#800080',
        "fluoroquinolone resistant gyrA": '#1E90FF',
        "23S rRNA with mutation conferring resistance to macrolide antibiotics": '#6e2c00',
        "streptogramin vat acetyltransferase": '#8A2BE2',
        "msr-type ABC-F protein": '#9370DB'
    },
    'ARGs': {
        "APH(2'')-If": '#FFA07A',
        "AAC(6')-Ie-APH(2'')-Ia bifunctional protein": '#FF8C00',
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
        "Campylobacter coli chloramphenicol acetyltransferase": '#008000',
        "Limosilactobacillus reuteri cat-TC": '#00FA9A',
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

# ------------------------------------------------------------------------------
# 4. Statistical comparison functions
# ------------------------------------------------------------------------------
def compare_resistance_across_countries(df, countries, species, categories):
    all_results = []
    for category_name, antibiotics in categories.items():
        for antibiotic in antibiotics:
            for i, country1 in enumerate(countries):
                for country2 in countries[i + 1:]:
                    for specie in species:
                        data1 = df[(df['Country'] == country1) & (df['Species_Name'] == specie)][antibiotic]
                        data2 = df[(df['Country'] == country2) & (df['Species_Name'] == specie)][antibiotic]

                        if data1.empty or data2.empty:
                            all_results.append({
                                'Category': category_name,
                                'Species': specie,
                                'Antibiotic': antibiotic,
                                'Country1': country1,
                                'Country2': country2,
                                'P-Value': np.nan,
                                'Test Used': 'Missing Data'
                            })
                            continue

                        table = [
                            [data1.sum(), max(0, data1.count() - data1.sum())],
                            [data2.sum(), max(0, data2.count() - data2.sum())]
                        ]

                        # Fisher for ARGs, Chi-square for classes (fallback Fisher)
                        if category_name == 'ARGs':
                            try:
                                _, p_value = fisher_exact(table, alternative='two-sided')
                                test_used = 'Fisher'
                            except Exception as e:
                                p_value = np.nan
                                test_used = 'Error in Fisher: ' + str(e)
                        else:
                            table_array = np.array(table)
                            if all(x > 0 for row_ in table_array for x in row_):
                                try:
                                    _, p_value, _, _ = chi2_contingency(table_array)
                                    test_used = 'Chi-squared'
                                except Exception as e:
                                    p_value = np.nan
                                    test_used = 'Error in Chi-squared: ' + str(e)
                            else:
                                try:
                                    _, p_value = fisher_exact(table_array, alternative='two-sided')
                                    test_used = 'Fisher (Fallback)'
                                except Exception as e:
                                    p_value = np.nan
                                    test_used = 'Error in Fisher (Fallback): ' + str(e)

                        all_results.append({
                            'Category': category_name,
                            'Species': specie,
                            'Antibiotic': antibiotic,
                            'Country1': country1,
                            'Country2': country2,
                            'P-Value': p_value,
                            'Test Used': test_used
                        })

    results_df = pd.DataFrame(all_results)
    if not results_df.empty:
        valid_p_values = results_df['P-Value'].notna()
        corrected_p_values = np.full(results_df.shape[0], np.nan)
        _, corrected_p_values[valid_p_values], _, _ = multipletests(
            results_df.loc[valid_p_values, 'P-Value'], method='fdr_bh'
        )
        results_df['Corrected P-Value'] = corrected_p_values
    else:
        results_df['Corrected P-Value'] = []

    return results_df


def compare_resistance_within_countries(df, countries, species, categories):
    all_results = []
    for category_name, antibiotics in categories.items():
        for antibiotic in antibiotics:
            for country in countries:
                for i, specie1 in enumerate(species):
                    for specie2 in species[i + 1:]:
                        data1 = df[(df['Country'] == country) & (df['Species_Name'] == specie1)][antibiotic]
                        data2 = df[(df['Country'] == country) & (df['Species_Name'] == specie2)][antibiotic]

                        if data1.empty or data2.empty:
                            all_results.append({
                                'Category': category_name,
                                'Country': country,
                                'Antibiotic': antibiotic,
                                'Species1': specie1,
                                'Species2': specie2,
                                'P-Value': np.nan,
                                'Test Used': 'Missing Data'
                            })
                            continue

                        table = [
                            [data1.sum(), max(0, data1.count() - data1.sum())],
                            [data2.sum(), max(0, data2.count() - data2.sum())]
                        ]

                        if category_name == 'ARGs':
                            try:
                                _, p_value = fisher_exact(table, alternative='two-sided')
                                test_used = 'Fisher'
                            except Exception as e:
                                p_value = np.nan
                                test_used = 'Error in Fisher: ' + str(e)
                        else:
                            table_array = np.array(table)
                            if all(x > 0 for row_ in table_array for x in row_):
                                try:
                                    _, p_value, _, _ = chi2_contingency(table_array)
                                    test_used = 'Chi-squared'
                                except Exception as e:
                                    p_value = np.nan
                                    test_used = 'Error in Chi-squared: ' + str(e)
                            else:
                                try:
                                    _, p_value = fisher_exact(table_array, alternative='two-sided')
                                    test_used = 'Fisher (Fallback)'
                                except Exception as e:
                                    p_value = np.nan
                                    test_used = 'Error in Fisher (Fallback): ' + str(e)

                        all_results.append({
                            'Category': category_name,
                            'Country': country,
                            'Antibiotic': antibiotic,
                            'Species1': specie1,
                            'Species2': specie2,
                            'P-Value': p_value,
                            'Test Used': test_used
                        })

    results_df = pd.DataFrame(all_results)
    if not results_df.empty:
        valid_p_values = results_df['P-Value'].notna()
        corrected_p_values = np.full(results_df.shape[0], np.nan)
        _, corrected_p_values[valid_p_values], _, _ = multipletests(
            results_df.loc[valid_p_values, 'P-Value'], method='fdr_bh'
        )
        results_df['Corrected P-Value'] = corrected_p_values
    else:
        results_df['Corrected P-Value'] = []

    return results_df


# ------------------------------------------------------------------------------
# 5. Run the comparisons
# ------------------------------------------------------------------------------
results_across_countries = compare_resistance_across_countries(
    data_df, countries_of_interest, species_of_interest, antibiotic_categories
)
results_within_countries = compare_resistance_within_countries(
    data_df, countries_of_interest, species_of_interest, antibiotic_categories
)

significant_results_across_countries = results_across_countries[results_across_countries['Corrected P-Value'] < 0.05]
significant_results_within_countries = results_within_countries[results_within_countries['Corrected P-Value'] < 0.05]

# ------------------------------------------------------------------------------
# 6. Calculate percentages and raw counts
# ------------------------------------------------------------------------------
percentages_dict = {}
raw_counts_dict = {}

for country in countries_of_interest:
    for specie in species_of_interest:
        for category_name, antibiotics in antibiotic_categories.items():
            subset_data = data_df[
                (data_df['Country'] == country) & 
                (data_df['Species_Name'] == specie)
            ]
            if subset_data.empty:
                continue

            for antibiotic in antibiotics:
                if antibiotic not in subset_data.columns:
                    continue

                if category_name == 'ARGs':
                    presence = (subset_data[antibiotic] > 0).sum()
                    total = subset_data[antibiotic].count()
                    percentage = (presence / total) * 100 if total > 0 else np.nan
                    raw_counts_dict[(country, specie, category_name, antibiotic)] = (presence, total)
                else:
                    total_usage = subset_data[antibiotics].sum().sum()
                    antibiotic_usage = subset_data[antibiotic].sum()
                    percentage = (antibiotic_usage / total_usage) * 100 if total_usage > 0 else np.nan
                    raw_counts_dict[(country, specie, category_name, antibiotic)] = (antibiotic_usage, total_usage)

                percentages_dict[(country, specie, category_name, antibiotic)] = percentage

# ------------------------------------------------------------------------------
# 7. Helper to attach percentages/raw counts to the significant results
# ------------------------------------------------------------------------------
def add_percentages(significant_df, percentages_d, raw_counts_d, cross_country=True):
    percentage1_list = []
    percentage2_list = []
    counts1_list = []
    counts2_list = []

    for _, row in significant_df.iterrows():
        category = row['Category']
        antibiotic = row['Antibiotic']

        if cross_country:
            country1 = row['Country1']
            country2 = row['Country2']
            species = row['Species']

            pct1 = percentages_d.get((country1, species, category, antibiotic), np.nan)
            pct2 = percentages_d.get((country2, species, category, antibiotic), np.nan)

            cnt1 = raw_counts_d.get((country1, species, category, antibiotic), (np.nan, np.nan))
            cnt2 = raw_counts_d.get((country2, species, category, antibiotic), (np.nan, np.nan))

        else:
            country = row['Country']
            specie1 = row['Species1']
            specie2 = row['Species2']

            pct1 = percentages_d.get((country, specie1, category, antibiotic), np.nan)
            pct2 = percentages_d.get((country, specie2, category, antibiotic), np.nan)

            cnt1 = raw_counts_d.get((country, specie1, category, antibiotic), (np.nan, np.nan))
            cnt2 = raw_counts_d.get((country, specie2, category, antibiotic), (np.nan, np.nan))

        percentage1_list.append(pct1)
        percentage2_list.append(pct2)
        counts1_list.append(cnt1)
        counts2_list.append(cnt2)

    significant_df['Percentage1'] = percentage1_list
    significant_df['Percentage2'] = percentage2_list

    # Unpack raw counts
    significant_df['Count1_Numerator'] = [c[0] for c in counts1_list]
    significant_df['Count1_Denominator'] = [c[1] for c in counts1_list]
    significant_df['Count2_Numerator'] = [c[0] for c in counts2_list]
    significant_df['Count2_Denominator'] = [c[1] for c in counts2_list]

    # Add "n=..." columns
    significant_df['n1'] = significant_df.apply(
        lambda r: f"n={int(r['Count1_Numerator'])}/{int(r['Count1_Denominator'])}" 
                  if not np.isnan(r['Count1_Numerator']) and not np.isnan(r['Count1_Denominator']) 
                  else "n=NA", axis=1
    )
    significant_df['n2'] = significant_df.apply(
        lambda r: f"n={int(r['Count2_Numerator'])}/{int(r['Count2_Denominator'])}" 
                  if not np.isnan(r['Count2_Numerator']) and not np.isnan(r['Count2_Denominator']) 
                  else "n=NA", axis=1
    )

    return significant_df

significant_results_across_countries = add_percentages(
    significant_results_across_countries, percentages_dict, raw_counts_dict, cross_country=True
)
significant_results_within_countries = add_percentages(
    significant_results_within_countries, percentages_dict, raw_counts_dict, cross_country=False
)

significant_results_across_countries['Comparison'] = 'Cross-Country'
significant_results_within_countries['Comparison'] = 'Within-Country'

combined_significant_results = pd.concat(
    [significant_results_across_countries, significant_results_within_countries],
    axis=0, ignore_index=True
)

# ------------------------------------------------------------------------------
# 8. Save ALL results (significant and non-significant) + significant only
# ------------------------------------------------------------------------------
output_dir = r"C:\Users\BLehr\OneDrive - London School of Hygiene and Tropical Medicine\Documents\OHPH_documents\Hub_summmary\meta_data_extra\group_useful\selected_contig\other_outputs\merge\outputs\Country\updated_2"
os.makedirs(output_dir, exist_ok=True)

# Save ALL p-values (both significant and non-significant) to Excel
excel_file_path_all_cross_country = os.path.join(output_dir, 'all_p_values_cross_country.xlsx')
excel_file_path_all_within_country = os.path.join(output_dir, 'all_p_values_within_country.xlsx')
excel_file_path_significant_cross_country = os.path.join(output_dir, 'significant_p_values_cross_country.xlsx')
excel_file_path_significant_within_country = os.path.join(output_dir, 'significant_p_values_within_country.xlsx')

results_across_countries.to_excel(excel_file_path_all_cross_country, index=False, sheet_name='AllCrossCountry')
results_within_countries.to_excel(excel_file_path_all_within_country, index=False, sheet_name='AllWithinCountry')
significant_results_across_countries.to_excel(excel_file_path_significant_cross_country, index=False, sheet_name='SignificantCrossCountry')
significant_results_within_countries.to_excel(excel_file_path_significant_within_country, index=False, sheet_name='SignificantWithinCountry')

print(f'All cross-country p-values saved to {excel_file_path_all_cross_country}')
print(f'All within-country p-values saved to {excel_file_path_all_within_country}')
print(f'Significant cross-country p-values saved to {excel_file_path_significant_cross_country}')
print(f'Significant within-country p-values saved to {excel_file_path_significant_within_country}')

# Save the combined significant results to CSV
output_path_combined = os.path.join(output_dir, 'significant_results_combined.csv')
combined_significant_results.to_csv(output_path_combined, index=False)
print(f"Significant results with percentages saved to {output_path_combined}")

# Also save per species
for sp in species_of_interest:
    output_path_species = os.path.join(output_dir, f'significant_results_{sp}.csv')
    subset = combined_significant_results[
        (combined_significant_results['Species'] == sp) |
        (combined_significant_results['Species1'] == sp) |
        (combined_significant_results['Species2'] == sp)
    ]
    subset.to_csv(output_path_species, index=False)
    print(f"Significant results for {sp} saved to {output_path_species}")

# ------------------------------------------------------------------------------
# 9. Create a DataFrame from percentages_dict (including raw counts)
# ------------------------------------------------------------------------------
percentages_list = []
for k, v in percentages_dict.items():
    country, species, category, antibiotic = k
    counts = raw_counts_dict.get(k, (np.nan, np.nan))
    percentages_list.append({
        'Country': country,
        'Species': species,
        'Category': category,
        'Antibiotic': antibiotic,
        'Percentage': v,
        'Numerator': counts[0],
        'Denominator': counts[1],
        'n': (f"n={int(counts[0])}/{int(counts[1])}" 
              if not np.isnan(counts[0]) and not np.isnan(counts[1]) else "n=NA")
    })

percentages_df = pd.DataFrame(percentages_list)

# Save percentages to CSV files for each species
for sp in species_of_interest:
    percentage_output_path = os.path.join(output_dir, f'percentages_{sp}.csv')
    subset = percentages_df[percentages_df['Species'] == sp]
    if not subset.empty:
        subset.to_csv(percentage_output_path, index=False)
        print(f"Percentages + raw counts for {sp} saved to {percentage_output_path}")
    else:
        print(f"No data available for species: {sp}")

# ------------------------------------------------------------------------------
# 10. Create the Manhattan Plot (ARGs grouped by drug class)
# ------------------------------------------------------------------------------
# Mapping ARG -> Drug Class
arg_to_drug_class = {
    "AAC(6')-Ie-APH(2'')-Ia bifunctional protein": 'AminoAB',
    "APH(2'')-If": 'AminoAB',
    "APH(2'')-Ig": 'AminoAB',
    "APH(3')-IIIa": 'AminoAB',
    "APH(3')-VIIa": 'AminoAB',
    "ANT(4')-Ia": 'AminoAB',
    "aad(6)": 'AminoAB',
    "ErmA": 'StreptoAB MacroAB LincosAB',
    "ErmB": 'StreptoAB MacroAB LincosAB',
    "OXA-184": 'Ceph Penam Carba',
    "OXA-185": 'Ceph Penam Carba',
    "OXA-193": 'Ceph Penam Carba',
    "OXA-447": 'Ceph Penam Carba',
    "OXA-449": 'Ceph Penam Carba',
    "OXA-465": 'Ceph Penam Carba',
    "OXA-576": 'Ceph Penam Carba',
    "OXA-605": 'Ceph Penam Carba',
    "OXA-617": 'Ceph Penam Carba',
    "OXA-622": 'Ceph Penam Carba',
    "OXA-623": 'Ceph Penam Carba',
    "OXA-625": 'Ceph Penam Carba',
    "OXA-631": 'Ceph Penam Carba',
    "OXA-632": 'Ceph Penam Carba',
    "OXA-635": 'Ceph Penam Carba',
    "OXA-638": 'Ceph Penam Carba',
    "cmeA": 'Ceph FusiAB MacroAB FluoroAB',
    "cmeB": 'Ceph FusiAB MacroAB FluoroAB',
    "cmeC": 'Ceph FusiAB MacroAB FluoroAB',
    "cmeR": 'Ceph FusiAB MacroAB FluoroAB',
    "tet(L)": 'TetraAB',
    "tet(O)": 'TetraAB',
    "tet(O/M/O)": 'TetraAB',
    "tet(W)": 'TetraAB',
    "tet(M)": 'TetraAB',
    "Campylobacter coli chloramphenicol acetyltransferase": 'PhenAB',
    "Limosilactobacillus reuteri cat-TC": 'PhenAB',
    "fexA": 'PhenAB',
    "lnuC": 'LincosAB',
    "lnuP": 'LincosAB',
    "SAT-4": 'NucleoAB',
    "optrA": 'PhenAB OxazoAB',
    "vatE": 'StreptoAB',
    "vatH": 'StreptoAB',
    "vgaD": 'StreptoAB',
    "mel": 'StreptoAB MacroAB',
    "T86I": 'FluoroAB',
    "A2075G": 'MacroAB',
}

drug_classes_order = [
    'AminoAB', 'Ceph Penam Carba', 'Ceph FusiAB MacroAB FluoroAB',
    'TetraAB', 'LincosAB', 'NucleoAB', 'StreptoAB', 'StreptoAB MacroAB', 
    'StreptoAB MacroAB LincosAB','PhenAB', 'PhenAB OxazoAB', 'FluoroAB', 'MacroAB',
]

# Build groupings
drug_class_to_args = {}
for arg in antibiotic_categories['ARGs']:
    dc = arg_to_drug_class.get(arg, 'Unknown')
    drug_class_to_args.setdefault(dc, []).append(arg)

args_grouped = []
drug_class_positions = []
current_pos = 0
for dc in drug_classes_order:
    these_args = drug_class_to_args.get(dc, [])
    args_grouped.extend(these_args)
    drug_class_positions.append((dc, current_pos, current_pos + len(these_args)))
    current_pos += len(these_args)

n_genes = len(args_grouped)

# Markers for species
markers = {
    'Campylobacter jejuni': 's',
    'Campylobacter coli': '^',
}

# Significance thresholds
significance_threshold_raw = 0.05
significance_threshold_corrected = 0.05

# Short rename dictionary for ARGs
rename_dict = {
    "APH(2'')-If": "aph(2'')-If",
    "AAC(6')-Ie-APH(2'')-Ia bifunctional protein": "aac(6')-Ie-aph(2'')-Ia",
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
    "Campylobacter coli chloramphenicol acetyltransferase": "Camp cat",
    "Limosilactobacillus reuteri cat-TC": "L. reuteri cat-TC",
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

# ------------------------------------------------------------------------------
# 11. Enhanced Manhattan plots function for publication quality
# ------------------------------------------------------------------------------
def create_manhattan_plots(
    show_x_title=True,
    y_axis_font_size=20,
    legend_font_size=24,
    subplot_label_size=32,
    show_sample_sizes=False,
    x_label_font_size=26,
    marker_size=350,
    legend_markerscale=2.0
):
    """
    marker_size: scatter s=â€¦ value for points
    legend_markerscale: how much to scale markers in the legend
    """

    fig = plt.figure(figsize=(50, 32))
    grid = plt.GridSpec(2, 2, figure=fig, hspace=0.5, wspace=0.25)

    country_positions = {
        'Bangladesh': (0, 0, 'A'),
        'India':      (0, 1, 'B'),
        'Vietnam':    (1, 0, 'C'),
        'Sri Lanka':  (1, 1, 'D')
    }

    for country, (r, c, label) in country_positions.items():
        ax = fig.add_subplot(grid[r, c])
        ax.set_facecolor('whitesmoke')

        # 1) Plot each species
        for sp in species_of_interest:
            df_sp = data_df[(data_df['Country']==country)&(data_df['Species_Name']==sp)]
            presence = [
                ((df_sp[gene]>0).sum()/df_sp[gene].count()*100) if gene in df_sp else 0
                for gene in args_grouped
            ]

            ax.scatter(
                np.arange(len(args_grouped)),
                presence,
                color=[custom_colors['ARGs'].get(g,'#CCCCCC') for g in args_grouped],
                marker=markers[sp],
                s=marker_size,          # â† adjustable
                label=rename_species.get(sp,sp) + (f" (n={len(df_sp)})" if show_sample_sizes else ""),
                edgecolors='black',
                linewidth=0.8
            )

        # 2) Stars for significance only if BOTH raw *and* BHâ€corrected p < 0.05
        for i, gene in enumerate(args_grouped):
            rows = results_within_countries[
                (results_within_countries['Country']==country)&
                (results_within_countries['Antibiotic']==gene)
            ]
            if any((rows['P-Value'] < significance_threshold_raw) &
                   (rows['Corrected P-Value'] < significance_threshold_corrected)):
                ax.scatter(i, 106,                    # â† moved from 102 to 106
                           marker='*',
                           s=marker_size*0.8,         # â† increased from 0.5 to 0.8
                           color='black',
                           zorder=10)

        # 3) Subplot label
        ax.text(0.5, 1.12, f'{label}) {country}',
                transform=ax.transAxes,
                ha='center', fontsize=subplot_label_size, fontweight='bold')

        # 4) Axes
        ax.set_xticks(np.arange(len(args_grouped)))
        ax.set_xticklabels([rename_dict.get(g,g) for g in args_grouped],
                           rotation=90, fontsize=x_label_font_size, fontweight='bold')
        for dc, start, end in drug_class_positions:
            ax.axvline(x=start-0.5, color='grey', linestyle='--', linewidth=0.5)
        ax.axhline(100, color='grey', linestyle='--', linewidth=0.5)

        ax.set_ylim(0, 110)
        ax.set_yticks(np.arange(0, 101, 10))
        if show_x_title:
            ax.set_xlabel('AMR Genes grouped by Drug Class',
                          fontsize=x_label_font_size, fontweight='bold')
        ax.set_ylabel('Percentage Presence',
                      fontsize=y_axis_font_size, fontweight='bold')
        ax.tick_params(axis='y', labelsize=y_axis_font_size)
        ax.tick_params(axis='x', labelsize=x_label_font_size)

    # 5) Legend with larger markers
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels,
               loc='lower center',
               bbox_to_anchor=(0.5, 0.02),
               ncol=len(handles),
               fontsize=legend_font_size,
               markerscale=legend_markerscale,  # â† adjustable
               frameon=True, fancybox=True, shadow=True)

    # 6) Layout & save
    plt.subplots_adjust(top=0.90, bottom=0.15,
                        left=0.06, right=0.98,
                        hspace=0.35, wspace=0.15)

    out = os.path.join(output_dir, 'combined_manhattan_plots_publication.png')
    plt.savefig(out, dpi=1200, bbox_inches='tight', facecolor='white')
    plt.show()
    
# ------------------------------------------------------------------------------
# 12. Enhanced Stacked Bar Charts
# ------------------------------------------------------------------------------
def create_stacked_bar_charts(
    legend_font_size=16, 
    show_percentages=False, 
    show_sample_sizes=True,
    subplot_label_size=24,
    y_axis_font_size=16,
    percentage_font_size=10,        # NEW: Control percentage text size
    min_percentage_to_show=5.0,     # NEW: Only show percentages >= this value
    round_percentages=True          # NEW: Round to nearest whole number
):
    fig = plt.figure(figsize=(24, 24))
    grid = plt.GridSpec(2, 2, figure=fig, hspace=0.5, wspace=0.5)

    country_positions = {
        'Bangladesh': (0, 0, 'E'),
        'India': (0, 1, 'F'),
        'Vietnam': (1, 0, 'G'),
        'Sri Lanka': (1, 1, 'H')
    }

    for country, (row, col, label) in country_positions.items():
        ax = fig.add_subplot(grid[row, col])
        
        # Summation of each drug class
        percentages_dc = pd.DataFrame(0, index=species_of_interest, columns=antibiotic_categories['Drug Classes'])
        sample_sizes = {}

        for sp in species_of_interest:
            subset_df = data_df[(data_df['Country'] == country) & (data_df['Species_Name'] == sp)]
            sample_sizes[sp] = len(subset_df)
            if subset_df.empty:
                percentages_dc.drop(sp, inplace=True)
                continue

            for dc in antibiotic_categories['Drug Classes']:
                if dc in subset_df.columns:
                    percentages_dc.loc[sp, dc] = subset_df[dc].sum()
                else:
                    percentages_dc.loc[sp, dc] = 0

        total_counts_dc = percentages_dc.sum(axis=1)
        percentages_dc = percentages_dc.div(total_counts_dc, axis=0).fillna(0) * 100

        # Build y-labels with a line break (species + n=xx)
        y_labels = []
        for sp in percentages_dc.index:
            short_sp = rename_species.get(sp, sp)
            if show_sample_sizes:
                y_labels.append(f"{short_sp}\n(n={sample_sizes[sp]})")
            else:
                y_labels.append(short_sp)

        bottom_dc = np.zeros(len(percentages_dc))
        for dc in antibiotic_categories['Drug Classes']:
            vals = percentages_dc[dc].values
            color_ = custom_colors['Drug Classes'].get(dc, '#CCCCCC')
            bars = ax.barh(np.arange(len(percentages_dc)), vals, left=bottom_dc, color=color_, label=dc)
            
            if show_percentages:
                for bar, val in zip(bars, vals):
                    # Only show percentage if it meets minimum threshold
                    if val >= min_percentage_to_show:
                        # Round to nearest whole number if requested
                        if round_percentages:
                            display_val = f"{round(val)}%"
                        else:
                            display_val = f"{val:.1f}%"
                        
                        ax.text(bar.get_x() + bar.get_width() / 2,
                                bar.get_y() + bar.get_height() / 2,
                                display_val,
                                ha='center', va='center', 
                                fontsize=percentage_font_size, 
                                color='black', 
                                fontweight='bold')
            bottom_dc += vals

        ax.set_yticks(np.arange(len(percentages_dc)))
        ax.set_yticklabels(y_labels, fontsize=y_axis_font_size, fontweight='bold')

        ax.text(-0.1, 1.15, f'{label}) {country}',
                transform=ax.transAxes,
                fontsize=subplot_label_size,
                fontweight='bold')
        
        ax.set_title(f'Stacked Bar Chart for Drug Classes in {country}', 
                     fontsize=18, fontweight='bold', pad=20)
        ax.set_xlabel('Total Percentage', fontsize=16, fontweight='bold')
        ax.tick_params(axis='y', labelsize=y_axis_font_size, pad=15)
        ax.tick_params(axis='x', labelsize=14)

    handles, labels_ = ax.get_legend_handles_labels()
    ncols = min(6, len(handles))
    fig.legend(handles, labels_,
               loc='center',
               bbox_to_anchor=(0.5, 0.02),
               ncol=ncols,
               fontsize=legend_font_size,
               frameon=True,
               fancybox=True,
               shadow=True)

    plt.subplots_adjust(top=0.92, bottom=0.18, left=0.12, right=0.95)

    bar_out = os.path.join(output_dir, 'combined_stacked_bar_charts_publication.png')
    plt.savefig(bar_out, dpi=1200, bbox_inches='tight', facecolor='white')
    plt.show()

# ------------------------------------------------------------------------------
# 13. Generate publication-quality plots
# ------------------------------------------------------------------------------
create_manhattan_plots(
    show_x_title=False,
    y_axis_font_size=20,
    legend_font_size=24,
    subplot_label_size=32,
    show_sample_sizes=False,
    x_label_font_size=26,
    marker_size=400,           # â† increase point size here
    legend_markerscale=2.5     # â† increase legend-symbol size here
)

# Example 1: Show percentages rounded to whole numbers, minimum 5%, font size 12
create_stacked_bar_charts(
    legend_font_size=18,
    show_percentages=True,           # Turn on percentages
    show_sample_sizes=True,
    subplot_label_size=28,
    y_axis_font_size=18,
    percentage_font_size=12,         # Size of percentage text
    min_percentage_to_show=5.0,      # Only show if >= 5%
    round_percentages=True           # Round 8.7% to 9%
)

# Example 2: Show percentages with decimals, minimum 3%, smaller font
# create_stacked_bar_charts(
#     legend_font_size=18,
#     show_percentages=True,
#     show_sample_sizes=True,
#     subplot_label_size=28,
#     y_axis_font_size=18,
#     percentage_font_size=10,         # Smaller percentage text
#     min_percentage_to_show=3.0,      # Show if >= 3%
#     round_percentages=False          # Keep decimals: 8.7%
# )

# Example 3: No percentages shown (clean chart)
# create_stacked_bar_charts(
#     legend_font_size=18,
#     show_percentages=False,          # No percentages
#     show_sample_sizes=True,
#     subplot_label_size=28,
#     y_axis_font_size=18
# )

print("All done! CSVs and Excel files saved, publication-quality Manhattan plots and stacked bar charts generated.")

#updated country 08-07-25
# Version 502 - Complete Adjusted Code country and species SIZE ADJUSTED
