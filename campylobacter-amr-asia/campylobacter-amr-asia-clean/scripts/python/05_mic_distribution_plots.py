"""
MIC Distribution Plots
======================
Distribution analysis of MIC values across antibiotics

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
###5000 distribution
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import surpyval as surv
import os
from matplotlib.lines import Line2D
from matplotlib import patches

class MICAnalyzer:
    def __init__(self, data_path):
        """Initialize the MIC analyzer with data and cutoff values"""
        self.data = pd.read_csv(data_path)
        self.interpretive_standards = {

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

            'Azithromycin':    {'C. jejuni': {'S': 0.25, 'R': 0.5}, 'C. coli': {'S': 0.5, 'R': 1}},
            'Ciprofloxacin':   {'C. jejuni': {'S': 0.5,  'R': 1},   'C. coli': {'S': 0.5,  'R': 1}},
            'Clindamycin':     {'C. jejuni': {'S': 0.5,  'R': 1},   'C. coli': {'S': 1,    'R': 2}},
            'Erythromycin':    {'C. jejuni': {'S': 4,    'R': 8},   'C. coli': {'S': 8,    'R': 16}},
            'Florfenicol':     {'C. jejuni': {'S': 4,    'R': 8},   'C. coli': {'S': 4,    'R': 8}},
            'Tetracycline':    {'C. jejuni': {'S': 1,    'R': 2},   'C. coli': {'S': 2,    'R': 4}},
            'Gentamicin':      {'C. jejuni': {'S': 2,    'R': 4},   'C. coli': {'S': 2,    'R': 4}},
            'Nalidixic Acid':  {'C. jejuni': {'S': 16,   'R': 32},  'C. coli': {'S': 16,   'R': 32}}
        }
        
        if 'Country' not in self.data.columns or 'Species_Name' not in self.data.columns:
            raise ValueError("Data must contain 'Country' and 'Species_Name' columns.")

    def get_breakpoints(self, antibiotic, species):
        """Get species-specific breakpoints for an antibiotic"""
        species_short = 'C. jejuni' if 'jejuni' in species else 'C. coli'
        return self.interpretive_standards[antibiotic][species_short]

    def prepare_interval_censored_data(self, antibiotic, country=None, species=None):
        """Prepare data for interval-censored survival analysis"""
        data = self.data.copy()
        if country:
            data = data[data['Country'] == country]
        if species:
            data = data[data['Species_Name'] == species]
        
        lower_bounds = []
        upper_bounds = []
        
        for val in data[antibiotic].dropna():
            val = str(val).strip()
            if val.startswith('<'):
                lower_bounds.append(-np.inf)
                upper_bounds.append(float(val[1:]))
            elif val.startswith('>'):
                lower_bounds.append(float(val[1:]))
                upper_bounds.append(np.inf)
            else:
                try:
                    mic = float(val)
                    lower_bounds.append(mic)
                    upper_bounds.append(mic)
                except ValueError:
                    continue
        
        return np.array(lower_bounds), np.array(upper_bounds)
    
    def analyze_interval_censored(self, antibiotic, output_dir=None):
        """Perform interval-censored survival analysis for one antibiotic"""
        print(f"\nInterval-Censored Survival Analysis for {antibiotic}")
        print("-" * 50)

        for species in self.data['Species_Name'].unique():
            for country in self.data['Country'].unique():
                lower_bounds, upper_bounds = self.prepare_interval_censored_data(antibiotic, country, species)
                n_samples = len(lower_bounds)
                if n_samples > 0:
                    try:
                        model = surv.NonParametric('Turnbull')
                        model.fit_intervals(lower_bounds, upper_bounds)

                        x = np.unique(np.concatenate([
                            lower_bounds[~np.isinf(lower_bounds)], 
                            upper_bounds[~np.isinf(upper_bounds)]
                        ]))
                        x.sort()
                        sf = model.sf(x)

                        median = model.median
                        print(f"\n{species}, {country}:")
                        print(f"Total samples: {n_samples}")
                        print(f"Median MIC: {median}")

                        # Get species-specific breakpoints
                        breakpoints = self.get_breakpoints(antibiotic, species)
                        s_breakpoint = breakpoints['S']
                        r_breakpoint = breakpoints['R']

                        plt.figure(figsize=(10, 6))
                        plt.step(x, sf, where='post')
                        plt.xscale('log', base=2)
                        plt.xlabel(f'{antibiotic} MIC (mg/L)')
                        plt.ylabel('Survival Probability')
                        plt.title(f'{antibiotic} - {species}, {country}\n(n={n_samples})')

                        plt.axvline(x=s_breakpoint, color='blue', linestyle='--', 
                                    alpha=0.6, label='Susceptible cutoff')
                        plt.axvline(x=r_breakpoint, color='red', linestyle='--', 
                                    alpha=0.6, label='Resistant cutoff')
                        plt.legend()

                        mic_values = sorted(set(
                            lower_bounds[~np.isinf(lower_bounds)].tolist() + 
                            upper_bounds[~np.isinf(upper_bounds)].tolist()
                        ))
                        n = max(1, len(mic_values) // 10)
                        mic_values_to_plot = mic_values[::n]
                        plt.xticks(mic_values_to_plot, 
                                   labels=[str(mic) for mic in mic_values_to_plot], 
                                   rotation=45)

                        plt.tight_layout()

                        if output_dir:
                            plot_filename = f"{antibiotic}_{species}_{country}_interval_censored.png"
                            plt.savefig(os.path.join(output_dir, plot_filename), 
                                        bbox_inches='tight', dpi=800)
                        plt.close()
                    except Exception as e:
                        print(f"Error analyzing {antibiotic} for {species}, {country}: {e}")
                else:
                    print(f"\n{species}, {country}: No data available.")

    def calculate_resistance_proportions(self, antibiotic):
        """Calculate proportions of susceptible, intermediate, and resistant isolates"""
        data = self.data.copy()
        data['Species_Name'] = data['Species_Name'].replace({
            'Campylobacter jejuni': 'C. jejuni',
            'Campylobacter coli':   'C. coli'
        })
        
        def categorize_mic(row):
            try:
                val = str(row[antibiotic]).strip()
                species = row['Species_Name']
                breakpoints = self.get_breakpoints(antibiotic, species)
                s_breakpoint = breakpoints['S']
                r_breakpoint = breakpoints['R']
                
                if val.startswith('<'):
                    mic = float(val[1:])
                    if mic <= s_breakpoint:
                        return 'Susceptible'
                    else:
                        return 'Intermediate'
                elif val.startswith('>'):
                    mic = float(val[1:])
                    if mic >= r_breakpoint:
                        return 'Resistant'
                    else:
                        return 'Intermediate'
                else:
                    mic = float(val)
                    if mic <= s_breakpoint:
                        return 'Susceptible'
                    elif mic >= r_breakpoint:
                        return 'Resistant'
                    else:
                        return 'Intermediate'
            except:
                return 'Unknown'
        
        data['Category'] = data.apply(categorize_mic, axis=1)
        
        proportions = data.groupby(['Country', 'Species_Name', 'Category']).size().unstack(fill_value=0)
        proportions['Total'] = proportions.sum(axis=1)
        proportions_percent = proportions.div(proportions['Total'], axis=0) * 100
        
        return proportions_percent.reset_index()
    
    def plot_resistance_proportions(self, antibiotic, output_dir=None):
        """Plot proportions of resistance categories"""
        proportions = self.calculate_resistance_proportions(antibiotic)
        
        proportions_melted = proportions.melt(
            id_vars=['Country', 'Species_Name', 'Total'], 
            value_vars=['Susceptible', 'Intermediate', 'Resistant'],
            var_name='Category', value_name='Proportion'
        )
        
        plt.figure(figsize=(12, 8))
        sns.barplot(data=proportions_melted, x='Country', y='Proportion', 
                    hue='Category', ci=None, estimator=sum)
        plt.title(f'Proportions of Resistance Categories for {antibiotic}')
        plt.ylabel('Proportion (%)')
        plt.xlabel('Country')
        plt.legend(title='Category')
        plt.tight_layout()
        
        for index, row in proportions.iterrows():
            plt.text(index, 100, f"n={int(row['Total'])}", ha='center', va='bottom')

        if output_dir:
            plot_filename = f"{antibiotic}_resistance_proportions.png"
            plt.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight', dpi=800)
        plt.close()

    def plot_mic_distribution(self, antibiotic, output_dir=None):
        """Plot the distribution of MIC concentrations as a bar chart stratified by species and country"""
        data = self.data.copy()
        
        data['Species_Name'] = data['Species_Name'].replace({
            'Campylobacter jejuni': 'C. jejuni',
            'Campylobacter coli':   'C. coli'
        })
        
        plot_data = []
        less_than_values = set()
        greater_than_values = set()
        exact_values = set()
        
        for _, row in data.iterrows():
            val = row.get(antibiotic, np.nan)
            if pd.isna(val):
                continue
            val_str = str(val).strip()
            try:
                if val_str.startswith('<'):
                    mic = float(val_str[1:])
                    mic_label = f"<{mic}"
                    less_than_values.add(mic_label)
                elif val_str.startswith('>'):
                    mic = float(val_str[1:])
                    mic_label = f">{mic}"
                    greater_than_values.add(mic_label)
                else:
                    mic = float(val_str)
                    mic_label = str(mic)
                    exact_values.add(mic_label)
                    
                plot_data.append({
                    'MIC': mic,
                    'MIC_Label': mic_label,
                    'Country': row['Country'],
                    'Species': row['Species_Name'],
                    'ST': str(row['ST']) if pd.notna(row['ST']) else 'Unknown'
                })
            except:
                continue

        plot_data = pd.DataFrame(plot_data)
        
        less_than_sorted    = sorted(less_than_values,    key=lambda x: float(x.replace('<','')))
        exact_sorted        = sorted(exact_values,        key=float)
        greater_than_sorted = sorted(greater_than_values, key=lambda x: float(x.replace('>','')))
        
        all_mic_values = less_than_sorted + exact_sorted + greater_than_sorted

        species_palette = {
            'C. jejuni': '#66CCFF',
            'C. coli':   '#FF9999'
        }
        
        plt.figure(figsize=(15, 12))
        
        n_countries = len(plot_data['Country'].unique())
        n_cols = 2
        n_rows = (n_countries + 1) // 2
        
        for i, country in enumerate(sorted(plot_data['Country'].unique())):
            ax = plt.subplot(n_rows, n_cols, i + 1)
            country_data = plot_data[plot_data['Country'] == country]
            
            bar_width = 0.35
            species_offsets = {'C. jejuni': -bar_width/2, 'C. coli': bar_width/2}
            
            for species in species_palette.keys():
                sp_data = country_data[country_data['Species'] == species]
                if sp_data.empty:
                    continue
                breakpoints = self.get_breakpoints(antibiotic, species)
                high_mic_threshold = breakpoints['R']
                
                mic_counts = pd.Series(0, index=all_mic_values, dtype=int)
                sp_counts  = sp_data.groupby('MIC_Label').size()
                mic_counts.update(sp_counts)
                
                x_positions = np.arange(len(all_mic_values)) + species_offsets[species]
                ax.bar(
                    x_positions,
                    mic_counts,
                    width=bar_width,
                    color=species_palette[species],
                    edgecolor='black',
                    linewidth=0.5,
                    label=species,
                    zorder=2
                )
                
                sp_high_mic = sp_data[sp_data['MIC'] >= high_mic_threshold]
                st_counts = sp_high_mic[sp_high_mic['ST'] != '-']['ST'].value_counts()
                sig_sts = st_counts[st_counts >= 5].index.tolist()
                sig_st_data = sp_data[sp_data['ST'].isin(sig_sts)]
                if not sig_st_data.empty:
                    overlay_color = '#003366' if species == 'C. jejuni' else '#993333'
                    st_mic_counts = pd.Series(0, index=all_mic_values, dtype=int)
                    st_temp = sig_st_data.groupby('MIC_Label').size()
                    st_mic_counts.update(st_temp)
                    
                    ax.bar(
                        x_positions,
                        st_mic_counts,
                        width=bar_width,
                        color=overlay_color,
                        alpha=0.5,
                        edgecolor='black',
                        linewidth=0.5,
                        zorder=3
                    )
                
                base_color = '#003366' if species == 'C. jejuni' else '#993333'
                
                def find_position(value):
                    for idx_mic, v in enumerate(all_mic_values):
                        if float(v.replace('<','').replace('>','')) >= value:
                            return idx_mic
                    return len(all_mic_values)
                
                s_pos = find_position(breakpoints['S'])
                r_pos = find_position(breakpoints['R'])
                
                ax.axvline(
                    x=s_pos + species_offsets[species],
                    color=base_color,
                    linestyle='-',
                    alpha=0.8,
                    linewidth=2,
                    zorder=4
                )
                ax.axvline(
                    x=r_pos + species_offsets[species],
                    color=base_color,
                    linestyle='--',
                    alpha=0.8,
                    linewidth=2,
                    zorder=4
                )
            
            ax.set_xticks(range(len(all_mic_values)))
            ax.set_xticklabels(all_mic_values, rotation=45)
            ax.set_ylabel('Number of Isolates')
            ax.set_xlabel(f'{antibiotic} MIC (mg/L)')
            ax.grid(True, axis='y', linestyle='-', alpha=0.2, color='gray')
            
            total_n = len(country_data)
            species_counts = country_data['Species'].value_counts().to_dict()
            
            species_stats = []
            for sp in species_palette.keys():
                sp_data = country_data[country_data['Species'] == sp]
                if sp_data.empty:
                    continue
                brk = self.get_breakpoints(antibiotic, sp)
                sp_r = sp_data[sp_data['MIC'] >= brk['R']]
                st_info = sp_r[sp_r['ST'] != '-']['ST'].value_counts()
                stats = f"{sp}: n={species_counts.get(sp, 0)}\n"
                stats += f"Sâ‰¤{brk['S']}, Râ‰¥{brk['R']}"
                if not st_info.empty:
                    dom_sts = [
                        f"ST{st}: n={st_info[st]}"
                        for st in st_info.index[:10]
                        if st_info[st] >= 5
                    ]
                    if dom_sts:
                        stats += f"\nDominant STs (R):\n" + "\n".join(dom_sts)
                species_stats.append(stats)
            
            species_n_text = "\n\n".join(species_stats)
            ax.set_title(f'{country}\nTotal n={total_n}')
            
            ax.text(
                1.02,
                0.98,
                species_n_text,
                transform=ax.transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
            )

        legend_fig = plt.figure(figsize=(12, 2))
        legend_ax = legend_fig.add_subplot(111)
        legend_ax.axis('off')
        
        legend_handles = []
        for sp in species_palette.keys():
            base_color = '#003366' if sp == 'C. jejuni' else '#993333'
            legend_handles.extend([
                Line2D([0], [0], color=species_palette[sp], lw=2, label=f'{sp}'),
                Line2D([0], [0], color=base_color, linestyle='-', 
                       alpha=0.8, linewidth=1, label=f'{sp} S breakpoint'),
                Line2D([0], [0], color=base_color, linestyle='--', 
                       alpha=0.8, linewidth=1, label=f'{sp} R breakpoint'),
                Line2D([0], [0], color=base_color,
                       linestyle='--', alpha=0.5, label=f'Dominant ST ({sp})')
            ])
        
        legend_ax.legend(
            handles=legend_handles,
            loc='center',
            ncol=8,
            frameon=True,
            title='Species and Cutoffs'
        )
        
        plt.figure(1)
        plt.tight_layout()
        
        # Save
        if output_dir:
            plot_filename = f"{antibiotic}_mic_distribution_with_st.png"
            plt.savefig(
                os.path.join(output_dir, plot_filename),
                bbox_inches='tight',
                dpi=800
            )
            
            legend_filename = f"{antibiotic}_mic_distribution_legend_with_st.png"
            legend_fig.savefig(
                os.path.join(output_dir, legend_filename),
                bbox_inches='tight',
                dpi=800
            )
        
        plt.close('all')

    # ---------------------------------------------------------------------
    # Combined Single-Figure (A-H) version
    # ---------------------------------------------------------------------
    def plot_combined_mic_distribution(
        self, 
        font_sizes={
            'title': 14, 
            'axis': 12, 
            'tick': 10, 
            'legend': 8, 
            'stats': 8
        }, 
        output_dir=None,
        bar_spacing=3.0
    ):
        """
        Your existing single combined approach. Omitted here for brevity,
        or you can keep your original code. 
        """
        pass

    # ---------------------------------------------------------------------
    # Combined Split Version (Aâ€“D, Eâ€“H) with shading
    # ---------------------------------------------------------------------
    def plot_combined_mic_distribution_split(
        self,
        font_sizes={
            'title': 22, 
            'axis': 18, 
            'tick': 16, 
            'legend': 20, 
            'stats': 18
        }, 
        output_dir=None,
        bar_spacing=3.5,    # Moderately increased bar spacing
        cat_spacing=2.5,    # Moderately increased category spacing
        cat_margin=0.35,    # Moderately increased margin
        figsize=(45, 40),   # Moderately increased figure size
        dpi=800,            
        bottom_adjust=0.12    
    ):
        """
        Splits the "combined" MIC distributions into two separate figures:
        Figure 1: Antibiotics A-D
        Figure 2: Antibiotics E-H
        """
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib import patches

        # A-D, E-H
        antibiotic_sets = [
            [
                ('Azithromycin', 'A'),
                ('Erythromycin', 'B'),
                ('Ciprofloxacin', 'C'),
                ('Nalidixic Acid', 'D')
            ],
            [
                ('Tetracycline', 'E'),
                ('Florfenicol', 'F'),
                ('Clindamycin', 'G'),
                ('Gentamicin', 'H')
            ]
        ]

        # Patterns, palette
        species_patterns = {
            'C. jejuni': '*',
            'C. coli': '-'
        }
        country_palette = {
            'Vietnam':    '#2ecc71',
            'India':      '#3498db',
            'Bangladesh': '#e74c3c',
            'Sri Lanka':  '#9b59b6'
        }

        def find_position(value, all_mic_values):
            for idx_mic, v in enumerate(all_mic_values):
                numeric_part = v.replace('<','').replace('>','')
                try:
                    cat_val = float(numeric_part)
                    if cat_val >= value:
                        return idx_mic
                except:
                    pass
            return len(all_mic_values) - 1

        # Build a sorted list of all possible (country, species) combos
        all_countries = sorted(list(self.data['Country'].unique()))
        all_species   = ['C. jejuni','C. coli']
        combos_in_data = []
        for c in all_countries:
            for s in all_species:
                sel = self.data[
                    (self.data['Country'] == c) &
                    (self.data['Species_Name'].str.contains(s.replace('C. ','')))
                ]
                if len(sel) > 0:
                    combos_in_data.append((c, s))
        combos_in_data = sorted(combos_in_data, key=lambda x: (x[0], x[1]))

        # Now build the plots
        for fig_index, subset in enumerate(antibiotic_sets, start=1):
            fig = plt.figure(figsize=figsize, constrained_layout=False)
            
            # Adjust subplot parameters
            plt.subplots_adjust(
                left=0.09,      # Slightly increased left margin
                right=0.96,     # Slightly reduced right margin
                top=0.96,       # Slightly reduced top margin
                bottom=bottom_adjust,
                wspace=0.20,    # Moderate spacing between plots horizontally
                hspace=0.25     # Moderate spacing between plots vertically
            )

            for subplot_idx, (antibiotic, label) in enumerate(subset, 1):
                ax = fig.add_subplot(2, 2, subplot_idx)

                data = self.data.copy()
                data['Species_Name'] = data['Species_Name'].replace({
                    'Campylobacter jejuni': 'C. jejuni',
                    'Campylobacter coli':   'C. coli'
                })

                # Gather MIC labels
                plot_data = []
                less_than_values = set()
                greater_than_values = set()
                exact_values = set()

                for _, row in data.iterrows():
                    val = row.get(antibiotic, None)
                    if pd.isna(val):
                        continue
                    val_str = str(val).strip()
                    try:
                        if val_str.startswith('<'):
                            mic = float(val_str[1:])
                            mic_label = f"<{mic}"
                            less_than_values.add(mic_label)
                        elif val_str.startswith('>'):
                            mic = float(val_str[1:])
                            mic_label = f">{mic}"
                            greater_than_values.add(mic_label)
                        else:
                            mic = float(val_str)
                            mic_label = str(mic)
                            exact_values.add(mic_label)
                        
                        plot_data.append({
                            'MIC': mic,
                            'MIC_Label': mic_label,
                            'Country':  row['Country'],
                            'Species':  row['Species_Name'],
                            'ST': str(row['ST']) if pd.notna(row['ST']) else 'Unknown'
                        })
                    except:
                        continue

                plot_data = pd.DataFrame(plot_data)
                
                less_than_sorted    = sorted(less_than_values,    key=lambda x: float(x.replace('<','')))
                exact_sorted        = sorted(exact_values,        key=float)
                greater_than_sorted = sorted(greater_than_values, key=lambda x: float(x.replace('>','')))
                all_mic_values = less_than_sorted + exact_sorted + greater_than_sorted

                # Define sub-slots
                combos = combos_in_data
                num_combos = len(combos)
                cluster_width = (cat_spacing - 2 * cat_margin)
                if num_combos > 0:
                    subslot_width = cluster_width / num_combos
                else:
                    subslot_width = cluster_width

                # Store data counts
                from collections import defaultdict
                counts_dict   = defaultdict(lambda: pd.Series(0, index=all_mic_values, dtype=int))
                overlay_dict  = defaultdict(lambda: pd.Series(0, index=all_mic_values, dtype=int))

                # Fill in counts
                for idx_row, row in plot_data.iterrows():
                    c  = row['Country']
                    sp = row['Species']
                    ml = row['MIC_Label']
                    if (c, sp) in combos:
                        counts_dict[(c, sp)][ml] += 1
                
                # For dominant ST
                for c, sp in combos:
                    brk = self.get_breakpoints(antibiotic, sp)
                    r_cut = brk['R']
                    sp_data = plot_data[
                        (plot_data['Country'] == c) &
                        (plot_data['Species'] == sp)
                    ]
                    high_data = sp_data[sp_data['MIC'] >= r_cut]
                    st_counts = high_data[high_data['ST'] != 'Unknown']['ST'].value_counts()
                    sig_sts   = st_counts[st_counts >= 5].index.tolist()
                    if len(sig_sts) > 0:
                        sig_data = high_data[high_data['ST'].isin(sig_sts)]
                        for ml, ct in sig_data.groupby('MIC_Label').size().items():
                            overlay_dict[(c, sp)][ml] += ct

                # Plot each category
                for i_cat, mic_label in enumerate(all_mic_values):
                    cat_left = i_cat * cat_spacing

                    for j, (country, species) in enumerate(combos):
                        x_left = cat_left + cat_margin + j * subslot_width
                        bar_w = subslot_width * 0.9

                        cnt  = counts_dict[(country, species)].get(mic_label, 0)
                        ovr  = overlay_dict[(country, species)].get(mic_label, 0)

                        color_main = country_palette.get(country, 'gray')
                        hatch_pat  = species_patterns.get(species, '/')
                        
                        # Main bar
                        ax.bar(
                            x_left,
                            cnt,
                            width=bar_w,
                            color=color_main,
                            edgecolor='black',
                            hatch=hatch_pat,
                            alpha=0.4,
                            linewidth=0.5,
                            zorder=2
                        )
                        # Overlay bar
                        ax.bar(
                            x_left,
                            ovr,
                            width=bar_w,
                            color=color_main,
                            edgecolor='black',
                            hatch=hatch_pat,
                            alpha=1.0,
                            linewidth=0.5,
                            zorder=3
                        )

                # Breakpoint lines
                j_brk = self.get_breakpoints(antibiotic, 'C. jejuni')
                c_brk = self.get_breakpoints(antibiotic, 'C. coli')
                s_j = find_position(j_brk['S'], all_mic_values)
                r_j = find_position(j_brk['R'], all_mic_values)
                s_c = find_position(c_brk['S'], all_mic_values)
                r_c = find_position(c_brk['R'], all_mic_values)

                s_j_x = s_j * cat_spacing
                r_j_x = r_j * cat_spacing
                s_c_x = s_c * cat_spacing
                r_c_x = r_c * cat_spacing

                if j_brk != c_brk:
                    ax.axvline(s_j_x, color='black', linestyle='-',  alpha=0.8, linewidth=1, zorder=4)
                    ax.axvline(r_j_x, color='black', linestyle='--', alpha=0.8, linewidth=1, zorder=4)
                    ax.axvline(s_c_x, color='black', linestyle=':',  alpha=0.8, linewidth=1, zorder=4)
                    ax.axvline(r_c_x, color='black', linestyle='-.', alpha=0.8, linewidth=1, zorder=4)
                else:
                    ax.axvline(s_j_x, color='black', linestyle='-',  alpha=0.8, linewidth=1, zorder=4)
                    ax.axvline(r_j_x, color='black', linestyle='--', alpha=0.8, linewidth=1, zorder=4)

                # Alternating shading
                ymin, ymax = ax.get_ylim()
                num_cats = len(all_mic_values)
                for i_cat2 in range(num_cats):
                    if i_cat2 % 2 == 0:
                        x_left = i_cat2 * cat_spacing
                        rect = patches.Rectangle(
                            (x_left, ymin),
                            cat_spacing,
                            ymax - ymin,
                            facecolor='lightgrey',
                            alpha=0.2,
                            zorder=0
                        )
                        ax.add_patch(rect)

                # Set ticks and labels
                cat_positions = (np.arange(num_cats) + 0.5) * cat_spacing
                ax.set_xticks(cat_positions)
                ax.set_xticklabels(all_mic_values, rotation=45, fontsize=font_sizes['tick'])
                ax.tick_params(axis='y', labelsize=font_sizes['tick'])
                ax.set_ylabel('Number of Isolates', fontsize=font_sizes['axis'])
                ax.set_xlabel('MIC (mg/L)', fontsize=font_sizes['axis'])

                ax.set_title(
                    f'{label}) {antibiotic}',
                    fontsize=font_sizes['title'],
                    pad=30,
                    loc='left',
                    weight='bold'
                )
                
                # Increase the size of the plot
                ax.set_position([
                    ax.get_position().x0 * 0.98,     # Slight adjustment
                    ax.get_position().y0 * 0.98,     # Slight adjustment
                    ax.get_position().width * 1.05,  # Increase width by 5%
                    ax.get_position().height * 1.05  # Increase height by 5%
                ])

                ax.grid(True, axis='y', linestyle='-', alpha=0.2, color='gray')

            # Build legend
            legend_handles = []
            for (country, species) in combos_in_data:
                patch_main = patches.Patch(
                    facecolor=country_palette.get(country, 'gray'),
                    edgecolor='black',
                    label=f'{country} - {species}',
                    hatch=species_patterns.get(species, '/'),
                    alpha=0.4
                )
                legend_handles.append(patch_main)
                
                patch_st = patches.Patch(
                    facecolor=country_palette.get(country, 'gray'),
                    edgecolor='black',
                    label=f'{country} - {species} (Dominant ST)',
                    hatch=species_patterns.get(species, '/'),
                    alpha=1.0
                )
                legend_handles.append(patch_st)

            # Add breakpoint lines to legend
            ab_label = subset[0][0]
            j_br = self.get_breakpoints(ab_label, 'C. jejuni')
            c_br = self.get_breakpoints(ab_label, 'C. coli')
            
            if j_br != c_br:
                legend_handles.extend([
                    Line2D([0], [0], color='black', linestyle='-',  
                        alpha=0.8, linewidth=1, label='C. jejuni S breakpoint'),
                    Line2D([0], [0], color='black', linestyle='--', 
                        alpha=0.8, linewidth=1, label='C. jejuni R breakpoint'),
                    Line2D([0], [0], color='black', linestyle=':',
                        alpha=0.8, linewidth=1, label='C. coli S breakpoint'),
                    Line2D([0], [0], color='black', linestyle='-.',
                        alpha=0.8, linewidth=1, label='C. coli R breakpoint')
                ])
            else:
                legend_handles.extend([
                    Line2D([0], [0], color='black', linestyle='-',
                        alpha=0.8, linewidth=1, label='Both species S breakpoint'),
                    Line2D([0], [0], color='black', linestyle='--',
                        alpha=0.8, linewidth=1, label='Both species R breakpoint')
                ])

            # Create and position legend
            legend_ax = fig.add_axes([0.1, 0.02, 0.8, 0.05])
            legend_ax.axis('off')
            legend_ax.legend(
                handles=legend_handles,
                loc='center',
                ncol=4,
                fontsize=font_sizes['legend'],
                frameon=True
            )

            # Save figure
            if output_dir:
                plot_filename = f"combined_mic_distribution_{fig_index}.png"
                plt.savefig(
                    os.path.join(output_dir, plot_filename),
                    bbox_inches='tight',
                    dpi=dpi
                )
            
            plt.close()

    def export_mic_data(self, antibiotic, output_dir=None):
        """Export MIC data to a CSV file for AI analysis"""
        data = self.data.copy()
        
        mic_values = []
        countries = []
        species_list = []
        other_columns = data.columns.difference([antibiotic, 'Country', 'Species_Name']).tolist()
        other_data = {col: [] for col in other_columns}
        
        for _, row in data.iterrows():
            val = row.get(antibiotic, np.nan)
            if pd.isna(val):
                continue
            val_str = str(val).strip()
            try:
                if val_str.startswith('<'):
                    mic = float(val_str[1:])
                elif val_str.startswith('>'):
                    mic = float(val_str[1:])
                else:
                    mic = float(val_str)
                mic_values.append(mic)
                countries.append(row['Country'])
                species_list.append(row['Species_Name'])
                for col in other_columns:
                    other_data[col].append(row[col])
            except:
                continue
        
        export_data = pd.DataFrame({
            'MIC': mic_values,
            'Country': countries,
            'Species': species_list,
            **other_data
        })
        
        if output_dir:
            csv_filename = f"{antibiotic}_mic_data.csv"
            export_data.to_csv(os.path.join(output_dir, csv_filename), index=False)
            print(f"MIC data exported to {csv_filename}")
        else:
            print("Output directory not specified. MIC data not exported.")

    def analyze_antibiotic(self, antibiotic, output_dir=None):
        """Perform complete analysis for one antibiotic"""
        print(f"\nAnalyzing {antibiotic}")
        print("=" * 50)
        
        # Interval-Censored Survival Analysis
        self.analyze_interval_censored(antibiotic, output_dir)
        
        # Plot Resistance Proportions
        self.plot_resistance_proportions(antibiotic, output_dir)
        
        # Plot MIC Distribution
        self.plot_mic_distribution(antibiotic, output_dir)
        
        # Export MIC Data
        self.export_mic_data(antibiotic, output_dir)
    
    def run_full_analysis(self, output_dir=None):
        """Run full analysis for all antibiotics"""
        for antibiotic in self.interpretive_standards.keys():
            self.analyze_antibiotic(antibiotic, output_dir)


# ------------- Example usage -------------
if __name__ == "__main__":
    # 1) Path to your MIC data CSV file
    data_path = r"C:\Users\BLehr\Documents\Paper_Backup\OHPH_documents\Hub_summary\meta_data_extra\group_useful\selected_contig\other_outputs\merge\for_mic\MIC_Exp_processing_v2.csv"
    #"C:\Users\lshml2\OneDrive - London School of Hygiene and Tropical Medicine\Documents\OHPH_documents\Hub_summmary\meta_data_extra\group_useful\selected_contig\other_outputs\merge\for_mic\MIC_Exp_processing.csv"
    
    # 2) Create output directory once
    output_dir = r"C:\Users\BLehr\Documents\Paper_Backup\OHPH_documents\Hub_summary\meta_data_extra\group_useful\selected_contig\other_outputs\merge\for_mic\distriubutionv2"
    #"C:\Users\lshml2\OneDrive - London School of Hygiene and Tropical Medicine\Documents\OHPH_documents\Hub_summmary\meta_data_extra\group_useful\selected_contig\other_outputs\merge\outputs\MIC\Country"
    os.makedirs(output_dir, exist_ok=True)
    
    analyzer = MICAnalyzer(data_path)

    font_sizes = {
        'title': 22, 
        'axis': 18, 
        'tick': 18, 
        'legend': 22,
        'stats': 18
    }

    # Single combined figure Aâ€“H
    analyzer.plot_combined_mic_distribution(
        font_sizes=font_sizes,
        output_dir=output_dir
    )
    
    # Split figure: Aâ€“D in one, Eâ€“H in another, with shading
    analyzer.plot_combined_mic_distribution_split(
        font_sizes=font_sizes,
        output_dir=output_dir,
        bar_spacing=3.0,
        cat_spacing=2.0
    )

    # Full analysis
    analyzer.run_full_analysis(output_dir)

#adjusted v5002
