"""
Geographic Map Visualization
============================
Study site map with sample points and regional coloring (Figure 1)

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
"""
Publication-Quality Map for AMR Campylobacter Paper (v3)
Features: Granular sample points, distinct regional coloring, global summary stats.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
import os

# --- OPTIONAL: Basemap ---
try:
    import contextily as ctx
    HAS_CTX = True
except ImportError:
    HAS_CTX = False
    print("Note: pip install contextily for basemap tiles")

# --- CONFIGURATION ---

# 1. Colors (User Requested)
# Vietnam: Greenish
# Bangladesh: Reddish
# Gujarat (India): Blueish
# Sri Lanka: Pink/Purplish
# Tamil Nadu: Orange/Gold (Distinct from others)
REGION_COLORS = {
    'Vietnam':    '#009E73', # Green
    'Bangladesh': '#D55E00', # Red
    'Gujarat':    '#0072B2', # Blue
    'Sri Lanka':  '#CC79A7', # Pink/Purple
    'Tamil Nadu': '#E69F00', # Orange/Gold
}

# 2. Region Mapping
REGION_MAPPING = {
    'Ahmedabad': 'Gujarat', 'Bharuch': 'Gujarat', 'Bhuj': 'Gujarat',
    'Godhra': 'Gujarat', 'Morbi': 'Gujarat', 'Panchmahal': 'Gujarat',
    'Rajkot': 'Gujarat', 'Surat': 'Gujarat', 'Vadodara': 'Gujarat',
    'Tamil Nadu': 'Tamil Nadu',
}

# 3. Text/Label Positions (Anchor Points)
# These are the "centers" of the regions
ANCHOR_POSITIONS = {
    'Gujarat':    {'lat': 22.5, 'lon': 71.5},
    'Tamil Nadu': {'lat': 11.0, 'lon': 78.5},
    'Bangladesh': {'lat': 24.0, 'lon': 90.0},
    'Sri Lanka':  {'lat': 7.5,  'lon': 80.8},
    'Vietnam':    {'lat': 21.0, 'lon': 105.8},
}

# 4. Annotation Box Offsets
# Carefully placed to push text into the oceans/empty space
# Format: (x_offset, y_offset, ha, va) relative to anchor
ANNOTATION_OFFSETS = {
    'Bangladesh': (2, -4.5, 'left', 'top'),   # Push South into Bay of Bengal
    'Vietnam':    (3, -2, 'left', 'top'),     # Push East/South into South China Sea
    'Gujarat':    (-3, -2, 'right', 'top'),   # Push West into Arabian Sea
    'Tamil Nadu': (3.5, 0, 'left', 'center'), # Push East into Bay of Bengal
    'Sri Lanka':  (3.5, -2, 'left', 'top'),   # Push South-East
}

# 5. Manual MIC Counts (if not in metadata)
MIC_COUNTS = {
    'Bangladesh': 105, 'Gujarat': 16, 'Tamil Nadu': 6,
    'Sri Lanka': 33, 'Vietnam': 183,
}

# --- HELPER FUNCTIONS ---

def add_jitter(arr, amount=0.08):
    """Adds slight random noise to coordinates to prevent perfect overlap"""
    return arr + np.random.uniform(-amount, amount, size=len(arr))

def draw_map_elements(ax, lon_min, lat_min, lat_max):
    """Draws Scale Bar and North Arrow"""
    # Scale bar (approx 500km) - Moved to avoid legend overlap
    bar_x, bar_y = lon_min + 5.0, lat_min + 2.0  # Moved right and up
    ax.plot([bar_x, bar_x + 4.5], [bar_y, bar_y], 'k-', lw=2, zorder=30)
    ax.text(bar_x + 2.25, bar_y + 0.3, '500 km', ha='center', va='bottom', 
            fontsize=10, fontweight='bold', zorder=30).set_path_effects([path_effects.withStroke(linewidth=2, foreground='white')])

    # North Arrow - Fixed overlapping issue
    arrow_x, arrow_y = lon_min + 5.0, lat_max - 3.0  # Moved right and down
    
    # Draw arrow first
    ax.annotate('', xy=(arrow_x, arrow_y + 1.5), xytext=(arrow_x, arrow_y),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=2, mutation_scale=20),
                zorder=30)
    
    # Then place "N" above the arrow (not at base)
    ax.text(arrow_x, arrow_y + 1.7, 'N', ha='center', va='bottom', 
            fontweight='bold', fontsize=14, zorder=30)

def load_data(gps_path, meta_path):
    """Merges GPS and Metadata"""
    # Load GPS
    gps_df = pd.read_csv(gps_path)
    gps_df.columns = [c.strip() for c in gps_df.columns]

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

    
    # Load Metadata
    meta_df = pd.read_excel(meta_path, header=1)
    # Adjust column names based on your file structure
    meta_df.columns = ['Strain_Name', 'Continent', 'Country', 'city', 'Site_description', 
                       'Poultry_Type', 'site_id', 'Date', 'Species_Name', 'ST'] + list(meta_df.columns[10:])
    
    # Merge
    merged = pd.merge(meta_df, gps_df, left_on='Strain_Name', right_on='Key', how='left')
    
    # Clean & Group
    merged['region'] = merged.apply(lambda x: REGION_MAPPING.get(x['city'], x['Country']) 
                                    if x['Country'] == 'India' else x['Country'], axis=1)
    
    # Standardize Terms
    merged['Site_description'] = merged['Site_description'].replace(
        {'Meat processing facility': 'Slaughter', 'Market': 'Market', 'Farm': 'Farm'})
    merged['Poultry_Type'] = merged['Poultry_Type'].replace(
        {'Hybrid': 'Hybrid', 'Coloured Feather Hybrid': 'CF Hybrid', 'Exotic broiler': 'Broiler'})

    return merged

def get_region_stats(df, region_name):
    """Calculates statistics for the annotation box"""
    subset = df[df['region'] == region_name]
    
    n_wgs = len(subset)
    n_mic = MIC_COUNTS.get(region_name, 0)
    
    # Diversity
    valid_sts = subset['ST'].dropna().astype(str)
    valid_sts = valid_sts[~valid_sts.isin(['-', 'nan', 'None', ''])]
    diversity_count = valid_sts.nunique()
    
    # Top STs
    if len(valid_sts) > 0:
        top_sts = valid_sts.value_counts().head(3)
        # Format: "ST-XX (n=Y)"
        top_st_str = "\n".join([f"ST-{st} ({n})" for st, n in top_sts.items()])
    else:
        top_st_str = "None"
        
    return {
        'wgs': n_wgs, 'mic': n_mic, 
        'diversity': diversity_count, 'top_sts': top_st_str,
        'lats': subset['latitude_x'].dropna().values,
        'lons': subset['longitude_x'].dropna().values
    }

def create_publication_map(gps_path, meta_path, output_dir):
    
    # 1. Prepare Data
    df = load_data(gps_path, meta_path)
    
    # Global Totals calculation
    total_wgs = len(df)
    total_mic = sum(MIC_COUNTS.values())
    
    # 2. Setup Figure
    fig, ax = plt.subplots(figsize=(16, 10), dpi=300) # Standard landscape ratio
    
    # Map Extents (South/SE Asia)
    ax.set_xlim(65, 112)
    ax.set_ylim(4, 28)
    
    # Basemap
    if HAS_CTX:
        try:
            # Positron is clean and grey/white, good for overlays
            ctx.add_basemap(ax, crs='EPSG:4326', 
                           source=ctx.providers.CartoDB.PositronNoLabels, 
                           zoom=6, alpha=0.9)
        except:
            ax.set_facecolor('#F4F6F7')
    else:
        ax.set_facecolor('#F4F6F7')

    # 3. Plotting Loop
    legend_handles = []
    
    for region, color in REGION_COLORS.items():
        if region not in df['region'].unique() and region != "Tamil Nadu": 
            # (Tamil Nadu check handles cases where it might be named differently in `region` col but key exists)
            continue
            
        stats = get_region_stats(df, region)
        pos = ANCHOR_POSITIONS.get(region, {'lat':0, 'lon':0})
        
        # A. The Scope Cloud (Individual Samples)
        # Smaller dots (s=25), semi-transparent to show density
        if len(stats['lats']) > 0:
            j_lats = add_jitter(stats['lats'], amount=0.12)
            j_lons = add_jitter(stats['lons'], amount=0.12)
            
            ax.scatter(j_lons, j_lats, 
                       c=color, s=25, alpha=0.6, 
                       edgecolors='white', linewidth=0.3, 
                       zorder=10)

        # B. The Anchor (Centroid)
        # Distinct ring to mark the region "center"
        ax.scatter(pos['lon'], pos['lat'], s=150, facecolors='none', edgecolors=color, linewidth=2, zorder=15)
        ax.scatter(pos['lon'], pos['lat'], s=40, c=color, zorder=15)

        # C. The Annotation Box (Data Card)
        ox, oy, ha, va = ANNOTATION_OFFSETS.get(region, (3, 0, 'left', 'center'))
        
        # Connection line
        ax.annotate("", xy=(pos['lon'], pos['lat']), xytext=(pos['lon']+ox, pos['lat']+oy),
                    arrowprops=dict(arrowstyle="-", color=color, lw=1.5, alpha=0.8), zorder=19)
        
        # Box Content
        box_text = (
            f"$\\mathbf{{{region}}}$\n"
            f"WGS: {stats['wgs']} | MIC: {stats['mic']}\n"
            f"ST Diversity: {stats['diversity']}\n"
            f"$\\it{{Dominant:}}$\n{stats['top_sts']}"
        )
        
        ax.text(pos['lon'] + ox, pos['lat'] + oy, box_text,
                ha=ha, va=va, fontsize=9, family='sans-serif',
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=color, lw=1.5, alpha=0.9),
                zorder=20)
        
        legend_handles.append(mpatches.Patch(color=color, label=region))

    # 4. Global Summary Box (Top Right Corner - usually empty in this map projection)
    summary_text = (
        f"$\\mathbf{{Study\\ Scope}}$\n"
        f"Total Sequenced: {total_wgs}\n"
        f"Total Phenotypic (MIC): {total_mic}\n"
        f"Species: $\\it{{C. jejuni}}$ & $\\it{{C. coli}}$"
    )
    
    # Placed in the top right (China area), usually empty or irrelevant for this specific dataset
    ax.text(111, 27, summary_text, ha='right', va='top', fontsize=11,
            bbox=dict(boxstyle="square,pad=0.6", fc="white", ec="black", lw=2),
            zorder=25)

    # 5. Main Title
    ax.set_title("Geographic Scope and Genetic Diversity of $\\it{Campylobacter}$ Isolates", 
                 fontsize=16, fontweight='bold', pad=15)

    # 6. Map Essentials
    draw_map_elements(ax, 65, 4, 28)
    
    # Legend (Bottom Left) - Moved further up to avoid scale bar
    ax.legend(handles=legend_handles, loc='lower left', 
              bbox_to_anchor=(0.02, 0.10),  # Changed from (0.02, 0.02) to (0.02, 0.10)
              title="Sampling Regions", title_fontsize=10, fontsize=9,
              frameon=True, fancybox=False, edgecolor='black')

    # Frame cleanup
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)

    # 7. Save
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "Figure1_Publication_Map_v3.png")
    
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.savefig(out_path.replace('.png', '.pdf'), bbox_inches='tight') # PDF for vector editing
    plt.savefig(out_path.replace('.png', '.tiff'), dpi=300, bbox_inches='tight') # TIFF for submission
    
    print(f"Map saved to: {out_path}")
    plt.close()

if __name__ == "__main__":
    # Update paths to your files
    gps_csv = r"E:\Dell_Windows_March_25\Paper_Backup\OHPH_documents\Hub_summary\WRITE_UP\Split_papers\Tables\Figures_to_send\V2\AfterV27\Samples_country_GPS.csv"
    metadata_xlsx = r"E:\Dell_Windows_March_25\Paper_Backup\OHPH_documents\Hub_summary\WRITE_UP\Split_papers\Tables\Figures_to_send\V2\AfterV27\Re_ OHPH papers\Table_S1_metadata.xlsx"
    output_dir = r"E:\Dell_Windows_March_25\Paper_Backup\OHPH_documents\Hub_summary\WRITE_UP\Split_papers\Tables\Figures_to_send\V2\AfterV27\Arg_map_output"
    
    create_publication_map(gps_csv, metadata_xlsx, output_dir)