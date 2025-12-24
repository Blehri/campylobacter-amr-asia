# Analysis Scripts Overview

Complete guide to all analysis scripts in this repository.

---

## Quick Reference

| # | Script | Figures/Tables | Description |
|---|--------|----------------|-------------|
| 01 | arg_resistance_visualizations.py | Fig 3-4 | ARG prevalence charts |
| 02 | sequence_type_analysis.py | Fig 5, Table S10 | ST-specific patterns |
| 03 | mic_heatmap_analysis.py | Fig 6, Tables S2-S3 | MIC heatmaps |
| 04 | mdr_prevalence_analysis.py | Fig 7, Table S12 | MDR analysis |
| 05 | mic_distribution_plots.py | Supplementary | MIC distributions |
| 06 | geographic_map_visualization.py | Fig 1 | Study site map |
| 07 | bird_type_comparisons.py | - | Poultry type stats |
| 08 | site_type_comparisons.py | Table S7 | Site type stats |

---

## Detailed Descriptions

### 01_arg_resistance_visualizations.py
**Purpose**: Generate stacked bar charts showing antimicrobial resistance gene (ARG) prevalence

**Outputs**:
- Figures 3-4 in main manuscript
- By country (for each species)
- By species (for each country)
- Three category levels: Drug Classes, ARG Classes, Individual ARGs

**Key Features**:
- Proportional stacking (100% scale)
- Custom color schemes for each category
- Sample sizes displayed
- Publication-quality at 300 DPI

**Usage**:
```python
# Update paths in script header
DATA_FILE = "data/My_merged_output_V5.csv"
OUTPUT_DIR = "outputs"

# Run
python scripts/python/01_arg_resistance_visualizations.py
```

---

### 02_sequence_type_analysis.py
**Purpose**: Analyze resistance patterns by sequence type (ST)

**Outputs**:
- Figure 5: ST-specific resistance distributions
- Table S10: ST comparisons by country
- Statistical tests (Fisher's exact, Chi-square with FDR)
- Excel files with multiple sheets

**Key Features**:
- Minimum 5 samples per ST group
- Separate analysis for C. jejuni and C. coli
- Identifies high-risk STs (e.g., ST-3236, ST-829)
- Country and species comparisons

**Usage**:
```python
python scripts/python/02_sequence_type_analysis.py
```

**Note**: Large script (~1269 lines) - may take several minutes to run

---

### 03_mic_heatmap_analysis.py
**Purpose**: Create interactive heatmaps of MIC (Minimum Inhibitory Concentration) testing results

**Outputs**:
- Figure 6: Interactive MIC heatmap (HTML)
- Table S2: ECOFF breakpoints
- Table S3: MIC concentration ranges
- Excel summary with statistics

**Key Features**:
- Color-coded by resistance level
- ECOFF comparison
- Interactive Plotly visualization
- Fisher's exact tests for comparisons

**Usage**:
```python
python scripts/python/03_mic_heatmap_analysis.py
```

**Output Format**: HTML file (open in browser), Excel spreadsheet

---

### 04_mdr_prevalence_analysis.py
**Purpose**: Analyze multi-drug resistance (MDR) patterns

**Outputs**:
- Figure 7: MDR prevalence by country and species
- Table S12: MDR classification summary
- JSON summary statistics

**Key Features**:
- MDR classification (0 to 5+ drug classes)
- High-MDR identification (≥5 classes)
- Prevalence calculations
- Visual breakdown by resistance level

**Usage**:
```python
python scripts/python/04_mdr_prevalence_analysis.py
```

---

### 05_mic_distribution_plots.py
**Purpose**: Generate distribution plots for MIC values across antibiotics

**Outputs**:
- Distribution plots for each antibiotic
- Kaplan-Meier style survival curves
- Statistical summaries

**Key Features**:
- Uses surpyval library for distribution fitting
- Censored data handling
- Multiple antibiotic comparisons

**Usage**:
```python
python scripts/python/05_mic_distribution_plots.py
```

**Dependencies**: Requires `surpyval` package (optional, listed in requirements.txt)

---

### 06_geographic_map_visualization.py
**Purpose**: Create publication-quality geographic map showing study sites

**Outputs**:
- Figure 1: Study site map with sample points
- Regional coloring
- Summary statistics overlay

**Key Features**:
- Granular sample points
- Country-specific coloring
- GPS coordinate plotting
- Professional cartographic design

**Usage**:
```python
python scripts/python/06_geographic_map_visualization.py
```

**Data Requirements**: Latitude and Longitude columns in dataset

---

### 07_bird_type_comparisons.py
**Purpose**: Statistical comparisons between different poultry types (chicken, duck, etc.)

**Outputs**:
- Excel file with statistical test results
- Across-country comparisons by bird type
- Within-country comparisons by bird type
- FDR-corrected p-values

**Key Features**:
- Fisher's exact or Chi-square tests (automatic selection)
- Multiple testing correction (FDR)
- Detailed contingency tables

**Usage**:
```python
python scripts/python/07_bird_type_comparisons.py
```

---

### 08_site_type_comparisons.py
**Purpose**: Statistical comparisons between production sites (farms, markets, slaughtering facilities)

**Outputs**:
- **Table S7**: Within-country comparisons by site type
- Excel file with all statistical tests
- FDR-corrected p-values

**Key Features**:
- Standardized terminology (Slaughtering facilities, Live bird markets, Farms)
- Within-country pairwise comparisons
- Species-stratified analysis
- Comprehensive statistical testing

**Usage**:
```python
python scripts/python/08_site_type_comparisons.py
```

**Important**: This script automatically renames site types to match published terminology

---

## Data Requirements

All scripts require a CSV file with these columns:

### Essential Metadata
- `Country`: Sample country
- `Species_Name`: Campylobacter species
- `Site_description`: Collection site type
- `ST`: Sequence type

### ARG Data (Binary 0/1 columns)
- Individual genes: `tet(O)`, `ermB`, `cmeB`, etc.
- Gene classes: `APH(2'')`, `OXA beta-lactamase`, etc.
- Drug classes: `AminoAB`, `TetraAB`, `MacroAB`, etc.

### MIC Data (Numeric columns, optional)
- `CIP_MIC`, `ERY_MIC`, `TET_MIC`, `AZM_MIC`, etc.

### Geographic Data (for map, optional)
- `Latitude`, `Longitude`

### Additional
- `Poultry_type`: For bird type comparisons

---

## Running All Scripts

To generate all outputs:

```bash
#!/bin/bash
# Run all analyses in order

echo "Starting Campylobacter AMR analysis pipeline..."

# Visualizations
python scripts/python/01_arg_resistance_visualizations.py
python scripts/python/06_geographic_map_visualization.py

# ST Analysis (may take several minutes)
python scripts/python/02_sequence_type_analysis.py

# MIC Analysis
python scripts/python/03_mic_heatmap_analysis.py
python scripts/python/05_mic_distribution_plots.py

# MDR Analysis
python scripts/python/04_mdr_prevalence_analysis.py

# Statistical Comparisons
python scripts/python/07_bird_type_comparisons.py
python scripts/python/08_site_type_comparisons.py

echo "Analysis complete! Check outputs/ directory"
```

---

## Customization

### Modifying Color Schemes

Edit the `custom_colors` dictionary in each visualization script:

```python
custom_colors = {
    'Drug Classes': {
        'AminoAB': '#FFA07A',    # Change this color
        'TetraAB': '#FF0000',     # Or this one
        # ...
    }
}
```

### Adding New Comparisons

In statistical scripts (07, 08):

1. Locate the `antibiotic_categories` dictionary
2. Add new antibiotics/genes to relevant category
3. Re-run the script

### Changing Output Format

Modify output functions to export as PDF instead of PNG:

```python
# Change from:
plt.savefig(filepath, dpi=300, bbox_inches='tight')

# To:
plt.savefig(filepath.replace('.png', '.pdf'), format='pdf', bbox_inches='tight')
```

---

## Troubleshooting

### Common Issues

**Issue**: "FileNotFoundError: data/My_merged_output_V5.csv"
- **Solution**: Place your data file in the `data/` directory OR update `DATA_FILE` path in script

**Issue**: "KeyError: 'Site_description'"
- **Solution**: Check your data file has this column, or adjust column name in script

**Issue**: "No module named 'plotly'"
- **Solution**: Install dependencies: `pip install -r requirements.txt`

**Issue**: Scripts run but no outputs
- **Solution**: Check `outputs/` directory exists and is writable

**Issue**: "TypeError: unsupported operand type(s)"
- **Solution**: Verify ARG columns contain binary (0/1) values, not strings

---

## Performance Notes

**Fast** (<1 minute):
- 01, 03, 04, 05, 06, 07, 08

**Moderate** (1-5 minutes):
- 02 (depends on number of STs)

**Memory Usage**:
- Most scripts: <1 GB RAM
- Script 02: Up to 2 GB with large datasets

---

## Output File Sizes

Typical output sizes (for 609 isolates):

- PNG figures: 0.5-2 MB each
- Excel statistical files: 100-500 KB each
- HTML heatmaps: 500 KB - 1 MB
- Total outputs: ~50-100 MB

---

## Citation

When using these scripts, cite both the repository and the publication:

```bibtex
@software{lehri2025_scripts,
  author = {Lehri, Burhan},
  title = {Campylobacter AMR Analysis Scripts},
  year = {2025},
  url = {https://github.com/[YOUR_USERNAME]/campylobacter-amr-asia}
}
```

---

**Last Updated**: December 2024  
**Version**: 1.0.0
