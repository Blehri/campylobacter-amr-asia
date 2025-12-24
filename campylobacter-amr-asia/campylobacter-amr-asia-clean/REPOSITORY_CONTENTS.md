# Repository Contents

Complete inventory of all files in this repository.

## Directory Structure

```
campylobacter-amr-asia/
│
├── README.md                           # Main repository documentation
├── SETUP_GUIDE.md                      # Step-by-step GitHub setup instructions
├── SCRIPTS_OVERVIEW.md                 # Detailed guide to all analysis scripts
├── REPOSITORY_CONTENTS.md              # This file
├── LICENSE                             # MIT License
├── requirements.txt                    # Python dependencies (pip)
├── environment.yml                     # Conda environment specification
├── .gitignore                          # Git ignore rules
│
├── config/
│   └── paths_template.yaml             # Configuration template
│
├── data/
│   └── README.md                       # Data requirements documentation
│
├── scripts/
│   ├── python/
│   │   ├── 01_arg_resistance_visualizations.py
│   │   ├── 02_sequence_type_analysis.py
│   │   ├── 03_mic_heatmap_analysis.py
│   │   ├── 04_mdr_prevalence_analysis.py
│   │   ├── 05_mic_distribution_plots.py
│   │   ├── 06_geographic_map_visualization.py
│   │   ├── 07_bird_type_comparisons.py
│   │   └── 08_site_type_comparisons.py
│   └── R/
│       └── README.md
│
└── outputs/
    └── .gitkeep                        # Preserves empty directory in Git
```

## File Descriptions

### Root Level Files

**README.md** (Main Documentation)
- Repository overview
- Installation instructions
- Usage examples
- Citation information
- Quick start guide

**SETUP_GUIDE.md** (GitHub Setup)
- Step-by-step GitHub repository creation
- Upload instructions (web and command line)
- Zenodo DOI setup
- Manuscript integration guidance
- Troubleshooting tips

**SCRIPTS_OVERVIEW.md** (Analysis Guide)
- Detailed description of each script
- Input/output specifications
- Usage examples
- Customization options
- Performance notes

**LICENSE** (MIT License)
- Open source license
- Usage permissions
- Liability disclaimers

**requirements.txt** (Python Dependencies)
- Core packages: pandas, numpy, scipy
- Statistical tools: statsmodels
- Visualization: matplotlib, seaborn, plotly
- Utilities: openpyxl, xlsxwriter
- Optional: surpyval, geopandas

**environment.yml** (Conda Environment)
- Complete conda environment specification
- Includes all dependencies
- Python version specification (≥3.9)

**.gitignore** (Git Ignore Rules)
- Excludes data files
- Excludes outputs
- Excludes Python cache
- Excludes IDE files

### Configuration

**config/paths_template.yaml**
- Template for local configuration
- Specifies data file paths
- Specifies output directories
- User copies to `paths.yaml` (not tracked by Git)

### Data

**data/README.md**
- Required data file specification
- Column descriptions
- Data format examples
- Privacy and sharing notes

### Python Scripts

**01_arg_resistance_visualizations.py** (~900 lines)
- Generates Figures 3-4
- Stacked bar charts
- By country and species
- Three category levels

**02_sequence_type_analysis.py** (~1269 lines)
- Generates Figure 5, Table S10
- ST-specific resistance patterns
- Statistical comparisons
- Identifies high-risk lineages

**03_mic_heatmap_analysis.py** (~477 lines)
- Generates Figure 6, Tables S2-S3
- Interactive MIC heatmaps
- ECOFF comparisons
- Phenotypic resistance visualization

**04_mdr_prevalence_analysis.py** (~300 lines)
- Generates Figure 7, Table S12
- MDR classification
- Prevalence analysis
- Summary statistics

**05_mic_distribution_plots.py** (~912 lines)
- MIC distribution analysis
- Distribution fitting
- Survival curves
- Statistical summaries

**06_geographic_map_visualization.py** (~279 lines)
- Generates Figure 1
- Study site mapping
- Regional visualization
- GPS coordinate plotting

**07_bird_type_comparisons.py** (~214 lines)
- Poultry type statistical tests
- Fisher's exact / Chi-square
- FDR correction
- Excel output

**08_site_type_comparisons.py** (~216 lines)
- Generates Table S7
- Site type comparisons
- Within-country analysis
- Terminology standardization

### R Scripts

**scripts/R/README.md**
- Placeholder for R analyses
- Instructions for adding R scripts
- Potential R implementations

### Outputs

**outputs/.gitkeep**
- Preserves directory in Git
- All generated files excluded by .gitignore

## Total Repository Size

**Without data/outputs**: ~2-3 MB
- Documentation: ~100 KB
- Python scripts: ~500 KB (combined)
- Configuration: ~5 KB

**With typical outputs**: ~50-100 MB
- Figures: ~20-30 MB
- Statistical files: ~5-10 MB
- Interactive visualizations: ~5-10 MB

## Version History

**v1.0.0** (Initial Release)
- All 8 analysis scripts
- Complete documentation
- Setup guides
- Ready for publication

## Next Steps After Download

1. Extract the ZIP file
2. Read `SETUP_GUIDE.md`
3. Follow GitHub upload instructions
4. Update placeholders with your information
5. Test scripts with your data
6. Create GitHub release
7. Obtain Zenodo DOI
8. Add to manuscript

## Support

For questions or issues:
- GitHub Issues: https://github.com/[YOUR_USERNAME]/campylobacter-amr-asia/issues
- Email: [your.email@institution.ac.uk]

---

**Created**: December 2024  
**For**: Microbial Genomics submission
**Study**: Campylobacter AMR in South and Southeast Asia
