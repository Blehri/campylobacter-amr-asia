# Campylobacter AMR Analysis - South and Southeast Asia

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Analysis scripts for: **"Geographic, species and lineage driven antimicrobial resistance patterns in poultry associated *Campylobacter* across South and Southeast Asia"**

*NA* (NA) | Burhan Lehri, et al.

---

## Quick Start

```bash
# Clone repository
git clone https://github.com/Blehri/campylobacter-amr-asia.git
cd campylobacter-amr-asia

# Install dependencies
pip install -r requirements.txt

# Place your data file in data/My_merged_output_V5.csv

# Run analyses
python scripts/python/01_arg_resistance_visualizations.py
python scripts/python/02_sequence_type_analysis.py
# ... (see Usage section for all scripts)
```

## Repository Contents

- **8 Python scripts** for comprehensive AMR analysis
- **Publication-quality visualizations** (Figures 1, 3-7)
- **Statistical comparisons** with FDR correction (Tables S5-S7)
- **Interactive MIC heatmaps**
- **MDR prevalence analysis**

## Installation

### Using pip:
```bash
pip install -r requirements.txt
```

### Using conda:
```bash
conda env create -f environment.yml
conda activate campy-amr
```

## Usage

See detailed instructions in README.md after cloning.

## Citation

```bibtex
@article{lehri2025campylobacter,
  title={Geographic, species and lineage driven antimicrobial resistance patterns},
  author={Lehri, Burhan and et al.},
  journal={NA},
  year={NA}
}
```

## License

MIT License - See LICENSE file

## Contact

**Burhan Lehri** | b.lehri@hotmail.co.uk
