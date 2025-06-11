# Post-Infectious CFS Molecular Signature Analysis

## Overview

This repository contains R code for analyzing gene expression differences between post-infectious Chronic Fatigue Syndrome (CFS) patients and healthy controls. The analysis supports research into the **Pathogen Effector Convergence Theory (PECT)** - the hypothesis that different pathogens target similar cellular pathways, explaining why diverse infections can lead to similar post-infectious syndromes.

## Key Findings

- **729 significantly altered genes** (3.27% of genome) in CFS patients
- **Strong upregulation bias**: 89% of significant changes are increases
- **Key affected pathways**: Immune dysfunction, DNA repair, cellular stress response, neurological function
- **Potential biomarkers**: CXCR4, ATRX, VCAN, SEC63

## Dataset

**Source**: [GSE14577](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14577) from Gene Expression Omnibus (GEO)
- **Platform**: Affymetrix Human Genome U95Av2 Array
- **Samples**: 8 post-infectious CFS patients vs 7 healthy controls
- **Publication**: Light AR, White AT, Hughen RW, Light KC. "Moderate exercise increases expression for sensory, adrenergic, and immune genes in chronic fatigue syndrome patients but not in normal subjects." J Pain. 2009;10(10):1099-112.

## Requirements

### R Packages
```r
# Bioconductor packages
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma"))

# CRAN packages  
install.packages("ggplot2")
```

### System Requirements
- R (≥ 4.0.0)
- Internet connection for data download
- ~2GB free disk space

## Quick Start

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/cfs-molecular-analysis.git
   cd cfs-molecular-analysis
   ```

2. **Run the analysis**:
   ```r
   source("cfs_analysis.R")
   ```

3. **View results**:
   - Check console output for summary statistics
   - Open generated CSV files for detailed results
   - View `CFS_volcano_plot.png` for visualization

## Output Files

| File | Description |
|------|-------------|
| `CFS_vs_Control_complete_results.csv` | Full differential expression results |
| `CFS_upregulated_genes.csv` | List of upregulated gene symbols |
| `CFS_downregulated_genes.csv` | List of downregulated gene symbols |
| `CFS_top20_significant_genes.csv` | Top 20 most significant genes with symbols |
| `CFS_analysis_summary.csv` | Summary statistics |
| `CFS_volcano_plot.png` | Volcano plot visualization |

## Key Results

### Top Upregulated Genes
1. **ASAP1** - Cell migration and invasion
2. **ATRX** - Chromatin remodeling, DNA repair
3. **VCAN** - Extracellular matrix, inflammation
4. **CXCR4** - Immune cell trafficking
5. **SEC63** - ER stress response

### Top Downregulated Genes
1. **SLC10A1** - Bile acid transport
2. **SCN10A** - Sodium channel, pain sensation
3. **OPHN1** - Neuronal development
4. **CHRNB2** - Neurotransmitter receptor

### Biological Pathways Affected
- **Immune system dysfunction** (CXCR4, ZMIZ1)
- **DNA repair defects** (ATRX, UBR5)
- **Cellular stress response** (SEC63, RNF10)
- **Neurological dysfunction** (SCN10A, OPHN1, CHRNB2)
- **RNA processing defects** (SRSF4, SRRM1)

## Applications for PECT Research

This gene signature can be compared with other post-infectious conditions to identify:

1. **Common molecular targets** across different pathogens
2. **Convergent biological pathways** 
3. **Potential therapeutic targets**
4. **Diagnostic biomarkers**

### Suggested Comparisons
- COVID-19 infection datasets
- Post-viral fatigue studies
- Long COVID molecular signatures
- Other post-infectious syndromes

## Limitations

- Small sample size (n=15 total)
- Single time point analysis
- Microarray technology (limited to known sequences)
- No clinical severity correlation

## Future Directions

1. **Validation studies** with larger cohorts
2. **RNA-seq analysis** for broader coverage
3. **Longitudinal tracking** of gene expression changes
4. **Single-cell analysis** to identify affected cell types
5. **Functional validation** of key pathways
6. **Drug screening** against identified targets

## Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## Citation

If you use this code or findings, please cite:

```
[Your Name] (2025). Post-Infectious CFS Molecular Signature Analysis. 
GitHub: https://github.com/yourusername/cfs-molecular-analysis
```

**Original dataset citation**:
Light AR, White AT, Hughen RW, Light KC. Moderate exercise increases expression for sensory, adrenergic, and immune genes in chronic fatigue syndrome patients but not in normal subjects. J Pain. 2009;10(10):1099-112.

## License

MIT License - see LICENSE file for details.

## Contact

- **Author**: [Your Name]
- **Email**: [your.email@domain.com]
- **Institution**: [Your Institution]
- **ORCID**: [Your ORCID ID]

## Acknowledgments

- Light et al. for making the original dataset publicly available
- Gene Expression Omnibus (GEO) for data hosting
- Bioconductor community for analysis tools

---

**Keywords**: Chronic Fatigue Syndrome, CFS, ME/CFS, post-infectious, gene expression, microarray, PECT, pathogen effector convergence, biomarkers, molecular signature