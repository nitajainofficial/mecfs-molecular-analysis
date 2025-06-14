# Post-Infectious ME/CFS Molecular Signature Analysis

## Overview

This repository contains R code for analyzing gene expression differences between post-infectious Myalgic Encephalomyelitis/Chronic Fatigue Syndrome (ME/CFS) patients and healthy controls. The analysis supports research into the hypothesis that different pathogens target similar cellular pathways, explaining why diverse infections can lead to similar post-infectious syndromes.

## Key Findings

- **729 significantly altered genes** (3.27% of genome) in ME/CFS patients
- **Strong upregulation bias**: 89% of significant changes are increases
- **Key affected pathways**: Immune dysfunction, DNA repair, cellular stress response, neurological function
- **Potential biomarkers**: CXCR4, ATRX, VCAN, SEC63

## Dataset

**Source**: [GSE14577](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14577) from Gene Expression Omnibus (GEO)
- **Platform**: Affymetrix Human Genome U95Av2 Array
- **Samples**: 8 post-infectious ME/CFS patients vs 7 healthy controls
- **Publication**: Gow JW, Hagan S, Herzyk P, Cannon C, Behan PO, Chaudhuri A. "A gene signature for post-infectious chronic fatigue syndrome." BMC Med Genomics. 2009 Jun 25;2:38.

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
   git clone https://github.com/nitajainofficial/mecfs-molecular-analysis.git
   cd mecfs-molecular-analysis
   ```

2. **Run the analysis**:
   ```r
   source("mecfs_analysis.R")
   ```

3. **View results**:
   - Check console output for summary statistics
   - Open generated CSV files for detailed results
   - View `MECFS_volcano_plot.png` for visualization

## Output Files

| File | Description |
|------|-------------|
| `MECFS_vs_Control_complete_results.csv` | Full differential expression results |
| `MECFS_upregulated_genes.csv` | List of upregulated gene symbols |
| `MECFS_downregulated_genes.csv` | List of downregulated gene symbols |
| `MECFS_top20_significant_genes.csv` | Top 20 most significant genes with symbols |
| `MECFS_analysis_summary.csv` | Summary statistics |
| `MECFS_volcano_plot.png` | Volcano plot visualization |

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
Jain, Nita (2025). Post-Infectious ME/CFS Molecular Signature Analysis. 
GitHub: https://github.com/nitajainofficial/mecfs-molecular-analysis
```

**Original dataset citation**:
Gow JW, Hagan S, Herzyk P, Cannon C, Behan PO, Chaudhuri A. A gene signature for post-infectious chronic fatigue syndrome. BMC Med Genomics. 2009 Jun 25;2:38. doi: 10.1186/1755-8794-2-38. PMID: 19555476; PMCID: PMC2716361.

## License

MIT License - see LICENSE file for details.

## Contact

- **Author**: Nita Jain
- **Email**: nitajain@timelessbiosciences.com
- **Institution**: Timeless Biosciences
- **ORCID**: https://orcid.org/0000-0002-4197-937X

## Acknowledgments

- Gow et al. for making the original dataset publicly available
- Gene Expression Omnibus (GEO) for data hosting
- Bioconductor community for analysis tools

---

**Keywords**: ME/CFS, post-infectious syndrome, gene expression, microarray, PECT, pathogen effector convergence, biomarkers, molecular signature
