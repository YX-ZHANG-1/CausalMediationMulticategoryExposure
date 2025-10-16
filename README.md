# Causal mediation analysis of the effect of depression on Alzheimer's disease risk in older adults
<p style="font-size:18px;">
  
# Description
  <p style="font-size:18px;">
This repository provides a reproducible <b>R code</b> for preprocessing DNA methylation data, screening ultra-high-dimensional CpG sites, and conducting causal mediation analysis with multi-category exposures.  
The code is designed for large-scale ADNI study and can be adapted to other datasets.
</p>

# R Scripts

- `preprocess_normalize_DNAm.R`: Preprocesses and normalizes DNA methylation data.
- `ballcor_screening_CpGs.R`: Screens CpGs using ball correlation.
- `medDML_multicategory_exposure.R`: Implements double machine learning (DML) to estimate total, direct, and indirect effects with multi-category exposures, using cross-fitting.
- `causal_mediation_analysis.R`: Runs full causal mediation analysis and summarizes results.
