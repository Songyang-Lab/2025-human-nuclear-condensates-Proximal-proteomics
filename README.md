# 2025-human-nuclear-condensates-Proximal-proteomics
This repository contains R scripts and validation/test files to reproduce the analyses and visualizations for the following figures in the manuscript:  
  - Fig 1d  
  - This script uses the GO_AllLists.csv file—downloaded from Metascape and containing GO analysis results of a mass spectrometry (MS) proteome list derived from a bait protein or nuclear condensate (NC)—to generate a barcode plot. GO terms with q-values ≤ 1×10⁻² were manually curated and grouped into distinct biological processes to highlight the functional characteristics of the bait or NC.

  - Fig 2a, Fig 2c  
  - This script calculates a Relevance-Associated (RA) matrix at the bait level. Marker proteins of 18 NCs were analyzed using PhastID. The normalized abundance values of enriched interactors (documented in NC_marker_matrix.txt) were used to compute pairwise RA scores between baits. The procedure consists of the following steps: Step 1: Abundance transformation: Reduces the influence of high-abundance variance and enhances sensitivity to low-abundance signals. Step 2: Correlation calculation: Computes Pearson correlation coefficients between bait pairs using transformed abundances. Negative correlations are set to zero. Step 3: Association normalization: Normalizes association strengths to the [0,1] range while preserving relative proportions using the formula: RA_ij = (r_ij/Σr_ik) * min(1/RA_ij), followed by a non-linear enhancement.

  - Fig 3a, Fig 3c, Extended Data Fig 3c
  - RA curves are used to predict protein localization and identify components of potential novel NC. The similarity between an unknown protein and NC markers is measured by Pearson correlation (r), along with its statistical significance (p-value).

  - NOVA map for Fig 4a, 4b, 5d
  - A similar RA-based approach is applied to over 300 published BioID datasets using the 18 NC markers. The resulting RA score matrix is visualized using UMAP to construct the NOVA map.
