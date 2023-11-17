# Eigenmodes of the Brain: Revisiting Connectomics and Geometry

This repository contains code and data supporting our recent manuscript titled "Eigenmodes of the Brain: Revisiting Connectomics and Geometry." Our study explores specific aspects discussed in the article ["Geometric Constraints on Human Brain Function" (Pang et al. 2023)](https://doi.org/10.1038/s41586-023-06098-1).

---

## Repository Contents

The repository includes Jupyter notebooks (located in the `notebooks` directory) that comprehensively describe the analytical procedures necessary to reproduce the results discussed in our manuscript. Additionally, several scripts (located in the `scripts` directory) utilized in our analyses are provided. Supplementary data required to replicate our findings is also included.

### Reproducible Notebooks

The `notebooks` directory contains three notebooks that collectively detail the steps to reproduce our results:

1. **Notebook 1: Alternative Brain Eigenmodes**: Scripts generating all brain eigenmodes utilized in this study.
2. **Notebook 2: Reconstruction Accuracy**: Scripts computing reconstruction accuracy from HCP's rest and task data.
3. **Notebook 3: Reproducible Results**: Scripts generating all figures presented in our manuscript.

### Accompanying Scripts

#### Automated Systematic Evaluations

We conducted a systematic analysis in our supplementary evaluations to assess the influence of various connectome mapping parameters and decisions on reconstruction accuracy. Bash and Python scripts for these evaluations (executed on a high-performance computing cluster) are available in the `scripts` directory.

#### Procrustes Analysis

As part of our supplementary evaluations, we conducted a Procrustes analysis using several Matlab scripts provided in the `scripts` directory.

### Provided Data

We've included minimal data within this repository to facilitate reproducibility of our findings:

- Our group average connectome (from the new tractography pipeline, without any smoothing) is available in the `data/connectomes` directory.
- Eigenmodes from our connectome reconstruction pipeline (using weighted, smoothed connectivity) are provided in the `data/eigenmodes` directory for easy access without the need to execute the scripts in Notebook 1.
- Additionally, several template files used in our analyses are included in the `data/templates` directory.

---
