# 📊 VarClustering - Advanced Variable Clustering in R

<div align="center">

![R Version](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/version-1.0.0-brightgreen)
![Build](https://img.shields.io/badge/build-passing-brightgreen)

*A comprehensive R package for clustering variables using hierarchical and partitioning methods*

[Installation](#-installation) • [Features](#-features) • [Quick Start](#-quick-start) • [Documentation](#-documentation)

</div>

---

## 📖 Table of Contents

- [Overview](#-overview)
- [Features](#-features)
- [Installation](#-installation)
- [Algorithms](#-algorithms)
  - [Hierarchical Clustering (HClustVar)](#1️⃣-hierarchical-clustering-hclustvar)
  - [K-means Variable Clustering (KmeansVariables)](#2️⃣-k-means-variable-clustering-kmeansvariables)
  - [Modality Clustering (ModCluster)](#3️⃣-modality-clustering-modcluster)
- [Interactive Shiny App](#-interactive-shiny-app)
- [Documentation](#-documentation)
- [Examples](#-examples)
- [Development](#-development)
- [Contributing](#-contributing)
- [License](#-license)
- [References](#-references)

---

## 🎯 Overview

**VarClustering** is an advanced R package designed for clustering **variables** rather than observations. It provides multiple algorithms to identify groups of correlated or associated variables, supporting:

- ✅ **Quantitative variables** (continuous data)
- ✅ **Qualitative variables** (categorical data)  
- ✅ **Mixed datasets** (quantitative + qualitative)

The package is ideal for:
- 🔬 **Dimension reduction** before regression or classification
- 📈 **Data exploration** to understand variable relationships
- 🧬 **Feature engineering** in machine learning pipelines

---

## ✨ Features

### 🔧 Core Functionality

| Feature | HClustVar | KmeansVariables | ModCluster |
|---------|:---------:|:---------------:|:----------:|
| **Quantitative variables** | ✅ | ✅ | ❌ |
| **Qualitative variables** | ✅ | ❌ | ✅ |
| **Mixed data** | ✅ | ❌ | ❌ |
| **Automatic type detection** | ✅ | ✅ | ✅ |
| **Multiple distance metrics** | ✅ | ✅ | ✅ |
| **Prediction on new observations** | ✅ | ✅ | ✅ |
| **Interactive visualizations** | ✅ | ✅ | ✅ |

### 📊 Visualization Tools

- 🌳 **Dendrograms** (hierarchical structure)
- 📉 **Elbow plots** (optimal K selection)
- 🎨 **Silhouette diagrams** (cluster quality)
- 🗺️ **MDS/PCA projections** (2D variable maps)
- 📈 **Contribution plots** (key variables per cluster)

### 🎮 Interactive Interface

- 🖥️ **Shiny web application** for exploratory analysis
- 📤 **Export results** (CSV, Excel)
- 🎯 **Real-time parameter tuning**
- 📊 **Dynamic visualizations**

---

## 📦 Installation

### From GitHub (Recommended)

This package is currently available on GitHub and can be installed using `devtools` or `remotes`:

```r
# Option 1: Using devtools
if (!require("devtools")) install.packages("devtools")
devtools::install_github("4xel-C/SISE_var_clustering")

# Option 2: Using remotes (lighter dependency)
if (!require("remotes")) install.packages("remotes")
remotes::install_github("4xel-C/SISE_var_clustering")
```

### Install with Vignettes

To build and install the package with vignettes (recommended for learning):

```r
devtools::install_github("4xel-C/SISE_var_clustering", 
                         build_vignettes = TRUE,
                         dependencies = TRUE)
```

### Development Version

For the latest development features:

```r
devtools::install_github("4xel-C/SISE_var_clustering", 
                         ref = "dev")
```

### Dependencies

The package automatically installs required dependencies. If needed manually:

```r
# Core dependencies (required)
install.packages(c(
  "R6",           # Object-oriented programming
  "FactoMineR",   # Factorial analysis (PCA, MCA, FAMD)
  "shiny",        # Interactive web application
  "DT"            # Interactive tables
))

# Optional dependencies (recommended)
install.packages(c(
  "readxl",       # Excel file support
  "ggplot2"      # Enhanced visualizations
))
```

### System Requirements

- **R version**: ≥ 4.0.0
- **Operating System**: Windows, macOS, Linux

---

## 🧮 Algorithms

### 1️⃣ Hierarchical Clustering (`HClustVar`)

**Agglomerative Hierarchical Clustering** adapted for variables.

#### 📌 Key Characteristics

- **Distance metrics**: 
  - `"rsquare"` - R² based (default, treats ±correlations equally)
  - `"r"` - Signed correlation (distinguishes positive/negative)
- **Linkage methods**: Ward.D, Ward.D2, Single, Complete, Average, McQuitty, Median, Centroid
- **Variable types**: Quantitative, Qualitative, Mixed

#### 💡 When to use

- ✅ You want to explore hierarchical variable relationships
- ✅ The number of clusters is unknown (use dendrogram)
- ✅ You need stable, reproducible results
- ✅ You have mixed data types
- ✅ Dataset not too large (due to HClust complexity)

#### 📝 Basic Example

```r
library(VarClustering)

# Load sample data
data(autos)

# Initialize model
hc <- HClustVar$new(
  vartype = "auto",
  dist.metric = "rsquare",
  cah.method = "ward.D"
)

# Fit on data
hc$fit(autos)

# Visualize dendrogram
hc$plot_dendrogram()

# Select optimal K
hc$plot_agg_levels()

# Cut tree
hc$cut_tree(k = 4)

# Analyze results
results <- hc$summary()
hc$plot_silhouette()
hc$mds_projection()
```

---

### 2️⃣ K-means Variable Clustering (`KmeansVariables`)

**Partitioning algorithm** optimized for variable clustering via iterative reallocation.

#### 📌 Key Characteristics

- **Distance metrics**:
  - `"r_squared"` - Groups variables with high |correlation| (default)
  - `"r_signed"` - Groups only positively correlated variables
- **Automatic K selection**: `n_clusters = "auto"` (via silhouette/Calinski-Harabasz). ⚠️ Compute time can be long because of the numerous iterations needed for this method!
- **Multiple initializations**: `n_init = 10` (avoids local minima)
- **Variable types**: Quantitative only

#### 💡 When to use

- ✅ You know the number of clusters (or use automatic selection)
- ✅ You want fast computation on large datasets
- ✅ You need to identify strongly correlated variable groups
- ✅ Your data is purely quantitative

#### 📝 Basic Example

```r
library(VarClustering)

# Load sample data
data(autos)

# Option 1: Manual K selection
km <- KmeansVariables$new(
  n_clusters = 5,
  n_init = 2,
  random_state = 42
)
km$fit(autos)

# Option 2: Automatic K selection
km_auto <- KmeansVariables$new(
  n_clusters = "auto",
  k_range = 2:10,
  selection_method = "silhouette"
)
km_auto$fit(autos)

# Analyze results
km$summary(print_summary = TRUE)
km$plot_elbow(k_range = 2:10)
km$plot_projection(show_labels = TRUE)
km$plot_contributions(top_n = 5)
```

---

### 3️⃣ Modality Clustering (`ModCluster`)

**Hierarchical clustering of modalities** (categories) using Dice distance for categorical data analysis.

#### 📌 Key Characteristics

- **Distance metric**: Dice coefficient (measures co-occurrence of categories)
- **Linkage methods**: Complete, Average, Single
- **Clustering target**: Modalities (categories) rather than variables
- **Visualization**: MCA projection for cluster interpretation
- **Variable types**: Qualitative only

#### 💡 When to use

- ✅ Pure categorical datasets (surveys, questionnaires)
- ✅ You want to cluster categories, not variables
- ✅ You need to identify which modalities co-occur
- ✅ You want to predict cluster membership for new observations

#### 📝 Basic Example

```r
library(VarClustering)

# Load sample data
data(vote)

# Initialize model
model <- ModCluster$new(
  method = "hclust",
  n_dimensions = 2,
  hclust_method = "average"
)

# Fit on data
model$fit(vote)

# Determine optimal K
model$plot_heights()

# Cut tree
model$cut_tree(3)

# Analyze results
summary(model)

# Visualizations
model$plot_dendrogram()
model$plot_mca()
```

---

## 🎮 Interactive Shiny App

Launch the interactive web interface for exploratory analysis:

```r
library(VarClustering)

# Launch Shiny app
varclust_gui()
```
<img width="1430" height="744" alt="Capture d&#39;écran 2025-11-22 174456" src="https://github.com/user-attachments/assets/23d951ba-7b77-4d53-b417-278b8983b0c2" />

### App Features

- 📁 **Data import**: CSV, Excel
- 🔧 **Variable configuration**: Set types (quant/qual), roles (active/illustrative)
- 🎯 **Algorithm selection**: HClustVar, KmeansVariables, Modality clustering
- 📊 **Live visualizations**: All plots update in real-time
- 📥 **Prediction**: Predict illustrative variables

---

## 📚 Documentation

### Vignettes

Access comprehensive tutorials:

```r
# View all available vignettes
browseVignettes("VarClustering_Package")
```

### Function Help

```r
# Class documentation
?HClustVar
?KmeansVariables
?ModCluster

# Method documentation
?HClustVar$fit
?KmeansVariables$predict
?ModCluster$plot_mca

# Dataset documentation
?autos
?players
?canines
?vote
```
---

## 🛠️ Development

This package was developed using `devtools` workflow for R package development.

### Local Development Setup

```r
# Clone the repository
git clone https://github.com/4xel-C/VarClustering.git
cd VarClustering

# Open in RStudio (or your favorite R IDE)
# Install development dependencies
devtools::install_dev_deps()

# Load package for interactive development
devtools::load_all()

# Run tests
devtools::test()

# Build documentation
devtools::document()

# Build vignettes
devtools::build_vignettes()

# Install locally
devtools::install()
```

### Package Structure

```
VarClustering/
├── R/                      # R source code
│   ├── ClusteringBase.R   # Parent class
│   ├── HClustVar.R        # Hierarchical clustering
│   ├── KmeansVariables.R  # K-means algorithm
│   ├── ModClust.R         # Modality clustering
│   ├── app_ui.R           # Shiny UI
│   ├── app_server.R       # Shiny server
│   └── data.R             # Dataset documentation
├── data/                   # Built-in datasets
│   ├── autos.rda
│   ├── players.rda
│   ├── canines.rda
│   └── vote.rda
├── man/                    # Documentation (auto-generated)
├── vignettes/             # Tutorials
│   └── Variable_clustering.Rmd
├── tests/                 # Unit tests
│   └── testthat/
├── notebooks/             # Notebooks of complexity testing  
├── doc                    # Vignette and application guide
├── DESCRIPTION           # Package metadata
├── NAMESPACE            # Exports (auto-generated)
└── README.md           # This file
```

### Development Tools Used

- **devtools**: Package development workflow
- **roxygen2**: Documentation generation
- **testthat**: Unit testing
- **usethis**: Package setup utilities
- **R6**: Object-oriented programming

### Building from Source

```r
# Build source package
devtools::build()

# Build binary package
devtools::build(binary = TRUE)
```
---

## 🤝 Contributing

Contributions are welcome! Please follow these guidelines:

### How to Contribute

1. 🍴 **Fork** the repository
2. 🌿 **Create** a feature branch
   ```bash
   git checkout -b feature/amazing-feature
   ```
3. ✍️ **Commit** your changes
   ```bash
   git commit -m 'Add amazing feature'
   ```
4. 📤 **Push** to the branch
   ```bash
   git push origin feature/amazing-feature
   ```
5. 🔀 **Open** a Pull Request

### Development Guidelines

- ✅ Follow [R package development best practices](https://r-pkgs.org/)
- ✅ Write unit tests for new features
- ✅ Document all functions using roxygen2
- ✅ Run `devtools::check()` before submitting PR

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) 2025 VarClustering Contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
```

---

## 📖 References

1. **Rakotomalala, R. (2022).** *Clustering de variables.* Support de cours, Université Lyon 2. https://eric.univ-lyon2.fr/ricco/cours/slides/classification_de_variables.pdf
  
2. **Rakotomalala, R. (2008).** *Classification de variables.* Tutoriels Tanagra pour le Data Mining. http://tutoriels-data-mining.blogspot.com/2008/03/classification-de-variables.html

4. **Vigneau, E., & Qannari, E. M. (2003).** *Clustering of variables around latent components.* Communications in Statistics-Simulation and Computation, 32(4), 1131-1150.  
   [DOI: 10.1081/SAC-120023882](https://doi.org/10.1081/SAC-120023882)

5. **Chavent, M., Kuentz-Simonet, V., Labenne, A., & Saracco, J. (2012).** *ClustOfVar: An R Package for the Clustering of Variables.* Journal of Statistical Software, 50(13), 1-16.  
   [DOI: 10.18637/jss.v050.i13](https://doi.org/10.18637/jss.v050.i13)

6. **Ward, J. H. (1963).** *Hierarchical Grouping to Optimize an Objective Function.* Journal of the American Statistical Association, 58(301), 236-244.

7. **Lebart, L., Morineau, A., & Piron, M. (2006).** *Statistique exploratoire multidimensionnelle : Visualisation et inférences en fouille de données (4e éd.).* Dunod.


---

<div align="center">

**⭐ If you find this package useful, please consider giving it a star on GitHub! ⭐**

Made with ❤️ using [devtools](https://devtools.r-lib.org/) and [R6](https://r6.r-lib.org/)

[⬆ Back to top](#-VarClustering---advanced-variable-clustering-in-r)

</div>
