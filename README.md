# 📊 VarClustering - Advanced Variable Clustering in R

<div align="center">

![R Version](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-in%20development-yellow)
![Build](https://img.shields.io/badge/build-passing-brightgreen)

*A comprehensive R package for clustering variables using hierarchical and partitioning methods*

[Installation](#-installation) • [Features](#-features) • [Quick Start](#-quick-start) • [Documentation](#-documentation) • [Examples](#-examples)

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
- [Quick Start](#-quick-start)
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
- 🔥 **Heatmaps** (correlation matrices)

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
devtools::install_github("4xel-C/VarClustering")

# Option 2: Using remotes (lighter dependency)
if (!require("remotes")) install.packages("remotes")
remotes::install_github("4xel-C/VarClustering")
```

### Install with Vignettes

To build and install the package with vignettes (recommended for learning):

```r
devtools::install_github("4xel-C/VarClustering", 
                         build_vignettes = TRUE,
                         dependencies = TRUE)
```

### Development Version

For the latest development features:

```r
devtools::install_github("4xel-C/VarClustering", 
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
  "ggplot2",      # Enhanced visualizations
  "plotly"        # Interactive plots
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
- **Automatic K selection**: `n_clusters = "auto"` (via silhouette/Calinski-Harabasz)
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
data(players)

# Option 1: Manual K selection
km <- KmeansVariables$new(
  n_clusters = 5,
  n_init = 20,
  random_state = 42
)
km$fit(players)

# Option 2: Automatic K selection
km_auto <- KmeansVariables$new(
  n_clusters = "auto",
  k_range = 2:10,
  selection_method = "silhouette"
)
km_auto$fit(players)

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
mc <- ModCluster$new(
  method = "hclust",
  n_dimensions = 5,
  hclust_method = "ward.D2"
)

# Fit on data
mc$fit(vote)

# Visualize dendrogram
mc$plot_dendrogram()

# Determine optimal K
mc$plot_heights()

# Cut tree
mc$cut_tree(k = 4)

# Analyze results
results <- mc$summary()
mc$plot_silhouette()
mc$plot_mca(axes = c(1, 2))

# Predict new observations
predictions <- mc$predict(new_data)
```

---

## 🚀 Quick Start

### 30-Second Example

```r
library(VarClustering)

# Load data
data(autos)

# Cluster variables
hc <- HClustVar$new()
hc$fit(autos)
hc$cut_tree(k = 3)

# View results
hc$summary()
```

### 5-Minute Tutorial

```r
# 1. Load package and data
library(VarClustering)
data(players)  # FIFA player stats

# 2. Hierarchical clustering
hc <- HClustVar$new(vartype = "quant", dist.metric = "rsquare")
hc$fit(players)
hc$plot_dendrogram()
hc$plot_agg_levels()  # Choose optimal K
hc$cut_tree(k = 5)

# 3. Analyze clusters
summary_data <- hc$summary()
print(summary_data$clust_summary)
print(summary_data$clust_members)

# 4. Visualize
hc$plot_silhouette()
hc$mds_projection()

# 5. Compare with K-means
km <- KmeansVariables$new(n_clusters = 5)
km$fit(players)
km$plot_projection()
```

---

## 🎮 Interactive Shiny App

Launch the interactive web interface for exploratory analysis:

```r
library(VarClustering)

# Launch Shiny app
varclust_gui()
```

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
browseVignettes("VarClustering")

# Introduction to variable clustering
vignette("introduction", package = "VarClustering")

# Hierarchical clustering guide
vignette("hierarchical-clustering", package = "VarClustering")

# K-means variable clustering
vignette("kmeans-variables", package = "VarClustering")

# Working with mixed data
vignette("mixed-data", package = "VarClustering")
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

# Check package (CRAN standards)
devtools::check()

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

# Build and check
devtools::build()
R CMD check VarClustering_*.tar.gz
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

### Scientific Publications

1. **Rakotomalala, R. (2020).** *Clustering de variables.* Support de cours, Université Lyon 2.

2. **Vigneau, E., & Qannari, E. M. (2003).** *Clustering of variables around latent components.* Communications in Statistics-Simulation and Computation, 32(4), 1131-1150.  
   [DOI: 10.1081/SAC-120023882](https://doi.org/10.1081/SAC-120023882)

3. **Chavent, M., Kuentz-Simonet, V., Labenne, A., & Saracco, J. (2012).** *ClustOfVar: An R Package for the Clustering of Variables.* Journal of Statistical Software, 50(13), 1-16.  
   [DOI: 10.18637/jss.v050.i13](https://doi.org/10.18637/jss.v050.i13)

4. **Ward, J. H. (1963).** *Hierarchical Grouping to Optimize an Objective Function.* Journal of the American Statistical Association, 58(301), 236-244.

5. **Dice, L. R. (1945).** *Measures of the amount of ecologic association between species.* Ecology, 26(3), 297-302.
   [DOI: 10.2307/1932409](https://doi.org/10.2307/1932409)

6. **Lebart, L., Morineau, A., & Piron, M. (2006).** *Statistique exploratoire multidimensionnelle : Visualisation et inférences en fouille de données (4e éd.).* Dunod.


---

<div align="center">

**⭐ If you find this package useful, please consider giving it a star on GitHub! ⭐**

Made with ❤️ using [devtools](https://devtools.r-lib.org/) and [R6](https://r6.r-lib.org/)

[⬆ Back to top](#-VarClustering---advanced-variable-clustering-in-r)

</div>
