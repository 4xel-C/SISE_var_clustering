# VarClustering Shiny Application Guide

A comprehensive guide to the interactive Shiny application for variable clustering.

---

## Table of Contents

1. [Launching the Application](#launching-the-application)
2. [Application Overview](#application-overview)
3. [User Workflow](#user-workflow)
   - [Step 1: Import Data](#step-1-import-data)
   - [Step 2: Configure Variables](#step-2-configure-variables)
   - [Step 3: Set Clustering Parameters](#step-3-set-clustering-parameters)
   - [Step 4: Run Clustering](#step-4-run-clustering)
   - [Step 5: Analyze Results](#step-5-analyze-results)
   - [Step 6: Export Results](#step-6-export-results)
4. [Available Algorithms](#available-algorithms)
5. [Output Tabs](#output-tabs)
6. [Troubleshooting](#troubleshooting)

---

## Launching the Application

### From R Console

```r
# Load the package
library(VarClustering)

# Launch the Shiny application
varclust_gui()
```

### Required Dependencies

The application requires these packages (automatically checked at launch):

- `shiny` - Web framework
- `DT` - Interactive tables
- `readxl` - Excel file import

If any package is missing, the app will display an error with installation instructions.

---

## Application Overview

The VarClustering Shiny application provides an interactive interface for performing variable clustering analysis without writing code. It supports both **Hierarchical Variable Clustering (HClustVar)** and **K-means Variable Clustering**.

### Interface Layout

The application consists of two main areas:

| Area | Description |
|------|-------------|
| **Sidebar (Left)** | Data import, variable configuration, clustering parameters, export |
| **Main Panel (Right)** | Dynamic tabs with visualizations and results |

<!-- 📸 SCREENSHOT: Full application interface overview -->
<!-- Place screenshot here: screenshots/app_overview.png -->
![Application Overview](screenshots/app_overview.png)

---

## User Workflow

### Step 1: Import Data

The application accepts **CSV** and **Excel** files (.xlsx, .xls).

#### Import Options

| Option | Description |
|--------|-------------|
| **First row as column names** | Check if your file has a header row |
| **Excel sheet number** | For Excel files with multiple sheets |

#### Supported Data Types

- **Quantitative variables**: Numeric columns
- **Qualitative variables**: Factor/character columns
- **Mixed**: Combination of both types

<!-- 📸 SCREENSHOT: Data import panel -->
<!-- Place screenshot here: screenshots/data_import.png -->
![Data Import](screenshots/data_import.png)

#### After Import

Once your file is uploaded:
1. The **Data Preview** tab shows your imported data
2. The **Configure Variables** button becomes active
3. Data summary displays number of observations and variables

<!-- 📸 SCREENSHOT: Data preview tab with imported data -->
<!-- Place screenshot here: screenshots/data_preview.png -->
![Data Preview](screenshots/data_preview.png)

---

### Step 2: Configure Variables

Click **"⚙️ Configure Variables"** to open the configuration modal.

#### Variable Roles

| Role | Description | Usage |
|------|-------------|-------|
| **Active** | Used in clustering analysis | Core variables for clustering |
| **Illustrative** | For interpretation only | Predicted to clusters after fitting |
| **Excluded** | Completely ignored | ID columns, irrelevant variables |

<!-- 📸 SCREENSHOT: Variable configuration modal -->
<!-- Place screenshot here: screenshots/variable_config_modal.png -->
![Variable Configuration Modal](screenshots/variable_config_modal.png)

#### Configuration Summary

After configuration, the **Variable Configuration** tab displays:
- ✅ Active variables (green panel)
- ℹ️ Illustrative variables (orange panel)
- ❌ Excluded variables (red panel)

<!-- 📸 SCREENSHOT: Variable configuration summary with three panels -->
<!-- Place screenshot here: screenshots/variable_config_summary.png -->
![Variable Configuration Summary](screenshots/variable_config_summary.png)

---

### Step 3: Set Clustering Parameters

Configure the clustering algorithm and its parameters.

#### Algorithm Selection

| Algorithm | Data Type | Description |
|-----------|-----------|-------------|
| **CAH Variables clustering** | Quant / Qual / Mixed | Hierarchical clustering with dendrogram |
| **K-means Variables clustering** | Quantitative only | Iterative reallocation with PCA centroids |

<!-- 📸 SCREENSHOT: Algorithm selection dropdown -->
<!-- Place screenshot here: screenshots/algorithm_selection.png -->
![Algorithm Selection](screenshots/algorithm_selection.png)

#### HClustVar Parameters

| Parameter | Options | Description |
|-----------|---------|-------------|
| **Variable type** | Auto-detect, Quantitative, Qualitative, Mixed | Type of variables to cluster |
| **Distance metric** | R², r | For quantitative variables |
| **Clustering method** | Ward D, Ward D2, Complete, Average, Single, etc. | Linkage method for CAH |
| **Number of clusters** | 2-10 | Target number of clusters |

<!-- 📸 SCREENSHOT: HClustVar parameters panel -->
<!-- Place screenshot here: screenshots/hclust_params.png -->
![HClustVar Parameters](screenshots/hclust_params.png)

#### K-means Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Distance metric** | R² | Correlation-based distance |
| **Max iterations** | 100 | Maximum iterations per run |
| **n_init (restarts)** | 10 | Number of random initializations |
| **Random seed** | 0 (random) | For reproducibility |
| **Number of clusters** | 3 | Target K value |

<!-- 📸 SCREENSHOT: K-means parameters panel -->
<!-- Place screenshot here: screenshots/kmeans_params.png -->
![K-means Parameters](screenshots/kmeans_params.png)

---

### Step 4: Run Clustering

Click **"🚀 Run Clustering"** to execute the analysis.

The application will:
1. Extract active variables from your data
2. Fit the selected clustering model
3. Generate all visualizations and summaries
4. Predict clusters for illustrative variables (if configured)

<!-- 📸 SCREENSHOT: Run Clustering button highlighted -->
<!-- Place screenshot here: screenshots/run_clustering.png -->
![Run Clustering](screenshots/run_clustering.png)

---

### Step 5: Analyze Results

After clustering completes, multiple tabs become available with results.

#### Summary Tab

The **Summary** tab provides comprehensive clustering results:

##### Quality Metrics

| Metric | Description | Interpretation |
|--------|-------------|----------------|
| **Average Silhouette** | Clustering quality score | > 0.5 is good, > 0.7 is excellent |
| **Total Variance Explained** | Proportion of variance captured | Higher is better |

<!-- 📸 SCREENSHOT: Quality metrics cards (silhouette and variance) -->
<!-- Place screenshot here: screenshots/quality_metrics.png -->
![Quality Metrics](screenshots/quality_metrics.png)

##### Cluster Summary Table

Shows for each cluster:
- Number of members
- Variation explained
- Proportion explained

<!-- 📸 SCREENSHOT: Cluster summary table -->
<!-- Place screenshot here: screenshots/cluster_summary_table.png -->
![Cluster Summary Table](screenshots/cluster_summary_table.png)

##### Variable Assignments Table

Details for each variable:
- Assigned cluster
- R² with own cluster
- R² with next closest cluster
- Ratio (quality indicator)

<!-- 📸 SCREENSHOT: Variable assignments table -->
<!-- Place screenshot here: screenshots/variable_assignments.png -->
![Variable Assignments](screenshots/variable_assignments.png)

##### Centroids Correlations

Correlation matrix between cluster centroids (lower correlations = better separation).

<!-- 📸 SCREENSHOT: Centroids correlation matrix -->
<!-- Place screenshot here: screenshots/centroids_correlations.png -->
![Centroids Correlations](screenshots/centroids_correlations.png)

##### R² Matrix

Complete matrix of R² values between all variables and all cluster centroids.

<!-- 📸 SCREENSHOT: Full R² matrix -->
<!-- Place screenshot here: screenshots/r2_matrix.png -->
![R² Matrix](screenshots/r2_matrix.png)

---

#### Visualization Tabs

##### Elbow Method

Plot of aggregation levels (HClustVar) or inertia (K-means) to help determine optimal K.

<!-- 📸 SCREENSHOT: Elbow plot -->
<!-- Place screenshot here: screenshots/elbow_plot.png -->
![Elbow Plot](screenshots/elbow_plot.png)

##### Factorial Projection

Variables projected in reduced dimensional space (MDS for HClustVar, PCA for K-means).

<!-- 📸 SCREENSHOT: Factorial projection plot with colored clusters -->
<!-- Place screenshot here: screenshots/factorial_projection.png -->
![Factorial Projection](screenshots/factorial_projection.png)

##### Dendrogram (HClustVar only)

Hierarchical tree structure with cluster coloring.

<!-- 📸 SCREENSHOT: Dendrogram with colored clusters -->
<!-- Place screenshot here: screenshots/dendrogram.png -->
![Dendrogram](screenshots/dendrogram.png)

##### Silhouette Plot (HClustVar only)

Silhouette coefficients for each variable showing cluster fit quality.

<!-- 📸 SCREENSHOT: Silhouette plot -->
<!-- Place screenshot here: screenshots/silhouette_plot.png -->
![Silhouette Plot](screenshots/silhouette_plot.png)

---

#### Predict Tab

If illustrative variables were configured, this tab shows their predicted cluster assignments.

##### Predicted Clusters

Table with cluster assignment for each illustrative variable.

<!-- 📸 SCREENSHOT: Prediction results table -->
<!-- Place screenshot here: screenshots/prediction_results.png -->
![Prediction Results](screenshots/prediction_results.png)

##### Proximity Matrix

R² values between illustrative variables and all cluster centroids.

<!-- 📸 SCREENSHOT: Prediction proximities matrix -->
<!-- Place screenshot here: screenshots/prediction_proximities.png -->
![Prediction Proximities](screenshots/prediction_proximities.png)

---

### Step 6: Export Results

Click **"Download Summary (CSV)"** in the sidebar to export clustering results.

The exported file contains:
- Variable names
- Cluster assignments
- Quality metrics

<!-- 📸 SCREENSHOT: Export button in sidebar -->
<!-- Place screenshot here: screenshots/export_button.png -->
![Export Button](screenshots/export_button.png)

---

## Available Algorithms

### HClustVar (Hierarchical Variable Clustering)

**Best for**: Exploratory analysis, mixed data types, when you want to see the hierarchical structure.

**Features**:
- Supports quantitative, qualitative, and mixed variables
- Multiple linkage methods (Ward, Complete, Average, etc.)
- Dendrogram visualization
- Cut tree at any level

**Parameters**:
```
Variable type: auto / quant / qual / mixed
Distance metric: R² / r (for quantitative)
Clustering method: ward.D / ward.D2 / complete / average / single
Number of clusters: 2-10
```

### KmeansVariables (K-means Variable Clustering)

**Best for**: Large quantitative datasets, when K is known or can be estimated.

**Features**:
- Fast computation with multiple initializations
- PCA-based centroids
- Automatic K selection (via elbow method)
- Reproducible results with random seed

**Parameters**:
```
Distance metric: R² / r
Max iterations: 100 (default)
n_init: 10 (default)
Random seed: 0 = random
Number of clusters: 2-10
```

---

## Output Tabs

| Tab | HClustVar | K-means | Description |
|-----|-----------|---------|-------------|
| 📊 Data Preview | ✅ | ✅ | View imported data |
| 🔍 Variable Configuration | ✅ | ✅ | Configure variable roles |
| 📋 Summary | ✅ | ✅ | Quality metrics and cluster details |
| 📈 Elbow Method | ✅ | ✅ | Optimal K determination |
| 🗺️ Factorial Projection | ✅ | ✅ | Variable visualization |
| ➡️ Predict | ✅ | ✅ | Illustrative variable predictions |
| 🌳 Dendrogram | ✅ | ❌ | Hierarchical tree |
| 📊 Silhouette Plot | ✅ | ❌ | Cluster quality per variable |

---

## Troubleshooting

### Common Issues

#### "No file selected"

**Solution**: Click "Browse" and select a CSV or Excel file.

#### "Error: requires numeric data"

**Cause**: K-means only supports quantitative variables.

**Solution**:
- Use HClustVar for qualitative/mixed data
- Or configure only numeric variables as "Active"

#### Empty visualizations

**Cause**: Clustering not yet executed.

**Solution**: Click "🚀 Run Clustering" button.

#### Poor silhouette scores

**Solutions**:
- Try different number of clusters
- Change clustering method (Ward usually works best)
- Check for problematic variables in Variable Configuration

#### Illustrative variables not predicted

**Cause**: No illustrative variables configured.

**Solution**: Open Variable Configuration and set some variables as "Illustrative".

---

## Best Practices

### 1. Data Preparation

- Remove ID columns (set as Excluded)
- Handle missing values before import
- Ensure correct variable types in your file

### 2. Variable Configuration

- Start with all variables as Active
- Move irrelevant variables to Excluded
- Use Illustrative for variables you want to predict

### 3. Parameter Selection

- Start with **Ward D2** method (usually best)
- Use **R²** distance metric for most cases
- Check **Elbow Method** tab to choose K

### 4. Result Interpretation

- Check **Average Silhouette** (> 0.5 is good)
- Examine **Factorial Projection** for cluster separation
- Review **Variable Assignments** for cluster composition

---

## Example Workflow

```r
# 1. Launch the application
library(VarClustering)
varclust_gui()

# 2. In the Shiny interface:
#    - Import your data file (CSV/Excel)
#    - Configure variables (Active/Illustrative/Excluded)
#    - Select algorithm (HClustVar or K-means)
#    - Set parameters (method, distance, K)
#    - Click "Run Clustering"
#    - Analyze results in the tabs
#    - Export summary if needed
```

---

## Screenshots Directory

To complete this documentation, add screenshots in a `screenshots/` folder:

```
vignettes/
├── Shiny_Application_Guide.md
└── screenshots/
    ├── app_overview.png
    ├── data_import.png
    ├── data_preview.png
    ├── variable_config_modal.png
    ├── variable_config_summary.png
    ├── algorithm_selection.png
    ├── hclust_params.png
    ├── kmeans_params.png
    ├── run_clustering.png
    ├── quality_metrics.png
    ├── cluster_summary_table.png
    ├── variable_assignments.png
    ├── centroids_correlations.png
    ├── r2_matrix.png
    ├── elbow_plot.png
    ├── factorial_projection.png
    ├── dendrogram.png
    ├── silhouette_plot.png
    ├── prediction_results.png
    ├── prediction_proximities.png
    └── export_button.png
```

---

## Session Info

This application was developed using:

- R version 4.x
- Shiny 1.x
- VarClustering package

---

*VarClustering Package - Eugénie Barlet, Modou Mboup, Axel Cano*
