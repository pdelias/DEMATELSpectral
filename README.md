# DEMATELSpectral

**Spectral Analysis and Sensitivity Analysis Extensions for DEMATEL Method**

This R package enhances the Decision-Making Trial and Evaluation Laboratory (DEMATEL) method with both spectral analysis capabilities and comprehensive sensitivity analysis, providing deeper insights into system dynamics, stability characteristics, and intervention potential through eigenvalue decomposition and sensitivity matrices.

## New in Version 0.2.0 🎉

- **Comprehensive Sensitivity Analysis**: Analyze how changes in direct influence relationships affect the system's dominant eigenvalue
- **Critical Relationship Identification**: Automatically identify the most sensitive relationships in your system
- **Intervention Analysis**: Find optimal leverage points for system modification
- **Advanced Visualizations**: Create heatmaps, network diagrams, and distribution plots
- **Enhanced Integration**: Combine spectral and sensitivity analysis in unified workflows

## Features

### Original Spectral Analysis
- Direct and indirect influence relationships
- Centrality measures (prominence and net effect)
- Impact-relation maps
- **System Dynamics**: Maximum influence propagation potential
- **Stability Analysis**: Convergence characteristics and condition numbers
- **Amplification Effects**: How initial changes amplify through the system
- **Structural Insights**: Dominant influence patterns and eigenvector analysis

### New Sensitivity Analysis
- **Sensitivity Matrix Computation**: Both numerical and analytical methods
- **Critical Relationships**: Identify relationships with highest impact on system stability
- **Intervention Planning**: Optimal leverage points for system modification
- **Amplification vs. Dampening**: Classify relationships by their effect type
- **Comprehensive Visualizations**: Heatmaps, networks, and statistical distributions
- **Integrated Reporting**: Combined spectral and sensitivity analysis reports

## Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("pdelias/DEMATELSpectral")

# Load the package
library(DEMATELSpectral)
```

## Quick Start

### Basic Spectral Analysis (Original Functionality)
```r
# Analyze DEMATEL matrices from CSV files
results <- analyze_dematel_files("your_D_matrix.csv", "your_T_matrix.csv")

# Extract key metrics
metrics <- extract_metrics(results)
```

### New: Complete Analysis from Original Matrix (Recommended)
```r
# Create test data with original A matrix
create_example_A_matrix(n = 4, case_name = "my_study")

# Perform complete analysis (most mathematically accurate)
results <- analyze_dematel_from_A("A_my_study.csv", 
                                 case_name = "My Research",
                                 include_sensitivity = TRUE)

# View comprehensive results  
print(results)
summary(results)

# Extract all metrics
metrics <- extract_enhanced_metrics(results)

# Create visualizations
plots <- visualize_sensitivity(results$sensitivity)

# Generate comprehensive report
report <- create_comprehensive_report(results, save_report = TRUE)
```

### Enhanced Integrated Analysis
```r
# Combined spectral and sensitivity analysis
results <- analyze_dematel_files_enhanced("D_matrix.csv", "T_matrix.csv",
                                         case_name = "My Study",
                                         include_sensitivity = TRUE)

# View comprehensive results
print(results)
summary(results)

# Extract enhanced metrics for research
metrics <- extract_enhanced_metrics(results)

# Generate comprehensive report
report <- create_comprehensive_report(results, save_report = TRUE)
```

### Core Functions

#### Recommended Workflow Functions
- `analyze_dematel_from_A()`: **Complete analysis starting from original A matrix (most accurate)**
- `create_example_A_matrix()`: Generate example data with A, D, and T matrices
- `extract_enhanced_metrics()`: Extract comprehensive metrics from complete analysis
- `create_comprehensive_report()`: Generate detailed analysis reports

#### Sensitivity Analysis Functions
- `DEMATEL_Sensitivity()`: Create sensitivity analysis object
- `compute_sensitivity_numerical()`: Numerical sensitivity computation
- `compute_sensitivity_analytical()`: Analytical sensitivity computation (faster)
- `identify_critical_relationships()`: Find most sensitive relationships
- `intervention_analysis()`: Analyze potential interventions
- `visualize_sensitivity()`: Create comprehensive visualizations
- `plot_sensitivity_network()`: Network visualization of critical relationships

### Enhanced Integration Functions
- `analyze_dematel_files_enhanced()`: Combined spectral and sensitivity analysis
- `extract_enhanced_metrics()`: Extract comprehensive metrics
- `create_comprehensive_report()`: Generate detailed analysis reports

### Original Functions (Still Available)
- `analyze_dematel_files()`: Core spectral analysis workflow
- `dematel_spectral_analysis()`: Spectral analysis function
- `extract_metrics()`: Extract spectral metrics
- `load_dematel_matrices()`: Load matrices from CSV
- `create_example_dematel()`: Generate example data

## Key Metrics

### Spectral Analysis Metrics
| Metric | Description | Interpretation |
|--------|-------------|----------------|
| λmax | Maximum eigenvalue | Maximum long-term influence propagation potential |
| Spectral Radius | Largest absolute eigenvalue | System's overall influence magnitude |
| Condition Number | λmax/λmin ratio | Numerical sensitivity |
| Amplification Factor | 1/(1-ρ(D)) | How much initial effects can amplify |
| Convergence Rate | -ln(λ2/λmax) | Speed of pattern stabilization |

### New Sensitivity Analysis Metrics
| Metric | Description | Interpretation |
|--------|-------------|----------------|
| Sensitivity Matrix | ∂λmax/∂aij | How changes in each relationship affect system stability |
| Critical Relationships | High-sensitivity relationships | Key leverage points for intervention |
| Amplifying vs. Dampening | Positive vs. negative sensitivity | Relationship effect classification |
| Mean Absolute Sensitivity | Average impact magnitude | Overall system sensitivity level |

## Example Workflow

```r
library(DEMATELSpectral)

# 1. Create example data (creates A, D, and T matrices)
create_example_A_matrix(n = 5, case_name = "food_safety")

# 2. Perform complete analysis (recommended approach)
results <- analyze_dematel_from_A("A_food_safety.csv",
                                 case_name = "Food Safety Study",
                                 include_sensitivity = TRUE)

# 3. Examine results
print(results)
summary(results)

# 4. Focus on sensitivity analysis
sens_analysis <- results$sensitivity
critical_relationships <- identify_critical_relationships(sens_analysis)
print(critical_relationships)

# 5. Create visualizations
plots <- visualize_sensitivity(sens_analysis, save_plots = TRUE)

# 6. Analyze intervention options
interventions <- intervention_analysis(sens_analysis, target_lambda_change = -0.1)
feasible_interventions <- interventions[interventions$feasible, ]
print(head(feasible_interventions, 5))

# 7. Generate comprehensive report
report <- create_comprehensive_report(results, save_report = TRUE)

# 8. Extract metrics for research collection
metrics <- extract_enhanced_metrics(results)
write.csv(metrics, "analysis_metrics.csv", row.names = FALSE)
```

## Multiple Case Analysis

For analyzing multiple studies (common in research):

```r
# Initialize results collector
all_metrics <- data.frame()

# Process multiple studies
case_list <- list(
  list(d = "study1_D.csv", t = "study1_T.csv", name = "Study 1"),
  list(d = "study2_D.csv", t = "study2_T.csv", name = "Study 2"),
  list(d = "study3_D.csv", t = "study3_T.csv", name = "Study 3")
)

# Analyze each case with sensitivity analysis
for (case in case_list) {
  if (file.exists(case$d) && file.exists(case$t)) {
    result <- analyze_dematel_files_enhanced(case$d, case$t, case$name, 
                                           include_sensitivity = TRUE, 
                                           verbose = FALSE)
    metrics <- extract_enhanced_metrics(result)
    all_metrics <- rbind(all_metrics, metrics)
  }
}

# Save consolidated results
write.csv(all_metrics, "comprehensive_analysis_results.csv", row.names = FALSE)

# Analyze patterns across studies
summary(all_metrics[c("lambda_max", "sensitivity_mean_abs", "n_critical_90th")])
```

## Data Format

Your CSV files should contain square numeric matrices representing:
- **D matrix**: Normalized direct influence matrix (values typically 0-1)
- **T matrix**: Total relations matrix T = D(I-D)⁻¹
- No headers by default (use `header = TRUE` if your files have column names)
- Square matrices (same number of rows and columns)
- Numeric values only

Example D matrix (3×3):
```
0.000,0.150,0.120
0.200,0.000,0.180
0.100,0.130,0.000
```

## Mathematical Foundation

### Spectral Decomposition
The package implements spectral decomposition of the total relations matrix T:
```
T = PΛP⁻¹
```
Where:
- Λ = diagonal matrix of eigenvalues λ₁, λ₂, ..., λₙ
- P = matrix of corresponding eigenvectors
- λmax = dominant eigenvalue revealing system characteristics

### Sensitivity Analysis
Sensitivity analysis computes how the dominant eigenvalue changes with respect to changes in the direct influence matrix:
```
∂λmax/∂aij = v^T (∂T/∂aij) u
```
Where v and u are the left and right eigenvectors of the dominant eigenvalue.

## Troubleshooting

### Common Issues

**"Matrix is not square" error:**
- Check if your CSV has headers when you specified `header = FALSE`
- Ensure equal number of rows and columns

**"Incomplete final line" warning:**
- Normal for some CSV files, doesn't affect results

**High computation time for sensitivity analysis:**
- Use `sensitivity_method = "analytical"` for faster computation
- Reduce `sensitivity_epsilon` for numerical method
- Consider analyzing smaller subsets first

## Research Applications

This package is particularly useful for:
- **Decision Support Systems**: Identifying critical decision factors
- **Risk Management**: Finding key risk propagation pathways
- **Policy Analysis**: Understanding policy intervention points
- **Supply Chain Management**: Analyzing critical relationships
- **Social Network Analysis**: Understanding influence patterns
- **System Dynamics**: Modeling complex adaptive systems

## Citation

If you use this package in your research, please cite:

```
Delias, P., & Kalkitsa, K. (2024). Towards Elevating DEMATEL with Spectral Analysis.
11th International Conference on Decision Support System Technology.
```

For the sensitivity analysis features:
```
Delias, P. (2025). DEMATELSpectral: Enhanced DEMATEL Analysis with Spectral 
and Sensitivity Analysis. R package version 0.2.0.
```

## Contributing

Issues and improvements welcome! Please report bugs or suggest enhancements through [GitHub issues](https://github.com/pdelias/DEMATELSpectral/issues).

## License

This package is released under the MIT License.

## Contact

**Package Author**: Pavlos Delias  
**GitHub**: [https://github.com/pdelias/DEMATELSpectral](https://github.com/pdelias/DEMATELSpectral)  
**Research Focus**: Spectral analysis extensions for decision support systems

---

*Enhanced with comprehensive sensitivity analysis for deeper system insights*
