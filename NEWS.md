# DEMATELSpectral 0.2.0

## Major New Features

### Comprehensive Sensitivity Analysis
- **New `DEMATEL_Sensitivity` class** for sensitivity analysis of DEMATEL systems
- **Numerical and analytical methods** for computing sensitivity matrices (∂λmax/∂aij)
- **Critical relationship identification** to find most sensitive system relationships
- **Intervention analysis** for optimal leverage points and system modification
- **Advanced visualizations** including heatmaps, network diagrams, and statistical distributions

### Enhanced Analysis Workflows
- **`analyze_dematel_from_A()`**: Complete analysis starting from original direct influence matrix A (recommended approach)
- **`analyze_dematel_files_enhanced()`**: Extended version of original function with optional sensitivity analysis
- **`create_example_A_matrix()`**: Generate example data with A, D, and T matrices
- **Mathematical accuracy**: Uses original A matrix for sensitivity analysis, avoiding reconstruction errors

### Comprehensive Reporting and Visualization
- **`visualize_sensitivity()`**: Create multiple visualization types (heatmaps, distributions, critical relationships)
- **`plot_sensitivity_network()`**: Network visualization of critical relationships
- **`create_comprehensive_report()`**: Generate detailed analysis reports combining spectral and sensitivity results
- **`extract_enhanced_metrics()`**: Extract comprehensive metrics for research applications

## New Functions

### Core Analysis Functions
- `DEMATEL_Sensitivity()`: Create sensitivity analysis object
- `compute_sensitivity_numerical()`: Numerical sensitivity computation
- `compute_sensitivity_analytical()`: Analytical sensitivity computation (faster)
- `analyze_dematel_from_A()`: Complete analysis from original A matrix
- `create_example_A_matrix()`: Generate example data

### Analysis and Reporting Functions
- `identify_critical_relationships()`: Find most sensitive relationships
- `intervention_analysis()`: Analyze potential interventions
- `get_sensitivity_stats()`: Compute sensitivity summary statistics
- `extract_enhanced_metrics()`: Extract comprehensive metrics
- `create_comprehensive_report()`: Generate detailed reports

### Visualization Functions
- `visualize_sensitivity()`: Comprehensive sensitivity visualizations
- `plot_sensitivity_network()`: Network visualization of critical relationships

## Key Metrics Added

### Sensitivity Analysis Metrics
- **Sensitivity Matrix**: ∂λmax/∂aij for each relationship
- **Critical Relationships**: High-sensitivity relationships above specified percentiles
- **Amplifying vs. Dampening**: Classification of relationship effects
- **Mean Absolute Sensitivity**: Overall system sensitivity level
- **Intervention Efficiency**: Optimal leverage points for system changes

## Dependencies Added
- `Matrix`: For efficient matrix operations
- `ggplot2`: Enhanced visualizations
- `reshape2`: Data manipulation for plotting
- `viridis`: Advanced color scales
- `corrplot`: Correlation visualizations
- `dplyr`: Data manipulation
- `gridExtra`: Multiple plot arrangements

## Enhanced Integration
- **Backward compatibility**: All original functions remain unchanged
- **Flexible workflows**: Choose between D/T file analysis or complete A matrix analysis
- **Research-ready outputs**: Standardized metrics extraction for multi-case studies
- **Comprehensive documentation**: Detailed examples and mathematical foundations

## Mathematical Improvements
- **Accurate sensitivity computation**: Uses original A matrix instead of reconstructed approximations
- **Robust numerical methods**: Improved stability and error handling
- **Multiple computation methods**: Choice between numerical and analytical approaches
- **Validation features**: Built-in checks for matrix consistency and mathematical accuracy

## Usage Examples

### Quick Start with New Approach
```r
# Create example data
create_example_A_matrix(n = 4, case_name = "study1")

# Complete analysis
results <- analyze_dematel_from_A("A_study1.csv", include_sensitivity = TRUE)

# View results
summary(results)
plots <- visualize_sensitivity(results$sensitivity)
```

### Research Workflow
```r
# Extract metrics for multiple studies
metrics <- extract_enhanced_metrics(results)
write.csv(metrics, "study_results.csv")

# Generate comprehensive report
report <- create_comprehensive_report(results, save_report = TRUE)
```

## Breaking Changes
- None. All original functions remain fully compatible.

## Bug Fixes
- Improved matrix loading and validation
- Enhanced error handling for edge cases
- Better handling of CSV file formats

---

# DEMATELSpectral 0.1.0

## Initial Release

### Spectral Analysis Features
- Spectral decomposition of DEMATEL total relations matrix
- Eigenvalue analysis for system dynamics
- Stability and convergence analysis
- Amplification factor computation
- Condition number analysis

### Core Functions
- `analyze_dematel_files()`: Main analysis workflow
- `dematel_spectral_analysis()`: Core spectral analysis
- `extract_metrics()`: Metrics extraction
- `load_dematel_matrices()`: Matrix loading utilities
- `create_example_dematel()`: Example data generation

### Research Applications
- System dynamics analysis
- Stability characteristics
- Intervention potential assessment
- Decision support system enhancement
