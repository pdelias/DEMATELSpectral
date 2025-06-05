# DEMATELSpectral

**Spectral Analysis Extension for DEMATEL Method**

This R package enhances the Decision-Making Trial and Evaluation Laboratory (DEMATEL) method with spectral analysis capabilities, providing deeper insights into system dynamics, stability characteristics, and intervention potential through eigenvalue decomposition.

Based on the research paper: *"Towards Elevating DEMATEL with Spectral Analysis"* by Delias & Kalkitsa.

## 🚀 Quick Start

### Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("pdelias/DEMATELSpectral")

# Load the package
library(DEMATELSpectral)
```

### Basic Usage

```r
# Analyze DEMATEL matrices from CSV files (no headers)
results <- analyze_dematel_files("your_D_matrix.csv", "your_T_matrix.csv")

# If your CSV files have headers
results <- analyze_dematel_files("your_D_matrix.csv", "your_T_matrix.csv", header = TRUE)

# Extract key metrics for data collection
metrics <- extract_metrics(results)
```

## 📊 What This Package Provides

### Traditional DEMATEL gives you:
- Direct and indirect influence relationships
- Centrality measures (prominence and net effect)
- Impact-relation maps

### **DEMATELSpectral adds:**
- **System Dynamics**: Maximum influence propagation potential
- **Stability Analysis**: Convergence characteristics and condition numbers  
- **Amplification Effects**: How initial changes amplify through the system
- **Structural Insights**: Dominant influence patterns and eigenvector analysis
- **Intervention Planning**: Optimal leverage points for system modification

## 🔧 Core Functions

### `analyze_dematel_files(d_file, t_file, case_name, header = FALSE)`
Complete analysis workflow from CSV files to spectral results.

**Parameters:**
- `d_file`: Path to normalized direct influence matrix (D) CSV file
- `t_file`: Path to total relations matrix (T) CSV file  
- `case_name`: Optional identifier for your analysis
- `header`: Whether CSV files contain headers (default: FALSE)

### `dematel_spectral_analysis(D, T, case_name)`
Core spectral analysis function for matrices already loaded in R.

### `extract_metrics(results)`
Extract key metrics as a data frame row for collecting results across multiple studies.

### `load_dematel_matrices(d_file, t_file, header = FALSE)`
Load D and T matrices from CSV files with proper validation.

## 📈 Key Outputs

| Metric | Description | Interpretation |
|--------|-------------|----------------|
| **λmax** | Maximum eigenvalue | Maximum long-term influence propagation potential |
| **Spectral Radius** | Largest absolute eigenvalue | System's overall influence magnitude |
| **Condition Number** | λmax/λmin ratio | Numerical sensitivity (>100 = high sensitivity) |
| **Amplification Factor** | 1/(1-ρ(D)) | How much initial effects can amplify |
| **Convergence Rate** | -ln(λ2/λmax) | Speed of pattern stabilization |
| **Concentration Ratio** | λmax/Σλi | Dominance of primary influence mode |
| **Eigenvector Stats** | SD and range of dominant eigenvector | Influence distribution patterns |

## 💡 Example Workflow

```r
library(DEMATELSpectral)

# Create example data for testing
create_example_dematel(n = 5)

# Perform spectral analysis
results <- analyze_dematel_files("D_example.csv", "T_example.csv", "Food Safety Study")

# View detailed results
print(results)

# Extract metrics for data collection
metrics <- extract_metrics(results)
print(metrics)

# Access specific values
cat("Maximum eigenvalue:", results$lambda_max)
cat("System is diagonalizable:", results$is_diagonalizable)
cat("Amplification factor:", results$amplification_factor)
```

## 📁 Input File Format

Your CSV files should contain **square numeric matrices** representing:

- **D matrix**: Normalized direct influence matrix (values typically 0-1)
- **T matrix**: Total relations matrix T = D(I-D)⁻¹

### CSV Format Requirements:
- **No headers by default** (use `header = TRUE` if your files have column names)
- **Square matrices** (same number of rows and columns)
- **Numeric values only**
- **No row names in the first column**

**Example D matrix (3×3):**
```
0.000,0.150,0.120
0.200,0.000,0.180
0.100,0.130,0.000
```

## 🔬 Research Applications

This package is designed for researchers analyzing complex systems such as:

- **Supply chain networks**
- **Organizational structures** 
- **Sustainability assessments**
- **Policy intervention systems**
- **Risk management frameworks**
- **Social network analysis**

## 📚 Multiple Case Analysis

For analyzing multiple cases (common in research):

```r
# Initialize results collector
all_metrics <- data.frame()

# Process multiple studies
case_list <- list(
  list(d = "study1_D.csv", t = "study1_T.csv", name = "Study 1"),
  list(d = "study2_D.csv", t = "study2_T.csv", name = "Study 2"),
  list(d = "study3_D.csv", t = "study3_T.csv", name = "Study 3")
)

# Analyze each case
for (case in case_list) {
  if (file.exists(case$d) && file.exists(case$t)) {
    result <- analyze_dematel_files(case$d, case$t, case$name, verbose = FALSE)
    metrics <- extract_metrics(result)
    all_metrics <- rbind(all_metrics, metrics)
  }
}

# Save consolidated results
write.csv(all_metrics, "spectral_analysis_results.csv", row.names = FALSE)

# Analyze patterns across studies
summary(all_metrics[c("lambda_max", "spectral_radius", "condition_number")])
```

## 🧮 Mathematical Background

The package implements spectral decomposition of the total relations matrix T:

**T = PΛP⁻¹**

Where:
- **Λ** = diagonal matrix of eigenvalues λ₁, λ₂, ..., λₙ
- **P** = matrix of corresponding eigenvectors
- **λmax** = dominant eigenvalue revealing system characteristics

Key mathematical relationships:
- **Diagonalizability**: Verified through eigenvector linear independence
- **Perron-Frobenius Theory**: Applied to ensure meaningful interpretation
- **Convergence Properties**: Derived from eigenvalue ratios

## 🐛 Troubleshooting

### Common Issues:

**"Matrix is not square" error:**
- Check if your CSV has headers when you specified `header = FALSE`
- Ensure equal number of rows and columns

**"Incomplete final line" warning:**
- Normal for some CSV files, doesn't affect results
- Use the updated `create_example_dematel()` to generate properly formatted files

**"Bad credentials" error:**
- Related to GitHub authentication, doesn't affect package usage
- Package functions will work normally

## 📖 Citation

If you use this package in your research, please cite:

```
Delias, P., & Kalkitsa, K. (2024). Towards Elevating DEMATEL with Spectral Analysis. 
11th International Conference on Decision Support System Technology.
```

## 🤝 Contributing

Issues and improvements welcome! Please report bugs or suggest enhancements through GitHub issues.

## 📄 License

This package is released under the MIT License.

---

**Package Author**: Pavlos Delias  
**GitHub**: https://github.com/pdelias/DEMATELSpectral  
**Research**: Spectral analysis extensions for decision support systems
