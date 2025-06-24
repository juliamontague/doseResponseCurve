# doseResponseCurve

**doseResponseCurve** is a minimal R package providing a single function to fit interactive dose-response curves and calculate IC50 values using a 4-parameter logistic model.

## Installation

Install the package from GitHub using devtools:

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("juliamontague/doseResponseCurve")
```

##  Usage

Load the package and the example data set, then call `fit_ic50_curve`.

```r
library(doseResponseCurve)

data(testdata)
fit_ic50_curve(testdata, molarity = "nm")
```

## License

This package is licensed under the MIT License. 

