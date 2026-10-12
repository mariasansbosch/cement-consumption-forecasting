# Install the R packages needed to reproduce the analysis
pkgs <- c("rmarkdown", "knitr", "nortest", "tseries", "lmtest")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing)
