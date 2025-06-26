required_packages <- c("readr", "readxl", "xlsx",
					   "msaenet", "caret",
					   "dplyr","purrr",
					   "ggplot2","ggpubr",
					   "doParallel","doRNG","parallel",
					   "MASS","splitstackshape")

# Install missing packages
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed_packages) {
    install.packages(pkg)
  }
}