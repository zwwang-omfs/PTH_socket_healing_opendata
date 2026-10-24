# Install the R packages used by lmm_section_staining.R into the default library.
# Versions used for the manuscript are listed in
# results/statistics/lmm_section_staining/fixed_effects/R_session_info.txt.

required_packages <- c("lme4", "lmerTest", "pbkrtest", "emmeans", "mvtnorm")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
}

for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required package was not installed: ", package)
  }
}

message("All required R packages are available.")
