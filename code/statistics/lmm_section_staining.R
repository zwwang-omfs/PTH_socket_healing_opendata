# Section-staining endpoints (Fig. 1F/H/J, Fig. 3B/C/E/F, Fig. 4C-F; Table S2).
#
# 2 x 2 factorial linear mixed-effects model with mouse as a random intercept
# (root regions are repeated measurements within mouse), fitted by REML with
# Kenward-Roger inference. Four prespecified contrasts are adjusted jointly by
# the multivariate-t method; Cy+iPTH vs Control is reported as exploratory.
#
# Usage (from any directory):
#   Rscript code/statistics/install_r_packages.R   # once
#   Rscript code/statistics/lmm_section_staining.R
#
# Input:  data/statistics/lmm_section_staining/*.csv (group, mouse_id, root, value)
# Output: results/statistics/lmm_section_staining/

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(pbkrtest)
  library(emmeans)
  library(mvtnorm)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(sub("^--file=", "", script_arg[1])) else "."
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/")

template_dir <- file.path(repo_root, "data", "statistics", "lmm_section_staining")
output_dir <- file.path(repo_root, "results", "statistics", "lmm_section_staining")
if (nzchar(Sys.getenv("LMM_OUTPUT_DIR"))) {
  output_dir <- Sys.getenv("LMM_OUTPUT_DIR")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fixed_effect_output_dir <- file.path(output_dir, "fixed_effects")
dir.create(fixed_effect_output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  "branch_length_d3.csv",
  "branch_length_d7.csv",
  "branch_roi_area_d3.csv",
  "branch_roi_area_d7.csv",
  "co_area_emcn_area_d3.csv",
  "co_area_emcn_area_d7.csv",
  "co_area_roi_area_d3.csv",
  "co_area_roi_area_d7.csv",
  "emcn_area_roi_area_d3.csv",
  "emcn_area_roi_area_d7.csv",
  "length_roi_area_d3.csv",
  "length_roi_area_d7.csv",
  "TRAP_d3.csv",
  "TRAP_d7.csv",
  "He_woven_bone_area_d7.csv"
)

group_levels <- c("C", "Cy", "Cy + iPTH", "iPTH")
root_levels <- c("d", "m", "p")

# emmeans orders the 2 x 2 grid as:
# C, Cy, iPTH, Cy+iPTH.
primary_contrasts <- list(
  "Cy - Control" = c(-1, 1, 0, 0),
  "iPTH - Control" = c(-1, 0, 1, 0),
  "Cy+iPTH - Cy" = c(0, -1, 0, 1),
  "Cy+iPTH - iPTH" = c(0, 0, -1, 1)
)
exploratory_contrast <- list(
  "Cy+iPTH - Control" = c(-1, 0, 0, 1)
)

symbol_from_p <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(
      p < 0.0001, "#",
      ifelse(
        p < 0.001, "***",
        ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))
      )
    )
  )
}

data_structure_note <- function(data) {
  pair_counts <- aggregate(
    value ~ subject_id + root,
    data = data,
    FUN = length
  )
  repeated <- pair_counts[pair_counts$value > 1, , drop = FALSE]
  if (nrow(repeated) == 0) {
    return("")
  }
  pairs <- paste0(repeated$subject_id, ":", repeated$root)
  paste0(
    "Repeated observations from the same mouse/root were retained as recorded: ",
    paste(pairs, collapse = ", "),
    "."
  )
}

prepare_data <- function(path) {
  data <- read.csv(path, check.names = FALSE)
  required <- c("group", "mouse_id", "root", "value")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(basename(path), " is missing columns: ", paste(missing, collapse = ", "))
  }

  data$group <- factor(trimws(data$group), levels = group_levels)
  data$root <- factor(tolower(trimws(data$root)), levels = root_levels)
  data$value <- as.numeric(data$value)
  data$cy <- factor(
    ifelse(data$group %in% c("Cy", "Cy + iPTH"), "Yes", "No"),
    levels = c("No", "Yes")
  )
  data$ipth <- factor(
    ifelse(data$group %in% c("iPTH", "Cy + iPTH"), "Yes", "No"),
    levels = c("No", "Yes")
  )
  contrasts(data$cy) <- contr.sum(2)
  contrasts(data$ipth) <- contr.sum(2)
  data$subject_id <- interaction(
    data$group,
    trimws(as.character(data$mouse_id)),
    drop = TRUE,
    sep = "_M"
  )

  if (anyNA(data[, required])) {
    stop(basename(path), " contains missing or invalid required values.")
  }
  if (length(unique(data$group)) != 4) {
    stop(basename(path), " does not contain all four treatment groups.")
  }
  data
}

kr_factorial_tests <- function(model) {
  table <- as.data.frame(
    anova(model, ddf = "Kenward-Roger", type = 3)
  )
  table$term <- rownames(table)
  rownames(table) <- NULL
  table <- table[table$term %in% c("cy", "ipth", "cy:ipth", "root"), , drop = FALSE]

  labels <- c(
    cy = "Overall Cy main effect",
    ipth = "Overall iPTH main effect",
    `cy:ipth` = "Cy x iPTH interaction",
    root = "Overall root effect"
  )

  data.frame(
    result_type = "factorial_test",
    inference_method = "kenward_roger",
    comparison_family = "not_applicable",
    comparison = unname(labels[table$term]),
    estimate = NA_real_,
    std_error = NA_real_,
    df = table$DenDF,
    ci_95_low = NA_real_,
    ci_95_high = NA_real_,
    test_statistic = table$`F value`,
    p_value_raw = table$`Pr(>F)`,
    p_value_reported = table$`Pr(>F)`,
    significance = symbol_from_p(table$`Pr(>F)`),
    multiplicity_method = "none; prespecified factorial model tests",
    stringsAsFactors = FALSE
  )
}

kr_fixed_effect_coefficients <- function(
  model,
  filename,
  data,
  note,
  random_variance,
  residual_variance
) {
  coefficient_table <- as.data.frame(
    coef(summary(model, ddf = "Kenward-Roger"))
  )
  coefficient_table$fixed_effect <- rownames(coefficient_table)
  rownames(coefficient_table) <- NULL

  coefficient_definitions <- c(
    `(Intercept)` = paste(
      "Model intercept: treatment-cell grand mean at root d under sum coding",
      "for Cy and iPTH."
    ),
    cy1 = paste(
      "Cy sum-contrast coefficient (No coded +1, Yes coded -1);",
      "positive values indicate a higher outcome without Cy, averaged over iPTH."
    ),
    ipth1 = paste(
      "iPTH sum-contrast coefficient (No coded +1, Yes coded -1);",
      "positive values indicate a higher outcome without iPTH, averaged over Cy."
    ),
    rootm = "Root m coefficient relative to reference root d.",
    rootp = "Root p coefficient relative to reference root d.",
    `cy1:ipth1` = paste(
      "Cy-by-iPTH interaction coefficient under sum coding;",
      "four times this coefficient equals the difference-in-differences."
    )
  )

  df <- coefficient_table$df
  critical_t <- qt(0.975, df = df)
  singular_fit <- isSingular(model, tol = 1e-4)

  data.frame(
    endpoint_file = filename,
    endpoint = tools::file_path_sans_ext(filename),
    fixed_effect = coefficient_table$fixed_effect,
    coefficient_definition = unname(
      coefficient_definitions[coefficient_table$fixed_effect]
    ),
    estimate = coefficient_table$Estimate,
    std_error = coefficient_table$`Std. Error`,
    df = df,
    ci_95_low = coefficient_table$Estimate - critical_t * coefficient_table$`Std. Error`,
    ci_95_high = coefficient_table$Estimate + critical_t * coefficient_table$`Std. Error`,
    test_statistic = coefficient_table$`t value`,
    exact_p_value = coefficient_table$`Pr(>|t|)`,
    inference_method = "Kenward-Roger coefficient t-test",
    multiplicity_method = "none; model coefficient tests",
    model = "value ~ Cy * iPTH + root + (1 | mouse)",
    cy_coding = "sum coding: No=+1, Yes=-1",
    ipth_coding = "sum coding: No=+1, Yes=-1",
    root_reference = "d",
    observations = nrow(data),
    mice = length(unique(data$subject_id)),
    random_intercept_variance = random_variance,
    residual_variance = residual_variance,
    icc = random_variance / (random_variance + residual_variance),
    singular_fit = singular_fit,
    data_note = note,
    stringsAsFactors = FALSE
  )
}

standardize_contrast_columns <- function(result) {
  names(result)[names(result) == "contrast"] <- "comparison"
  names(result)[names(result) == "SE"] <- "std_error"
  names(result)[names(result) == "lower.CL"] <- "ci_95_low"
  names(result)[names(result) == "upper.CL"] <- "ci_95_high"
  names(result)[names(result) == "t.ratio"] <- "test_statistic"
  names(result)[names(result) == "p.value"] <- "p_value_reported"
  result
}

kr_primary_and_exploratory_contrasts <- function(model) {
  estimated_means <- emmeans(
    model,
    ~ cy * ipth,
    lmer.df = "kenward-roger"
  )

  expected_grid <- data.frame(
    cy = factor(c("No", "Yes", "No", "Yes"), levels = c("No", "Yes")),
    ipth = factor(c("No", "No", "Yes", "Yes"), levels = c("No", "Yes"))
  )
  actual_grid <- as.data.frame(estimated_means)[, c("cy", "ipth")]
  if (!identical(actual_grid, expected_grid)) {
    stop("Unexpected emmeans grid order; contrast coefficients were not applied.")
  }

  primary_raw <- as.data.frame(
    summary(
      contrast(estimated_means, method = primary_contrasts, adjust = "none"),
      infer = c(TRUE, TRUE)
    )
  )
  primary_adjusted <- as.data.frame(
    summary(
      contrast(estimated_means, method = primary_contrasts, adjust = "mvt"),
      infer = c(TRUE, TRUE)
    )
  )
  primary_raw <- standardize_contrast_columns(primary_raw)
  primary_adjusted <- standardize_contrast_columns(primary_adjusted)
  primary_adjusted$p_value_raw <- primary_raw$p_value_reported
  primary_adjusted$comparison_family <- "four_primary"
  primary_adjusted$multiplicity_method <- paste(
    "multivariate-t across four prespecified primary comparisons"
  )

  exploratory <- as.data.frame(
    summary(
      contrast(
        estimated_means,
        method = exploratory_contrast,
        adjust = "none"
      ),
      infer = c(TRUE, TRUE)
    )
  )
  exploratory <- standardize_contrast_columns(exploratory)
  exploratory$p_value_raw <- exploratory$p_value_reported
  exploratory$comparison_family <- "exploratory"
  exploratory$multiplicity_method <- "none; exploratory nominal p value"

  result <- rbind(primary_adjusted, exploratory)
  result$significance <- symbol_from_p(result$p_value_reported)
  result$result_type <- "planned_contrast"
  result$inference_method <- "kenward_roger"

  result[, c(
    "result_type",
    "inference_method",
    "comparison_family",
    "comparison",
    "estimate",
    "std_error",
    "df",
    "ci_95_low",
    "ci_95_high",
    "test_statistic",
    "p_value_raw",
    "p_value_reported",
    "significance",
    "multiplicity_method"
  )]
}

analyse_file <- function(filename) {
  data <- prepare_data(file.path(template_dir, filename))
  note <- data_structure_note(data)

  model <- lmer(
    value ~ cy * ipth + root + (1 | subject_id),
    data = data,
    REML = TRUE
  )

  output <- rbind(
    kr_factorial_tests(model),
    kr_primary_and_exploratory_contrasts(model)
  )

  random_variance <- as.numeric(VarCorr(model)$subject_id[1])
  residual_variance <- sigma(model)^2
  output$endpoint_file <- filename
  output$model <- "value ~ Cy * iPTH + root + (1 | mouse)"
  output$observations <- nrow(data)
  output$mice <- length(unique(data$subject_id))
  output$random_intercept_variance <- random_variance
  output$residual_variance <- residual_variance
  output$icc <- random_variance / (random_variance + residual_variance)
  output$singular_fit <- isSingular(model, tol = 1e-4)
  output$data_note <- note
  output$analysis_rule <- paste(
    "2 x 2 factorial LMM with Kenward-Roger inference;",
    "Cy x iPTH interaction as the primary rescue-related factorial test;",
    "multivariate-t adjustment across four prespecified primary comparisons;",
    "Cy+iPTH - Control reported separately as an unadjusted exploratory comparison."
  )
  output$symbol_definition <- paste(
    "*=p<0.05; **=p<0.01; ***=p<0.001;",
    "#=p<0.0001; ns=p>=0.05"
  )

  output_path <- file.path(
    output_dir,
    paste0(
      tools::file_path_sans_ext(filename),
      "_factorial_LMM_KR_primary4_MVT_exploratory_results.csv"
    )
  )
  write.csv(
    output,
    output_path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  message("Wrote: ", output_path)

  fixed_effect_output <- kr_fixed_effect_coefficients(
    model = model,
    filename = filename,
    data = data,
    note = note,
    random_variance = random_variance,
    residual_variance = residual_variance
  )
  fixed_effect_output_path <- file.path(
    fixed_effect_output_dir,
    paste0(
      tools::file_path_sans_ext(filename),
      "_LMM_KR_fixed_effect_coefficients.csv"
    )
  )
  write.csv(
    fixed_effect_output,
    fixed_effect_output_path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  message("Wrote: ", fixed_effect_output_path)

  fixed_effect_output
}

set.seed(20260607)

fixed_effect_tables <- lapply(input_files, analyse_file)
fixed_effect_master <- do.call(rbind, fixed_effect_tables)
master_output_path <- file.path(
  fixed_effect_output_dir,
  "Supplementary_Table_LMM_fixed_effect_coefficients_master.csv"
)
write.csv(
  fixed_effect_master,
  master_output_path,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
message("Wrote: ", master_output_path)

publication_effect_labels <- c(
  `(Intercept)` = "Intercept",
  cy1 = "Cy",
  ipth1 = "iPTH",
  rootm = "Root m vs d",
  rootp = "Root p vs d",
  `cy1:ipth1` = "Cy x iPTH"
)
publication_table <- data.frame(
  Outcome = gsub("_", " ", fixed_effect_master$endpoint),
  `Fixed effect` = unname(
    publication_effect_labels[fixed_effect_master$fixed_effect]
  ),
  Estimate = fixed_effect_master$estimate,
  SE = fixed_effect_master$std_error,
  df = fixed_effect_master$df,
  `t statistic` = fixed_effect_master$test_statistic,
  `Exact p-value` = fixed_effect_master$exact_p_value,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
publication_output_path <- file.path(
  fixed_effect_output_dir,
  "Supplementary_Table_LMM_fixed_effects_for_manuscript.csv"
)
write.csv(
  publication_table,
  publication_output_path,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
message("Wrote: ", publication_output_path)

table_notes <- c(
  "# Supplementary Table: LMM fixed-effect coefficients",
  "",
  "Suggested title: Fixed-effect estimates from the linear mixed-effects models.",
  "",
  paste0(
    "Model: value ~ Cy * iPTH + root + (1 | mouse). ",
    "Models were fitted by restricted maximum likelihood."
  ),
  paste0(
    "Degrees of freedom and coefficient t-tests used the Kenward-Roger ",
    "approximation."
  ),
  paste0(
    "Cy and iPTH used sum coding (No = +1, Yes = -1); root d was the ",
    "reference level for the root-region coefficients."
  ),
  paste0(
    "The coefficient p-values are nominal model-coefficient tests and were ",
    "not multiplicity-adjusted. Planned contrast p-values are reported ",
    "separately using the prespecified MVT adjustment."
  ),
  paste0(
    "Under the stated sum coding, positive Cy or iPTH coefficients indicate ",
    "higher outcomes in the corresponding No condition, averaged over the ",
    "other treatment factor."
  )
)
notes_output_path <- file.path(
  fixed_effect_output_dir,
  "Supplementary_Table_LMM_fixed_effects_notes.md"
)
writeLines(table_notes, notes_output_path, useBytes = TRUE)
message("Wrote: ", notes_output_path)

session_output_path <- file.path(
  fixed_effect_output_dir,
  "R_session_info.txt"
)
writeLines(capture.output(sessionInfo()), session_output_path, useBytes = TRUE)
message("Wrote: ", session_output_path)
