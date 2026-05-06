library(data.table)
library(rpart)
library(rpart.plot)
library(dplyr)
library(tidyr)
library(ggplot2)

cell_lines <- c("H2009","HCC1395","HCC1937","HCC1954","H1437","Hs578T","HG008","COLO829")

# adjust to wherever your per-cell-line SV TSVs live
input_dir <- "../SVs/SV_tables"
setwd("../SVs")
dir.create("Trees",        showWarnings = FALSE)
dir.create("Trees/Plots",  showWarnings = FALSE)
dir.create("Trees/Rules",  showWarnings = FALSE)


###########################################################################################################################################################################
########################################################################### Helper functions #############################################################################
###########################################################################################################################################################################

# load one cell line's SV TSV and apply the minimal cleaning needed for rpart
load_cell_line <- function(cell_line) {
  variants <- as.data.frame(fread(file.path(input_dir, paste0(cell_line, "_SV_variants_table_annotated.tsv"))))
  # PR/SR breaks R as a column name; rename it
  variants <- variants %>% dplyr::rename(PR_SR_ratio = `PR/SR`)
  # SV_type as factor for native categorical splits
  # in_bed as logical (comes in as "True"/"False" strings from Python)
  # is_correct as logical
  variants <- variants %>%
    mutate(SV_type    = factor(SV_type),
           in_bed     = as.logical(in_bed),
           is_correct = as.logical(is_correct))
  variants
}

# load the pooled file containing variants from all cell lines
load_all_file <- function() {
  variants <- as.data.frame(fread(file.path(input_dir, "ALL_cell_lines_SV_variants_table_annotated.tsv")))
  variants <- variants %>% dplyr::rename(PR_SR_ratio = `PR/SR`)
  variants <- variants %>%
    mutate(SV_type    = factor(SV_type),
           in_bed     = as.logical(in_bed),
           is_correct = as.logical(is_correct))
  variants
}

# extract the rules of a (pruned) tree as a tidy data frame for inspection
tree_rules_df <- function(tree, cell_line) {
  if (nrow(tree$frame) == 1) {
    return(data.frame(Cell_line = cell_line, Rule = "ROOT (no splits)", Predicted_class = NA, Node_n = nrow(tree$frame)))
  }
  rules <- rpart.rules(tree, cover = TRUE, nn = TRUE, style = "tall")
  rules$Cell_line <- cell_line
  rules
}

# fit, prune to cp.min and 1-SE, return everything in one list
fit_tree <- function(df, cell_line) {
  set.seed(42)

  # grow a deep tree (small cp, small minsplit) — pruning will trim it.
  # parms$prior = c(0.5, 0.5) treats the two classes as equally weighted regardless
  # of the actual ~87% FALSE / ~13% TRUE class imbalance. Without this, the tree may
  # predict FALSE everywhere and stop, since that's already 87% accurate.
  # in_bed added to the feature set alongside the original predictors.
  tree <- rpart(is_correct ~ SV_type + SS + SR + PR + PR_SR_ratio + ND + ND2 + in_bed,
                data    = df,
                method  = "class",
                parms   = list(prior = c(0.5, 0.5)),
                control = rpart.control(cp = 0.01, minsplit = 10, minbucket = 3, xval = 10))

  # cp table: each row is a candidate pruning level with its CV error
  cp_table <- as.data.frame(tree$cptable)

  # cp.min = the cp that minimises cross-validated error (smallest CV-error tree)
  cp_min   <- cp_table$CP[which.min(cp_table$xerror)]

  # 1-SE rule: the simplest tree whose CV error is within 1 SE of the minimum
  # (smaller, more parsimonious — usually preferred for interpretability)
  min_idx  <- which.min(cp_table$xerror)
  threshold_1se <- cp_table$xerror[min_idx] + cp_table$xstd[min_idx]
  cp_1se   <- cp_table$CP[which(cp_table$xerror <= threshold_1se)[1]]

  list(
    full       = tree,
    pruned_min = prune(tree, cp = cp_min),
    pruned_1se = prune(tree, cp = cp_1se),
    cp_table   = cp_table,
    cp_min     = cp_min,
    cp_1se     = cp_1se,
    cell_line  = cell_line
  )
}

# write tree plot, rules, and cp table to disk for one cell line
save_tree_outputs <- function(fit, which = c("pruned_1se", "pruned_min", "full")) {
  which     <- match.arg(which)
  tree      <- fit[[which]]
  cell_line <- fit$cell_line

  # plot
  png(file.path("Trees/Plots", paste0(cell_line, "_", which, "_tree.png")),
      width = 12, height = 8, units = "in", res = 300)
  rpart.plot(tree,
             type        = 4,           # split labels above, factor levels shown explicitly
             extra       = 104,         # show fitted prob and percent of obs at each node
             box.palette = "RdYlGn",    # red = FALSE, green = TRUE
             branch      = 0.3,
             tweak       = 1.1,
             main        = paste0(cell_line, " — ", which, " (cp = ",
                                  round(if (which == "pruned_min") fit$cp_min
                                        else if (which == "pruned_1se") fit$cp_1se
                                        else min(fit$cp_table$CP), 4), ")"))
  dev.off()

  # rules — the human-readable thresholds your boss wants to see
  rules_file <- file.path("Trees/Rules", paste0(cell_line, "_", which, "_rules.txt"))
  sink(rules_file)
  cat("Cell line:", cell_line, "\n")
  cat("Tree:     ", which, "\n")
  cat("Root n:   ", nrow(tree$frame), "nodes\n\n")
  cat("=== Decision rules ===\n")
  print(rpart.rules(tree, cover = TRUE, nn = TRUE, style = "tall"))
  cat("\n=== Variable importance ===\n")
  if (!is.null(tree$variable.importance)) {
    print(round(tree$variable.importance, 2))
  } else {
    cat("(no splits — root-only tree)\n")
  }
  cat("\n=== CP table (cross-validated error per pruning level) ===\n")
  print(fit$cp_table)
  cat("\nSelected cp_min =", round(fit$cp_min, 4), "(minimum CV error)\n")
  cat("Selected cp_1se =", round(fit$cp_1se, 4), "(1-SE rule, simpler tree)\n")
  sink()

  invisible(NULL)
}

# --- evaluate a tree on a data set and return confusion-matrix metrics ---
eval_tree <- function(tree, test_df) {
  feature_cols <- c("SV_type","SS","SR","PR","PR_SR_ratio","ND","ND2","in_bed")
  pred <- predict(tree, test_df[, feature_cols], type = "class")
  TP <- sum(pred == TRUE  & test_df$is_correct == TRUE)
  FP <- sum(pred == TRUE  & test_df$is_correct == FALSE)
  FN <- sum(pred == FALSE & test_df$is_correct == TRUE)
  TN <- sum(pred == FALSE & test_df$is_correct == FALSE)
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  f1        <- if (!is.na(precision) & !is.na(recall) & (precision + recall) > 0)
                 2 * precision * recall / (precision + recall) else NA
  list(TP = TP, FP = FP, FN = FN, TN = TN,
       Precision = round(precision, 3),
       Recall    = round(recall,    3),
       F1        = round(f1,        3))
}

# --- extract root split feature and threshold from a fitted tree ---
root_split <- function(tree) {
  if (nrow(tree$frame) == 1) return(list(var = "(no splits)", threshold = NA))
  splits   <- as.data.frame(tree$splits)
  root_var <- rownames(splits)[1]

  # for numeric features, splits$index is the threshold — directly readable
  # for factor features, splits$index is an internal lookup index, not meaningful;
  # read the actual level partition from tree$csplit instead
  if (!is.null(tree$csplit) && root_var %in% names(tree$xlevels)) {
    levels    <- tree$xlevels[[root_var]]
    csplit_row <- tree$csplit[1, ]          # 1 = goes left, 3 = goes right
    left_levels  <- levels[csplit_row == 1]
    right_levels <- levels[csplit_row == 3]
    threshold <- paste0(paste(left_levels,  collapse = ","), " | ",
                        paste(right_levels, collapse = ","))
  } else {
    threshold <- splits$index[1]
  }
  list(var = root_var, threshold = threshold)
}


###########################################################################################################################################################################
############################################################################# Fit per cell line ##########################################################################
###########################################################################################################################################################################

all_fits     <- list()
cell_results <- data.frame()   # collects one row of metrics per cell line

for (cell_line in cell_lines) {
  message("Fitting tree for ", cell_line)
  df <- load_cell_line(cell_line)

  # quick label-balance report
  tab <- table(df$is_correct)
  message("  ", nrow(df), " variants: ",
          tab["TRUE"], " TRUE / ", tab["FALSE"], " FALSE")

  fit <- fit_tree(df, cell_line)
  all_fits[[cell_line]] <- fit

  # write both pruning choices so you can compare
  save_tree_outputs(fit, "pruned_1se")
  save_tree_outputs(fit, "pruned_min")

  # evaluate pruned_1se tree on its own training data and collect metrics
  # note: evaluated on training data → optimistic; use LOCO results for honest estimates
  rs   <- root_split(fit$pruned_1se)
  perf <- eval_tree(fit$pruned_1se, df)
  cell_results <- rbind(cell_results, data.frame(
    Cell_line           = cell_line,
    Top_split_var       = rs$var,
    Top_split_threshold = rs$threshold,
    TP        = perf$TP,
    FP        = perf$FP,
    FN        = perf$FN,
    TN        = perf$TN,
    Precision = perf$Precision,
    Recall    = perf$Recall,
    F1        = perf$F1
  ))
}

# MEAN and SD rows, then write the summary TSV
cell_summary <- rbind(
  cell_results,
  data.frame(Cell_line = "MEAN", Top_split_var = NA, Top_split_threshold = NA,
             TP = NA, FP = NA, FN = NA, TN = NA,
             Precision = round(mean(cell_results$Precision, na.rm = TRUE), 3),
             Recall    = round(mean(cell_results$Recall,    na.rm = TRUE), 3),
             F1        = round(mean(cell_results$F1,        na.rm = TRUE), 3)),
  data.frame(Cell_line = "SD",   Top_split_var = NA, Top_split_threshold = NA,
             TP = NA, FP = NA, FN = NA, TN = NA,
             Precision = round(sd(cell_results$Precision, na.rm = TRUE), 3),
             Recall    = round(sd(cell_results$Recall,    na.rm = TRUE), 3),
             F1        = round(sd(cell_results$F1,        na.rm = TRUE), 3))
)
write.table(cell_summary, "Trees/Rules/_per_cell_line_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# also fit one tree on the pooled ALL file (all cell lines combined)
message("Fitting tree on ALL cell lines combined")
df_all <- load_all_file()
tab    <- table(df_all$is_correct)
message("  ", nrow(df_all), " variants: ",
        tab["TRUE"], " TRUE / ", tab["FALSE"], " FALSE")
fit_all <- fit_tree(df_all, "ALL")
all_fits[["ALL"]] <- fit_all
save_tree_outputs(fit_all, "pruned_1se")
save_tree_outputs(fit_all, "pruned_min")

# the summary blocks below loop over this vector, which now also includes ALL
entries <- c(cell_lines, "ALL")


###########################################################################################################################################################################
############################################################################# Cross-cell-line summary ####################################################################
###########################################################################################################################################################################

# variable importance side-by-side across cell lines + ALL (1-SE trees)
imp_long <- bind_rows(lapply(entries, function(cl) {
  tree <- all_fits[[cl]]$pruned_1se
  if (is.null(tree$variable.importance)) return(NULL)
  data.frame(Cell_line = cl,
             Feature   = names(tree$variable.importance),
             Importance = unname(tree$variable.importance))
}))

if (nrow(imp_long) > 0) {
  imp_long %>%
    mutate(Cell_line = factor(Cell_line, levels = entries)) %>%   # keep ALL last in legend
    ggplot(aes(x = Feature, y = Importance, fill = Cell_line)) +
    geom_col(position = "dodge") +
    coord_flip() +
    labs(title = "Decision tree variable importance per cell line + ALL (1-SE pruning)",
         x = "Feature", y = "rpart variable importance") +
    theme_minimal()
  ggsave("Trees/Plots/_variable_importance_across_cell_lines.png",
         plot = last_plot(), width = 12, height = 7, units = "in", dpi = 300)
}

# extract the top-level split from each tree — i.e. the single most important
# threshold the data picks first, per cell line (and for ALL). This is often the headline result.
top_splits <- bind_rows(lapply(entries, function(cl) {
  tree <- all_fits[[cl]]$pruned_1se
  if (nrow(tree$frame) == 1) {
    return(data.frame(Cell_line = cl, Top_split_var = "(no splits)", Top_split_threshold = NA))
  }
  splits <- as.data.frame(tree$splits)
  # the first row of $splits corresponds to the root split
  data.frame(Cell_line           = cl,
             Top_split_var       = rownames(splits)[1],
             Top_split_threshold = splits$index[1])
}))
write.table(top_splits, "Trees/Rules/_top_splits_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

message("\nDone. Outputs in:")
message("  Trees/Plots/   tree plots + cross-cell-line importance plot")
message("  Trees/Rules/   per-cell-line rule listings + _top_splits_summary.tsv")
message("  Trees/Rules/   _per_cell_line_summary.tsv (training-set evaluation — optimistic)")