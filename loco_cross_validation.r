library(data.table)
library(rpart)
library(rpart.plot)
library(dplyr)
library(tidyr)
library(ggplot2)

cell_lines   <- c("H2009","HCC1395","HCC1937","HCC1954","H1437","Hs578T","HG008","COLO829")
feature_cols <- c("SV_type","SS","SR","PR","PR_SR_ratio","ND","ND2","in_bed")

# adjust these two paths to your setup
input_dir <- "../SVs/SV_tables_ND2_filters"
setwd("../SVs")

dir.create("LOCO",        showWarnings = FALSE)
dir.create("LOCO/Plots",  showWarnings = FALSE)
dir.create("LOCO/Rules",  showWarnings = FALSE)


###########################################################################################################################################################################
########################################################################### Helper functions #############################################################################
###########################################################################################################################################################################

clean_variants <- function(df) {
  df %>%
    dplyr::rename(PR_SR_ratio = `PR/SR`) %>%
    mutate(SV_type    = factor(SV_type),
           in_bed     = as.logical(in_bed),
           is_correct = as.logical(is_correct))
}

load_cell_line <- function(cell_line) {
  as.data.frame(fread(file.path(input_dir,
    paste0(cell_line, "_SV_variants_table_annotated.tsv")))) %>% clean_variants()
}

fit_tree <- function(df, label) {
  set.seed(42)
  tree <- rpart(is_correct ~ .,
                data    = df[, c("is_correct", feature_cols)],
                method  = "class",
                parms   = list(prior = c(0.5, 0.5)),
                control = rpart.control(cp        = 0.01,
                                        minsplit  = 10,
                                        minbucket = 3,
                                        xval      = 10))
  cp_table <- as.data.frame(tree$cptable)
  min_idx  <- which.min(cp_table$xerror)
  cp_min   <- cp_table$CP[min_idx]
  cp_1se   <- cp_table$CP[which(cp_table$xerror <=
                 cp_table$xerror[min_idx] + cp_table$xstd[min_idx])[1]]
  list(full       = tree,
       pruned_min = prune(tree, cp = cp_min),
       pruned_1se = prune(tree, cp = cp_1se),
       cp_table   = cp_table,
       cp_min     = cp_min,
       cp_1se     = cp_1se,
       label      = label)
}

eval_tree <- function(tree, test_df) {
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

root_split <- function(tree) {
  if (nrow(tree$frame) == 1) return(list(var = "(no splits)", threshold = NA))
  splits <- as.data.frame(tree$splits)
  list(var = rownames(splits)[1], threshold = splits$index[1])
}

save_outputs <- function(fit, which_tree, test_df, held_out) {
  tree   <- fit[[which_tree]]
  label  <- fit$label
  cp_val <- if (which_tree == "pruned_1se") fit$cp_1se else fit$cp_min

  # plot
  png(file.path("LOCO/Plots", paste0(label, "_", which_tree, "_tree.png")),
      width = 12, height = 8, units = "in", res = 300)
  rpart.plot(tree,
             type        = 4,
             extra       = 104,
             box.palette = "RdYlGn",
             branch      = 0.3,
             tweak       = 1.1,
             main        = paste0(label, " — ", which_tree,
                                  " (cp = ", round(cp_val, 4), ")"))
  dev.off()

  # rules file
  sink(file.path("LOCO/Rules", paste0(label, "_", which_tree, "_rules.txt")))
  cat("Held-out cell line:", held_out, "\n")
  cat("Tree:               ", which_tree, "\n\n")
  cat("=== Decision rules ===\n")
  if (nrow(tree$frame) == 1) {
    cat("(no splits — root-only tree)\n")
  } else {
    print(rpart.rules(tree, cover = TRUE, nn = TRUE, style = "tall"))
  }
  cat("\n=== Variable importance ===\n")
  if (!is.null(tree$variable.importance)) {
    print(round(tree$variable.importance, 2))
  } else {
    cat("(none — root-only tree)\n")
  }
  cat("\n=== Performance on held-out cell line (", held_out, ") ===\n")
  perf <- eval_tree(tree, test_df)
  cat(sprintf("  TP=%d  FP=%d  FN=%d  TN=%d\n",
              perf$TP, perf$FP, perf$FN, perf$TN))
  cat(sprintf("  Precision=%.3f  Recall=%.3f  F1=%.3f\n",
              perf$Precision, perf$Recall, perf$F1))
  cat("\n=== CP table ===\n")
  print(fit$cp_table)
  cat("\ncp_min =", round(fit$cp_min, 4), "(minimum CV error)\n")
  cat("cp_1se =", round(fit$cp_1se, 4), "(1-SE rule — simpler tree)\n")
  sink()
}


###########################################################################################################################################################################
########################################################################### LOCO cross-validation ########################################################################
###########################################################################################################################################################################
# For each of the 8 cell lines:
#   - Hold it out completely (it is never seen during training)
#   - Train a tree on the other 7 cell lines
#   - Evaluate that tree on the held-out cell line
#   - Record Precision / Recall / F1
#
# Average the 8 scores across folds → honest performance estimate
# (what to expect when applying the model to a new, unseen cell line)
#
# Comparing root splits across folds → stability check
# (if all 8 folds pick the same feature/threshold, a universal rule is defensible)
###########################################################################################################################################################################

message("Loading per-cell-line data...")
combined <- bind_rows(lapply(cell_lines, function(cl) {
  load_cell_line(cl) %>% mutate(Cell_line = cl)
}))
message("  Total: ", nrow(combined), " variants across ", length(cell_lines), " cell lines")

loco_results <- data.frame()

for (held_out in cell_lines) {
  message("\nFold: holding out ", held_out)

  train_df <- combined %>% filter(Cell_line != held_out)
  test_df  <- combined %>% filter(Cell_line == held_out)

  tab_tr <- table(train_df$is_correct)
  tab_te <- table(test_df$is_correct)
  message("  train : ", nrow(train_df),
          " variants  (TRUE: ", tab_tr["TRUE"], " / FALSE: ", tab_tr["FALSE"], ")")
  message("  test  : ", nrow(test_df),
          " variants  (TRUE: ", tab_te["TRUE"], " / FALSE: ", tab_te["FALSE"], ")")

  fold_label <- paste0("LOCO_holdout-", held_out)
  fit <- fit_tree(train_df, fold_label)

  save_outputs(fit, "pruned_1se", test_df, held_out)
  save_outputs(fit, "pruned_min", test_df, held_out)

  rs   <- root_split(fit$pruned_1se)
  perf <- eval_tree(fit$pruned_1se, test_df)

  loco_results <- rbind(loco_results, data.frame(
    Held_out            = held_out,
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

  message("  Precision=", perf$Precision,
          "  Recall=",    perf$Recall,
          "  F1=",        perf$F1)
}


###########################################################################################################################################################################
########################################################################### Summary ######################################################################################
###########################################################################################################################################################################

mean_precision <- round(mean(loco_results$Precision, na.rm = TRUE), 3)
mean_recall    <- round(mean(loco_results$Recall,    na.rm = TRUE), 3)
mean_f1        <- round(mean(loco_results$F1,        na.rm = TRUE), 3)
sd_precision   <- round(sd(loco_results$Precision,   na.rm = TRUE), 3)
sd_recall      <- round(sd(loco_results$Recall,      na.rm = TRUE), 3)
sd_f1          <- round(sd(loco_results$F1,          na.rm = TRUE), 3)

# save summary table with a MEAN and SD row at the bottom
loco_summary <- rbind(
  loco_results,
  data.frame(Held_out = "MEAN", Top_split_var = NA, Top_split_threshold = NA,
             TP = NA, FP = NA, FN = NA, TN = NA,
             Precision = mean_precision, Recall = mean_recall, F1 = mean_f1),
  data.frame(Held_out = "SD",   Top_split_var = NA, Top_split_threshold = NA,
             TP = NA, FP = NA, FN = NA, TN = NA,
             Precision = sd_precision,   Recall = sd_recall,   F1 = sd_f1)
)
write.table(loco_summary, "LOCO/Rules/_loco_stability_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

message("\n=== LOCO performance estimate ===")
message("  Precision : ", mean_precision, " ± ", sd_precision)
message("  Recall    : ", mean_recall,    " ± ", sd_recall)
message("  Mean F1   : ", mean_f1,        " ± ", sd_f1)

# plot 1: Precision / Recall / F1 per fold with mean dashed lines
loco_results %>%
  pivot_longer(cols = c(Precision, Recall, F1),
               names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Held_out, y = Value, colour = Metric, group = Metric)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_hline(data = data.frame(
               Metric = c("Precision","Recall","F1"),
               Mean   = c(mean_precision, mean_recall, mean_f1)),
             aes(yintercept = Mean, colour = Metric),
             linetype = "dashed", linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "LOCO cross-validation: performance on held-out cell line (1-SE tree)",
       subtitle = paste0("Dashed lines = fold means  |  ",
                         "Mean F1 = ", mean_f1, " ± ", sd_f1),
       x = "Held-out cell line", y = "Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("LOCO/Plots/_loco_performance.png",
       plot = last_plot(), width = 10, height = 6, units = "in", dpi = 300)

# plot 2: which feature was the root split in each fold
loco_results %>%
  ggplot(aes(x = Held_out, y = Top_split_var, colour = Top_split_var)) +
  geom_point(size = 5) +
  labs(title = "LOCO stability: root split feature per fold (1-SE tree)",
       subtitle = "All dots on the same row = consistent root split across cell lines",
       x = "Held-out cell line", y = "Root split feature") +
  theme_minimal() +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave("LOCO/Plots/_loco_root_splits.png",
       plot = last_plot(), width = 10, height = 5, units = "in", dpi = 300)

message("\nDone. Outputs in:")
message("  LOCO/Plots/   one tree PNG per fold + _loco_performance.png + _loco_root_splits.png")
message("  LOCO/Rules/   one rule file per fold + _loco_stability_summary.tsv")