## =================================================================================================
## Plot SVD results
## =================================================================================================
set.seed(51681684)

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggsignif)
library(ModelMetrics)

source(file.path("Utils", "rundata_utils.R"))
source(file.path("Utils", "plot_utils.R"))


## Constants
# Overall theme
PLOT_THEME <- theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        strip.text = element_text(size = 6, margin = margin(2.5, 2.5, 2.5, 2.5)),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.6, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.margin = margin(t = 2.5, r = 5.5, b = 2.5, l = 0),
        panel.grid = element_blank(),
        plot.background = element_blank())

theme_set(PLOT_THEME)

# Order here determines plotting order:
RUN_NAMES <- c(AllGenomeFeatures_SVD = "Genome composition",
							 AllGenomeFeatures_and_SVD_clover_rank12 = "Genome composition\n+\nObserved network",
							 AllGenomeFeatures_and_SVD_trefle_rank12 = "Genome composition\n+\nImputed network")
RUN_IDS <- names(RUN_NAMES)

BEST_RUN_ID <- "SVD_trefle_rank12"
FULL_RUN_ID <- "AllGenomeFeatures_and_SVD_trefle_rank12"

## Raw clover data:
clover <- read_csv(file.path("ExternalData", "clover.csv"))

## Read predictions:
predictions_test_boot <- read_multiple_runs(RUN_IDS, "_Predictions.rds") %>% 
	filter(.data$Dataset == "test")

predictions_bagged <- lapply(RUN_IDS, function(x) 
  readRDS(file.path("RunData", "BaggedModels", sprintf("%s_Bagged_predictions.rds", x))))

cutoffs_bagged <- sapply(RUN_IDS, function(x) 
  readRDS(file.path("RunData", "BaggedModels", sprintf("%s_Bagging_cutoff.rds", x))))

names(predictions_bagged) <- RUN_IDS
names(cutoffs_bagged) <- RUN_IDS

## Read combined model (used for variable importance comparison)
model_fits <- read_run_rds(FULL_RUN_ID, "_ModelFits.rds")
iteration_data <- read_run_rds(FULL_RUN_ID, "_CalculatedData.rds")


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Distribution of AUC values -----------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Clean up and merge predictions
combined_test_boot <- predictions_test_boot %>% 
	filter(.data$Dataset == "test") %>% 
	mutate(InfectsHumans = .data$InfectsHumans == "True",
				 RunName = RUN_NAMES[.data$RunID],
				 RunName = factor(.data$RunName, levels = RUN_NAMES))

predictions_bagged <- predictions_bagged %>% 
  bind_rows(.id = "RunID") %>% 
  mutate(InfectsHumans = .data$InfectsHumans == "True",
         RunName = RUN_NAMES[.data$RunID],
         RunName = factor(.data$RunName, levels = RUN_NAMES))


# Calculate AUC
auc_test_boot <- combined_test_boot %>% 
	group_by(.data$RunName, .data$Iteration) %>% 
	summarise(AUC = auc(actual = .data$InfectsHumans, predicted = .data$RawScore))

auc_bagged <- predictions_bagged %>% 
  group_by(.data$RunName) %>% 
  summarise(AUC = auc(actual = .data$InfectsHumans, predicted = .data$BagScore))


# Plot
auc_plot <- ggplot(auc_test_boot, aes(x = RunName, y = AUC)) +
  geom_violin(colour = NA, fill = "grey90") +
	geom_boxplot(colour = "grey40", width = 0.4, fill = NA) +
  
  geom_point(shape = 23, size = 2, colour = "grey40", fill = "#4477AA", data = auc_bagged) +
	
	scale_y_continuous(breaks = seq(0, 1, by = 0.05), limits = c(0.5, NA),
	                   expand = expansion(add = c(0.0015, 0.03))) +
	xlab("Feature set") +
	theme(plot.margin = margin(t = 14, r = 5.5, b = 5.5, l = 5.5),
				panel.border = element_rect(fill = NA, size = 0.3)) # Turning off clipping makes this appear thicker than in other plots


# Add p-values for all comparisons:
P_VAL_COMPARISONS <- list(
  Level1 = list(c("Genome composition", "Genome composition\n+\nObserved network")),
  Level2 = list(c("Genome composition", "Genome composition\n+\nImputed network")),
  Level3 = list(c("Genome composition\n+\nObserved network", "Genome composition\n+\nImputed network"))
)

P_VAL_HEIGHTS = c(Level1 = 0.96, Level2 = 0.99, Level3 = 1.02)


for (lvl in names(P_VAL_COMPARISONS)) {
  p_vals <- c()
  
  for (comp in P_VAL_COMPARISONS[[lvl]]) {
    dta <- auc_test_boot %>% 
      filter(.data$RunName %in% comp)
    
    test_res <- kruskal.test(x = dta$AUC, g = dta$RunName)
    p_val <- test_res$p.value
    p_val <- if_else(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
    
    auc_plot <- auc_plot +
      geom_signif(comparisons = list(comp), 
                  y_position = P_VAL_HEIGHTS[[lvl]],
                  annotations = p_val,
                  tip_length = 0.01, textsize = 2.5)
  }
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Top features, by feature set ---------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Get variable importance
#  - Variable importance / effect on train set predictions in each run (using training set, so these are 
#    conceptually similar to coefficients in a glm)
#  - Reporting the variation across iterations, not across viruses
#  			- Thus, for each level of summary reported (e.g. variable type, variable cluster, etc.),
#  			  calculate average importance across all train viruses in a given iteration
#  			- Then report mean + sd of variation of this total importance across runs...
train_data_splits <- lapply(iteration_data, function(x) x$train)
used_features <- colnames(train_data_splits$`1`)
used_features <- used_features[used_features != "InfectsHumans"]

# Get per-virus variable importance:
varimp_raw <- mapply(FUN = get_shap_contributions, 
										 trainedModel = model_fits, newdata = train_data_splits,
										 SIMPLIFY = FALSE)


## Calculate importance across all models
#  - Since we are primarily interested in the bagged model, calculate the mean effect size
#    across all iterations
featureset_importance <- bind_rows(varimp_raw) %>% 
	id_variable_types("Feature") %>% 
	group_by(.data$Feature, .data$VariableType, .data$Gene) %>% 
	summarise(mean_importance = mean(abs(.data$SHAP))) %>% 
	
	ungroup() %>% 
	mutate(Rank = rank(-.data$mean_importance, ties.method = "random", na.last = "keep"),
				 Rank = ifelse(is.na(.data$Rank), Inf, .data$Rank))


## Clean up names for plot:
featureset_importance <- featureset_importance %>% 
	ungroup() %>% 
	id_feature_set() %>% 
  
  mutate(SetName = case_when(.data$VariableType == "svd" ~ "Imputed\nnetwork",
                             .data$VariableType %in% c("GenomicDensity", "VirusDirect") ~ "Genome\ncomposition",
                             TRUE ~ .data$VariableType)) %>% 
	mutate(SetName = factor(.data$SetName, levels = c("Genome\ncomposition", "Imputed\nnetwork")),
				 SetName_numeric = as.numeric(.data$SetName),
				 SetName_numeric = .data$SetName_numeric - min(.data$SetName_numeric) + 1)


featureset_imp_plot <- ggplot(featureset_importance, aes(x = SetName, y = Rank, colour = SetName)) +
	geom_blank() +  # stops next layer from converting x-axis to numeric
	geom_segment(aes(x = SetName_numeric - 0.5, xend = SetName_numeric + 0.5, yend = Rank), colour = "grey80") +
	
	geom_boxplot(outlier.shape = NA, size = 1, fill = "white", alpha = 0.3) +
	
	scale_colour_brewer(palette = "Set2", guide = FALSE) +
	scale_y_reverse() +
	labs(x = "Feature type", y = "Feature rank")


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Compare predictions ------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
best_detection <- clover %>% 
  filter(.data$Host == "Homo sapiens") %>% 
  group_by(.data$Virus) %>% 
  summarise(DetectionMethod = case_when(any(.data$Detection_Isolation) ~ "Isolation",
                                        any(.data$Detection_Genetic) ~ "Genetic",
                                        any(.data$Detection_Serology) ~ "Serology",
                                        any(.data$Detection_NotSpecified) ~ "Not specified"))

# Plot virus ranks
bagged_final <- predictions_bagged %>% 
  filter(.data$RunID == FULL_RUN_ID) %>% 
  left_join(best_detection, by = c("LatestSppName" = "Virus")) %>% 
  mutate(DetectionMethod = if_else(is.na(.data$DetectionMethod), "No detection", .data$DetectionMethod),
         DetectionMethod = factor(.data$DetectionMethod, 
                                  levels = rev(c("Isolation", "Genetic", "Serology", 
                                                 "Not specified", "No detection")))) %>% 
  arrange(.data$BagScore) %>% 
  mutate(LatestSppName = factor(.data$LatestSppName, levels = .data$LatestSppName))


y_cutoff <- cutoffs_bagged[FULL_RUN_ID]

x_cutoff <- bagged_final %>% 
  filter(.data$BagScore < y_cutoff) %>% 
  top_n(1, wt = .data$BagScore) %>% 
  pull(.data$LatestSppName) %>% 
  as.numeric() + 0.5

p_ranks <- ggplot(bagged_final, aes(x = LatestSppName, y = BagScore, colour = InfectsHumans)) +
  geom_linerange(aes(ymin = BagScore_Lower, ymax = BagScore_Upper)) +
  geom_step(group = 1, colour = "grey20", size = 1) +
  
  geom_hline(yintercept = y_cutoff, linetype = 2, colour = "grey20") +
  geom_vline(xintercept = x_cutoff, linetype = 2, colour = "grey20") +
  
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_colour_manual(breaks = c("TRUE", "FALSE"),
                      values = c("TRUE" = "#EE6677", "FALSE" = "#4477AA"),
                      labels = c("TRUE" = "True", "FALSE" = "False")) +
  labs(x = NULL, y = "Predicted probability", colour = "Strong evidence\nof human infection") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = -5.5, l = 5.5))

p_detection <- ggplot(bagged_final, aes(x = LatestSppName, y = DetectionMethod, fill = InfectsHumans)) +
  geom_tile() +
  geom_vline(xintercept = x_cutoff, linetype = 2, colour = "grey20") +
  scale_fill_manual(breaks = c("TRUE", "FALSE"),
                    values = c("TRUE" = "#EE6677", "FALSE" = "#4477AA"),
                    labels = c("TRUE" = "True", "FALSE" = "False"),
                    guide = FALSE) +
  labs(x = "Virus species", y = "Detection method") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_combined_ranks <- plot_grid(p_ranks, p_detection, nrow = 2, align = "v", axis = "lr",
                              rel_heights = c(2, 1))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Combine and save ---------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
main_plot <- plot_grid(auc_plot, featureset_imp_plot,
                       ncol = 2, rel_widths = c(3, 1),
                       align = "h", axis = "tb",
                       labels = c("A", "B"))


dir.create("Plots")
ggsave2("Plots/human_models_main.pdf", main_plot, width = 7.5, height = 3.5)
ggsave2("Plots/human_models_virus_ranks.pdf", p_combined_ranks, width = 7.5, height = 3.8)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Saved cleaned table of bagged predictions --------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Add: zoonotic/human/other, svd name, prediction, (category?)
name_data <- readRDS(file.path("CalculatedData", "SplitData_Training.rds")) %>% 
	distinct(.data$UniversalName, .data$LatestSppName)

svd_names <- read_csv(file.path("InternalData", "svd_name_matching.csv"))

clean_bagged <- predictions_bagged %>% 
	left_join(name_data, by = "UniversalName") %>% 
	left_join(svd_names, by = "LatestSppName") %>% 
	arrange(-.data$BagScore)

stopifnot(!any(is.na(clean_bagged$svd_name)))

clean_bagged %>% 
	mutate(observed = .data$InfectsHumans == "True",
				 predicted = .data$BagScore > 0.5) %>% 
	select(virus = .data$svd_name, 
				 .data$observed,
				 .data$predicted,
				 probability = .data$BagScore,
				 probability_lower = .data$BagScore_Lower,
				 probability_upper = .data$BagScore_Upper) %>% 
	write_excel_csv("svd_zoonotic_rank_predictions.csv")


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Print AUC values ---------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cat("\nBootstrap AUCs:\n")
auc_test_boot %>% 
	group_by(.data$RunName) %>% 
	summarise(mean = mean(.data$AUC),
						sd = sd(.data$AUC)) %>% 
	print()

cat("\nBagged AUC for best model:\n")
predictions_bagged %>% 
	summarise(AUC = auc(actual = .data$InfectsHumans == "True", predicted = .data$BagScore)) %>% 
	print()

cat("\nSensitivity/specificity for best model (cutoff = 0.5):\n")
predictions_bagged %>% 
	mutate(Prediction = if_else(.data$BagScore > 0.5, "True", "False")) %>% 
	group_by(.data$InfectsHumans) %>% 
	summarise(proportion_accurate = sum(.data$InfectsHumans == .data$Prediction)/n()) %>% 
	print()
