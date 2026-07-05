
install.packages("kknn")
library(kknn)
library(tidyverse)
library(mice)
library(kknn)
install.packages("VIM")
library(VIM)


# 0. LOAD DATA

# Adjust paths if needed
scores_3_neg <- read.csv("scores_3_neg.csv")
scores_3_pos <- read.csv("scores_3_pos.csv")
scores_4_neg <- read.csv("scores_4_neg.csv")
scores_4_pos <- read.csv("scores_4_pos.csv")

# Working on all four datasets with the same pipeline.

datasets <- list(
  scores_3_neg = scores_3_neg,
  scores_3_pos = scores_3_pos,
  scores_4_neg = scores_4_neg,
  scores_4_pos = scores_4_pos
)

# Helper: separate id column from score columns
scores_only <- function(df) df %>% select(-id)
candidates   <- names(scores_only(scores_3_neg))  # same for all datasets


# 1. EXPLORATORY: NA OVERVIEW


na_summary <- function(df, name) {
  s <- scores_only(df)
  tibble(
    dataset    = name,
    n_voters   = nrow(s),
    total_cells = nrow(s) * ncol(s),
    n_na       = sum(is.na(s)),
    pct_na     = round(100 * sum(is.na(s)) / (nrow(s) * ncol(s)), 1)
  )
}

na_overview <- imap_dfr(datasets, na_summary)
print(na_overview)

# Per-candidate NA counts (using scores_3_neg as example)
candidate_na <- scores_only(scores_3_neg) %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "candidate", values_to = "n_na") %>%
  arrange(desc(n_na))
print(candidate_na)


# ============================================================
# 2. IMPUTATION FUNCTIONS
# ============================================================

# Each function takes a score-only data frame (no id column)
# and returns a completed data frame with NAs filled.

# --- 2.1 Lowest-Rank Imputation ---
# Assign minimum value of the scale to all NAs
impute_lowest <- function(df) {
  min_val <- min(df, na.rm = TRUE)
  df %>% mutate(across(everything(), ~replace_na(., min_val)))
}

# --- 2.2 Complete-Case Deletion ---
# Remove all rows with any NA
impute_complete_case <- function(df) {
  df %>% drop_na()
}

# --- 2.3 Mean Imputation ---
# Replace each NA with the column (candidate) mean
impute_mean <- function(df) {
  df %>% mutate(across(everything(), ~replace_na(., mean(., na.rm = TRUE))))
}

# --- 2.4 Cluster-Based Imputation ---
# Group voters by their overall response pattern (k-means),
# then impute within each cluster using the cluster mean.
# Note: requires voters.csv for demographic clusters if available.
# Here we use response-pattern clustering as a proxy.
impute_cluster <- function(df, k = 5) {

  # Nur numerische Daten
  df_num <- as.data.frame(lapply(df, as.numeric))

  # Vollständige Fälle für Clustering
  complete_idx <- complete.cases(df_num)
  df_complete <- df_num[complete_idx, , drop = FALSE]

  # Sicherheitscheck
  if (nrow(df_complete) < k) {
    stop("Weniger vollständige Fälle als Cluster vorhanden.")
  }

  set.seed(42)

  km <- kmeans(
    df_complete,
    centers = k,
    nstart = 25
  )

  df_imputed <- df_num

  incomplete_idx <- which(!complete_idx)

  for (i in incomplete_idx) {

    row <- as.numeric(df_num[i, ])

    obs_cols <- which(!is.na(row))
    na_cols  <- which(is.na(row))

    if (length(na_cols) == 0) next

    # Falls eine Zeile komplett leer ist
    if (length(obs_cols) == 0) {

      df_imputed[i, na_cols] <-
        colMeans(df_complete, na.rm = TRUE)[na_cols]

    } else {

      centers_obs <- km$centers[, obs_cols, drop = FALSE]

      dists <- apply(
        centers_obs,
        1,
        function(center) {
          sum((row[obs_cols] - center)^2)
        }
      )

      nearest <- which.min(dists)

      replacement <-
        as.numeric(
          km$centers[nearest, na_cols, drop = FALSE]
        )

      df_imputed[i, na_cols] <- replacement
    }
  }

  return(df_imputed)
}
# --- 2.5 kNN Imputation ---
# Uses kknn package to impute each column via k nearest neighbours
impute_knn <- function(df, k = 5) {

  imp <- VIM::kNN(
    df,
    k = k,
    imp_var = FALSE
  )

  as_tibble(imp)
}
s <- scores_only(scores_3_neg)

test_knn <- impute_knn(s)

print(class(test_knn))
print(dim(test_knn))
print(sum(is.na(test_knn)))

# --- 2.6 Multiple Imputation (mice) ---
# Generates m=5 imputed datasets and pools candidate scores
impute_multiple <- function(df, m = 5) {
  imp <- mice(df, m = m, method = "pmm", seed = 42, printFlag = FALSE)

  # Pool: average imputed values across m datasets
  completed_list <- lapply(1:m, function(i) complete(imp, i))
  pooled <- Reduce("+", completed_list) / m
  as_tibble(pooled)
}
s <- scores_only(scores_3_neg)

test_knn <- impute_knn(s)

dim(test_knn)
sum(is.na(test_knn))

# ============================================================
# 3. APPLY ALL 6 STRATEGIES TO ALL 4 DATASETS
# ============================================================

run_all_imputations <- function(df, dataset_name) {
  s <- scores_only(df)

  message("Processing: ", dataset_name)

  list(
    dataset       = dataset_name,
    lowest_rank   = impute_lowest(s),
    complete_case = impute_complete_case(s),
    mean_imp      = impute_mean(s),
    cluster_imp   = impute_cluster(s),
    knn_imp       = impute_knn(s),
    multiple_imp  = impute_multiple(s)
  )
}

results <- imap(datasets, run_all_imputations)

# ============================================================
# 4. COMPUTE CANDIDATE RANKINGS
# ============================================================

# Range Voting aggregation: sum of scores per candidate
# (missing values already handled by imputation)

compute_ranking <- function(imputed_df) {
  imputed_df %>%
    summarise(across(everything(), sum, na.rm = TRUE)) %>%
    pivot_longer(everything(), names_to = "candidate", values_to = "total_score") %>%
    arrange(desc(total_score)) %>%
    mutate(rank = row_number())
}

# Extract rankings for each dataset × strategy
rankings <- map(results, function(res) {
  strategies <- c(
  "lowest_rank",
  "complete_case",
  "mean_imp",
  "cluster_imp",
  "knn_imp",
  "multiple_imp"
)
  map_dfr(strategies, function(strat) {
    compute_ranking(res[[strat]]) %>%
      mutate(strategy = strat, dataset = res$dataset)
  })
})

rankings_df <- bind_rows(rankings)


# 4b. Baseline Systmens: Aproval Voting and Borda Count


approval <- read.csv("approval.csv")
borda    <- read.csv("borda.csv")

compute_baseline_ranking <- function(df, system_name) {
  df %>%
    select(-id) %>%
    summarise(across(everything(), sum, na.rm = TRUE)) %>%
    pivot_longer(everything(), names_to = "candidate", values_to = "total_score") %>%
    arrange(desc(total_score)) %>%
    mutate(
      rank = row_number(),
      system = system_name
    )
}

approval_ranking <- compute_baseline_ranking(approval, "Approval Voting")
borda_ranking    <- compute_baseline_ranking(borda, "Borda Count")

baseline_rankings <- bind_rows(
  approval_ranking,
  borda_ranking
)

baseline_winners <- baseline_rankings %>%
  filter(rank == 1) %>%
  select(system, winner = candidate, total_score)

print(baseline_winners)

write_csv(baseline_rankings, "baseline_rankings.csv")
write_csv(baseline_winners, "baseline_winners.csv")

# ============================================================
# 5. COMPARE RANKINGS ACROSS STRATEGIES
# ============================================================

# 5.1 Winner per strategy × dataset
winners <- rankings_df %>%
  filter(rank == 1) %>%
  select(dataset, strategy, winner = candidate, total_score)

print(winners)

# 5.2 Rank shifts: compare each strategy to lowest_rank as reference
rank_comparison <- rankings_df %>%
  select(dataset, strategy, candidate, rank) %>%
  pivot_wider(names_from = strategy, values_from = rank)

print(rank_comparison)

# 5.3 Maximum rank shift per candidate across strategies
rank_sensitivity <- rankings_df %>%
  group_by(dataset, candidate) %>%
  summarise(
    min_rank  = min(rank),
    max_rank  = max(rank),
    rank_range = max(rank) - min(rank),
    .groups = "drop"
  ) %>%
  arrange(desc(rank_range))

print(rank_sensitivity)


# ============================================================
# 6. VISUALISATIONS
# ============================================================

# 6.1 Ranking heatmap per dataset
plot_ranking_heatmap <- function(dataset_name) {
  rankings_df %>%
    filter(dataset == dataset_name) %>%
    mutate(
      candidate = candidate %>%
        str_replace_all("\\.", " ") %>%
        str_replace("Jean Luc Mélenchon", "Mélenchon") %>%
        str_replace("Emmanuel Macron", "Macron") %>%
        str_replace("Marine Le Pen", "Le Pen"),

      strategy = factor(
        strategy,
        levels = c(
          "lowest_rank",
          "complete_case",
          "mean_imp",
          "cluster_imp",
          "knn_imp",
          "multiple_imp"
        )
      )
    ) %>%
    ggplot(aes(x = strategy, y = reorder(candidate, -rank), fill = rank)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = rank), size = 3.5, colour = "black") +
    scale_fill_gradient(
      low = "#F6F1E9",
      high = "#8C7B75",
      name = "Rank"
    ) +
    labs(
      title    = paste("Candidate Rankings across Imputation Methods –", dataset_name),
      subtitle = "Lower rank = higher score",
      x        = "Imputation Method",
      y        = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title  = element_text(face = "bold"),
      panel.grid  = element_blank()
    )
}

# Plot for each dataset
walk(names(datasets), function(d) {
  p <- plot_ranking_heatmap(d)
  ggsave(paste0("ranking_heatmap_", d, ".png"), p, width = 10, height = 6, dpi = 150)
  message("Saved: ranking_heatmap_", d, ".png")
})

# 6.2 Rank sensitivity plot (range of ranks per candidate)
p_sensitivity <- rank_sensitivity %>%
  filter(dataset == "scores_3_neg") %>%
  mutate(candidate = str_trunc(candidate, 15)) %>%
  ggplot(aes(x = reorder(candidate, rank_range), y = rank_range)) +
  geom_col(fill = "#38BDF8", alpha = 0.85) +
  geom_text(aes(label = rank_range), hjust = -0.2, size = 3.5) +
  coord_flip() +
  labs(
    title    = "Rank Sensitivity by Candidate (scores_3_neg)",
    subtitle = "How much does a candidate's rank vary across imputation strategies?",
    x        = NULL,
    y        = "Rank range (max – min rank)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("rank_sensitivity.png", p_sensitivity, width = 9, height = 6, dpi = 150)
message("Saved: rank_sensitivity.png")

# 6.3 Winner stability table
p_winners <- winners %>%
  mutate(strategy = factor(strategy, levels = c("lowest_rank", "complete_case",
                                                 "mean_imp", "cluster_imp",
                                                 "knn_imp", "multiple_imp"))) %>%
  ggplot(aes(x = strategy, y = dataset, fill = winner)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(
  aes(label = gsub("Jean.Luc.Mélenchon", "Mélenchon", winner)),
  size = 3.5,
  colour = "black"
) +
  scale_fill_manual(
  values = c("Jean.Luc.Mélenchon" = "#CFC6B8"),
  name = "Winner"
) +
  labs(
    title = "Winner by Imputation Strategy and Dataset",
    x     = "Imputation Strategy",
    y     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title  = element_text(face = "bold"),
    panel.grid  = element_blank()
  )

ggsave("winner_stability.png", p_winners, width = 9, height = 4, dpi = 150)
message("Saved: winner_stability.png")


# ============================================================
# 7. EXPORT RESULTS
# ============================================================

write_csv(rankings_df,      "rankings_all_strategies.csv")
write_csv(winners,          "winners_by_strategy.csv")
write_csv(rank_sensitivity, "rank_sensitivity.csv")

message("Done. All results saved.")
