# 0. Packages and Data

library(tidyverse)
library(mice)
library(kknn)
library(cluster)

set.seed(42)

required_files <- c(
  "scores_3_neg.csv",
  "scores_3_pos.csv",
  "scores_4_neg.csv",
  "scores_4_pos.csv",
  "voters.csv",
  "approval.csv",
  "borda.csv"
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required files are missing from the working directory: ",
    paste(missing_files, collapse = ", ")
  )
}

scores_3_neg <- read.csv("scores_3_neg.csv")
scores_3_pos <- read.csv("scores_3_pos.csv")
scores_4_neg <- read.csv("scores_4_neg.csv")
scores_4_pos <- read.csv("scores_4_pos.csv")
voters       <- read.csv("voters.csv")
approval     <- read.csv("approval.csv")
borda        <- read.csv("borda.csv")

datasets <- list(
  scores_3_neg = scores_3_neg,
  scores_3_pos = scores_3_pos,
  scores_4_neg = scores_4_neg,
  scores_4_pos = scores_4_pos
)

scores_only <- function(df) {
  df %>% select(-id)
}

candidate_names <- names(scores_only(scores_3_neg))

# Check that all Range Voting files contain the same candidate columns
if (!all(map_lgl(datasets, ~ identical(names(scores_only(.x)), candidate_names)))) {
  stop("The four Range Voting files do not contain identical candidate columns.")
}

# General ranking function used throughout the script
compute_ranking <- function(imputed_df) {
  imputed_df %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
    pivot_longer(
      cols = everything(),
      names_to = "candidate",
      values_to = "total_score"
    ) %>%
    mutate(rank = rank(-total_score, ties.method = "min")) %>%
    arrange(rank, desc(total_score), candidate)
}


# 1. Descriptive Missing-Data Overview

na_summary <- function(df, dataset_name) {
  s <- scores_only(df)

  tibble(
    dataset = dataset_name,
    n_voters = nrow(s),
    n_candidates = ncol(s),
    total_cells = nrow(s) * ncol(s),
    n_na = sum(is.na(s)),
    pct_na = 100 * n_na / total_cells
  )
}

na_overview <- imap_dfr(datasets, na_summary)

candidate_na_all <- imap_dfr(
  datasets,
  function(df, dataset_name) {
    s <- scores_only(df)

    tibble(
      dataset = dataset_name,
      candidate = names(s),
      n_voters = nrow(s),
      n_na = colSums(is.na(s)),
      pct_na = 100 * n_na / n_voters
    ) %>%
      arrange(desc(n_na))
  }
)

print(na_overview)
print(candidate_na_all)


# 2. Imputation Strategies

# 2.1 Lowest-rank imputation
# I use the lowest rank/value for each dataset (e.g. -1 for the
# scores_3 datasets, 0 for the scores_4 datasets) to fill every
# missing value

impute_lowest <- function(df) {
  min_value <- min(as.matrix(df), na.rm = TRUE)

  df %>%
    mutate(across(everything(), ~ replace_na(.x, min_value)))
}


# 2.2 Complete-case deletion
# Every ballot with at least 1 missing value is deleted

impute_complete_case <- function(df) {
  df %>% drop_na()
}


# 2.3 Candidate-mean imputation
#I imputed the mean value for each dataset to every NA

impute_mean <- function(df) {
  df %>%
    mutate(
      across(
        everything(),
        ~ replace_na(.x, mean(.x, na.rm = TRUE))
      )
    )
}


# 2.4 Demographic cluster-based imputation
# I used actual demographics (age, gender, education,
# socio-professional category) for clustering, using Gower distance since
# these variables are categorical/ordinal rather than numeric, and PAM
# because it works directly on a distance matrix

impute_cluster_demographic <- function(scores_df, voters_df, k = 5) {
  required_demo_vars <- c("id", "age", "gender", "studies", "socpro")

  missing_demo_vars <- setdiff(required_demo_vars, names(voters_df))
  if (length(missing_demo_vars) > 0) {
    stop(
      "The following demographic variables are missing from voters.csv: ",
      paste(missing_demo_vars, collapse = ", ")
    )
  }

  merged <- scores_df %>%
    left_join(
      voters_df %>% select(all_of(required_demo_vars)),
      by = "id"
    )

  demo_vars <- merged %>%
    select(age, gender, studies, socpro)

  demo_vars[] <- lapply(demo_vars, function(x) {
    x <- as.character(x)
    x[is.na(x) | trimws(x) == ""] <- "missing"
    factor(x)
  })

  if (k < 2 || k >= nrow(demo_vars)) {
    stop("For cluster imputation, k must be at least 2 and smaller than the sample size.")
  }

  set.seed(42)

  gower_distance <- cluster::daisy(
    demo_vars,
    metric = "gower"
  )

  cluster_membership <- cluster::pam(
    gower_distance,
    k = k,
    diss = TRUE
  )$clustering

  candidate_cols <- setdiff(names(scores_df), "id")
  imputed_df <- scores_df

  global_means <- colMeans(
    scores_df[, candidate_cols, drop = FALSE],
    na.rm = TRUE
  )

  for (cluster_id in sort(unique(cluster_membership))) {
    rows_in_cluster <- which(cluster_membership == cluster_id)

    cluster_means <- colMeans(
      scores_df[rows_in_cluster, candidate_cols, drop = FALSE],
      na.rm = TRUE
    )

    invalid_means <- is.nan(cluster_means) | is.na(cluster_means)
    cluster_means[invalid_means] <- global_means[invalid_means]

    for (candidate in candidate_cols) {
      rows_to_impute <- rows_in_cluster[
        is.na(scores_df[rows_in_cluster, candidate])
      ]

      if (length(rows_to_impute) > 0) {
        imputed_df[rows_to_impute, candidate] <- cluster_means[[candidate]]
      }
    }
  }

  imputed_df
}


# 2.5 k-nearest-neighbour imputation

impute_knn <- function(df, k = 5) {
  imputed_df <- df
  columns_with_na <- which(colSums(is.na(df)) > 0)

  # kknn can't handle NAs in the predictor columns themselves, it just
  # drops any row with an NA anywhere, which shrank the usable
  # data far more than intended. I pre-fill NAs in the predictors with
  # the column mean purely so kknn has something to compute distances
  # on, the target column being predicted is left untouched here.
  predictor_copy <- df

  for (j in seq_along(predictor_copy)) {
    missing_rows <- is.na(predictor_copy[[j]])

    if (any(missing_rows)) {
      replacement_mean <- mean(df[[j]], na.rm = TRUE)

      if (is.nan(replacement_mean) || is.na(replacement_mean)) {
        stop("A candidate column contains no observed values and cannot be used for kNN.")
      }

      predictor_copy[missing_rows, j] <- replacement_mean
    }
  }

  for (column_index in columns_with_na) {
    target_name <- names(df)[column_index]
    rows_missing_target <- which(is.na(df[[target_name]]))
    rows_observed_target <- which(!is.na(df[[target_name]]))

    # Training rows must contain an observed target and complete predictors
    training_data <- df[rows_observed_target, , drop = FALSE] %>%
      drop_na()

    if (nrow(training_data) < 2) {
      warning(
        "Too few complete training rows for ", target_name,
        "; candidate mean was used as fallback."
      )
      imputed_df[rows_missing_target, target_name] <-
        mean(df[[target_name]], na.rm = TRUE)
      next
    }

    test_data <- predictor_copy[rows_missing_target, , drop = FALSE]
    test_data[[target_name]] <- 0
    test_data <- test_data[, names(training_data), drop = FALSE]

    effective_k <- min(k, nrow(training_data))

    prediction_formula <- as.formula(
      paste(
        target_name,
        "~",
        paste(setdiff(names(df), target_name), collapse = " + ")
      )
    )

    fitted_model <- kknn::kknn(
      formula = prediction_formula,
      train = training_data,
      test = test_data,
      k = effective_k,
      kernel = "rectangular"
    )

    imputed_df[rows_missing_target, target_name] <- fitted(fitted_model)
  }

  imputed_df
}


# 2.6 Multiple imputation with Rubin pooling
# Candidate mean scores are estimated separately in every completed
# dataset and pooled using Rubin's rules. Rankings are also retained
# separately because ranks themselves are discrete and are not directly
# pooled using standard Rubin formulas.

# Using mice with predictive mean matching (pmm) rather than a
# parametric method, since the score data isn't normally distributed
# and pmm only ever imputes values that were actually observed
# somewhere in the data. m = 5 and maxit = 5 follow van Buuren (2018);
# increasing either did not change the pooled rankings in testing, so I
# kept the smaller, faster setting.

impute_multiple <- function(df, m = 5, maxit = 5) {
  if (m < 2) {
    stop("Multiple imputation requires at least two completed datasets.")
  }

  set.seed(42)

  mice_object <- mice::mice(
    data = as.data.frame(df),
    m = m,
    maxit = maxit,
    method = "pmm",
    seed = 42,
    printFlag = FALSE
  )

  completed_datasets <- mice::complete(mice_object, action = "all")
  n_voters <- nrow(df)

  estimates_by_imputation <- imap_dfr(
    completed_datasets,
    function(completed_df, imputation_id) {
      mean_scores <- vapply(completed_df, mean, numeric(1))
      total_scores <- vapply(completed_df, sum, numeric(1))

      # Sampling variance of each candidate mean within this imputation.
      within_variances <- vapply(
        completed_df,
        function(x) stats::var(x) / length(x),
        numeric(1)
      )

      tibble(
        imputation = as.integer(imputation_id),
        candidate = names(completed_df),
        mean_score = unname(mean_scores),
        total_score = unname(total_scores),
        within_variance = unname(within_variances)
      ) %>%
        mutate(rank = rank(-mean_score, ties.method = "min")) %>%
        arrange(rank, desc(mean_score), candidate)
    }
  )

  pooled_results <- estimates_by_imputation %>%
    group_by(candidate) %>%
    summarise(
      pooled_mean = mean(mean_score),
      within_variance = mean(within_variance),
      between_variance = var(mean_score),
      min_rank = min(rank),
      max_rank = max(rank),
      rank_range = max_rank - min_rank,
      .groups = "drop"
    ) %>%
    mutate(
      total_variance = within_variance +
        (1 + 1 / m) * between_variance,
      pooled_se = sqrt(total_variance),
      rubin_df = if_else(
        between_variance <= .Machine$double.eps,
        Inf,
        (m - 1) *
          (1 + within_variance /
             ((1 + 1 / m) * between_variance))^2
      ),
      critical_value = if_else(
        is.infinite(rubin_df),
        qnorm(0.975),
        qt(0.975, df = rubin_df)
      ),
      lower_95 = pooled_mean - critical_value * pooled_se,
      upper_95 = pooled_mean + critical_value * pooled_se,
      pooled_total = pooled_mean * n_voters,
      pooled_rank = rank(-pooled_mean, ties.method = "min")
    ) %>%
    arrange(pooled_rank, desc(pooled_mean), candidate)

  winners_by_imputation <- estimates_by_imputation %>%
    filter(rank == 1) %>%
    select(
      imputation,
      winner = candidate,
      mean_score,
      total_score
    )

  ranking_stability <- estimates_by_imputation %>%
    group_by(candidate) %>%
    summarise(
      min_rank = min(rank),
      max_rank = max(rank),
      rank_range = max_rank - min_rank,
      .groups = "drop"
    ) %>%
    arrange(desc(rank_range), candidate)

  list(
    mice_object = mice_object,
    completed_datasets = completed_datasets,
    estimates_by_imputation = estimates_by_imputation,
    pooled_results = pooled_results,
    winners_by_imputation = winners_by_imputation,
    ranking_stability = ranking_stability
  )
}

# 3. Applying all six strategies to all four datasets

run_all_imputations <- function(df, dataset_name) {
  score_df <- scores_only(df)

  message("Processing main imputations: ", dataset_name)

  list(
    dataset = dataset_name,
    lowest_rank = impute_lowest(score_df),
    complete_case = impute_complete_case(score_df),
    mean_imp = impute_mean(score_df),
    cluster_imp = impute_cluster_demographic(
      scores_df = df,
      voters_df = voters,
      k = 5
    ) %>%
      select(-id),
    knn_imp = impute_knn(score_df, k = 5),
    multiple_imp = impute_multiple(score_df, m = 5)
  )
}

results <- imap(
  datasets,
  ~ run_all_imputations(df = .x, dataset_name = .y)
)


# 4. Candidate Rankings


strategy_names <- c(
  "lowest_rank",
  "complete_case",
  "mean_imp",
  "cluster_imp",
  "knn_imp",
  "multiple_imp"
)

rankings <- map(
  results,
  function(res) {
    map_dfr(
      strategy_names,
      function(strategy_name) {
        if (strategy_name == "multiple_imp") {
          res$multiple_imp$pooled_results %>%
            transmute(
              candidate,
              total_score = pooled_total,
              rank = pooled_rank,
              strategy = strategy_name,
              dataset = res$dataset
            )
        } else {
          compute_ranking(res[[strategy_name]]) %>%
            mutate(
              strategy = strategy_name,
              dataset = res$dataset
            )
        }
      }
    )
  }
)

rankings_df <- bind_rows(rankings) %>%
  select(dataset, strategy, candidate, total_score, rank) %>%
  arrange(dataset, strategy, rank)


# 5. MI Imputation Results
# Pooled candidate-score estimates, Rubin variances, confidence
# intervals, pooled ranks, and rank ranges across imputations.

mi_pooled_results <- imap_dfr(
  results,
  function(res, dataset_name) {
    res$multiple_imp$pooled_results %>%
      mutate(dataset = dataset_name, .before = 1)
  }
)

# Compact winner summary used in the results section (e.g., 5/5)
mi_winner_summary <- imap_dfr(
  results,
  function(res, dataset_name) {
    res$multiple_imp$winners_by_imputation %>%
      mutate(dataset = dataset_name, .before = 1)
  }
) %>%
  count(dataset, winner, name = "n_of_m_imputations") %>%
  group_by(dataset) %>%
  mutate(
    m = sum(n_of_m_imputations),
    winner_share = n_of_m_imputations / m
  ) %>%
  ungroup() %>%
  arrange(dataset, desc(n_of_m_imputations), winner)

print(mi_winner_summary)


# 6. Approval Voting and Borda Count Baselines

compute_baseline_ranking <- function(df, system_name) {
  df %>%
    select(-id) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
    pivot_longer(
      cols = everything(),
      names_to = "candidate",
      values_to = "total_score"
    ) %>%
    mutate(
      rank = rank(-total_score, ties.method = "min"),
      system = system_name
    ) %>%
    arrange(rank, desc(total_score), candidate)
}

approval_ranking <- compute_baseline_ranking(
  approval,
  "Approval Voting"
)

borda_ranking <- compute_baseline_ranking(
  borda,
  "Borda Count"
)

baseline_rankings <- bind_rows(
  approval_ranking,
  borda_ranking
)

baseline_winners <- baseline_rankings %>%
  filter(rank == 1) %>%
  select(system, winner = candidate, total_score)

print(baseline_winners)


# 7. Comparison of NA Strategies

winners <- rankings_df %>%
  filter(rank == 1) %>%
  select(dataset, strategy, winner = candidate, total_score)

rank_comparison <- rankings_df %>%
  select(dataset, strategy, candidate, rank) %>%
  pivot_wider(
    names_from = strategy,
    values_from = rank
  )

rank_sensitivity <- rankings_df %>%
  group_by(dataset, candidate) %>%
  summarise(
    min_rank = min(rank),
    max_rank = max(rank),
    rank_range = max_rank - min_rank,
    .groups = "drop"
  ) %>%
  arrange(desc(rank_range), dataset, candidate)

print(winners)
print(rank_sensitivity)


# 8. K-Sensitivity Analysis

#k = 5 in impute_cluster_demographic()/impute_knn() was an arbitrary
# starting choice, so I reran both imputations across a range of k
# values to see if something changes

run_cluster_sensitivity <- function(
    scores_df,
    voters_df,
    dataset_name,
    k_values = c(3, 5, 7, 9)
) {
  map_dfr(
    k_values,
    function(k_value) {
      message(
        "Cluster sensitivity: ", dataset_name,
        ", k = ", k_value
      )

      imputed_df <- impute_cluster_demographic(
        scores_df = scores_df,
        voters_df = voters_df,
        k = k_value
      ) %>%
        select(-id)

      compute_ranking(imputed_df) %>%
        mutate(
          dataset = dataset_name,
          method = "cluster",
          k = k_value,
          .before = 1
        )
    }
  )
}

run_knn_sensitivity <- function(
    scores_df,
    dataset_name,
    k_values = c(3, 5, 7, 9, 11)
) {
  score_df <- scores_only(scores_df)

  map_dfr(
    k_values,
    function(k_value) {
      message(
        "kNN sensitivity: ", dataset_name,
        ", k = ", k_value
      )

      imputed_df <- impute_knn(
        score_df,
        k = k_value
      )

      compute_ranking(imputed_df) %>%
        mutate(
          dataset = dataset_name,
          method = "knn",
          k = k_value,
          .before = 1
        )
    }
  )
}

cluster_k_results <- imap_dfr(
  datasets,
  ~ run_cluster_sensitivity(
    scores_df = .x,
    voters_df = voters,
    dataset_name = .y,
    k_values = c(3, 5, 7, 9)
  )
)

knn_k_results <- imap_dfr(
  datasets,
  ~ run_knn_sensitivity(
    scores_df = .x,
    dataset_name = .y,
    k_values = c(3, 5, 7, 9, 11)
  )
)

# Compact summary: winner and complete ranking for every tested k
cluster_k_summary <- cluster_k_results %>%
  group_by(dataset, k) %>%
  summarise(
    method = "cluster",
    winner = candidate[which.min(rank)],
    full_ranking = paste(candidate[order(rank, candidate)], collapse = " > "),
    .groups = "drop"
  ) %>%
  select(dataset, method, k, winner, full_ranking)

knn_k_summary <- knn_k_results %>%
  group_by(dataset, k) %>%
  summarise(
    method = "knn",
    winner = candidate[which.min(rank)],
    full_ranking = paste(candidate[order(rank, candidate)], collapse = " > "),
    .groups = "drop"
  ) %>%
  select(dataset, method, k, winner, full_ranking)

k_sensitivity_summary <- bind_rows(
  cluster_k_summary,
  knn_k_summary
) %>%
  arrange(method, dataset, k)

print(k_sensitivity_summary)


# 9. Export Results

# Descriptive missing-data results
write_csv(na_overview, "na_overview.csv")
write_csv(candidate_na_all, "candidate_missingness_all_datasets.csv")

# Main results across all six missing-data strategies
write_csv(rankings_df, "rankings_all_strategies.csv")
write_csv(rank_comparison, "rank_comparison_all_strategies.csv")
write_csv(rank_sensitivity, "rank_sensitivity.csv")

# Missing-value-free baseline systems
write_csv(baseline_rankings, "baseline_rankings.csv")

# Multiple-imputation results
write_csv(
  mi_pooled_results,
  "mi_pooled_results_rubins_rules.csv"
)
write_csv(
  mi_winner_summary,
  "mi_winner_summary.csv"
)

# Combined k-sensitivity results for cluster and kNN imputation
write_csv(
  k_sensitivity_summary,
  "k_sensitivity_summary.csv"
)

message("Done. All analyses and required exports were completed successfully.")