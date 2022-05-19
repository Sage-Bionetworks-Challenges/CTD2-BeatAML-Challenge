library(tidyverse)

data.dir <- "../Data/"

# Calculate three kernels between the following specimens
# - [train x train] 
# - [train x validation] 
# - [train x leaderboard]

# Jaccard similarity is calculated based on binary mutation profiles.
# Additional mutations from categorical clinical data added.

df_mutation_extra_vars <- 
  readRDS(
    "mutation_vars_from_categorical_clinical_data.rds"
  ) %>% 
  tidyr::gather(
    key = Hugo_Symbol,
    value = value,
    -.data$lab_id, -.data$dataset
  ) %>% 
  dplyr::filter(
    .data$value == 1
  ) %>% 
  dplyr::select(
    -.data$value
  )


# --------------------------------------------------------------------------------------------------------------
# Read training mutation data ----------------------------------------------------------------------------------
df_mutation_train <- 
  readr::read_csv(
    paste0(data.dir, "training/dnaseq.csv"),
    col_names = TRUE
  ) %>% 
  dplyr::bind_rows(
    df_mutation_extra_vars %>% 
      dplyr::filter(
        .data$dataset == "training"
      )
  )
names(df_mutation_train)

df_mutation_train$lab_id %>% n_distinct()       # 211
df_mutation_train$Hugo_Symbol %>% n_distinct()  # 245

# Specimens in the training AUC dataset (labels)
label_specimens_train <- 
  readr::read_csv(
    paste0(data.dir, "training/aucs.csv"),
    col_names = TRUE
  ) %>% 
  dplyr::pull(.data$lab_id) %>% 
  unique()
# overlap
intersect(label_specimens_train, df_mutation_train$lab_id) %>% length() == length(label_specimens_train) # FALSE
intersect(label_specimens_train, df_mutation_train$lab_id) %>% length()  # 211
# Considered mutations not present for 2 training specimens 


# Read validation mutation data ----------------------------------------------------------------------------
df_mutation_validation <- 
  readr::read_csv(
    paste0(data.dir, "validation/dnaseq.csv"),
    col_names = TRUE
  ) %>% 
  dplyr::bind_rows(
    df_mutation_extra_vars %>% 
      dplyr::filter(
        .data$dataset == "validation"
      )
  )
names(df_mutation_validation)

df_mutation_validation$lab_id %>% n_distinct()       # 63
df_mutation_validation$Hugo_Symbol %>% n_distinct()  # 65 

intersect(df_mutation_validation$Hugo_Symbol, df_mutation_train$Hugo_Symbol) %>% length()   # 48

# Specimens in the validation AUC dataset (labels)
label_specimens_validation <- 
  readr::read_csv(
    paste0(data.dir, "validation/aucs.csv"),
    col_names = TRUE
  ) %>% 
  dplyr::pull(.data$lab_id) %>% 
  unique()
# overlap
intersect(label_specimens_validation, df_mutation_validation$lab_id) %>% length() == length(label_specimens_train) # FALSE
intersect(label_specimens_validation, df_mutation_validation$lab_id) %>% length() # 63
# No mutations for 2 validation specimens


# Read leaderboard mutation data ----------------------------------------------------------------------------
df_mutation_leaderboard <- 
  readr::read_csv(
    paste0(data.dir, "leaderboard/dnaseq.csv"),
    col_names = TRUE
  ) %>% 
  dplyr::bind_rows(
    df_mutation_extra_vars %>% 
      dplyr::filter(
        .data$dataset == "leaderboard"
      )
  )
names(df_mutation_leaderboard)

df_mutation_leaderboard$lab_id %>% n_distinct()       # 79
df_mutation_leaderboard$Hugo_Symbol %>% n_distinct()  # 152

intersect(df_mutation_leaderboard$Hugo_Symbol, df_mutation_train$Hugo_Symbol) %>% length()        # 88
intersect(df_mutation_leaderboard$Hugo_Symbol, df_mutation_validation$Hugo_Symbol) %>% length()   # 41

# Specimens in the leaderboard AUC dataset (labels)
label_specimens_leaderboard <- 
  readr::read_csv(
    paste0(data.dir, "leaderboard/aucs.csv"),
    col_names = TRUE
  ) %>% 
  dplyr::pull(.data$lab_id) %>% 
  unique()
# overlap
intersect(label_specimens_leaderboard, df_mutation_leaderboard$lab_id) %>% length() == length(label_specimens_train) # FALSE
intersect(label_specimens_leaderboard, df_mutation_leaderboard$lab_id) %>% length() # 79
# No mutations for 1 leaderboard specimen


# Process to [specimen x mutation] format ------------------------------------------------------------------
# - Consider only mutations in the training data
# - Include specimens that are not present in the mutation data (with mutation status set to 0)

df_mutation_train <- 
  # Include missing specimens
  dplyr::left_join(
    tibble::tibble(
      specimen = label_specimens_train
    ),
    df_mutation_train %>% 
      dplyr::select(
        specimen = .data$lab_id,
        mutation = .data$Hugo_Symbol
      ) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(val = 1) %>% 
      tidyr::spread(
        key = .data$mutation,
        value = .data$val
      ),
    by = "specimen"
  ) %>% 
  dplyr::mutate_at(
    .vars = vars(-.data$specimen),
    .funs = function(.x) dplyr::if_else(condition = is.na(.x), true = 0, false = .x)
  )
  
df_mutation_validation <- 
  # Include missing specimens
  dplyr::left_join(
    tibble::tibble(
      specimen = label_specimens_validation
    ),
    df_mutation_validation %>% 
      dplyr::select(
        specimen = .data$lab_id,
        mutation = .data$Hugo_Symbol
      ) %>% 
      dplyr::distinct() %>% 
      # Consider only mutations in the training data
      dplyr::filter(
        .data$mutation %in% (df_mutation_train %>% dplyr::select(-.data$specimen) %>% names)
      ) %>% 
      dplyr::mutate(val = 1) %>% 
      tidyr::spread(
        key = .data$mutation,
        value = .data$val
      ),
    by = "specimen"
  ) %>% 
  # Replace NA with 0 - negative mutation status
  dplyr::mutate_at(
    .vars = vars(-.data$specimen),
    .funs = function(.x) dplyr::if_else(condition = is.na(.x), true = 0, false = .x)
  )

df_mutation_leaderboard <- 
  # Include missing specimens
  dplyr::left_join(
    tibble::tibble(
      specimen = label_specimens_leaderboard
    ),
    df_mutation_leaderboard %>% 
      dplyr::select(
        specimen = .data$lab_id,
        mutation = .data$Hugo_Symbol
      ) %>% 
      dplyr::distinct() %>% 
      # Consider only mutations in the training data
      dplyr::filter(
        .data$mutation %in% (df_mutation_train %>% dplyr::select(-.data$specimen) %>% names)
      ) %>% 
      dplyr::mutate(val = 1) %>% 
      tidyr::spread(
        key = .data$mutation,
        value = .data$val
      ),
    by = "specimen"
  ) %>% 
  # Replace NA with 0 - negative mutation status
  dplyr::mutate_at(
    .vars = vars(-.data$specimen),
    .funs = function(.x) dplyr::if_else(condition = is.na(.x), true = 0, false = .x)
  )

# Put together mutation data for training, validation and leaderboard specimens
df_mutation_all <- 
  dplyr::bind_rows(
    df_mutation_train,
    df_mutation_validation,
    df_mutation_leaderboard
  ) %>% 
  # Replace NA with 0 - negative mutation status 
  # (since not all the mutations in the training data are present in the validation and leaderboard data)
  dplyr::mutate_at(
    .vars = vars(-.data$specimen),
    .funs = function(.x) dplyr::if_else(condition = is.na(.x), true = 0, false = .x)
  )

# Ensure the lab ids are distinct
print(intersect(df_mutation_train$specimen, df_mutation_validation$specimen))
print(intersect(df_mutation_train$specimen, df_mutation_leaderboard$specimen))
print(intersect(df_mutation_leaderboard$specimen, df_mutation_validation$specimen))

# Subset the all matrix after replacing NAs with 0
df_mutation_train <- subset(df_mutation_all, specimen %in% df_mutation_train$specimen)
write.table(df_mutation_train, file = paste0(data.dir, "/training/mutation_recoded.csv"), sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

df_mutation_validation <- subset(df_mutation_all, specimen %in% df_mutation_validation$specimen)
write.table(df_mutation_validation, file = paste0(data.dir, "/validation/mutation_recoded.csv"), sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

df_mutation_leaderboard <- subset(df_mutation_all, specimen %in% df_mutation_leaderboard$specimen)
write.table(df_mutation_leaderboard, file = paste0(data.dir, "/leaderboard/mutation_recoded.csv"), sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)


df_mutation_all %>% dplyr::select(-.data$specimen) %>% unlist() %>% as.vector() %>% unique()
# 0 1

# # Save gene-based variables included in the mutation data to file
# df_mutation_all %>% 
#   dplyr::select(-.data$specimen) %>% 
#   names() %>% 
#   saveRDS(
#     file = "gene_symbols_in_mutation_data.rds"
#   )
  


# Matrix format
mat_mutation_all <- 
  df_mutation_all %>% 
  tibble::column_to_rownames(var = "specimen") %>% 
  as.matrix()
dim(mat_mutation_all) 

