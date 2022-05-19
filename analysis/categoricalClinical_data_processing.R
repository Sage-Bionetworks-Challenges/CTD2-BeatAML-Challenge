library(tidyverse)

data.dir <- "../Data/"

# Read training clinical data ----------------------------------------------------------------------------------
df_clinical_train <- 
  readr::read_csv(
    paste0(data.dir, "training/clinical_categorical.csv"),
    col_names = TRUE
  ) %>% 
  dplyr::select(
    .data$lab_id,
    #.data$priorMDS,
    .data$priorMDSMoreThanTwoMths,
    #.data$priorMDSMPN,
    .data$priorMDSMPNMoreThanTwoMths,
    #.data$priorMPN,
    .data$priorMPNMoreThanTwoMths,
    .data$`FAB/Blast.Morphology`,
    .data$Karyotype,
    .data$finalFusion
  ) %>% 
  # Karyotype: based on Jeff's advise consolidate 46,XX[20], 46,XY, 46,XY[19], 46,XY[20] 
  # (coding: 1, 2, 3, 4)
  dplyr::mutate(
    Karyotype = dplyr::case_when(
      Karyotype == 1 ~ 6,
      Karyotype == 2 ~ 6,
      Karyotype == 3 ~ 6,
      Karyotype == 4 ~ 6,
      TRUE           ~ Karyotype
    )
  )

# Specimens in the training AUC dataset (labels)
label_specimens_train <- 
  readr::read_csv(
    paste0(data.dir, "training/aucs.csv"),
    col_names = TRUE
  ) %>% 
  dplyr::pull(.data$lab_id) %>% 
  unique()
# overlap
intersect(label_specimens_train, df_clinical_train$lab_id) %>% length() == length(label_specimens_train) # TRUE

which(is.na(df_clinical_train)) %>% length()
# no missing values

# Variable overview
df_clinical_train %>% 
  dplyr::select(-.data$lab_id) %>% 
  purrr::map(
    ~unique(.x)
  )


# Convert categorical variables to binary ones,
# e.g., each "Final fusion" will become a separate binary indicator variable
df_clinical_train_recoded <-
  # Keep specimen order
  tibble::tibble(lab_id = df_clinical_train$lab_id) %>% 
  dplyr::left_join(
  # 
    df_clinical_train %>% 
    tidyr::gather(
      key = "variable",
      value = "value",
      -lab_id
    ) %>% 
    dplyr::left_join(
      readr::read_csv(
        paste0(data.dir, "training/clinical_categorical_legend.csv"),
        col_names = TRUE
      ) %>% 
      dplyr::rename(
        "variable" = "column",
        "value" = "enum",
        "variable_coding" = "value"
      ),
      by = c("variable", "value")
    ) %>% 
    dplyr::group_by(.data$variable) %>% 
    # Variable type
    dplyr::mutate(
      var_type = dplyr::if_else(
        condition = length(unique(variable_coding)) <= 2,
        true = "binary",
        false = "cat"
      )
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      variable_new = dplyr::if_else(
        condition = var_type == "binary",
        true = variable,
        # recode, new binary variable
        false = paste0(variable, "_", variable_coding)
      ),
      value_new = dplyr::if_else(
        condition = var_type == "binary",
        true = value,
        # new binary variable
        false = 1
      )
    ) %>% 
    dplyr::select(
      .data$lab_id,
      .data$variable_new,
      .data$value_new
    ) %>% 
    tidyr::spread(
      key = variable_new,
      value = value_new
    ),
    by = "lab_id"
  ) %>% 
  # Replace NA with 0 
  dplyr::mutate_at(
    .vars = vars(-.data$lab_id),
    .funs = function(.x) dplyr::if_else(condition = is.na(.x), true = 0, false = .x)
  )

# Quick check for the values in the data frame
df_clinical_train_recoded %>% 
  dplyr::select(-.data$lab_id) %>% 
  unlist() %>% 
  as.vector() %>% 
  unique()
# 1 0
  
# Variables in the validation and leaderboard data will be created based on the training data
# (some variables might not be present; potentially, there could be some new variables,
# not present in the training data so they will not be considered)
training_vars <- 
  df_clinical_train_recoded %>% 
  dplyr::select(-.data$lab_id) %>% 
  names()
 
# Original variable types based on training data
training_org_var_types <- 
  df_clinical_train %>% 
  tidyr::gather(
    key = "variable",
    value = "value",
    -lab_id
  ) %>% 
  dplyr::left_join(
    readr::read_csv(
      paste0(data.dir, "training/clinical_categorical_legend.csv"),
      col_names = TRUE
    ) %>% 
      dplyr::rename(
        "variable" = "column",
        "value" = "enum",
        "variable_coding" = "value"
      ),
    by = c("variable", "value")
  ) %>% 
  dplyr::group_by(.data$variable) %>% 
  # Variable type
  dplyr::mutate(
    var_type = dplyr::if_else(
      condition = length(unique(variable_coding)) <= 2,
      true = "binary",
      false = "cat"
    )
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(
    .data$variable,
    .data$var_type
  ) %>% 
  dplyr::distinct()
  

# --------------------------------------------------------------------------------------------------------------
# Read validation clinical data ----------------------------------------------------------------------------------
df_clinical_validation <- 
  readr::read_csv(
    paste0(data.dir, "validation/clinical_categorical.csv"),
    col_names = TRUE
  ) %>% 
  dplyr::select(
    .data$lab_id,
    #.data$priorMDS,
    .data$priorMDSMoreThanTwoMths,
    #.data$priorMDSMPN,
    .data$priorMDSMPNMoreThanTwoMths,
    #.data$priorMPN,
    .data$priorMPNMoreThanTwoMths,
    .data$`FAB/Blast.Morphology`,
    .data$Karyotype,
    .data$finalFusion
  ) %>% 
  # Karyotype: consolidate 46,XX[20], 46,XY, 46,XY[19], 46,XY[20] 
  # (coding: 1, 2, 3, 4)
  dplyr::mutate(
    Karyotype = dplyr::case_when(
      Karyotype == 1 ~ 6,
      Karyotype == 2 ~ 6,
      Karyotype == 3 ~ 6,
      Karyotype == 4 ~ 6,
      TRUE           ~ Karyotype
    )
  )

# Specimens in the validation AUC dataset (labels)
label_specimens_validation <- 
  readr::read_csv(
    paste0(data.dir, "validation/aucs.csv"),
    col_names = TRUE
  ) %>% 
  dplyr::pull(.data$lab_id) %>% 
  unique()
# overlap
intersect(label_specimens_validation, df_clinical_validation$lab_id) %>% length() == length(label_specimens_validation) # TRUE

which(is.na(df_clinical_validation)) %>% length()
# no missing values

# Variable overview
df_clinical_validation %>% 
  dplyr::select(-.data$lab_id) %>% 
  purrr::map(
    ~unique(.x)
  )

# Convert categorical variables to binary ones based on training data
df_clinical_validation_recoded <-
  # Keep specimen order
  tibble::tibble(lab_id = df_clinical_validation$lab_id) %>% 
  dplyr::left_join(
    # 
    df_clinical_validation %>% 
      tidyr::gather(
        key = "variable",
        value = "value",
        -lab_id
      ) %>% 
      dplyr::left_join(
        readr::read_csv(
          paste0(data.dir, "validation/clinical_categorical_legend.csv"),
          col_names = TRUE
        ) %>% 
          dplyr::rename(
            "variable" = "column",
            "value" = "enum",
            "variable_coding" = "value"
          ),
        by = c("variable", "value")
      ) %>% 
      # Variable types based on training data
      dplyr::left_join(
        training_org_var_types,
        by = "variable"
      ) %>% 
      dplyr::mutate(
        variable_new = dplyr::if_else(
          condition = var_type == "binary",
          true = variable,
          # recode, new binary variable
          false = paste0(variable, "_", variable_coding)
        ),
        value_new = dplyr::if_else(
          condition = var_type == "binary",
          true = value,
          # new binary variable
          false = 1
        )
      ) %>% 
      # Consider only variables in the training data
      dplyr::filter(
        variable_new %in% training_vars
      ) %>% 
      dplyr::select(
        .data$lab_id,
        .data$variable_new,
        .data$value_new
      ) %>% 
      tidyr::spread(
        key = variable_new,
        value = value_new
      ),
    by = "lab_id"
  ) %>% 
  # Replace NA with 0 
  dplyr::mutate_at(
    .vars = vars(-.data$lab_id),
    .funs = function(.x) dplyr::if_else(condition = is.na(.x), true = 0, false = .x)
  )

# Quick check for the values in the data frame
df_clinical_validation_recoded %>% 
  dplyr::select(-.data$lab_id) %>% 
  unlist() %>% 
  as.vector() %>% 
  unique()
# 0 1


# --------------------------------------------------------------------------------------------------------------
# Read leaderboard clinical data ----------------------------------------------------------------------------------
df_clinical_leaderboard <- 
  readr::read_csv(
    paste0(data.dir, "leaderboard/clinical_categorical.csv"),
    col_names = TRUE
  ) %>% 
  dplyr::select(
    .data$lab_id,
    #.data$priorMDS,
    .data$priorMDSMoreThanTwoMths,
    #.data$priorMDSMPN,
    .data$priorMDSMPNMoreThanTwoMths,
    #.data$priorMPN,
    .data$priorMPNMoreThanTwoMths,
    .data$`FAB/Blast.Morphology`,
    .data$Karyotype,
    .data$finalFusion
  ) %>% 
  # Karyotype: based on Jeff's advise consolidate 46,XX[20], 46,XY, 46,XY[19], 46,XY[20] 
  # (coding: 1, 2, 3, 4)
  dplyr::mutate(
    Karyotype = dplyr::case_when(
      Karyotype == 1 ~ 6,
      Karyotype == 2 ~ 6,
      Karyotype == 3 ~ 6,
      Karyotype == 4 ~ 6,
      TRUE           ~ Karyotype
    )
  )

# Specimens in the leaderboard AUC dataset (labels)
label_specimens_leaderboard <- 
  readr::read_csv(
    paste0(data.dir, "leaderboard/aucs.csv"),
    col_names = TRUE
  ) %>% 
  dplyr::pull(.data$lab_id) %>% 
  unique()
# overlap
intersect(label_specimens_leaderboard, df_clinical_leaderboard$lab_id) %>% length() == length(label_specimens_leaderboard) # TRUE

which(is.na(df_clinical_leaderboard)) %>% length()
# no missing values

# Variable overview
df_clinical_leaderboard %>% 
  dplyr::select(-.data$lab_id) %>% 
  purrr::map(
    ~unique(.x)
  )

# Convert categorical variables to binary ones based on training data
df_clinical_leaderboard_recoded <-
  # Keep specimen order
  tibble::tibble(lab_id = df_clinical_leaderboard$lab_id) %>% 
  dplyr::left_join(
    # 
    df_clinical_leaderboard %>% 
      tidyr::gather(
        key = "variable",
        value = "value",
        -lab_id
      ) %>% 
      dplyr::left_join(
        readr::read_csv(
          paste0(data.dir, "leaderboard/clinical_categorical_legend.csv"),
          col_names = TRUE
        ) %>% 
          dplyr::rename(
            "variable" = "column",
            "value" = "enum",
            "variable_coding" = "value"
          ),
        by = c("variable", "value")
      ) %>% 
      # Variable types based on training data
      dplyr::left_join(
        training_org_var_types,
        by = "variable"
      ) %>% 
      dplyr::mutate(
        variable_new = dplyr::if_else(
          condition = var_type == "binary",
          true = variable,
          # recode, new binary variable
          false = paste0(variable, "_", variable_coding)
        ),
        value_new = dplyr::if_else(
          condition = var_type == "binary",
          true = value,
          # new binary variable
          false = 1
        )
      ) %>% 
      # Consider only variables in the training data
      dplyr::filter(
        variable_new %in% training_vars
      ) %>% 
      dplyr::select(
        .data$lab_id,
        .data$variable_new,
        .data$value_new
      ) %>% 
      tidyr::spread(
        key = variable_new,
        value = value_new
      ),
    by = "lab_id"
  ) %>% 
  # Replace NA with 0 
  dplyr::mutate_at(
    .vars = vars(-.data$lab_id),
    .funs = function(.x) dplyr::if_else(condition = is.na(.x), true = 0, false = .x)
  )

# Quick check for the values in the data frame
df_clinical_leaderboard_recoded %>% 
  dplyr::select(-.data$lab_id) %>% 
  unlist() %>% 
  as.vector() %>% 
  unique()
# 0 1



# ----------------------------------------------------------------------------------------------------------------
# Put together clinical data for training, validation and leaderboard specimens
df_clinical_all <- 
  dplyr::bind_rows(
    df_clinical_train_recoded,
    df_clinical_validation_recoded,
    df_clinical_leaderboard_recoded
  ) %>% 
  # Replace NA with 0
  # (since not all the binary training variables are present in the leaderboard and validation data
  dplyr::mutate_at(
    .vars = vars(-.data$lab_id),
    .funs = function(.x) dplyr::if_else(condition = is.na(.x), true = 0, false = .x)
  )

# Ensure the lab ids are distinct
intersect(df_clinical_train_recoded$lab_id, df_clinical_validation_recoded$lab_id)
intersect(df_clinical_train_recoded$lab_id, df_clinical_leaderboard_recoded$lab_id)
intersect(df_clinical_leaderboard_recoded$lab_id, df_clinical_validation_recoded$lab_id)

# Subset the all matrix after replacing NAs with 0
df_clinical_train_recoded <- subset(df_clinical_all, lab_id %in% df_clinical_train_recoded$lab_id)
write.table(df_clinical_train_recoded, file = paste0(data.dir, "/training/clinical_categorical_recoded.csv"), sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

df_clinical_validation_recoded <- subset(df_clinical_all, lab_id %in% df_clinical_validation_recoded$lab_id)
write.table(df_clinical_validation_recoded, file = paste0(data.dir, "/validation/clinical_categorical_recoded.csv"), sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

df_clinical_leaderboard_recoded <- subset(df_clinical_all, lab_id %in% df_clinical_leaderboard_recoded$lab_id)
write.table(df_clinical_leaderboard_recoded, file = paste0(data.dir, "/leaderboard/clinical_categorical_recoded.csv"), sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

df_clinical_all %>% dplyr::select(-.data$lab_id) %>% unlist() %>% as.vector() %>% unique()
# 1 0

# Number of specimens with "positive stataus" per variable
# colSums(df_clinical_all %>% dplyr::select(-.data$lab_id)) %>% boxplot()

# Matrix format
mat_clinical_all <- 
  df_clinical_all %>% 
  dplyr::rename("specimen" = "lab_id") %>% 
  tibble::column_to_rownames(var = "specimen") %>% 
  as.matrix()
dim(mat_clinical_all)
