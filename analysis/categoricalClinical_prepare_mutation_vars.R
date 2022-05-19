library(tidyverse)

# Prepare categorical clinical variables to be moved to mutation data

data.dir <- "../Data/"


# Read data -------------------------------------------------------------------------------------------------
df_clinical_all <- 
  dplyr::bind_rows(
    # training
    readr::read_csv(
      paste0(data.dir, "training/clinical_categorical.csv"),
      col_names = TRUE
    ) %>% 
    dplyr::mutate(
      dataset = "training"
    ),
    # leaderboard
    readr::read_csv(
      paste0(data.dir, "leaderboard/clinical_categorical.csv"),
      col_names = TRUE
    ) %>% 
    dplyr::mutate(
      dataset = "leaderboard"
    ),
    # validation
    readr::read_csv(
      paste0(data.dir, "validation/clinical_categorical.csv"),
      col_names = TRUE
    ) %>% 
    dplyr::mutate(
      dataset = "validation"
    )
  )


# Convert categorical variables to binary ones -------------------------------------------------------------
df_clinical_all_recoded <-
  # Keep specimen order
  df_clinical_all %>% 
    dplyr::select(
      .data$lab_id,
      .data$dataset
    ) %>% 
  dplyr::left_join(
    # 
    df_clinical_all %>% 
      dplyr::select(-.data$dataset) %>% 
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
    .vars = vars(-.data$lab_id, -.data$dataset),
    .funs = function(.x) dplyr::if_else(condition = is.na(.x), true = 0, false = .x)
  )

df_clinical_all_recoded %>% dplyr::group_by(dataset) %>% dplyr::summarise(n())

# Quick check for the values in the data frame
df_clinical_all_recoded %>% 
  dplyr::select(-.data$lab_id, -.data$dataset) %>% 
  unlist() %>% 
  as.vector() %>% 
  unique()
# 1 0


# Prepare mutation variables ---------------------------------------------------------------------------

df_mutation <- 
  df_clinical_all_recoded %>% 
  dplyr::select(
    .data$lab_id,
    .data$dataset,
    #
    CBFB_MYH11_1 = .data$`specificDxAtAcquisition_AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22); CBFB-MYH11`,
    CBFB_MYH11_2 = .data$`specificDxAtAcquisition_AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22); CBFB-MYH11`,
    CBFB_MYH11_3 = .data$`finalFusion_CBFB-MYH11`,
    # 
    RPN1_EVI1_1 = .data$`specificDxAtInclusion_AML with inv(3)(q21q26.2) or t(3;3)(q21;q26.2); RPN1-EVI1`,
    RPN1_EVI1_2 = .data$`specificDxAtAcquisition_AML with inv(3)(q21q26.2) or t(3;3)(q21;q26.2); RPN1-EVI1`,
    #
    CEBPA_1 = .data$`specificDxAtInclusion_AML with mutated CEBPA`,
    CEBPA_2 = .data$`specificDxAtAcquisition_AML with mutated CEBPA`,
    #
    RUNX1_RUNX1T1_1 = .data$`specificDxAtInclusion_AML with t(8;21)(q22;q22); RUNX1-RUNX1T1`,
    RUNX1_RUNX1T1_2 = .data$`specificDxAtAcquisition_AML with t(8;21)(q22;q22); RUNX1-RUNX1T1`,
    RUNX1_RUNX1T1_3 = .data$`finalFusion_RUNX1-RUNX1T1`,
    #
    MLLT3_MLL_1 = .data$`specificDxAtInclusion_AML with t(9;11)(p22;q23); MLLT3-MLL`,
    MLLT3_MLL_2 = .data$`specificDxAtAcquisition_AML with t(9;11)(p22;q23); MLLT3-MLL`,
    MLLT3_MLL_3 = .data$`finalFusion_MLLT3-KMT2A`,
    #
    PML_RARA_1 = .data$`specificDxAtInclusion_Acute promyelocytic leukaemia with t(15;17)(q22;q12); PML-RARA`,
    PML_RARA_2 = .data$`specificDxAtAcquisition_Acute promyelocytic leukaemia with t(15;17)(q22;q12); PML-RARA`,
    PML_RARA_3 = .data$`FAB/Blast.Morphology_M3`,
    PML_RARA_4 = .data$`finalFusion_PML-RARA`,
    #
    FLT3_ITD = .data$`FLT3-ITD`,
    #
    NPM1_1 = .data$NPM1,
    NPM1_2 = .data$`specificDxAtInclusion_AML with mutated NPM1`,
    NPM1_3 = .data$`specificDxAtAcquisition_AML with mutated NPM1`
  ) %>% 
  #
  # Consolidate as advised by Jeff
  dplyr::mutate(
    CBFB_MYH11 = CBFB_MYH11_1 + CBFB_MYH11_2 + CBFB_MYH11_3,
    RPN1_EVI1 = RPN1_EVI1_1 + RPN1_EVI1_2,
    CEBPA = CEBPA_1 + CEBPA_2,
    RUNX1_RUNX1T1 = RUNX1_RUNX1T1_1 + RUNX1_RUNX1T1_2 + RUNX1_RUNX1T1_3,
    MLLT3_MLL = MLLT3_MLL_1 + MLLT3_MLL_2 + MLLT3_MLL_3,
    PML_RARA = PML_RARA_1 + PML_RARA_2 + PML_RARA_3 + PML_RARA_4,
    NPM1 = NPM1_1 + NPM1_2 + NPM1_3
  ) %>% 
  dplyr::select(
    .data$lab_id,
    .data$dataset,
    .data$CBFB_MYH11,
    .data$RPN1_EVI1,
    .data$CEBPA,
    .data$RUNX1_RUNX1T1,
    .data$MLLT3_MLL,
    .data$PML_RARA,
    .data$NPM1,
    .data$FLT3_ITD
  )

purrr::map(
  setdiff(names(df_mutation), c("lab_id", "dataset")),
  ~ df_mutation %>% 
    dplyr::group_by(!!rlang::sym(.x)) %>% 
    dplyr::summarise(n())
)


df_mutation <- 
  df_mutation %>% 
  dplyr::mutate(
    CBFB_MYH11 = dplyr::if_else(CBFB_MYH11 == 0, 0, 1),
    RPN1_EVI1 = dplyr::if_else(RPN1_EVI1 == 0, 0, 1),
    CEBPA = dplyr::if_else(CEBPA == 0, 0, 1),
    RUNX1_RUNX1T1 = dplyr::if_else(RUNX1_RUNX1T1 == 0, 0, 1),
    MLLT3_MLL = dplyr::if_else(MLLT3_MLL == 0, 0, 1),
    PML_RARA = dplyr::if_else(PML_RARA == 0, 0, 1),
    NPM1 = dplyr::if_else(NPM1 == 0, 0, 1),
    FLT3_ITD = dplyr::if_else(FLT3_ITD == 0, 0, 1)
  )
  
purrr::map(
  setdiff(names(df_mutation), c("lab_id", "dataset")),
  ~ df_mutation %>% 
    dplyr::group_by(!!rlang::sym(.x)) %>% 
    dplyr::summarise(n())
)

# Save to file
saveRDS(
  object = df_mutation,
  file = "mutation_vars_from_categorical_clinical_data.rds"
)
