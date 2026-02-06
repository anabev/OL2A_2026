# ------ Libraries ------
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)
library(nnet)
library(pscl)
library(caret)
library(car)

# ------ Data input ------
# Reading clinical and timeline data files from TCGA-BRCA cohort

clin_patient <- read.delim("C:\\Users\\ana_v\\OneDrive\\Área de Trabalho\\elixir\\brca_tcga_gdc\\data_clinical_patient.txt", stringsAsFactors = FALSE)
clin_sample  <- read.delim("C:\\Users\\ana_v\\OneDrive\\Área de Trabalho\\elixir\\brca_tcga_gdc\\data_clinical_sample.txt", stringsAsFactors = FALSE)
timeline_dx  <- read.delim("C:\\Users\\ana_v\\OneDrive\\Área de Trabalho\\elixir\\brca_tcga_gdc\\data_timeline_diagnosis.txt", stringsAsFactors = FALSE)
timeline_tx  <- read.delim("C:\\Users\\ana_v\\OneDrive\\Área de Trabalho\\elixir\\brca_tcga_gdc\\data_timeline_treatment.txt", stringsAsFactors = FALSE)


# ------ Sample selection ------
# Selection of samples used in this study
# Only samples with A2 code were retained, according to the study design

clin_patient <- clin_patient[23:124, ]
clin_sample  <- clin_sample[23:124, ]
timeline_dx  <- timeline_dx[27:159, ]
timeline_tx  <- timeline_tx[47:471, ]

# ------ Treatment count ------
# Counting the number of treatment types per patient

contagem_amostras <- timeline_tx %>%
  count(PATIENT_ID, name = "NUMBER_TREATMENT_TYPE")

# ------ Handling repeated measurements ------
# Some data frames contain multiple records per patient
# Only the first observation per PATIENT_ID is retained

clin_patient <- clin_patient %>%
  group_by(PATIENT_ID) %>%
  slice(1) %>%
  ungroup()

clin_sample <- clin_sample %>%
  group_by(PATIENT_ID) %>%
  slice(1) %>%
  ungroup()

timeline_dx <- timeline_dx %>%
  group_by(PATIENT_ID) %>%
  slice(1) %>%
  ungroup()

timeline_tx <- timeline_tx %>%
  group_by(PATIENT_ID) %>%
  slice(1) %>%
  ungroup()

# ------ Data integration ------
# Merging all datasets into a single analytical data frame

df_final <- clin_patient %>%
  left_join(clin_sample, by = "PATIENT_ID") %>%
  left_join(timeline_dx, by = "PATIENT_ID") %>%
  left_join(contagem_amostras, by = 'PATIENT_ID') %>%
  left_join(timeline_tx, by = "PATIENT_ID")

# ------ Molecular subtype information ------
# Adding molecular subtype annotations (TCGA classification)

subtypes <- read.csv(
  "C:/Users/ana_v/OneDrive/Área de Trabalho/elixir/subtipos TCGA-BRCA-A2.csv",
  stringsAsFactors = FALSE
)

subtypes <- subtypes %>%
  group_by(PATIENT_ID) %>%
  slice(1) %>%
  ungroup()

df_final <- df_final %>%
  left_join(
    subtypes %>% select(PATIENT_ID, SUBTYPE),
    by = "PATIENT_ID"
  )

# ------ Variable type inspection ------
# Identifying numeric and categorical variables

num_cols <- df_final %>%
  select(where(is.numeric)) %>%
  colnames()

cat_cols <- df_final %>%
  select(where(~ is.character(.) | is.factor(.))) %>%
  colnames()

num_cols
cat_cols

# ------ Manual type conversion ------
# All variables were initially read as categorical
# Manual conversion of selected variables to numeric type

cols_numericas <- c(
  "AGE",
  "DAYS_LAST_FOLLOWUP",
  "DAYS_TO_BIRTH",
  "DAYS_TO_DEATH",
  "DFS_MONTHS",
  "TMB_NONSYNONYMOUS",
  "START_DATE.x",
  "AGE_AT_DIAGNOSIS",
  "START_DATE.y",
  "NUMBER_TREATMENT_TYPE",
  "STOP_DATE.y",
  "NUMBER_OF_CYCLES"
)

df_final <- df_final %>%
  mutate(
    across(all_of(cols_numericas), as.numeric)
  )

# ------ Missing data imputation ------
# Simple imputation strategy:
# - Mean imputation for numerical variables
# - Mode imputation for categorical variables

moda <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

df_final <- df_final %>%
  mutate(
    across(
      where(~ is.character(.) | is.factor(.)),
      ~ na_if(trimws(as.character(.)), "")
    )
  )

df_final <- df_final %>%
  mutate(
    across(
      all_of(num_cols),
      ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
    )
  )

df_final <- df_final %>%
  mutate(
    across(
      all_of(cat_cols),
      ~ {
        x <- .
        x[is.na(x)] <- moda(x)
        x
      }
    )
  )

df_final <- df_final %>%
  mutate(
    across(
      all_of(cat_cols),
      as.factor
    )
  )

# ------ Variable exclusion ------
# Removing variables not relevant to the analytical objectives
# or redundant due to prior aggregation steps

cols_remover <- c(
  "BIOPSY_SITE",
  "ETHNICITY",
  "ICD_10",
  "OTHER_PATIENT_ID",
  "PRIMARY_SITE_PATIENT",
  "PRIOR_MALIGNANCY",
  "PRIOR_TREATMENT",
  "SEX",
  "PROJECT_ID",
  "PROJECT_NAME",
  "PROJECT_STATE",
  "SAMPLE_ID",
  "IS_FFPE",
  "OTHER_SAMPLE_ID",
  "SAMPLE_TYPE",
  "SAMPLE_TYPE_ID",
  "ONCOTREE_CODE",
  "CANCER_TYPE",
  "CANCER_TYPE_DETAILED",
  "START_DATE.x",
  "STOP_DATE.x",
  "EVENT_TYPE.x",
  "GRADE",
  "EVENT_TYPE.y",
  "NUMBER_OF_CYCLES",
  "THERAPEUTIC_AGENT",
  "TREATMENT_ANATOMIC_SITE",
  "TREATMENT_EFFECT",
  "TREATMENT_END_REASON",
  "TREATMENT_TYPE",
  "TREATMENT_OUTCOME"
)

df_final <- df_final %>%
  select(-all_of(cols_remover))

# ------ Treatment variable adjustment ------
# Since most patients received more than one treatment type,
# TREATMENT_TYPE was excluded and replaced by NUMBER_TREATMENT_TYPE
# The variable was coerced to integer format

df_final <- df_final %>%
  mutate(
    NUMBER_TREATMENT_TYPE = as.integer(NUMBER_TREATMENT_TYPE)
  )

# ------ Descriptive statistics and contingency tables ------
# Identification of numeric and categorical variables
num_cols <- df_final %>%
  select(where(is.numeric)) %>%
  colnames()

cat_cols <- df_final %>%
  select(where(is.factor) | where(is.character)) %>%
  colnames()

# Summary statistics for numerical variables
summary_numericas <- lapply(
  df_final[num_cols],
  function(x) {
    c(
      N      = sum(!is.na(x)),
      Mean   = mean(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE),
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE)
    )
  }
)

summary_numericas <- as.data.frame(do.call(rbind, summary_numericas))
summary_numericas

# Contingency tables of categorical variables by molecular subtype
tables_por_subtype <- lapply(
  setdiff(cat_cols, "SUBTYPE"),
  function(var) {
    table(df_final[[var]], df_final$SUBTYPE, useNA = "ifany")
  }
)

names(tables_por_subtype) <- setdiff(cat_cols, "SUBTYPE")
tables_por_subtype

# Summary of numerical variables stratified by molecular subtype
summary_num_por_subtype <- lapply(
  num_cols,
  function(var) {
    df_final %>%
      group_by(SUBTYPE) %>%
      summarise(
        N      = sum(!is.na(.data[[var]])),
        Mean   = mean(.data[[var]], na.rm = TRUE),
        SD     = sd(.data[[var]], na.rm = TRUE),
        Median = median(.data[[var]], na.rm = TRUE),
        Q1     = quantile(.data[[var]], 0.25, na.rm = TRUE),
        Q3     = quantile(.data[[var]], 0.75, na.rm = TRUE),
        .groups = "drop"
      )
  }
)

names(summary_num_por_subtype) <- num_cols
summary_num_por_subtype


# ------ Boxplot visualization of selected numerical variables ------
# Reshaping data to long format for visualization

df_long <- df_final %>%
  select(
    SUBTYPE,
    AGE,
    NUMBER_TREATMENT_TYPE,
    TMB_NONSYNONYMOUS
  ) %>%
  pivot_longer(
    cols = -SUBTYPE,
    names_to = "Variable",
    values_to = "Value"
  )

# Recoding variable names and defining subtype order
df_long <- df_long %>%
  mutate(
    Variable = recode(
      Variable,
      AGE = "Age",
      NUMBER_TREATMENT_TYPE = "Number of Treatment Types",
      TMB_NONSYNONYMOUS = "Tumor Mutational Burden"
    ),
    SUBTYPE = factor(
      SUBTYPE,
      levels = c("Basal", "Her2", "LumA", "LumB", "Normal")
    )
  )

# Boxplot of numerical variables stratified by molecular subtype
ggplot(df_long, aes(x = SUBTYPE, y = Value, fill = SUBTYPE)) +
  geom_boxplot(
    alpha = 0.85,
    outlier.size = 1.5,
    outlier.alpha = 0.8
  ) +
  facet_wrap(~ Variable, scales = "free_y", nrow = 1) +
  labs(
    title = "Distribution of Numerical Variables by Molecular Subtype",
    x = "Molecular Subtype",
    y = "Value"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )


# ------ Pathological T stage recoding ------
# Grouping pathological T stage categories

df_final <- df_final %>%
  mutate(
    PATH_T_STAGE_AGRUPADO = case_when(
      PATH_T_STAGE %in% c("T1", "T1b", "T1c") ~ "T1",
      PATH_T_STAGE == "T2" ~ "T2",
      PATH_T_STAGE == "T3" ~ "T3",
      PATH_T_STAGE == "T4b" ~ "T4",
      TRUE ~ as.character(PATH_T_STAGE)
    )
  )

table_path_stage <- table(df_final$PATH_T_STAGE_AGRUPADO, df_final$SUBTYPE)
table_path_stage

# ------ Pathological N stage recoding ------
# Grouping pathological N stage categories
df_final <- df_final %>%
  mutate(
    PATH_N_STAGE_AGRUPADO = case_when(
      PATH_N_STAGE %in% c("N0", "N0 (i-)", "N0 (i+)") ~ "N0",
      PATH_N_STAGE %in% c("N1", "N1a", "N1mi")      ~ "N1",
      PATH_N_STAGE %in% c("N2", "N2a")              ~ "N2",
      PATH_N_STAGE %in% c("N3", "N3a", "N3c")       ~ "N3",
      TRUE ~ as.character(PATH_N_STAGE)
    )
  )

table_path_n <- table(df_final$PATH_N_STAGE_AGRUPADO, df_final$SUBTYPE)
table_path_n

#-----------------------------------------------------------------------
# Setting Luminal A as the reference category for molecular subtype
df_final$SUBTYPE <- relevel(df_final$SUBTYPE, ref = "LumA")

# Multinomial logistic regression for pathological T stage
md1 <- multinom(PATH_T_STAGE_AGRUPADO ~ SUBTYPE, data = df_final)
summary(md1)

# Extraction of coefficients and standard errors
coefs <- summary(md1)$coefficients
SE <- summary(md1)$standard.errors

# Odds ratios and 95% confidence intervals
OR <- exp(coefs)
lower <- exp(coefs - 1.96 * SE)
upper <- exp(coefs + 1.96 * SE)

# Results table
resultado1 <- data.frame(
  Outcome = rep(rownames(coefs), each = ncol(coefs)),
  Comparacao = rep(colnames(coefs), times = nrow(coefs)),
  OR = as.vector(OR),
  Lower95 = as.vector(lower),
  Upper95 = as.vector(upper)
)

resultado1

# Tidy model output with confidence intervals
resultado <- tidy(md1, exponentiate = TRUE, conf.int = TRUE)
resultado

#------------------------------------------------

# Resetting reference category for molecular subtype
df_final$SUBTYPE <- relevel(df_final$SUBTYPE, ref = "LumA")

#-- Pathological T stage model
md_t <- multinom(PATH_T_STAGE_AGRUPADO ~ SUBTYPE, data = df_final)
summary(md_t)

coefs_t <- summary(md_t)$coefficients
SE_t <- summary(md_t)$standard.errors

OR_t <- exp(coefs)
lower_t <- exp(coefs_t - 1.96 * SE_t)
upper_t <- exp(coefs_t + 1.96 * SE_t)

resultado_t <- data.frame(
  Outcome_t = rep(rownames(coefs_t), each = ncol(coefs_t)),
  Comparacao_t = rep(colnames(coefs_t), times = nrow(coefs_t)),
  OR_t = as.vector(OR_t),
  Lower95_t = as.vector(lower_t),
  Upper95_t = as.vector(upper_t)
)

resultado_t

#-- Pathological N stage model
md_n <- multinom(PATH_N_STAGE_AGRUPADO ~ SUBTYPE, data = df_final)
summary(md_n)

coefs_n <- summary(md_n)$coefficients
SE_n <- summary(md_n)$standard.errors

OR_n <- exp(coefs_n)
lower_n <- exp(coefs_n - 1.96 * SE_n)
upper_n <- exp(coefs_n + 1.96 * SE_n)

resultado_n <- data.frame(
  Outcome_n = rep(rownames(coefs_n), each = ncol(coefs_n)),
  Comparacao_n = rep(colnames(coefs_n), times = nrow(coefs_n)),
  OR_n = as.vector(OR_n),
  Lower95_n = as.vector(lower_n),
  Upper95_n = as.vector(upper_n)
)

resultado_n

#-- Pathological M stage model
md_m <- multinom(PATH_M_STAGE ~ SUBTYPE, data = df_final)
summary(md_m)

coefs_m <- summary(md_m)$coefficients
SE_m    <- summary(md_m)$standard.errors

OR_m <- exp(coefs_m)
lower_m <- exp(coefs_m - 1.96 * SE_m)
upper_m <- exp(coefs_m + 1.96 * SE_m)

resultado_m <- data.frame(
  Comparacao_m = names(coefs_m),
  OR_m = OR_m,
  Lower95_m = lower_m,
  Upper95_m = upper_m
)

resultado_m

#-- Vital status model
md_vs <- multinom(VITAL_STATUS ~ SUBTYPE, data = df_final)
summary(md_vs)

coefs_vs <- summary(md_vs)$coefficients
SE_vs    <- summary(md_vs)$standard.errors

OR_vs <- exp(coefs_vs)
lower_vs <- exp(coefs_vs - 1.96 * SE_vs)
upper_vs <- exp(coefs_vs + 1.96 * SE_vs)

resultado_vs <- data.frame(
  Comparacao_vs = names(coefs_vs),
  OR_vs = OR_vs,
  Lower95_vs = lower_vs,
  Upper95_vs = upper_vs
)

resultado_vs

#-- Disease-free survival (DFS) status model
md_dfs <- multinom(DFS_STATUS ~ SUBTYPE, data = df_final)
summary(md_dfs)

coefs_dfs <- summary(md_dfs)$coefficients
SE_dfs    <- summary(md_dfs)$standard.errors

OR_dfs <- exp(coefs_dfs)
lower_dfs <- exp(coefs_dfs - 1.96 * SE_dfs)
upper_dfs <- exp(coefs_dfs + 1.96 * SE_dfs)

resultado_dfs <- data.frame(
  Comparacao_dfs = names(coefs_dfs),
  OR_dfs = OR_dfs,
  Lower95_dfs = lower_dfs,
  Upper95_dfs = upper_dfs
)

resultado_dfs

#-- Function to extract odds ratios from multinomial models
get_or_multinom <- function(model){
  
  coefs <- summary(model)$coefficients
  SE    <- summary(model)$standard.errors
  
  OR    <- exp(coefs)
  lower <- exp(coefs - 1.96 * SE)
  upper <- exp(coefs + 1.96 * SE)
  
  data.frame(
    Outcome   = rep(rownames(coefs), each = ncol(coefs)),
    Predictor = rep(colnames(coefs), times = nrow(coefs)),
    OR        = as.vector(OR),
    Lower95   = as.vector(lower),
    Upper95   = as.vector(upper)
  )
}

#-- Age model
modelo_age <- multinom(SUBTYPE ~ AGE, data = df_final)
resultado_age <- get_or_multinom(modelo_age)
resultado_age

#-- Tumor mutational burden model
modelo_TMB_NONSYNONYMOUS <- multinom(SUBTYPE ~ TMB_NONSYNONYMOUS, data = df_final)
resultado_TMB_NONSYNONYMOUS <- get_or_multinom(modelo_TMB_NONSYNONYMOUS)
resultado_TMB_NONSYNONYMOUS

#-- Number of treatment types model
modelo_NUMBER_TREATMENT_TYPE <- multinom(SUBTYPE ~ NUMBER_TREATMENT_TYPE, data = df_final)
resultado_NUMBER_TREATMENT_TYPE <- get_or_multinom(modelo_NUMBER_TREATMENT_TYPE)
resultado_NUMBER_TREATMENT_TYPE

#-- Adjusted multinomial logistic regression model
modelo_ajustado <- multinom(
  SUBTYPE ~ 
    PATH_T_STAGE_AGRUPADO +
    PATH_N_STAGE_AGRUPADO +
    PATH_M_STAGE +
    VITAL_STATUS +
    DFS_STATUS +
    AGE +
    TMB_NONSYNONYMOUS +
    NUMBER_TREATMENT_TYPE,
  data = df_final
)

# Adjusted odds ratios and confidence intervals
coefs_ajustado <- summary(modelo_ajustado)$coefficients
SE_ajustado    <- summary(modelo_ajustado)$standard.errors
OR_ajustado    <- exp(coefs_ajustado)
lower_ajustado <- exp(coefs_ajustado - 1.96 * SE_ajustado)
upper_ajustado <- exp(coefs_ajustado + 1.96 * SE_ajustado)

resultado_ajustado <- data.frame(
  Outcome_ajustado = rep(rownames(coefs_ajustado), each = ncol(coefs_ajustado)),
  Predictor_ajustado = rep(colnames(coefs_ajustado), times = nrow(coefs_ajustado)),
  OR_ajustado = as.vector(OR_ajustado),
  Lower95_ajustado = as.vector(lower_ajustado),
  Upper95_ajustado = as.vector(upper_ajustado)
)

resultado_ajustado

#-- Pearson and deviance residuals
res_pearson <- residuals(modelo_ajustado, type = "pearson")
summary(as.vector(res_pearson))

res_deviance <- residuals(modelo_ajustado, type = "deviance")
summary(as.vector(res_deviance))

par(mfrow = c(1,2))

hist(
  as.vector(res_pearson),
  main = "Pearson Residuals",
  xlab = "Value",
  col = "gray"
)

hist(
  as.vector(res_deviance),
  main = "Deviance Residuals",
  xlab = "Value",
  col = "gray"
)

par(mfrow = c(1,1))

#-- Model validation
#-- Likelihood ratio test
modelo_nulo <- multinom(SUBTYPE ~ 1, data = df_final)
anova(modelo_nulo, modelo_ajustado, test = "Chisq")

#-- Pseudo-R²
pR2(modelo_ajustado)

#-- Variance inflation factor (VIF)
vif_modelo <- vif(lm(
  as.numeric(SUBTYPE) ~ 
    PATH_T_STAGE_AGRUPADO +
    PATH_N_STAGE_AGRUPADO +
    PATH_M_STAGE +
    VITAL_STATUS +
    DFS_STATUS +
    AGE +
    TMB_NONSYNONYMOUS +
    NUMBER_TREATMENT_TYPE,
  data = df_final
))
vif_modelo

#-- Confusion matrix and accuracy
pred_classe <- predict(modelo_ajustado, type = "class")
tab_confusao <- table(
  Observed = df_final$SUBTYPE,
  Predicted = pred_classe
)
tab_confusao

acuracia <- sum(diag(tab_confusao)) / sum(tab_confusao)
acuracia

#-- Cross-validation
set.seed(123)

ctrl <- trainControl(
  method = "cv",
  number = 10
)

cv_model <- train(
  SUBTYPE ~ 
    PATH_T_STAGE_AGRUPADO +
    PATH_N_STAGE_AGRUPADO +
    PATH_M_STAGE +
    VITAL_STATUS +
    DFS_STATUS +
    AGE +
    TMB_NONSYNONYMOUS +
    NUMBER_TREATMENT_TYPE,
  data = df_final,
  method = "multinom",
  trControl = ctrl,
  trace = FALSE
)

cv_model

#-- Simulated envelope plot of PIT residuals
prob_pred <- predict(modelo_ajustado, type = "probs")

set.seed(123)

res_pit <- sapply(1:nrow(df_final), function(i) {
  classe_obs <- df_final$SUBTYPE[i]
  p <- prob_pred[i, classe_obs]
  runif(1, min = 0, max = p)
})

n_sim <- 1000
res_simulados <- matrix(NA, nrow = nrow(df_final), ncol = n_sim)

for (s in 1:n_sim) {
  y_sim <- apply(prob_pred, 1, function(p)
    sample(colnames(prob_pred), size = 1, prob = p)
  )
  
  for (i in 1:nrow(df_final)) {
    p <- prob_pred[i, y_sim[i]]
    res_simulados[i, s] <- runif(1, 0, p)
  }
}

lim_inf <- apply(res_simulados, 1, quantile, 0.025)
lim_sup <- apply(res_simulados, 1, quantile, 0.975)

plot(
  sort(res_pit),
  type = "p",
  pch = 16,
  cex = 0.6,
  ylab = "PIT Residuals",
  xlab = "Ordered observations"
)

lines(sort(lim_inf), col = "red", lty = 2)
lines(sort(lim_sup), col = "red", lty = 2)

#-- Hosmer-Lemeshow goodness-of-fit test
library(ResourceSelection)

prob_pred <- predict(modelo_ajustado, type = "probs")

y_luma <- ifelse(df_final$SUBTYPE == "LumA", 1, 0)
hl_luma <- hoslem.test(y_luma, prob_pred[, "LumA"], g = 10)
hl_luma

subtipos <- colnames(prob_pred)

hl_resultados <- lapply(subtipos, function(s) {
  y_bin <- ifelse(df_final$SUBTYPE == s, 1, 0)
  hoslem.test(y_bin, prob_pred[, s], g = 10)
})

names(hl_resultados) <- subtipos
hl_resultados