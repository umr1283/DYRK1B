
#### Envir =========================================================================================
output_directory <-"outputs/DYRK1B_RADIO_mist"
dir.create(path = output_directory, recursive = TRUE, mode = "0775")

library(readxl)
library(writexl)
source("mist.R")
library(rlang)
unnest_dt <- function(tbl, col) {
  tbl <- as.data.table(tbl)
  col <- ensyms(col)
  clnms <- syms(setdiff(colnames(tbl), as.character(col)))
  tbl <- as.data.table(tbl)
  tbl <- eval(
    expr(tbl[, as.character(unlist(!!!col)), by = list(!!!clnms)])
  )
  colnames(tbl) <- c(as.character(clnms), as.character(col))
  tbl
}
library(data.table)

#### Data ==========================================================================================

## to get IID and add PCS in pheno RADIO (for enfants)
PhenotypeDir <- "data/Phenotypes/"
Phenotype_file <- paste0(PhenotypeDir, "Phenotyes_latest.xlsx")
Phenotypes <- as.data.frame(readxl::read_excel(
  path = Phenotype_file, 
  col_names = TRUE,
  na = "NA",
  guess_max = 9365
))

# genetic data adult analyses
data_adultes <- readxl::read_excel("data/DYRK1B/DYRK1B-adultes.xlsx", guess_max = 7000)
str(data_adultes)

data_adultes$LOG_FI <- log(data_adultes$FI_D0)
data_adultes$LOG_triglycerides <- log(data_adultes$triglycerides)
data_adultes <- merge(data_adultes, Phenotypes[, c("ADN", "IID")], by = "ADN")

# genetic data children analyses
data_enfants <- readxl::read_excel("data/DYRK1B/DYRK1B-enfants.xlsx", guess_max = 5000)
str(data_enfants)
data_enfants <- merge(data_enfants, Phenotypes[, c("ADN", "IID")], by = "ADN")

covar_selected <- c("SEX", "AGE",  paste0("PC", 1:5))

samples_IID <- unique(c(data_enfants$IID, data_adultes$IID))
samples_ADN <- unique(c(data_enfants$ADN, data_adultes$ADN))

#### genotypes & annotation ========================================================================

annotation <- rbindlist(l = list(
  setDT(readxl::read_xlsx(path = "data/DYRK1B/DYRK1B-adultes.xlsx", guess_max = 7000))[
    , .SD, .SDcols = c("ADN",
      "rare-neutral-PLPall", "rare-neutral-PLPAB",
      "rare-neutral-PLPA",
      "rare-all", "P/LP-all", 
      "P/LP-AB",
      "P/LP-C", "P/LP-B",
      "P/LP-A"
    )
  ][, sheet := "adults"], 
  setDT(readxl::read_xlsx(path = "data/DYRK1B/DYRK1B-enfants.xlsx", guess_max = 5000))[
    , .SD, .SDcols = c("ADN",
      "rare-neutral-PLPall", "rare-neutral-PLPB",
      "rare-all", "P/LP-all",
      "P/LP-C", "P/LP-B"
    )
  ][, sheet := "children"]
), use.names = TRUE, fill = TRUE)
head(annotation)
str(annotation)
# tail(annotation)

annotation <- data.table::melt(
  annotation,
  id.vars = c("ADN", "sheet"),
  measure.vars = unique(c(
    "rare-neutral-PLPall", "rare-neutral-PLPAB",
    "rare-neutral-PLPA",
    "rare-all", "P/LP-all", 
    "P/LP-AB",
    "P/LP-C", "P/LP-B",
    "P/LP-A", 
    "rare-neutral-PLPall", "rare-neutral-PLPB",
    "rare-all", "P/LP-all",
    "P/LP-C", "P/LP-B"
  )), 
  variable.name = "selection_name", value.name = "ACMG"
)[
  , selection_name := as.character(selection_name)
][!is.na(ACMG), ]
# head(annotation)
# tail(annotation)
# str(annotation)
# table(annotation$ACMG, useNA = "ifany")

#### Build the matrix of genotypes ####
geno_mat <- dcast(
  data = unique(annotation[, .(ACMG, ADN)])[, geno := 1], ## toutes les mutations sont hétérozygotes 
  formula = ACMG ~ ADN, 
  value.var = c("geno")
)
#dim(geno_mat)
#geno_mat[1:5, 1:5]
geno_mat <- transpose(l = geno_mat, make.names = "ACMG", keep.names = "ADN")
#dim(geno_mat)
#geno_mat[1:5, 1:5]
geno_mat <- rbindlist(l = list(
  geno_mat,
  data.table(ADN = setdiff(x = samples_ADN, geno_mat$ADN)) # complete the genotype matrix with all WT
), use.names = TRUE, fill = TRUE)

geno_mat[is.na(geno_mat)] <- 0 # complete all non carriers are WT (DP checked before)

data_adultes <- merge(data_adultes, geno_mat, by = "ADN")
data_enfants <- merge(data_enfants, geno_mat, by = "ADN")

message("Save Genotypes")
fwrite(
  x = geno_mat, file = file.path(output_directory, "DYRK1B_geno_hg19.tsv.gz"), 
  row.names = FALSE, col.names = TRUE
)

#### Analyses adults  ==============================================================================
message("Adultes ...")
traits_adults <- data.table(
  data_from = "DYRK1B-adulte", 
  traits = c(
    "BMI", 
    "CC_T2D56", #with BMI
    "CC_T2D56",
    "CC_T2D61", #with BMI
    "CC_T2D61", 
    "CC_OBESITE", 
    "CC_OWT",
    "FG_D0", #with BMI
    "FG_D0", 
    "LOG_FI", #with BMI
    "LOG_FI", 
    "HDL", #with BMI
    "HDL", 
    "LOG_triglycerides",  #with BMI
    "LOG_triglycerides",
    "A", "B", "C", "D", "E", "MS",  #with BMI
    "A", "B", "C", "D", "E", "MS"
  ), 
  is_binary = c(
    FALSE, #"BMI", 
    TRUE, #"CC_T2D56", #with BMI
    TRUE, #"CC_T2D56",
    TRUE, #"CC_T2D61",#with BMI
    TRUE, #"CC_T2D61", 
    TRUE, #"CC_OBESITE", 
    TRUE, #"CC_OWT",
    FALSE, #"FG_D0", #with BMI
    FALSE, #"FG_D0", 
    FALSE, #"LOG_FI", #with BMI
    FALSE, #"LOG_FI", 
    FALSE, #"HDL",  #with BMI
    FALSE, #"HDL",
    FALSE, #"LOG_triglycerides",  #with BMI
    FALSE, #"LOG_triglycerides",
    TRUE, TRUE,TRUE,TRUE,TRUE, TRUE,#"A", "B", "C", "D", "E", "MS" #with BMI
    TRUE,TRUE,TRUE,TRUE,TRUE, TRUE#"A", "B", "C", "D", "E", "MS"
  ), 
  covars = list(
    list(covar_selected), # "BMI", 
    list(c("BMI", covar_selected)), #"CC_T2D56", #with BMI
    list(covar_selected), # "CC_T2D56",
    list(c("BMI", covar_selected)), #"CC_T2D61", #with BMI
    list(covar_selected), # "CC_T2D61",
    list(covar_selected), # "CC_OBESITE", 
    list(covar_selected), # "CC_OWT",
    list(c("BMI", covar_selected)), #"FG_D0", #with BMI
    list(covar_selected), # "FG_D0", 
    list(c("BMI", covar_selected)), #"LOG_FI", #with BMI
    list(covar_selected), # "LOG_FI", 
    list(c("BMI", covar_selected)), #"HDL", #with BMI
    list(covar_selected), # "HDL", 
    list(c("BMI", covar_selected)), #"LOG_triglycerides", #with BMI
    list(covar_selected), # "LOG_triglycerides", 
    list(c("BMI", covar_selected)), # A
    list(c("BMI", covar_selected)), 
    list(c("BMI", covar_selected)), 
    list(c("BMI", covar_selected)), 
    list(c("BMI", covar_selected)),
    list(c("BMI", covar_selected)), # "A", "B", "C", "D", "E", "MS"  #with BMI
    list(covar_selected), # A
    list(covar_selected), 
    list(covar_selected), 
    list(covar_selected), 
    list(covar_selected),
    list(covar_selected) # "A", "B", "C", "D", "E", "MS"
  )
)
#tail(traits_adults)

results_adults <- lapply(X = seq_len(nrow(traits_adults)), FUN = function(iline) {
  message("[results_adults] trait ", iline, " ", traits_adults$traits[iline])
  results_ad <- cbind(
    traits_adults[iline, ], 
    data.frame(
      analyses = c(
        "rare-neutral-PLPall", "rare-neutral-PLPAB",
        "rare-neutral-PLPA",
        "rare-all", "P/LP-all", 
        "P/LP-AB",
        "P/LP-C", "P/LP-B",
        "P/LP-A"
      )
    ))
  results_ad$covars = NULL
  results_ad$covars = paste0(unlist(traits_adults[["covars"]][iline]), collapse = ";")
  results_ad$cluster_size <- unlist(lapply(X = seq_len(nrow(results_ad)), FUN = function(ianalysis_line) {
    ianalysis <- results_ad$analyses[[ianalysis_line]]
    sum(data_adultes[!is.na(data_adultes[[traits_adults$traits[iline]]]), ][, c(ianalysis)] != 0, na.rm = T)
  }))
  results_ad$nb_variants <- unlist(lapply(X = seq_len(nrow(results_ad)), FUN = function(ianalysis_line) {
    ianalysis <- results_ad$analyses[[ianalysis_line]]
    length(na.omit(unique(data_adultes[!is.na(data_adultes[[traits_adults$traits[iline]]]), ][, c(ianalysis)])))
  }))
  results_ad$sample_size <- unlist(lapply(X = seq_len(nrow(results_ad)), FUN = function(ianalysis_line) {
    ianalysis <- results_ad$analyses[[ianalysis_line]]
    variant_analysis <- annotation$ACMG[annotation$sheet %in% "adults" & annotation$selection_name %in% ianalysis]
    nrow(na.omit(data_adultes[,
      c(
        traits_adults[["traits"]][iline],
        unlist(traits_adults[["covars"]][iline]),
        variant_analysis
      )
    ]))
  }))
  message("[results_adults] Run Mist ...")
  if (traits_adults[["is_binary"]][iline]) { 
    results_ad$mist_outputs <- lapply(X = results_ad$analyses, FUN = function(ianalysis) {
      message(ianalysis)
      variant_analysis <- annotation$ACMG[annotation$sheet %in% "adults" & annotation$selection_name %in% ianalysis]
      data_here <- na.omit(
        data_adultes[,
          c(
            traits_adults[["traits"]][iline],
            unlist(traits_adults[["covars"]][iline]),
            variant_analysis
          )
        ]
      )
      mist_logit(
        y = data_here[[traits_adults[["traits"]][iline]]], 
        X = data_here[, unlist(traits_adults[["covars"]][iline])], 
        G = data_here[, c(variant_analysis), drop = FALSE], 
        Z = rep(1, length(variant_analysis))
      )
    })
  } else {
    results_ad$mist_outputs <- lapply(X = results_ad$analyses, FUN = function(ianalysis) {
      message(ianalysis)
      variant_analysis <- annotation$ACMG[annotation$sheet %in% "adults" & annotation$selection_name %in% ianalysis]
      data_here <- na.omit(
        data_adultes[,
          c(
            traits_adults[["traits"]][iline],
            unlist(traits_adults[["covars"]][iline]),
            variant_analysis
          )
        ]
      )
      mist_linear(
        y = data_here[[traits_adults[["traits"]][iline]]], 
        X = data_here[, unlist(traits_adults[["covars"]][iline])], 
        G = data_here[, c(variant_analysis), drop = FALSE], 
        Z = rep(1, length(variant_analysis))
      )
    })
  }
  results_ad$mist_scores <- lapply(X = results_ad$mist_outputs, FUN = function(ianalysis) {
    as.data.frame(ianalysis$out_MiST)
  })
  results_ad$burden_res <- lapply(X = results_ad$mist_outputs, FUN = function(ianalysis) {
    as.data.frame(ianalysis$out_rare)
  })
  
  #### save res =====================================================================================
  res_final_ad <- do.call("cbind", list(
    results_ad[, c("data_from", "traits", "covars", "analyses", "sample_size", "cluster_size", "nb_variants")],
    rbindlist(results_ad$mist_scores),
    rbindlist(results_ad$burden_res)
  ))
  message("[results_adults] return")
  res_final_ad
})

#### save res adults  ================================================================================
res_final_adults_traits <- rbindlist(l = results_adults, use.names = TRUE, fill = TRUE)
res_final_adults_traits
writexl::write_xlsx(
  x = res_final_adults_traits, 
  path = file.path(output_directory, "DYRK1B_adultes_res.xlsx")
)

#### Analyses enfants ==============================================================================
message("Enfants ...")

results_enfants <- data.table(
  data_from = "DYRK1B-enfant",
  traits = "CC_Obesity_child",
  covars_list = list(
    list(covar_selected),
    list(covar_selected),
    list(covar_selected),
    list(covar_selected),
    list(covar_selected),
    list(covar_selected)
  ),
  analyses = c(
    "rare-neutral-PLPall", "rare-neutral-PLPB",
    "rare-all", "P/LP-all",
    "P/LP-C", "P/LP-B"
  )
)
results_enfants

results_enfants$covars <- unlist(lapply(X = seq_len(nrow(results_enfants)), FUN = function(ianalysis_line) {
  paste0(unlist(results_enfants[["covars_list"]][ianalysis_line]), collapse = ";")
}))

results_enfants$cluster_size <- unlist(lapply(X = seq_len(nrow(results_enfants)), FUN = function(ianalysis_line) {
  ianalysis <- results_enfants$analyses[[ianalysis_line]]
  sum(data_enfants[!is.na(data_enfants[[results_enfants$traits[ianalysis_line]]]), ][, c(ianalysis)] != 0, na.rm = T)
}))

results_enfants$nb_variants <- unlist(lapply(X = seq_len(nrow(results_enfants)), FUN = function(ianalysis_line) {
  ianalysis <- results_enfants$analyses[[ianalysis_line]]
  length(na.omit(unique(data_enfants[!is.na(data_enfants[[results_enfants$traits[ianalysis_line]]]), ][, c(ianalysis)])))
}))

results_enfants$sample_size <- unlist(lapply(X = seq_len(nrow(results_enfants)), FUN = function(ianalysis_line) {
  ianalysis <- results_enfants$analyses[[ianalysis_line]]
  variant_analysis <- annotation$ACMG[annotation$sheet %in% "children" & annotation$selection_name %in% ianalysis]
  nrow(na.omit(data_enfants[,
    c(
      results_enfants[["traits"]][ianalysis_line],
      unlist(results_enfants[["covars_list"]][ianalysis_line]),
      variant_analysis
    )
  ]))
}))

message("Run Mist ...")
results_enfants$mist_outputs <- lapply(X = results_enfants$analyses, FUN = function(ianalysis) {
  message(ianalysis)

  message(ianalysis)
  variant_analysis <- annotation$ACMG[annotation$sheet %in% "children" & annotation$selection_name %in% ianalysis]
  data_here <- na.omit(
    data_enfants[,
      c(
        "CC_Obesity_child", # only one trait here ! directly named
        covar_selected, # directly its covars
        variant_analysis
      )
    ]
  )
  mist_logit( # only a CC analysis
    y = data_here$CC_Obesity_child,
    X = data_here[, covar_selected], 
    G = data_here[, c(variant_analysis), drop = FALSE], 
    Z = rep(1, length(variant_analysis))
  )
})
results_enfants$mist_scores <- lapply(X = results_enfants$mist_outputs, FUN = function(ianalysis) {
  as.data.frame(ianalysis$out_MiST)
})
results_enfants$burden_res <- lapply(X = results_enfants$mist_outputs, FUN = function(ianalysis) {
  as.data.frame(ianalysis$out_rare)
})

#### save res child ================================================================================
res_final_enfants <- do.call("cbind", list(
  results_enfants[, c("data_from", "traits", "covars", "analyses", "sample_size", "cluster_size", "nb_variants")],
  rbindlist(results_enfants$mist_scores),
  rbindlist(results_enfants$burden_res)
))
res_final_enfants
writexl::write_xlsx(
  x = res_final_enfants, 
  path = file.path(output_directory, "DYRK1B_enfant_res.xlsx"), 
  col_names = TRUE
)

#### end ===========================================================================================
message("Success ! ")
