  
  
  # ####################################################################
  # Nutrition_MDD_Women: Modeling Script
  # Author: [Name]
  # Date: December 2024
  # ####################################################################
  
  # Load Required Libraries
  library(data.table)
  library(INLA)
  library(sf)
  library(glue)
  library(tictoc)
  library(dplyr)
  
  # Data Preparation ---------------------------------------------------
  
  # Define Input Variables
  country <- "Nepal"  # Replace with the country of interest
  survey_year <- 2022
  survey_type <- "DHS"
  indicator <- "nt_wm_mdd"
  
  # Load Survey Data
  data <- read_dta("path_to_survey_data.dta")  # Replace with dataset path
  gps_data <- read.dbf("path_to_gps_data.dbf")  # Replace with GPS file path
  
  # Data Cleaning and Merging
  setnames(data, c('v001'), c('cluster_id'))
  gps_vars <- c('dhsyear', 'dhsclust', 'latnum', 'longnum')
  gps_data <- gps_data[, ..gps_vars]
  setnames(gps_data, c('dhsyear', 'dhsclust', 'latnum', 'longnum'), c('year', 'cluster_id', 'latitude', 'longitude'))
  merged_data <- merge(data, gps_data, by = 'cluster_id')
  
  # Variable Selection
  selected_vars <- c('year', 'cluster_id', 'nt_wm_mdd', 'latitude', 'longitude', 'residence', 
                     'region', 'wmage', 'employment', 'hhwealth', 'edu', 'ph_num_members')
  df <- merged_data[, ..selected_vars]
  
  # Covariate Preparation ----------------------------------------------
  # Load country shapefile and define the polygon of interest
  country_shapefile <- st_read("path_to_shapefile.shp")  # Replace with shapefile path
  country_polygon <- country_shapefile[country_shapefile$ADM0_NAME == country, ]
  
  # Load Covariates
  covariates <- list(
    evi = "mean", worldpop = "mean", nyprotein = "mean", nyvitamin_a = "mean",
    motorizedhf = "mean", walkinghf = "mean", nightlights = "mean"
  )
  
  # Extract Covariates
  source("path_to_extraction_function.R")  # Replace with extraction function path
  covariate_layers <- load_and_crop_covariates(
    covariates, country_polygon, start_year = survey_year, end_year = survey_year
  )
  df <- extract_covariates(df, covariate_layers)
  
  # Variable Selection with VIF ----------------------------------------
  vif_cutoff <- 5
  complete_vars <- names(df)[sapply(names(df), function(v) mean(!is.na(df[[v]])) > 0.9)]
  selected_vars <- reduce_vif(df[, ..complete_vars], vif_cutoff)
  
  # Modeling -----------------------------------------------------------
  # Define Model Formula
  formula <- as.formula(glue("
    nt_wm_mdd ~ 1 + {paste(selected_vars, collapse = ' + ')} +
    f(stratum, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = c(1.5, 0.05)))) +
    f(psu, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = c(1.5, 0.05)))) +
    f(hh, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = c(1.5, 0.05))))
  "))
  
  # Fit the Model
  tic("Model fit timer")
  model <- INLA::inla(
    formula, family = "binomial", data = df,
    control.predictor = list(compute = TRUE),
    control.fixed = list(prec = list(prior = 'pc.prec', param = c(3, 0.05))),
    control.compute = list(dic = TRUE, waic = TRUE)
  )
  toc()
  
  # Summarize Results
  summary(model)
  
  # Save Odds Ratios
  odds_ratios <- model$summary.fixed[-1, c("mean", "0.025quant", "0.975quant")] %>%
    exp() %>%
    as.data.table() %>%
    setnames(c("mean", "lower", "upper")) %>%
    .[, significance := ifelse(lower > 1, "Positive", ifelse(upper < 1, "Negative", "Not Significant"))]
  fwrite(odds_ratios, file = "model_effects.csv")  # Replace with file path
  
  # Output Summary Metrics
  print(paste("DIC:", model$dic$dic))
  print(paste("WAIC:", model$waic$waic))
  
  # Random Effects Summary
  random_effects <- model$summary.hyperpar %>%
    .[, `:=`(
      SD = 1 / sqrt(mean),
      OR = exp(1 / sqrt(mean))
    )]
  print(random_effects)
