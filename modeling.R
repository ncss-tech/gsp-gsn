
library(terra)
library(caret)
library(ranger)
library(Boruta)
library(data.table)


# covariates ----
# fp <- "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data"
fp_cov100 <- file.path("./cov100/covariates")
fp <- "./data"

covs <- read.csv(file.path(fp, "covariates_SoilGrids_USA48_Covs100m.csv"))
covs <- read.csv(file.path(fp_cov100, "SoilGrids_USA48_Covs100m.csv"))


tiles <- list.files("./tiles", pattern = ".tif$", full.names = TRUE) |>
  cbind(fp = _) |>
  as.data.frame()
tiles$tile <- tiles$fp |>
  strsplit(split = "_|.tif|_pmkind-ssurgo.tif$") |>
  sapply(function(x) x[length(x)]) |>
  as.integer()
tiles$fn <- tiles$fp |>
  strsplit(split= "\\./tiles/|_[0-9]{1,3}.tif$|_[0-9]{1,3}_pmkind-ssurgo.tif$") |>
  sapply(function(x)  x[length(x)])
tiles <- poorman::arrange(tiles, tile)
tiles <- subset(tiles, !grepl("predict", fp))



# load point data snapshots ----
# fp <- "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data"

load(file = file.path(fp, "gsn_seg_t_df.RData"))
h <- gsn_seg_df; rm(gsn_seg_df)
vars <- names(h)[2:16]
# h_wi <- data.table::as.data.table(h)|>
#   data.table::dcast(formula = pedon_key ~ segment_id, value.var = vars) |>
#   as.data.table()
h_wi <- reshape(
  h,
  direction = "wide",
  idvar = "pedon_key",
  timevar = "segment_id",
  v.names = vars
)

scd_sf <- readRDS(file.path(fp, "scd_site_nasis_sf.rds"))
s <- cbind(scd_sf, sf::st_coordinates(scd_sf)) |>
  sf::st_drop_geometry() |>
  subset(valid_xy == TRUE) |>
  transform(
    decade = substr(year, 1, 3) |> as.integer() * 10,
    NPKyear = cut(
      year, 
      breaks = seq(1885, 2027, 6),
      labels = seq(1885, 2021, 6),
      right = FALSE
    ) |> as.character() |> as.integer(),
    CECyear = cut(
      year, 
      breaks = seq(1851, 2027, 43),
      labels = seq(1851, 2021, 43),
      right = FALSE
    ) |> as.character() |> as.integer(),
    Cyear = cut(
      year, 
      breaks = seq(1885, 2026, 23),
      labels = seq(1885, 2022, 23),
      right = FALSE
    ) |> as.character() |> as.integer()
  )

ex <- readRDS(file = file.path(fp, "scd_site_extract_sf.rds")) |>
  sf::st_drop_geometry() |>
  transform(
    NLCD116_crops = ifelse(as.character(NLCD116) %in% c("82"), 1, 0) |> as.factor()
  )
nm <- names(ex)
idx <- grepl("P[0-9]{2}CHL5|N[0-9]{2}PRI5|LNDCOV6" , nm)
tm <- nm[idx]
ex2 <- ex[, !nm %in% tm]



vars_s  <- c("site_key", "pedon_key", "pedlabsampnum", "year", "X", "Y", "x_precision", "y_precision", "dups", "n_coords", "n_year", "valid_xy", "NPKyear", "CECyear", "Cyear")
vars_ex <- names(ex2)[c(1, 24:ncol(ex2))]

dat <- merge(s[vars_s], h_wi, by = "pedon_key", all.x = TRUE, sort = FALSE) |>
  merge(ex2[vars_ex], by = "site_key", all.x = TRUE, sort = FALSE)




# feature selection (Boruta) ----
s_vars <- c("N_pt_ppm_log.000-030", "P_pt_ppm_log.000-030", "K_pt_ppm_log.000-030")
s_vars2 <- sapply(s_vars, function(x) strsplit(x, "_")[[1]][1])

fs_bor <- lapply(s_vars[3], function(x) {
  
  vars <- c(x, vars_ex)
  
  idx <- NA
  if (x ==  "N_pt_ppm_log.000-030") idx <- dat[gsub("_log", "", x)] < 20
  if (x ==  "P_pt_ppm_log.000-030") idx <- dat[gsub("_log", "", x)] < 3000
  if (x ==  "K_pt_ppm_log.000-030") idx <- dat[gsub("_log", "", x)] < 900
  if (all(is.na(idx))) idx <- rep(TRUE, nrow(dat))
  
  idx <- complete.cases(dat[vars]) & idx
  
  fs_bor <- Boruta(
    y = dat[idx, x],
    x = dat[idx, vars_ex[-1]], 
    maxRuns = 35, 
    doTrace = 1
  )
  
  saveRDS(fs_bor, file.path(fp, paste0("fs_bor_", strsplit(x, "_")[[1]][1], ".rds")))
}) 


fs_bor <- lapply(s_vars, function(x) {
  x2 <- strsplit(x, "_")[[1]][1]
  fs <- readRDS(file.path(fp, paste0("fs_bor_", x2, ".rds")))
})
names(fs_bor) <- s_vars


Boruta:::plot.Boruta(fs_bor[[1]], horizontal = TRUE, ylab = "", xlab = "Importance")
fs_vars <- getSelectedAttributes(fs_bor, withTentative = FALSE)


fs_vars <- lapply(fs_bor, function(x) {
  attStats(x) |>
    poorman::arrange(- meanImp)
})

fs_vars2 <- lapply(fs_vars, function(x) {
  x[x$meanImp > floor(x$min[1:10]), ] |> row.names()
})


# calibration ----

# 3 - QRF Model calibration with ranger =======================
## 3.1 - Set training parameters ------------------------------

sf <- function(data, lev = NULL, model = "ranger") {
  data2 <- exp(data) + 0.01
  defaultSummary(data2, lev = lev, model = model)
}

fitControl <- trainControl(
  summaryFunction = sf,
  method  = "repeatedcv",
  number  = 10,
  repeats = 10,
  savePredictions = TRUE
)


## 3.2 - Tune hyperparameters ---------------------------------

rf_par <- lapply(s_vars2, function(x) {
  
  pick_y <- function(x) {
    switch(
      x, 
      "N" = "N_pt_ppm_log.000-030",
      "P" = "P_pt_ppm_log.000-030",
      "K" = "K_pt_ppm_log.000-030"
    )}
  
  pick_year <- function(x) {
    switch(
      x,
      "N" = "NPKyear",
      "P" = "NPKyear",
      "K" = "NPKyear"
    )}
  
  y        <- pick_y(x)
  var_year <- pick_year(x)
  vars     <- c(var_year, fs_vars2[[y]])
  
  mtry <- round(length(vars) / 3)
  tuneGrid <-  expand.grid(
    mtry = abs(c(
      mtry-round(mtry / 2),
      mtry-round(mtry / 3), 
      mtry, 
      mtry + round(mtry / 3),
      mtry + round(mtry / 2)
    )),
    min.node.size = 5,
    splitrule = c("variance", "extratrees", "maxstat")
  )
  
  list(tuneGrid = tuneGrid, vars = vars, y = y)
})
names(rf_par) <- s_vars



## 3.3 - Calibrate the ranger model ---------------------------

rf_l <- lapply(s_vars, function(x) {
  
  tuneGrid <- rf_par[[x]]$tuneGrid
  y        <- rf_par[[x]]$y
  vars     <- rf_par[[x]]$vars
  
  idx <- NA
  if (x ==  "N_pt_ppm_log.000-030") idx <- dat[gsub("_log", "", x)] < 60000
  if (x ==  "P_pt_ppm_log.000-030") idx <- dat[gsub("_log", "", x)] < 900
  if (x ==  "K_pt_ppm_log.000-030") idx <- dat[gsub("_log", "", x)] < 8000
  if (all(is.na(idx))) idx <- rep(TRUE, nrow(dat))
  
  idx <- complete.cases(dat[c(y, vars)]) & idx
  cat(summary(idx), "\n")
  
  # ord <- tapply(dat[[x]], dat$GEOUSG6, mean, na.rm = TRUE)
  # ord[is.na(ord)] <- mean(ord, na.rm = TRUE)
  # lev <- names(sort(ord, decreasing = TRUE))
  # dat$GEOUSG6 <- factor(as.character(dat$GEOUSG6), levels = lev, ordered = TRUE)
  # 
  # ord <- tapply(dat[[x]], dat$mlra52_aea_100m, mean, na.rm = TRUE)
  # ord[is.na(ord)] <- mean(ord, na.rm = TRUE)
  # lev <- names(sort(ord, decreasing = TRUE))
  # dat$mlra52_aea_100m <- factor(as.character(dat$mlra52_aea_100m), levels = lev, ordered = TRUE)
  # 
  # ord <- tapply(dat[[x]], dat$DRNGSS7_f , mean, na.rm = TRUE)
  # ord[is.na(ord)] <- mean(ord, na.rm = TRUE)
  # lev <- names(sort(ord, decreasing = TRUE))
  # dat$DRNGSS7_f  <- factor(as.character(dat$DRNGSS7_f ), levels = lev, ordered = TRUE)
  # 
  # ord <- tapply(dat[[x]], dat$NLCD116, mean, na.rm = TRUE)
  # ord[is.na(ord)] <- mean(ord, na.rm = TRUE)
  # lev <- names(sort(ord, decreasing = TRUE))
  # dat$NLCD116 <- factor(as.character(dat$NLCD116), levels = lev, ordered = TRUE)
  # 
  # ord <- tapply(dat[[x]], dat$PVEGKT6, mean, na.rm = TRUE)
  # ord[is.na(ord)] <- mean(ord, na.rm = TRUE)
  # lev <- names(sort(ord, decreasing = TRUE))
  # dat$PVEGKT6 <- factor(as.character(dat$PVEGKT6), levels = lev, ordered = TRUE)
  # 
  
  model <- caret::train(
    y = dat[idx, y], x = dat[idx, vars],
    method = "ranger",
    quantreg = TRUE,
    importance = "permutation",
    respect.unordered.factors = TRUE,
    keep.inbag = TRUE,
    trControl = fitControl,
    # trControl = trainControl(method = "cv", number = 3, summaryFunction = sf), 
    verbose = TRUE,
    tuneGrid = tuneGrid,
    # tuneGrid = data.frame(mtry = 24, min.node.size = 5, splitrule = "variance")
  )
  
  saveRDS(model, file = file.path(fp, paste0("model_rf_", x, ".rds")))
})

rf_l <- lapply(s_vars, function(x) {
  model_rn2 <- readRDS(file.path(fp, paste0("model_rf_", x, ".rds")))
})
names(rf_l) <- s_vars

rf_l



# model ----
vars <- c("N_pt_log.000-030", "NPKyear", "NLCD116_crops", row.names(fs_vars[fs_vars$meanImp >= 15, ]))
train <- dat[vars]
test <- ranger(x = train[-1], y = train[, 1], quantreg = TRUE, importance = "permutation")
test

vip::vip(model_rn2$finalModel)
pdp::partial(model_rn2$finalModel, pred.var = "NPKyear", plot = TRUE, rug = TRUE, smooth = F, train = model_rn2$trainingData, title = "Partial Effect of the 6 Year Interval in the Random Forest", grid.resolution = 20)




# predict ----

pfun1 <- function(model, ...) exp(ranger:::predict.ranger(model, ..., na.exclude = TRUE)$predictions) - 0.01

pfun2 <- function(model, ...) {
  ranger:::predict.ranger(model, ...)$predictions |> 
    t() |>
    exp() - 0.01
}


cfun <- function(i, var, model, nt = NULL, nc = 1, cond = FALSE) {
  
  cat(i, as.character(Sys.time()), "\n")
  
  tiles2 <- tiles[tiles$tile == i, ]
  rs <- rast(tiles2$fp)
  
  names(rs)[names(rs) == "Layer_1"] <- "NLCD116"
  names(rs) <- gsub("-", ".", names(rs))
  
  # create missing rasters
  NLCD116_crops <- rs[["NLCD116"]]
  vals <- values(NLCD116_crops)
  values(NLCD116_crops) <- ifelse(!is.na(vals) & vals == 82, 1, 0)
  names(NLCD116_crops)  <- "NLCD116_crops"
  
  NPKyear <- rast(resolution = res(rs), extent = ext(rs), crs = crs(rs), vals = 2017, names = "NPKyear")
  
  rs <- c(rs, NPKyear, NLCD116_crops)
  
  idx <- names(rs) %in% names(model$trainingData)
  rs <- rs[[idx]]
  
  # some rasters are missing values while others aren't, this is causing an issues during predict, seems to be related to factor variables
  # I think terra is sampling the SpatRaster
  
  # find tiles with missing values
  idx2 <- as.int(!anyNA(rs))
  idx2[idx2 == 0]  <- NA
  
  
  # trim tiles with NA values
  if (all(minmax(idx2) != "NaN")) {
    idx2 <- trim(idx2)
  } 
  
  # crop the raster stack to match the factor rasters
  if (dim(rs)[1] != dim(idx2)[1]) {
    rs <- crop(rs, trim(idx2))
  }
  
  test <- idx2 |> minmax() |> sum() > 0
  if (is.na(test)) test <- FALSE
  
  if (test) {
    
    terra::predict(
      rs,
      model$finalModel,
      fun       = if (cond == TRUE) pfun2 else pfun1,
      na.rm     = TRUE,
      type      = ifelse(cond == TRUE, "quantiles", "response"),
      what      = if (cond == TRUE) mean else NULL,
      num.threads = nt,
      cores = nc,
      cpkgs = if (nc > 1) "ranger",
      progress  = "text",
      filename  = paste0(
        "./tiles/a_predict_", var, "_", if (cond == TRUE) "avg_" else "", i, ".tif"
      ),
      overwrite = TRUE
    )
    
    # terra::predict(
    #   rs,
    #   model$finalModel,
    #   fun       = if (cond == TRUE) pfun2 else pfun1,
    #   na.rm     = TRUE,
    #   type      = ifelse(cond == TRUE, "quantiles", "se"),
    #   se.method = "infjack",
    #   what      = if (cond == TRUE) sd else NULL,
    #   progress  = "text",
    #   filename  = paste0(
    #     "./tiles/a_predict_", var, "_", if (cond == TRUE) "sd_" else "", i, ".tif"
    #   ),
    #   overwrite = TRUE
    # )
  }
  gc(reset = TRUE)
}

# test2 <- list.files("./tiles", pattern = "a_predict_N_pt_ppm_log")
# idx <- (1L:104L)[! 1L:104L %in% unique(test2$tile)] |> dput()


future::plan(future::sequential)
future::plan(future::multisession, workers = 5)

system.time(future.apply::future_mapply(cfun, i = c(15, 19), MoreArgs = list(var = s_vars[1], model = rf_l$`N_pt_ppm_log.000-030`)))
future.apply::future_mapply(cfun, i = c(23, 26), MoreArgs = list(var = s_vars[2], model = rf_l$`P_pt_ppm_log.000-030`))
future.apply::future_mapply(cfun, i = c(1:23, 26:104), MoreArgs = list(var = s_vars[3], model = rf_l$`K_pt_ppm_log.000-030`))

system.time(lapply(39:43, function(i) cfun(i, s_vars[1], rf_l$`N_pt_ppm_log.000-030`)))


lapply(1:104, function(i) cfun(i, s_vars[1], rf_l$`N_pt_ppm_log.000-030`))


tiles_N <- list.files("./tiles", pattern = "a_predict_N_pt_ppm_log", full.names = TRUE) 
plot(vrt(tiles_N), col = viridis::viridis(10), breaks = quantile(dat$`N_pt_ppm.000-030`, probs = seq(0, 1, 0.1), na.rm = TRUE), main = "Total N ppm")

tiles_P <- list.files("./tiles", pattern = "a_predict_P_pt_ppm_log", full.names = TRUE) 
plot(vrt(tiles_P), col = viridis::viridis(10), breaks = quantile(dat$`P_pt_ppm.000-030`, probs = seq(0, 1, 0.1), na.rm = TRUE), main = "Available P ppm")

tiles_K <- list.files("./tiles", pattern = "a_predict_K_pt_ppm_log", full.names = TRUE) 
plot(vrt(tiles_K), col = viridis::viridis(10), breaks = quantile(dat$`K_pt_ppm.000-030`, probs = seq(0, 1, 0.1), na.rm = TRUE), main = "Available K ppm")

tiles_N_r <- tiles_N |> as.list() |> lapply(rast)
tiles_P_r <- tiles_P |> as.list() |> lapply(rast)
tiles_K_r <- tiles_K |> as.list() |> lapply(rast)
merge(sprc(tiles_N_r), filename = "a_predict_N.tif", overwrite = TRUE, progress = 1)
merge(sprc(tiles_P_r), filename = "a_predict_P.tif", overwrite = TRUE, progress = 1)
merge(sprc(tiles_K_r), filename = "a_predict_K.tif", overwrite = TRUE, progress = 1)


merge(sprc(tiles_N_r), filename = "a_predict_N_2.tif", overwrite = TRUE, progress = 1, gdal = c("COMPRESS=DEFLATE"))






