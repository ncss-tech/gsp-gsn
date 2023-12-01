
library(terra)
library(caret)
library(ranger)
library(Boruta)
library(data.table)


# covariates ----
fp <- "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data"
fp_cov100 <- file.path("./cov100/covariates")
fp <- "./data"

covs <- read.csv(file.path(fp, "covariates_SoilGrids_USA48_Covs100m.csv"))
covs <- read.csv(file.path(fp, "SoilGrids_USA48_Covs100m.csv"))


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


# Boruta ----
vars <- c("N_pt_log.000-030", vars_ex)
idx <- complete.cases(dat[vars]) & dat$`N_pt.000-030` < 20
fs_bor <- Boruta(
  y = dat$`N_pt_log.000-030`[idx],
  x = dat[idx, vars_ex[-1]], 
  maxRuns = 35, 
  doTrace = 1
)
# saveRDS(fs_bor, file.path(fp, "fs_bor_N.rds"))
fs_bor <- readRDS(file.path(fp, "fs_bor_N.rds"))
Boruta:::plot.Boruta(fs_bor, horizontal = TRUE, ylab = "", xlab = "Importance")

fs_vars <- getSelectedAttributes(fs_bor, withTentative = FALSE)
fs_vars <- attStats(fs_bor) |>
  poorman::arrange(- meanImp)



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
fs_vars2 <- c("NPKyear", "NLCD116_crops", row.names(fs_vars[fs_vars$meanImp >= 15, ]))

mtry <- round(length(fs_vars2) / 3)
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

## 3.3 - Calibrate the ranger model ---------------------------
soilatt <- "N_pt_log.000-030"
vars <- c(soilatt, fs_vars2)
idx <- complete.cases(dat[vars]) & dat$`N_pt.000-030` < 20

model_rn <- caret::train(
  y = dat[idx, soilatt], x = dat[idx, fs_vars2],
  method = "ranger",
  quantreg = TRUE,
  importance = "permutation",
  trControl = fitControl,
  verbose = TRUE,
  tuneGrid = tuneGrid
)
print(model_rn)
print(model_rn$bestTune)



# model ----
vars <- c("N_pt_log.000-030", "NPKyear", "NLCD116_crops", row.names(fs_vars[fs_vars$meanImp >= 15, ]))
train <- dat[vars]
test <- ranger(x = train[-1], y = train[, 1], quantreg = TRUE, importance = "permutation")
test

vip::vip(test)
pdp::partial(test, pred.var = "NPKyear", plot = TRUE, rug = TRUE, smooth = F, train = train, title = "Partial Effect of the 6 Year Interval in the Random Forest", grid.resolution = 20)




# predict ----

pfun <- function(model, ...) exp(ranger:::predict.ranger(model, ...)$predictions) - 0.01

cfun <- function(i) {
  rs <- rast(tiles$fp[tiles$tile == i])
  
  names(rs)[names(rs) == "Layer_1"] <- "NLCD116"
  
  NLCD116_crops <- rs[["NLCD116"]]
  vals <- values(NLCD116_crops)
  values(NLCD116_crops) <- ifelse(!is.na(vals) & vals == 82, 1, 0)
  names(NLCD116_crops) <- "NLCD116_crops"
  
  NPKyear <- rast(resolution = res(rs), extent = ext(rs), crs = crs(rs), vals = 2017, names = "NPKyear")
  rs <- c(rs, NPKyear, NLCD116_crops)
  terra::predict(rs, test, fun = pfun, index = 1, progress = "text", na.rm = TRUE, filename = paste0("./tiles/a_predict_N_", i, ".tif"), overwrite = TRUE)
}

future::plan(future::multisession, workers = 20) # 30)
future.apply::future_lapply(1:104, cfun)
lapply(102, cfun)


tiles_fp <- list.files("./tiles", pattern = "a_predict_N_[0-9]{1,3}.tif$", full.names = TRUE) 
tiles_r <- tiles_fp |> 
  as.list() |>
  lapply(rast)

plot(vrt(tiles_fp))
merge(sprc(p), filename = "a_predict_N.tif", overwrite = TRUE, progress = 1)








