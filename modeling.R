
library(terra)
library(caret)
library(ranger)
library(Boruta)
library(data.table)


# load snapshots ----
fp <- "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data"
covs <- read.csv(file.path(fp, "covariates_SoilGrids_USA48_Covs100m.csv"))


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
Boruta:::plot.Boruta(fs_bor, horizontal = TRUE, ylab = "", xlab = "Importance")
                     
fs_vars <- getSelectedAttributes(fs_bor, withTentative = FALSE)
fs_vars <- attStats(fs_bor) |>
  poorman::arrange(- meanImp)


test <- ranger(x = dat[idx, c("NPKyear", "NLCD116_crops", row.names(fs_vars[fs_vars$meanImp >= 15, ]))], y = dat$`N_pt_log.000-030`[idx], quantreg = TRUE, importance = "permutation")

