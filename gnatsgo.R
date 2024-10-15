
library(aqp)
library(soilDB)
library(dplyr)
library(sf)
library(terra)


# load data ----
ss <- get_soilseries_from_NASIS()


## gNATSGO ----
fp_nat <- file.path("D:/geodata/soils/gNATSGO_CONUS_Oct2023")
dsn <- file.path(fp_nat, "gNATSGO_CONUS.gdb")
gnatsgo <- fetchGDB(dsn, child = TRUE)
# save(gnatsgo, file = "D:/geodata/soils/gNATSGO_CONUS_Oct2023/gNATSGO_CONUS_Oct2023_fetchGDB.RData")
load(file = file.path(fp_nat, "gNATSGO_CONUS_Oct2023_fetchGDB.RData"))
mu <- get_mapunit_from_GDB(dsn, stats = TRUE) |>
  transform(mukey = as.integer(mukey)) |>
  subset(!duplicated(mukey))


# find missing pmkind
mukey_mis_pmkind <- gnatsgo |>
  site() |>
  subset(is.na(pmkind) & compkind != "Miscellaneous area", select = c(mukey, cokey))
n <- nrow(mukey_mis_pmkind)
brks <- classInt::classIntervals(1:n, n = 10)
mukey_mis_pmkind$idx <- 1:n |>
  cut(x = _, breaks = brks$brks)

test <- split(mukey_mis_pmkind, list(mukey_mis_pmkind$idx)) |>
  lapply(function(x) {
    fetchSDA(paste0("mukey IN ", soilDB::format_SQL_in_statement(x$mukey)), duplicates = TRUE)
  }) 
# save(test, file = file.path(fp_nat, "mukey_mis_pmkind_SDA.RData"))
test2 <- lapply(test, function(x) site(x)) |> do.call("rbind", args = _)
test2 <- test2[test2$cokey %in% mukey_mis_pmkind$cokey, ]
table(test2$pmkind)


# gNATSGO raster
gnatsgo_r <- rast(file.path(fp_nat, "FY2024_gNATSGO_mukey_grid.tif"))


# create tiles
tiles <- ext(gnatsgo_r) |> 
  # st_make_grid(cellsize = 15000) |> 
  st_make_grid(cellsize = 10000 * 30) |>
  # st_make_grid(n = c(10, 7)) |> 
  st_as_sf(crs = crs(gnatsgo_r, proj = TRUE)) |>
  transform(id = 1:length(x)) |>
  vect()
plot(tiles)

makeTiles(gnatsgo_r, tiles, filename = file.path(fp_nat, "gnatsgo_.tif"), overwrite = TRUE)
l_rat <- list.files(path = fp_nat, pattern = "gnatsgo_[0-9]{1,3}.tif$", full.names = TRUE)


future::plan(future::multisession, workers = 8)
idx <- future.apply::future_lapply(l_rat, function(x) {
  cat(x, as.character(Sys.time()), "\n")
  any(values(rast(x)) > 0)
  }) |>
  unlist()
id <- sapply(l_rat, function(x) {
  strsplit(x, "_|\\.")[[1]][4] |> as.integer()
})
df <- data.frame(id, idx, l_rat)
df <- df[order(df$id), ]

tiles <- tiles[df$idx, ]
tiles <- cbind(tiles, df[df$idx, -1])
# saveRDS(tiles, file = file.path(fp_nat, "tiles.rds"))
tiles_nona <- readRDS(file.path(fp_nat, "tiles.rds"))
write_sf(st_as_sf(tiles_nona), file.path(fp_nat, "tiles.shp"))
l_rat <- l_rat[l_rat %in% tiles_nona$l_rat]




## STATSGO ----
fp_sta <- "D:/geodata/soils/STATSGO2/wss_gsmsoil_US_[2016-10-13]"
statsgo_sf <- read_sf(file.path(fp_sta, "US/spatial"), layer = "soilmu_a_us") |> 
  transform(MUKEY = as.integer(MUKEY)) |>
  st_transform(crs = crs(gnatsgo_r, proj = TRUE))
statsgo_v  <- statsgo_sf |> vect()


statsgo_r <- rast(
  vals = NA_integer_, 
  extent = ext(gnatsgo_r), 
  crs = crs(gnatsgo_r),
  resolution = res(gnatsgo_r)
  )
statsgo_r <- rasterize(
  statsgo_v,
  statsgo_r,
  field = "MUKEY",
  filename = file.path(fp_sta, "statsgo.tif"),
  wopt = list(datatype = "INT4U"),
  # overwrite = TRUE
  )
statsgo_r <- rast(file.path(fp_sta, "statsgo.tif"))
NAflag(statsgo_r) <- 0


makeTiles(statsgo_r, tiles, filename = file.path(fp_sta, "statsgo_.tif"), overwrite = TRUE)

l_rat <- list.files(path = fp_sta, pattern = "statsgo_[0-9]{1,3}.tif$")
tiles_nona <- readRDS(file.path(fp_nat, "tiles.rds"))
tiles_nona$sta <- gsub(paste0(fp_nat, "/"), "", tiles_nona$l_rat)
tiles_nona$sta <- gsub("gnatsgo", "statsgo", tiles_nona$sta)

file.remove(file.path(fp_sta, l_rat[!l_rat %in% tiles_nona$sta]))


# harmonize ----
fp_nat <- file.path("D:/geodata/soils/gNATSGO_CONUS_Oct2023")
gnatsgo_l_top <- list.files(path = fp_nat, pattern = "gnatsgo_[0-9]{1,3}_pmkind_top2.tif$", full = TRUE)
gnatsgo_l_bot <- list.files(path = fp_nat, pattern = "gnatsgo_[0-9]{1,3}_pmkind_bot2.tif$")

fp_sta <- "D:/geodata/soils/STATSGO2/wss_gsmsoil_US_[2016-10-13]"
statsgo_l <- list.files(path = fp_sta, pattern = "statsgo_[0-9]{1,3}.tif$", full = TRUE)


hmz <- function(input, zonal, workers = 8) {
  
  f <- function(input, zonal) {
    cat("processing ", input, as.character(Sys.time()), "\n")
    ri <- rast(input) |> values()
    rz <- rast(zonal) |> values()
    
    r <- cbind(ri, rz) |> data.table::as.data.table()
    names(r) <- c("input", "zonal")
    r <- r[, .N, by = .(zonal, input)]
    r <- r[order(zonal, input, -N)] 
  }
  
  future::plan(future::multisession, workers = workers)
  future.apply::future_mapply(f, input, zonal, SIMPLIFY = FALSE)
}

test <- hmz(gnatsgo_l_top, statsgo_l, 6)
test <- do.call("rbind", test)
# saveRDS(test, file = file.path(fp_sta, "pmkind_top2_lookup.rds"))
pmkind_bot2_lookup <- readRDS(file = file.path(fp_sta, "pmkind_bot2_lookup.rds"))

statsgo_lookup <- pmkind_bot2_lookup[
  ,.(mukey = zonal, pmkind_bot2 = input, N)][
    !is.na(pmkind_bot2) & mukey > 0][
      , .(N = sum(N)), by = .(mukey, pmkind_bot2)][
        order(mukey, -N)][
          !duplicated(mukey)] |>
  as.data.frame()



# pmkind ----
s <- gnatsgo |>
  site() |>
  transform(mukey = as.integer(mukey))
s$pmkind_top <- sapply(s$pmkind, function(x) {
  strsplit(x, " over ") |> unlist() ->.;
  .[1]
})
s$pmkind_bot <- sapply(s$pmkind, function(x) {
  strsplit(x, " over ") |> unlist() ->.;
  .[length(.)]
})

gen_pmkind <- function(x) {

    pmkind2 = tolower(x)
    pmkind2 = ifelse(grepl("mine spoil|dredge|human",      pmkind2), "anthropogenic", pmkind2)
    pmkind2 = ifelse(grepl("organic|coprogenic",           pmkind2), "organic",       pmkind2)
    pmkind2 = ifelse(grepl("eolian|parna",                 pmkind2), "eolian",        pmkind2)
    pmkind2 = ifelse(grepl("^ash| ash|pumice|tuff|tephra", pmkind2), "ash",           pmkind2)
    pmkind2 = ifelse(grepl("loess",                        pmkind2), "loess",         pmkind2)
    pmkind2 = ifelse(grepl("glaciofluvial deposits",       pmkind2), "outwash",       pmkind2)
    pmkind2 = ifelse(grepl("till|drift|cryoturbate",       pmkind2), "till",          pmkind2)
    pmkind2 = ifelse(grepl("pedisediment|alluvium",        pmkind2), "alluvium",      pmkind2)
    pmkind2 = ifelse(grepl("marine|estuarine|marl|beach",  pmkind2), "marine",        pmkind2)
    pmkind2 = ifelse(grepl("lacustrine|backswamp|overbank|diatomaceou|diamicton|lagoonal", pmkind2), "lacustrine", pmkind2)
    pmkind2 = ifelse(grepl("colluvium|creep",              pmkind2), "colluvium", pmkind2)
    pmkind2 = ifelse(grepl("residuum|grus|siltstone|limestone|dolomite|igneous|quartzite|metamorphic|sandstone|shale|conglomerate|cinders|saprolite|greensands|scoria|mixed|breccia", pmkind2),  "residuum",  pmkind2)
    pmkind2 = ifelse(grepl("earth|debris|valley side|slope|talus|mass movement|spread|slump block|solifluction|lahar|pyroclastic|avalanche|scree|rockfall|topple|slide |fall |flow ", pmkind2), "mass movement", pmkind2)
  }


co_pmkind <- s |>
  left_join(mu[c("mukey", "muacres", "areasymbol")], by = "mukey") |>
  within({
    coacres = muacres * comppct_r / 100
    coacres = ifelse(is.na(coacres), 0, coacres)
    compname2 = tolower(compname)
    compname2 = sapply(compname2, function(x) strsplit(x, ",")[[1]][1])
    ss = compname2 %in% tolower(ss$soilseriesname)
    taxpartsize2 = taxpartsize
    taxpartsize2 = gsub("-skeletal", "", taxpartsize2)
    taxpartsize2 = gsub("sandy or sandy", "sandy", taxpartsize2)
  }) |>
  within({
    pmkind2     = gen_pmkind(pmkind)
    pmkind_top2 = gen_pmkind(pmkind_top)
    pmkind_bot2 = gen_pmkind(pmkind_bot)
  })


## impute missing ----
mis <- co_pmkind |>
  filter(ss == TRUE) |>
  group_by(compname2, pmkind2, pmkind_top2, pmkind_bot2) |>
  summarize(
    acres = sum(coacres, na.rm = TRUE),
    pmkind_pct   = round(acres / muacres[1] * 100)
    ) |>
  ungroup() |>
  arrange(compname2, - acres) |>
  filter(! duplicated(compname2)) |>
  select(compname2, pmkind2_mis = pmkind2, pmkind_top2_mis = pmkind_top2, pmkind_bot2_mis = pmkind_bot2)

co_pmkind2 <- co_pmkind |>
  left_join(mis, by = "compname2") |>
  mutate(
    pmkind2     = ifelse(!is.na(pmkind2),     pmkind2,     pmkind2_mis),
    pmkind_top2 = ifelse(!is.na(pmkind_top2), pmkind_top2, pmkind_top2_mis),
    pmkind_bot2 = ifelse(!is.na(pmkind_bot2), pmkind_bot2, pmkind_bot2_mis),
    pmkind2_mis     = NULL,
    pmkind_top2_mis = NULL,
    pmkind_bot2_mis = NULL
  )


## dominant condition ----
dcon_pmkind <- co_pmkind2 %>%
  group_by(mukey, pmkind2) %>%
  summarize(
    acres = sum(coacres, na.rm = TRUE),
    pmkind_pct   = round(acres / muacres[1] * 100)
  ) %>%
  ungroup() %>%
  arrange(mukey, - acres) %>%
  filter(! duplicated(mukey)) %>%
  select(mukey, pmkind_pct, pmkind2)


dcon_pmkind_top <- co_pmkind2 %>%
  group_by(mukey, pmkind_top2) %>%
  summarize(
    acres = sum(coacres, na.rm = TRUE),
    pmkind_top_pct   = round(acres / muacres[1] * 100)
  ) %>%
  ungroup() %>%
  arrange(mukey, - acres) %>%
  filter(! duplicated(mukey)) %>%
  select(mukey, pmkind_top_pct, pmkind_top2)


dcon_pmkind_bot <- co_pmkind2 %>%
  group_by(mukey, pmkind_bot2) %>%
  summarize(
    acres = sum(coacres, na.rm = TRUE),
    pmkind_bot_pct   = round(acres / muacres[1] * 100)
  ) %>%
  ungroup() %>%
  arrange(mukey, - acres) %>%
  filter(! duplicated(mukey)) %>%
  select(mukey, pmkind_bot_pct, pmkind_bot2)

dcon_pmkind_all <- mu["mukey"] |>
  merge(dcon_pmkind,     by = "mukey", all.x = TRUE) |>
  merge(dcon_pmkind_top, by = "mukey", all.x = TRUE) |>
  merge(dcon_pmkind_bot, by = "mukey", all.x = TRUE) |>
  as.data.frame()
# saveRDS(dcon_pmkind_all, file.path(fp_nat, "dcon_pmkind.rds"))
dcon_pmkind_all <- readRDS(file.path(fp_nat, "dcon_pmkind.rds"))



# deratify ----
tiles <- readRDS(file.path(fp_nat, "tiles.rds"))
l_rat <- tiles$l_rat[tiles$idx]

derat <- function(l, dat, vars, workers = 8) {
  
  future::plan(future::multisession, workers = workers)
  
  lapply(vars, function(var) {
    future.apply::future_lapply(l, function(x) {
      # lapply(l, function(x) {
      cat("processing ", var, x, as.character(Sys.time()), "\n")
      r <- rast(x)
      val <- values(r) |> as.integer() |> data.table::data.table(mukey = _)
      if (any(val$mukey %in% dat$mukey)) {
          val <- data.table::merge.data.table(
            val, 
            dat[c("mukey", var)], 
            by = "mukey", all.x = TRUE, 
            sort = FALSE
            ) |> as.data.frame()
          
            if (!is.numeric(val[[var]])) {
              lev <- dat[[var]] |> as.factor() |> levels()
              val[[var]] <- val[[var]] |> factor(levels = lev) |> as.integer()
            }
          values(r) <- val[[var]]
          writeRaster(r, file = gsub(".tif", paste0("_", var, ".tif"), x), overwrite = TRUE)
      }
    })})
}


## pmkind ----
derat(l_rat, dcon_pmkind_all, c("pmkind2", "pmkind_top2", "pmkind_bot2")[3])

var <- "pmkind_bot2"
l2 <- list.files(path = fp_nat, pattern = paste0("gnatsgo_[0-9]{1,3}_", var, ".tif$"), full.names = TRUE)
plot(vrt(l2), col = viridis::viridis(13)); plot(tiles_nona, add = T)
test <- sprc(l2)
fn <- file.path(fp_nat, paste0("gnatsgo_Oct22_", var, ".tif"))
merge(test, filename = fn, datatype = "INT1U", overwrite = TRUE)


## statsgo pmkind ----
fp_sta <- "D:/geodata/soils/STATSGO2/wss_gsmsoil_US_[2016-10-13]"
statsgo_l <- list.files(path = fp_sta, pattern = "statsgo_[0-9]{1,3}.tif$", full = TRUE)

pmkind_lookup <- readRDS(file = file.path(fp_sta, "pmkind_top2_lookup.rds"))
statsgo_lookup <- pmkind_lookup[
    !is.na(input) & zonal > 0][
  , .(N = sum(N)), by = .(zonal, input)][
    order(zonal, -N)][
      !duplicated(zonal)][
        ,.(
          mukey       = zonal, 
          pmkind_top2 = input, 
          N)] |>
  as.data.frame()

View(subset(statsgo_lookup, !is.na(pmkind_bot2) & mukey > 0 & mukey == 658396))
View(subset(pmkind_bot2_lookup, zonal == 660574))

derat(statsgo_l, statsgo_lookup, c("pmkind_top2"))

l2 <- list.files(path = fp_sta, pattern = "statsgo_[0-9]{1,3}_pmkind_top2.tif$", full.names = TRUE)
plot(vrt(l2), col = viridis::viridis(13)); plot(tiles_nona, add = T)

test <- sprc(l2)
fn <- file.path(fp_sta, "statsgo_pmkind_top.tif")
merge(test, filename = fn, datatype = "INT1U", overwrite = TRUE)



# merge, resample, crop ----
ex <- rast("D:/geodata/project_data/ADIUCL5.tif")

gnatsgo_pmkind_top <- file.path(fp_nat, "gnatsgo_Oct22_pmkind_top.tif") |> rast()
gnatsgo_pmkind_bot <- file.path(fp_nat, "gnatsgo_Oct22_pmkind_bot.tif") |> rast()

statsgo_pmkind_top <- file.path(fp_sta, "statsgo_pmkind_top.tif") |> rast()
statsgo_pmkind_bot <- file.path(fp_sta, "statsgo_pmkind_bot.tif") |> rast()


# top
crs(statsgo_pmkind_top) <- crs(ex)
statsgo_pmkind_top <- statsgo_pmkind_top |>
  merge(gnatsgo_pmkind_top, first = TRUE) |>
  resample(ex, method = "mode") |>
  crop(ex, filename = file.path(fp_sta, "statsgo_pmkind_top_100m_crop.tif"))


# bot
crs(statsgo_pmkind_bot) <- crs(ex)
statsgo_pmkind_bot <- statsgo_pmkind_bot |>
  merge(gnatsgo_pmkind_bot, first = TRUE) |>
  resample(ex, method = "mode") |>
  crop(ex, filename = file.path(fp_sta, "statsgo_pmkind_bot_100m_crop.tif"))


# look table
table(dcon_pmkind_all$pmkind_top2) |> 
  as.data.frame() |> 
  cbind(id = 1:13) |> 
  write.csv(x = _, file.path(fp_sta, "statsgo_pmkind_lookup.csv"), row.names = FALSE)







