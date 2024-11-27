
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


# create tiles ----
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

file.remove(l_rat[!l_rat %in% tiles_nona$l_rat])



## STATSGO ----
fp_sta <- "D:/geodata/soils/STATSGO2/wss_gsmsoil_US_[2016-10-13]"
statsgo_sf <- read_sf(file.path(fp_sta, "US/spatial"), layer = "soilmu_a_us") |>   transform(MUKEY = as.integer(MUKEY)) |>
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

l_rat_sta <- list.files(path = fp_sta, pattern = "statsgo_[0-9]{1,3}.tif$")
tiles_nona <- readRDS(file.path(fp_nat, "tiles.rds"))
tiles_nona$sta <- gsub(paste0(fp_nat, "/"), "", tiles_nona$l_rat)
tiles_nona$sta <- gsub("gnatsgo", "statsgo", tiles_nona$sta)

file.remove(file.path(fp_sta, l_rat_sta[!l_rat_sta %in% tiles_nona$sta]))


# harmonize ----
fp_nat <- file.path("D:/geodata/soils/gNATSGO_CONUS_Oct2023")
gnatsgo_l <- list.files(path = fp_nat, pattern = "gnatsgo_[0-9]{1,3}.tif$", full = TRUE)
fp_sta <- "D:/geodata/soils/STATSGO2/wss_gsmsoil_US_[2016-10-13]"
statsgo_l <- list.files(path = fp_sta, pattern = "statsgo_[0-9]{1,3}.tif$", full = TRUE)


hmz_fun <- function(input, zonal, workers = 8) {
  
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

test <- hmz_fun(gnatsgo_l, statsgo_l, 6)
test <- do.call("rbind", test)
lookup <- test[
  ,.(statsgo_mukey = zonal, gnatsgo_mukey = input, N)][
    !is.na(gnatsgo_mukey) & statsgo_mukey > 0][
      , .(N = sum(N)), by = .(statsgo_mukey, gnatsgo_mukey)][
        order(statsgo_mukey, -N)] |>
  as.data.frame()

# saveRDS(lookup, file = file.path(fp_sta, "statsgo_gnatsgo_lookup.rds"))
sgo_lookup <- readRDS(file = file.path(fp_sta, "statsgo_gnatsgo_lookup.rds"))



# gNATSGO derivatives ----
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
  left_join(mu[c("mukey", "muacres", "areasymbol", "musym")], by = "mukey") |>
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
  }) |>
  filter(musym != "NOTCOM")


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
  filter(! is.na(pmkind2)) %>%
  filter(! duplicated(compname2)) |>
  select(compname2, pmkind2_mis = pmkind2, pmkind_top2_mis = pmkind_top2, pmkind_bot2_mis = pmkind_bot2)

co_pmkind2 <- co_pmkind |>
  left_join(mis, by = "compname2") |>
  mutate(
    pmkind2     = ifelse(!is.na(pmkind2),     pmkind2,     pmkind2_mis),
    pmkind_top2 = ifelse(!is.na(pmkind_top2), pmkind_top2, pmkind_top2_mis),
    pmkind_bot2 = ifelse(!is.na(pmkind_bot2), pmkind_bot2, pmkind_bot2_mis),
    
    pmkind2     = ifelse(is.na(pmkind2)     & grepl("rock outcrop", compname2),     "residuum", pmkind2),
    pmkind_top2 = ifelse(is.na(pmkind_top2) & grepl("rock outcrop", compname2), "residuum", pmkind_top2),
    pmkind_bot2 = ifelse(is.na(pmkind_bot2) & grepl("rock outcrop", compname2), "residuum", pmkind_bot2),
    
    pmkind2     = ifelse(is.na(pmkind2)     & grepl("alluvi|^fan ", landform),     "alluvium", pmkind2),
    pmkind_top2 = ifelse(is.na(pmkind_top2) & grepl("alluvi|^fan ", landform), "alluvium", pmkind_top2),
    pmkind_bot2 = ifelse(is.na(pmkind_bot2) & grepl("alluvi|^fan ", landform), "alluvium", pmkind_bot2),
    
    pmkind2_mis     = NULL,
    pmkind_top2_mis = NULL,
    pmkind_bot2_mis = NULL
  )


## dominant condition ----
### pmkind ----
dcon_pmkind <- co_pmkind2 %>%
  group_by(mukey, pmkind2) %>%
  summarize(
    acres = sum(coacres, na.rm = TRUE),
    pmkind_pct   = round(sum(comppct_r, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  arrange(mukey, - acres, - pmkind_pct) %>%
  filter(! is.na(pmkind2)) %>%
  filter(! duplicated(mukey)) %>%
  select(mukey, pmkind_pct, pmkind2)


dcon_pmkind_top <- co_pmkind2 %>%
  group_by(mukey, pmkind_top2) %>%
  summarize(
    acres = sum(coacres, na.rm = TRUE),
    pmkind_top_pct   = round(sum(comppct_r, na.rm = TRUE)),
  ) %>%
  ungroup() %>%
  arrange(mukey, - acres, - pmkind_top_pct) %>%
  filter(! is.na(pmkind_top2)) %>%
  filter(! duplicated(mukey)) %>%
  select(mukey, pmkind_top_pct, pmkind_top2)


dcon_pmkind_bot <- co_pmkind2 %>%
  group_by(mukey, pmkind_bot2) %>%
  summarize(
    acres = sum(coacres, na.rm = TRUE),
    pmkind_bot_pct   = round(sum(comppct_r, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  arrange(mukey, - acres, - pmkind_bot_pct) %>%
  filter(! is.na(pmkind_bot2)) %>%
  filter(! duplicated(mukey)) %>%
  select(mukey, pmkind_bot_pct, pmkind_bot2)

### drainage ----
lev <- c("Excessively drained", "Somewhat excessively drained", "Well drained", "Moderately well drained", "Somewhat poorly drained", "Poorly drained", "Very poorly drained", "Subaqueous")

dcon_drainagecl_poorly_pct <- co_pmkind2 %>%
  mutate(
    drainagecl = factor(drainagecl, levels = lev, ordered = TRUE),
    drainagecl_poorly = as.integer(drainagecl >= "Somewhat poorly drained")
    ) |>
  group_by(mukey) %>%
  summarize(
    acres = sum(coacres, na.rm = TRUE),
    drainagecl_poorly_pct = round(sum(comppct_r[drainagecl_poorly == 1], na.rm = TRUE)),
    water = sum(grepl("^Water$|Water | water", compname[which.max(comppct_r)]), na.rm = TRUE) > 0
  ) %>%
  ungroup() %>%
  # arrange(mukey, -drainagecl_poorly_pct) %>%
  # filter(! duplicated(mukey)) %>%
  mutate(
    drainagecl_poorly_pct = ifelse(is.na(drainagecl_poorly_pct) & water != TRUE, 0, drainagecl_poorly_pct),
    drainagecl_poorly_pct = ifelse(water == TRUE, NA, drainagecl_poorly_pct)
    ) %>%
  select(mukey, drainagecl_poorly_pct)


dcon_drainagecl <- co_pmkind2 %>%
  group_by(mukey, drainagecl) %>%
  summarize(
    acres = sum(coacres, na.rm = TRUE),
    drainagecl_pct = round(sum(comppct_r, na.rm = TRUE)),
    water = any(sum(grepl("^Water$|Water | water", compname[which.max(comppct_r)]), na.rm = TRUE))
  ) %>%
  ungroup() %>%
  arrange(mukey, - acres, - drainagecl_pct) %>%
  # filter(! is.na(drainagecl)) %>%
  filter(! duplicated(mukey)) %>%
  mutate(
    drainagecl = ifelse(is.na(drainagecl) & water != TRUE, "Well drained", drainagecl),
    drainagecl = factor(drainagecl, levels = lev, ordered = TRUE),
    drainagecl = ifelse(water == TRUE, NA, drainagecl),
    drainagecl = factor(drainagecl, levels = 1:8, labels = lev, ordered = TRUE)
  ) %>%
  select(mukey, drainagecl_pct, drainagecl)


### combine ----
dcon_all <- mu["mukey"] |>
  merge(dcon_pmkind,                by = "mukey", all.x = TRUE) |>
  merge(dcon_pmkind_top,            by = "mukey", all.x = TRUE) |>
  merge(dcon_pmkind_bot,            by = "mukey", all.x = TRUE) |>
  merge(dcon_drainagecl_poorly_pct, by = "mukey", all.x = TRUE) |>
  merge(dcon_drainagecl,            by = "mukey", all.x = TRUE) |>
  as.data.frame()
# saveRDS(dcon_all, file.path(fp_nat, "dcon_all.rds"))
dcon_all <- readRDS(file.path(fp_nat, "dcon_all.rds"))



# deratify ----
## gNATSGO ----
derat_fun <- function(l, dat, vars, workers = 8) {
  # l = list of tiles
  # dat = data frame with variables to be deratified
  # vars = variables within the data frame
  
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
        
        if (is.character(val[[var]])) {
          lev <- dat[[var]] |> as.factor() |> levels()
          val[[var]] <- val[[var]] |> factor(levels = lev) |> as.integer()
        }
        if (is.factor(val[[var]])) {
          val[[var]] <- val[[var]] |> as.integer()
        }
        values(r) <- val[[var]]
        writeRaster(r, file = gsub(".tif", paste0("_", var, ".tif"), x), overwrite = TRUE)
      }
    })})
}


fp_nat <- file.path("D:/geodata/soils/gNATSGO_CONUS_Oct2023")
tiles <- readRDS(file.path(fp_nat, "tiles.rds"))
l_rat <- tiles$l_rat[tiles$idx]


### deratify ----
vars <- c("pmkind2", "pmkind_top2", "pmkind_bot2", "drainagecl_poorly_pct", "drainagecl")
derat_fun(l_rat, dcon_all, vars[c(1, 3)])


### mosaic ----
lapply(vars[c(1, 3)], function(x) {
  l2 <- list.files(
    path = fp_nat, 
    pattern = paste0("gnatsgo_[0-9]{1,3}_", x, ".tif$"), full.names = TRUE
    )
  n <- dcon_all[[x]] |> as.factor() |> as.integer() |> unique() |> length()
  plot(vrt(l2), col = viridis::turbo(n)[sample(1:n)], main = x)
  plot(tiles_nona, add = TRUE)
  test <- sprc(l2)
  fn <- file.path(fp_nat, paste0("gnatsgo_Oct22_", x, ".tif"))
  merge(test, filename = fn, datatype = "INT1U", overwrite = TRUE)
})



## STATSGO ----
fp_sta <- "D:/geodata/soils/STATSGO2/wss_gsmsoil_US_[2016-10-13]"
statsgo_l <- list.files(path = fp_sta, pattern = "statsgo_[0-9]{1,3}.tif$", full = TRUE)

lookup <- readRDS(file = file.path(fp_sta, "statsgo_gnatsgo_lookup.rds")) |>
  merge(dcon_all, by.x = "gnatsgo_mukey", by.y = "mukey", all.x = TRUE)
idx <- sapply(lookup, is.character)
lookup[idx] <- lapply(lookup[idx], as.factor)

### aggregate gnatsgo vs statsgo ----
vars <- c("pmkind2", "pmkind_top2", "pmkind_bot2", "drainagecl_poorly_pct", "drainagecl")
dtype <- c("d", "d", "d", "c", "d")
vd <- data.frame(vars, dtype)

vd <- split(vd, vd$"vars")
vd <- lapply(vd, function(x) {
  
  vars2 <- c("statsgo_mukey", x$vars, "N")
  vars3 <- vars2[1:2]
  
  if (x$dtype == "d") {
    y <- data.table::as.data.table(lookup)[
      , ..vars2][
        , .(N = sum(N, na.rm = TRUE)), by = vars3][
          order(statsgo_mukey, -N)][
           !duplicated(statsgo_mukey)][
             order(statsgo_mukey)
            ] |>
      as.data.frame()
  }
  if (x$dtype == "c") {
    y <- data.table::as.data.table(lookup)[
      , ..vars2][
        , lapply(.SD, mean, na.rm=TRUE), by = statsgo_mukey][
          order(statsgo_mukey)
        ] |>
      as.data.frame()
  }
  return(y)
})
# check for nonmatches
lapply(vd, function(x) x$statsgo_mukey) |> 
  do.call("cbind", args = _) |>  
  apply(X = _, 1, function(x) length(unique(x)) == 1) |> 
  all()
vd <- do.call("cbind", vd)
nm <- sapply(names(vd), function(x) strsplit(x, "\\.")[[1]][2]) |> unname()
idx <- !duplicated(nm)
vd <- vd[idx]
names(vd) <- nm[idx]
statsgo_lookup <- vd
names(statsgo_lookup)[1] <- "mukey"



### deratify ----
vars <- c("pmkind2", "pmkind_top2", "pmkind_bot2", "drainagecl_poorly_pct", "drainagecl")
derat_fun(statsgo_l, statsgo_lookup, vars)


### mosaic ----
lapply(vars, function(x) {
  l2 <- list.files(
    path = fp_sta, 
    pattern = paste0("statsgo_[0-9]{1,3}_", x, ".tif$"), full.names = TRUE
  )
  n <- dcon_all[[x]] |> as.factor() |> as.integer() |> unique() |> length()
  plot(vrt(l2), col = viridis::turbo(n)[sample(1:n)], main = x)
  plot(tiles_nona, add = TRUE)
  test <- sprc(l2)
  fn <- file.path(fp_sta, paste0("statsgo_", x, ".tif"))
  merge(test, filename = fn, datatype = "INT1U", overwrite = TRUE)
})



# merge, resample, crop ----
ex <- rast("D:/geodata/project_data/ADIUCL5.tif")

gnatsgo <- c(
  nat_pmkind_top = "gnatsgo_Oct22_pmkind_top2.tif", 
  nat_pmkind_bot = "gnatsgo_Oct22_pmkind_bot2.tif", 
  nat_drngcl_prp = "gnatsgo_Oct22_drainagecl_poorly_pct.tif",
  nat_drngcl     = "gnatsgo_Oct22_drainagecl.tif"
  )
gnatsgo_l <- lapply(
  gnatsgo, function(x) file.path(
    fp_nat, x) |> rast())
names(gnatsgo_l) <- names(gnatsgo)

statsgo <- c(
  sta_pmkind_top = "statsgo_pmkind_top2.tif", 
  sta_pmkind_bot = "statsgo_pmkind_bot2.tif", 
  sta_drngcl_prp = "statsgo_drainagecl_poorly_pct.tif",
  sta_drngcl     = "statsgo_drainagecl.tif"
)
statsgo_l <- lapply(
  statsgo, function(x) file.path(
    fp_sta, x) |> rast())
names(statsgo_l) <- names(statsgo)


mrc_fun <- function(x, y, fp, fn, method) {
  
  cat("processing", varnames(x), "using", method, as.character(Sys.time()), "\n")
  
  crs(x) <- crs(ex)
  crs(y) <- crs(ex)
  
  x <- x |>
    merge(y, first = TRUE) |>
    resample(ex, method = method) |>
    crop(ex, filename = file.path(
      fp,
      fn
      ),
      overwrite = TRUE
    )
  
  return(x)
}

fn <- gsub(".tif", "_100m_crop.tif", statsgo)
method <- rep("mode", 3)
method <- c(method[1:2], "average", method[3])

mapply(mrc_fun, statsgo_l, gnatsgo_l, fp_sta, fn, method)
fn <- gsub(".tif", "_100m_crop.tif", gnatsgo)
mapply(mrc_fun, gnatsgo_l, statsgo_l, fp_nat, fn, method)



# look table
table(dcon_all$pmkind_top2) |> 
  as.data.frame() |> 
  cbind(id = 1:13) |> 
  write.csv(x = _, file.path(fp_sta, "statsgo_gnatsgo_pmkind_lookup.csv"), row.names = FALSE)

table(dcon_all$drainagecl) |> 
  as.data.frame() |> 
  cbind(id = 1:8) |> 
  write.csv(x = _, file.path(fp_sta, "statsgo_gnatsgo_drainagecl_lookup.csv"), row.names = FALSE)






