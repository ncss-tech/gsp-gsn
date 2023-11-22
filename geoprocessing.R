
library(sf)
library(terra)


# covariates ----
fp_cov100 <- file.path("./cov100/covariates")
covs <- read.csv(file.path(fp_cov100, "SoilGrids_USA48_Covs100m.csv"))
bio <- subset(covs, SERIES_NAME == "Bioclimatic layers")
bio_r <- rast(file.path(fp_cov100, paste0(bio$WORLDGRIDS_CODE, ".tif")))


lf <- list.files(fp_cov100, pattern = ".tif$", full.names = TRUE)
lf2 <- list.files(file.path(fp_cov100, "gSSURGO_filled"), ".tif$", full.names = TRUE)
lf <- c(lf, lf2)
lf_ext <- sapply(lf, function(x) rast(x) |> ext() |> _[1:4]) |> t()
lf_crs <- sapply(lf, function(x) rast(x) |> crs(proj = TRUE))
lf <- cbind(fp = lf, ext = lf_ext, crs = lf_crs) |> data.frame(row.names = NULL)
lf$fn <- gsub("\\./cov100/covariates/|\\.tif$", "", lf$fp)


idx <- grepl("NAD83", lf$crs)
rs1 <- rast(lf$fp[!idx])
rs2 <- rast(lf$fp[ idx])
rs3 <- project(rs2, rs1, align = TRUE)
rs  <- c(rs1, rs3)
rs  <- rast(lf$fp)

bb_tiles  <- ext(rs) |> 
  vect(crs = st_crs(5070)$proj4string) |> 
  st_as_sf() |> 
  st_make_grid(cellsize = 4000 * 100) |>
  st_as_sf() # |>
# st_transform(crs = 4326)
bb_tiles$tile <- 1:nrow(bb_tiles)
# bb_tiles <- ext(rs1) |> 
#   as.matrix() |> 
#   t() |> 
#   as.data.frame() |> 
#   st_as_sf(coords = c("V1", "V2"), crs = 5070) |> 
#   st_make_grid() |>
#   st_transform(crs = 4326)

write_sf(bb_tiles, dsn = "solus.sqlite", layer = "bb")


# r <- rast(ncols = 13, nrows = 8, crs = crs(rs), extent = ext(rs))
# values(r) <- 1:ncell(r)
# r2 <- r[1,]

fun <- function(x) {
  r <- rast(x)
  makeTiles(r, 4000, filename = file.path(fp_cov100, "tiles", paste0(names(r), "_.tif")), gdal="COMPRESS=NONE")
}

# fun2 <- function(x) {
#   r <- rast(x)
#   
#   r2 <- rast()
#   ext(r2) <- c(1460000, 1860000, 1319000, 1719000)
#   crop(r, r2, filename = file.path(fp_cov100, "tiles", paste0(names(r), "_50.tif")), gdal="COMPRESS=NONE")
# }

future::plan(future::multisession, workers = 120)
future.apply::future_lapply(lf$fp[idx], fun)



# copy to VM ----
lf_tiles <- list.files(file.path(fp_cov100, "tiles")) |>
  cbind(fp = _) |>
  as.data.frame()
lf_tiles$tile <- strsplit(lf_tiles$fp, split = "_|.tif") |>
  sapply(function(x) x[length(x)]) |>
  as.integer()
lf_tiles <- poorman::arrange(lf_tiles, tile)

idx <- grepl("mlra52|PMTGSS7_f|DRNGSS7_f", lf_tiles)
# file.copy(file.path(fp_cov100, "tiles", lf_tiles$fp[idx]), "./tiles")



# capture tile extents ----
tiles <- list.files("./tiles", pattern = ".tif$", full.names = TRUE) |>
  cbind(fp = _) |>
  as.data.frame()
tiles$tile <- tiles$fp |>
  strsplit(split = "_|.tif") |>
  sapply(function(x) x[length(x)]) |>
  as.integer()
tiles$fn <- tiles$fp |>
  strsplit(split= "\\./tiles/|_[0-9]{1,3}.tif$") |>
  sapply(function(x)  x[length(x)])
tiles <- poorman::arrange(tiles, tile)


tiles_bb <- tiles |>
  subset(! duplicated(tile)) |>
  _[["fp"]] |>
  lapply(function(x) {
    tmp <- rast(x) |>
      ext() |>
      st_bbox() |>
      st_as_sfc()
    return(tmp)
  }) |> 
  do.call("rbind", args = _) |>
  st_as_sfc() |>
  st_as_sf() |>
  st_set_crs(value = 5070) |>
  st_transform(crs = 4326) |>
  transform(tile = unique(tiles$tile))

mapview::mapview(tiles_bb)



# soil data ----
h <- readRDS("./data/gsn_seg_df.rds")
s <- readRDS("./data/scd_site_sf.rds")
st_crs(s) <- 4326

idx <- st_intersects(s, tiles_bb) |> 
  as.integer()
# sapply(function(x) if (length(x) > 0) x else NA)
s$tile <- idx
s2 <- s |> 
  st_transform(crs = 5070)



# extract ----
tiles_l <- tiles |>
  subset(tile %in% s2$tile)
tiles_l <- split(tiles_l, tiles_l$tile)

s2_l <- split(s2, s2$tile)


extract_fun <- function(tile, s) {
  r <- rast(tile$fp)
  extract(r, vect(s))
}

system.time(
  test2 <- extract_fun(tiles_l[[1]], s2_l[[1]])
)
cbind(s2_l[[1]], test2) |> st_transform(crs = 4326) |> mapview::mapview()


test2 <- mapply(extract_fun, tiles_l, s2_l, SIMPLIFY = FALSE)
test3 <- mapply(cbind, s2_l, test2, SIMPLIFY = FALSE) |>
  do.call("rbind", args = _)
names(test3)[names(test3) == "Layer_1"] <- "NLCD116"

vars <- c("PVEGKT6", "DRNGSS7", "PMTGSS7", "NLCD116", "LNDCOV6", "COUNTY6", "GESUSG6")
vars <- names(test3)[names(test3) %in% vars]

test3[vars] <- lapply(vars, function(x) {
  val <- rast(lf$fp[lf$fn == x]) |> unique() |> _[[1]]
  x2  <- test3[[x]] |> factor(levels = val)
  return(x2)
})


# saveRDS(test3, file = "./data/scd_site_extract_sf.rds")
test3 <- readRDS(file = "./data/scd_site_extract_sf.rds")








