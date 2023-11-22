
library(sf)
library(ggplot2)
library(dplyr)
library(data.table)


# load snapshots ----
fp <- "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data"
scd_l <- readRDS(file.path(fp, "ncss-scd_sda_20231102.rds"))
# scd_l2 <- readRDS(file = file.path(fp, "ncss_labdata.rds"))
f <- readRDS(file.path(fp, "fetchNASIS_spc_20230926.rds"))



# metadata
dm <- dm::dm_from_src(soilDB:::.openNASISchannel(), learn_keys = TRUE)
dm$system |> 
  dm::collect() |> 
  subset(sysver == "Lab SDA Data Mart 1.0") |> 
  as.data.frame()
md <- dm$attribute |> dm::collect() |> 
  subset(sysiidref == 41045) |>
  merge(
    dm$uom |> dm::collect() |> subset(select = c(uomiid, uomsym)), 
    by.x = "uomiidref", by.y = "uomiid", all.x = TRUE, sort = FALSE
    )
md <- md[c(ncol(md), 1:(ncol(md) - 1))]



# convert dates ----
s <- f@site
anyDuplicated(paste(s$peiid, s$obsdate))
s$peiid <- as.integer(s$peiid)


scd_nasis <- merge(
  scd_l$combine_nasis_ncss, 
  s[c("peiid", "obsdate", "obsdatekind", "geocoordsource")],
  by.x = "pedoniid", by.y = "peiid",
  all.x = TRUE,
  sort = FALSE
)



scd_nasis <- base::within(scd_nasis, {
  samp_classdate2 = strptime(samp_classdate, "%m/%e/%Y %H:%M:%S %p")
  corr_classdate2 = strptime(corr_classdate, "%m/%e/%Y %H:%M:%S %p")
  SSL_classdate2  = strptime(SSL_classdate,  "%m/%e/%Y %H:%M:%S %p")
  site_obsdate2   = strptime(site_obsdate,   "%m/%e/%Y %H:%M:%S %p")
  obsdate = as.Date(obsdate, "%Y-%m-%d")
  
  samp_year    = format(samp_classdate2, "%Y") |> as.integer()
  corr_year    = format(corr_classdate2, "%Y") |> as.integer()
  SSL_year     = format(SSL_classdate2,  "%Y") |> as.integer()
  site_obsyear = format(site_obsdate2, "%Y") |> as.integer()
  obsyear      = format(obsdate,         "%Y") |> as.integer()
  year         = apply(cbind(samp_year, corr_year, SSL_year, site_obsyear, obsyear), 1, min, na.rm = TRUE)
  year         = ifelse(is.infinite(year), as.integer(substr(upedonid, 2, 5)), year)
  year         = ifelse(is.na(year),       as.integer(substr(usiteid,  2, 5)), year)
  year         = ifelse(is.na(year),       as.integer(substr(pedlabsampnum, 1, 2)) + 2000L, year)
  # year      = ifelse(is.na(obsyear), year, obsyear)
})


# saveRDS(scd_nasis, file = "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data/scd_nasis.rds")
scd_nasis <- readRDS(file = "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data/scd_nasis.rds")



# coordinates ----
## rename columns ----
scd_s <- scd_l$site
nm <- names(scd_s) 
idx <- grepl("longitude", nm)
names(scd_s)[idx] <- sub("longitude", "lon", nm[idx])
idx <- grepl("latitude", nm)
names(scd_s)[idx] <- sub("latitude",  "lat", nm[idx])
names(scd_s)[nm == "horizontal_datum_name"] <- "datum"

# nm <- names(scd_s) 
# var <- "std_decimal_degrees"
# idx <- which(grepl(var, nm))
# names(scd_s)[idx] <- gsub(var, "dd", nm[idx])

table(dms = complete.cases(scd_s[3:11]), dd = complete.cases(scd_s[12:13]))


## fix datum ----
scd_s <- within(scd_s, {
  ellps = NA
  ellps[datum == "old hawaiian"] = "GRS80"  # EPSG = 6135
  ellps[datum == "Guam1963"]     = "clrk66" # EPSG = 4675
  datum[!is.na(ellps)]           = NA
  datum = sub("North American Datum of 1927", "NAD27",  datum)
  datum = sub("North American Datum of 1983", "NAD83",  datum)
  datum = sub("World Geodetic System 1984",   "WGS84",  datum)
  datum[datum == "NAD82"] <- "NAD83"
  datum[datum == "WGS85"] <- "WGS84"
  datum[datum == "WGS86"] <- "WGS84"
  
  lon_direction = tolower(lon_direction)
  lat_direction = tolower(lat_direction)
  
  lon_degrees = abs(lon_degrees)
  lon_minutes = abs(lon_minutes)
  lon_seconds = abs(lon_seconds)
  
  lat_degrees = abs(lat_degrees)
  lat_minutes = abs(lat_minutes)
  lat_seconds = abs(lat_seconds)
  
  x = (lon_degrees + lon_minutes / 60 + lon_seconds / 3600) * ifelse(lon_direction == "west",  -1, 1)
  y = (lat_degrees + lat_minutes / 60 + lat_seconds / 3600) * ifelse(lat_direction == "south", -1, 1)
})
table(datum = scd_s$datum, ellps =  scd_s$ellps, useNA = "always")



# estimate datum for missing
s2 <- aggregate(year ~ site_key, data = scd_nasis, min, na.rm = TRUE)
scd_s <- merge(scd_s, s2, by = "site_key", all.x = TRUE, sort = FALSE)

test <- scd_s |>
  transform(xy = complete.cases(x, y),
            year = as.integer(year)
            ) |>
  aggregate(xy ~ year + datum, data = _, sum, na.rm = TRUE)

subset(test, year %in% 1920:2023) |>
  ggplot(aes(x = year, y = xy, col = datum)) +
  geom_line(linewidth = 1, alpha = 0.7)


scd_s <- within(scd_s, {
  idx = (is.na(datum) | datum == "") & is.na(ellps)
  datum[idx  & year <= 2012]  = "NAD83"
  datum[idx & year  >   2012] = "WGS84"
  datum[is.na(datum) & is.na(ellps)] = "NAD83"
  idx = NULL
})
table(datum = scd_s$datum, ellps =  scd_s$ellps, useNA = "always")



## convert dms to dd ----
s2 <- {
  split(scd_s, ~ paste(datum, ellps), drop = TRUE) ->.;
  lapply(., function(x) {
    x2 <- subset(x, complete.cases(x, y))
    cat("converting", x2$datum[1], x2$ellps[1], "\n")
    proj4 = paste0(
      "+proj=longlat ", 
      ifelse(!is.na(x2$ellps[1]), "+ellps=",   "+datum="), 
      ifelse(!is.na(x2$ellps[1]), x2$ellps[1], x2$datum[1]),
      " +no_defs"
      )
    x2 <- st_as_sf(
      x2,
      coords = c("x", "y"),
      crs = proj4
    )
    x2 <- st_transform(x2, crs = 4326)
    x2 <- cbind(x2, st_coordinates(x2))
    st_drop_geometry(x2)
  }) ->.;
  do.call("rbind", .)
}

vars <- c("site_key", "X", "Y")
scd_site <- merge(scd_s, s2[vars], by = "site_key", all.x = TRUE, sort = FALSE)
scd_site <- within(scd_site, {
  lon_dd = ifelse(is.na(lon_std_decimal_degrees), X, lon_std_decimal_degrees)
  lat_dd = ifelse(is.na(lat_std_decimal_degrees), Y, lat_std_decimal_degrees) 
})


## coordinate precision ----
coord_precision <- function(df, digits) {
  test <- round(df, digits) == df
  apply(test, 1, all)
}

coord_prec2 <- function(x, y) {
  
  x_n1 <- regexpr("\\.", x)
  x_n2 <- nchar(x)
  
  y_n1 <- regexpr("\\.", y)
  y_n2 <- nchar(y)
  
  test <- cbind(
    x_precision = substr(x, x_n1 + 1, x_n2) |> nchar(), 
    y_precision = substr(y, y_n1 + 1, y_n2) |> nchar()
  )
}

# 5 digits = ~ 1m; 4 digits = ~ 11m; 3 digits = ~ 111m
test <- coord_prec2(scd_site$lon_dd, scd_site$lat_dd)
{apply(test, 1, max) >= 3} |> summary()

scd_site <- cbind(scd_site, test)



## project dd ----
scd_sf <- subset(scd_site, complete.cases(lon_dd, lat_dd))
# scd_sf <- subset(scd_site, complete.cases(lon_std_decimal_degrees, lat_std_decimal_degrees))
scd_sf <- st_as_sf(
  scd_sf,
  coords = c("lon_dd", "lat_dd"),
  # coords = c("longitude_std_decimal_degrees", "latitude_std_decimal_degrees"),
  crs = 4326
)
st_bbox(scd_sf)  



# intersect with the GADM reference
world_bndy <- geodata::world(path = getwd(), resolution = 1) |> 
  st_as_sf() |>
  st_make_valid() |>
  st_transform(crs = 4326)

idx <- st_intersects(scd_sf, world_bndy) |> 
  lapply(function(x) x[1]) |> 
  unlist()
scd_sf <- scd_sf |> cbind(st_drop_geometry(world_bndy[idx, 1:2]))

table(USA = scd_sf$GID_0 == "USA", Count = !is.na(scd_sf$site_key), useNA = "always")

scd_sf[names(scd_sf) %in% c("X", "Y", "year")] <- NULL

# saveRDS(scd_sf, file = "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data/scd_site_sf.rds")
scd_sf <- readRDS(file = "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data/scd_site_sf.rds")



# extract variables ----

pat_phys <- c("^sand|^silt|^clay|^texture|^particle_size|frag|bulk_density_third")
nms_phys  <- names(scd_l$physical_properties)
vars_phys <- nms_phys[
  grepl(pat_phys, nms_phys) 
  & !grepl("disp", nms_phys)
  ]

pat_chem <- c("total_nitrogen|^phosphorus_|^new_zealand_phos|^potassium_|^cec_nh4|ph_h2o|ph_cacl2|ph_kcl|total_carbon|organic_carbon|caco3_lt_2")
nms_chem  <- names(scd_l$chemical_properties)
vars_chem <- nms_chem[
  grepl(pat_chem, nms_chem) 
  & !grepl("disp", nms_chem)
  ]

gsn_df <- scd_l$layer |>
  merge(scd_l$physical_properties[c(1:3, match(vars_phys, nms_phys))], by = "labsampnum", all.x = TRUE, sort = FALSE) |>
  merge(scd_l$chemical_properties[c(1:3, match(vars_chem, nms_chem))], 
                                   by = "labsampnum", all.x = TRUE, sort = FALSE)


## uncode method codes ----
idx <- grep("method", names(gsn_df))
test <- lapply(idx, function(x) {
  factor(gsn_df[[x]], 
         levels = scd_l$method_code$proced_code, 
         labels = scd_l$method_code$proced_name) |> 
    droplevels()
}) |> 
  as.data.frame()
names(test) <- paste0(names(gsn_df)[idx], "_desc")
gsn_df <- cbind(gsn_df, test)



## summarize ----
test <- scd_l$chemical_properties %>%
  filter(!is.na(total_nitrogen_ncs)) %>%
  inner_join(scd_l$layer, by = "labsampnum") %>%
  inner_join(scd_nasis, by = "pedon_key")
test2 <- test %>% 
  group_by(pedon_key) %>% 
  summarize(n_N = sum(!is.na(total_nitrogen_ncs)) > 0) %>%
  inner_join(scd_nasis, by = "pedon_key") %>%
  select(pedon_key, n_N, samp_year, corr_year, SSL_year, year)


# tally the number of pedon per decade with 1 or more measurements for each soil property
test <- as.data.table(gsn_df)
vars <- c(vars_phys, vars_chem)
vars <- vars[!grepl("method", vars)]
test <- test[, lapply(.SD, function(x) any(!is.na(x))), .SDcols = vars, by = c("site_key", "pedon_key")]
vars <- c("site_key", "pedon_key", "year")
test <- merge(scd_nasis[vars], test, by = vars[1:2], all.x = TRUE, sort = FALSE)
test$decade <- substr(test$year, 1, 3) |> as.integer() * 10L
test2 <- aggregate(. ~ decade, data = test[- c(1:3)], sum, na.rm = TRUE)|>
  t() |>
  as.data.frame()
names(test2) <- as.character(test2[1, ])
test2 <- test2[-1, ]
test2 <- cbind(var = row.names(test2), test2)
row.names(test2) <- NULL
View(test2)



## nitrogen ----
nm  <- names(gsn_df)
idx <- which(grepl("nitrogen", nm) & !grepl("method", nm))
summary(gsn_df[idx])
aggregate(total_nitrogen_ncs ~ total_nitrogen_ncs_method_desc, data = gsn_df, function(x) {
  summary(x) |> round(2) |> t()
  })

table(gsn_df$total_nitrogen_ncs_method_desc, gsn_df$total_nitrogen_ncs_method)


### pedotransfer function from Tomer 2017 ----
gsn_df <- within(gsn_df, {
  N_pt = ifelse(
    !grepl("Kjeldahl|Unknown", total_nitrogen_ncs_method_desc), 
    0.032 + 1.008 * total_nitrogen_ncs, 
    total_nitrogen_ncs
    )
  })

# correlation
with(gsn_df, cor(total_nitrogen_ncs, N_pt, use = "complete.obs"))


# scatter plot
with(gsn_df, plot(total_nitrogen_ncs, N_pt))
abline(0, 1)



## phosphorus ----
nm  <- names(gsn_df)
idx <- which(grepl("phosphorus", nm) & !grepl("method", nm))
summary(gsn_df[idx])


### subset non-missing data ----
gsn_df[idx] > 0 ->.;
# .[rowSums(., na.rm = TRUE) > 0, ] |>
rowSums(., na.rm = TRUE) |> 
  table()
gsn_df_P <- gsn_df[rowSums(., na.rm = TRUE) > 0, c(1:3, idx)]

gsn_df_P_layer <- merge(
  scd_l$layer, 
  gsn_df_P, 
  by = "labsampnum", 
  all.y = TRUE, 
  sort = FALSE
  )
gsn_df_P_layer <- as.data.table(gsn_df_P_layer)
vars <- nm[idx]
gsn_df_P_agg <- gsn_df_P_layer[, lapply(.SD, function(x) any(!is.na(x))), .SDcols = vars, by = "site_key"]
# gsn_df_P_agg <- aggregate(cbind(phosphorus_bray1, phosphorus_olsen, phosphorus_mehlich_3, phosphorus_mehlich3_extractable) ~ as.character(site_key), data = gsn_df_P_agg, FUN = function(x) any(x > 0), na.action = na.omit)

scd_sf_P <- merge(scd_sf, gsn_df_P_agg, by = "site_key", all.x = TRUE, sort = FALSE)
mapview::mapview(scd_sf_P)


# summarize missing
summary(gsn_df_P)


# Bray 1 overlapping missing
gsn_df_P |>
  group_by(bray1 = !is.na(phosphorus_bray1)) |>
  summarize(
    meh3  = sum(!is.na(phosphorus_mehlich_3), na.rm = TRUE),
    meh3x = sum(!is.na(phosphorus_mehlich3_extractable), na.rm = TRUE)
    )

gsn_df_P |>
  filter(is.na(phosphorus_bray1)) |>
  group_by(meh3 = !is.na(phosphorus_mehlich_3)) |>
  summarize(
    meh3x = sum(!is.na(phosphorus_mehlich3_extractable), na.rm = TRUE)
  )


### pedotransfer function for missing Bray 1 ----
table(gsn_df_P$phosphorus_mehlich_3 > 1000)
filter(gsn_df_P, phosphorus_mehlich_3 < 1000) %>%
  ggplot(aes(x = phosphorus_bray1, y = phosphorus_mehlich_3)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline() +
  coord_fixed()

table(gsn_df_P$phosphorus_mehlich_3 > 500)
filter(gsn_df_P, phosphorus_mehlich_3 < 600) %>%
  ggplot(aes(x = phosphorus_mehlich3_extractable, y = phosphorus_mehlich_3)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline() +
  coord_fixed()

meh_bray_lm <- lm(phosphorus_mehlich_3 ~ phosphorus_bray1, data = gsn_df)
meh_meh_lm  <- lm(phosphorus_mehlich_3 ~ phosphorus_mehlich3_extractable, data = gsn_df)
summary(meh_bray_lm)
summary(meh_meh_lm)

gsn_df <- within(gsn_df, {
  P_pt = phosphorus_mehlich_3
  P_pt = ifelse(is.na(P_pt), predict(meh_bray_lm, data.frame(phosphorus_bray1)),                P_pt)
  P_pt = ifelse(is.na(P_pt), predict(meh_meh_lm,  data.frame(phosphorus_mehlich3_extractable)), P_pt)
  })


# saveRDS(gsn_df, file = file.path(fp, "gsn_df.rds"))
gsn_df <- readRDS(file = file.path(fp, "gsn_df.rds"))
# load(file = file.path(fp, "gsn_df.RData"))



# remove duplicates ----
scd_nasis <- readRDS(file = "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data/scd_nasis.rds")
scd_sf <- readRDS(file = "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data/scd_site_sf.rds")


scd_sf <- cbind(
  st_drop_geometry(scd_sf), 
  st_coordinates(scd_sf)
  ) |>
  merge(scd_nasis, by = "site_key", all.y = TRUE, sort = FALSE) |>
  within({
    coords = paste(round(X, 5), round(Y, 5))
    dups   = coords %in% coords[duplicated(coords)]
  }) |>
  # why does LDM have sites with no pedlabsampnum?
  subset(!is.na(pedlabsampnum))
  

sum(scd_sf$dups)


# test2 <- aggregate(site_key ~ coords, data = scd_sf, length)
# names(test2)[2] <- "n_coords"
test <- as.data.table(scd_sf)[, 
  .(n_coords = .N,
    n_year = names(sort(table(year), decreasing = TRUE))[1]
    ),
  by = "coords"
  ]
table(test$n_coords)
sum(test$n_coords)


scd_sf <- merge(scd_sf, test, by = "coords", all.x = TRUE, sort = FALSE) |>
  subset(complete.cases(X, Y))


scd_sf <- st_as_sf(
  scd_sf, 
  coords = c("X", "Y"), 
  crs = 4326
)


scd_sf$valid_xy <- with(
  scd_sf, 
  n_coords <= 5 
  & (x_precision >= 4 | y_precision >= 4)
  # & !(n_coords > 1 & geocoordsource %in% c("estimated from other source", "unknown") & (x_precision >= 6 | y_precision >= 6))
  & !(n_coords > 1 & year != n_year)
  # & in the same transect
  )
table(scd_sf$valid_xy)
table(substr(scd_sf$year, 1, 3) |> as.integer() * 10, scd_sf$valid_xy)
table(as.integer(substr(scd_sf$year, 1, 3)) * 10, scd_sf$n_coords)
View(scd_sf[idx, c("coords", "user_site_id", "pedlabsampnum", "n_coords", "year", "valid_xy")])


dsn <- file.path(fp, "scd_site_nasis_sf.sqlite")
file.remove(dsn)
write_sf(scd_sf, dsn)
# saveRDS(scd_sf, file = file.path(fp, "scd_site_nasis_sf.rds"))
scd_sf <- readRDS(file = "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data/scd_site_nasis_sf.rds")



vars <- c("layer", "physical_properties", "chemical_properties", "")
sapply(scd_l[vars], function(x) sum(duplicated(x$labsampnum), na.rm = TRUE))

aggregate(siteiid ~ site_key, data = scd_nasis, length) |> summary()
aggregate(pedoniid ~ pedon_key, data = scd_nasis, length) |> summary()
aggregate(pedon_key ~ site_key, data = scd_nasis, length) |> summary()
aggregate(pedlabsampnum ~ pedon_key, data = scd_nasis, length) |> summary()



# check horizon depths ----
source("C:/workspace2/github/ncss-tech/gsp-bs/validate.R")
test <- validate_depths(gsn_df, id = "pedon_key", top = "hzn_top", bot = "hzn_bot")
# test <- aqp::checkHzDepthLogic(gsn_df, hzdepths = c("hzn_top", "hzn_bot"), idname = "pedon_key", byhz = TRUE)

gsn_df2 <- gsn_df


# bad O horizons ----
aggregate(total_carbon_ncs ~ grepl("O", hzn_desgn), data = gsn_df2, quantile, probs = seq(0, 1, 0.1))

idx  <- with(gsn_df2, hzn_top > hzn_bot & (grepl("O", hzn_desgn) | (is.na(hzn_desgn) & total_carbon_ncs > 10)))
summary(idx)

vars <- c("hzn_top", "hzn_bot")
gsn_df2[vars] <- gsn_df[vars] * ifelse(idx, -1, 1)
summary(gsn_df2[vars])                       
  
  
idx_bad_o <- gsn_df2[idx, ]$pedon_key |> unique()

h2_bad_o_fixed <- gsn_df2 %>%
  filter(pedon_key %in% idx_bad_o) |>
  arrange(pedon_key, hzn_top) |>
  mutate(thk = abs(hzn_bot - hzn_top)) |>
  group_by(pedon_key) |>
  mutate(hzn_top = c(0, cumsum(thk)[-length(thk)]),
         hzn_bot = cumsum(thk)
  ) |>
  ungroup() |>
  mutate(thk = NULL)

h2_good_o <- gsn_df2[! gsn_df2$pedon_key %in% idx_bad_o, ]

gsn_df2 <- rbind(h2_good_o, h2_bad_o_fixed)

test2 <- validate_depths(gsn_df2, id = "pedon_key", top = "hzn_top", bot = "hzn_bot")
# test2 <- aqp::checkHzDepthLogic(gsn_df2, hzdepths = c("hzn_top", "hzn_bot"), idname = "pedon_key", byhz = TRUE)

summary(test)
summary(test2)


# # bad bottom depths
# h2_bad_bot_dep_fixed <- gsn_df2 |>
#   group_by(pedon_key) |>
#   filter(
#     (is.na(hzn_bot) | max(hzn_top) == max(hzn_bot)) 
#     & hzn_top == max(hzn_top)) |> 
#   ungroup() |>
#   mutate(hzn_bot = hzn_top + 1)
# 
# h2_good_bot_dep <- gsn_df2[! gsn_df$layer_key %in% h2_bad_bot_dep_fixed$layer_key, ]
# 
# gsn_df3 <- rbind(h2_good_bot_dep, h2_bad_bot_dep_fixed)
# 
# test3 <- aqp::checkHzDepthLogic(gsn_df3, hzdepths = c("hzn_top", "hzn_bot"), idname = "pedon_key", byhz = TRUE)

  
gsn_df3 <- cbind(gsn_df2, test2[1:4 * -1])
  
# save(gsn_df3, file = file.path(fp, "gsn_df3.RData"))
load(file = file.path(fp, "gsn_df3.RData"))
gsn_df3 <- gsn_df


# check units ----
gt0 <- function(x) x >= 0
bt0and100 <- function(x) x >= 0 & x <= 100
sum100 <- function(df) {
  as.integer(round(rowSums(df, na.rm = TRUE))) == 100L
  }

vars100 <- c("clay_total", "sand_total", "silt_total", "total_carbon_ncs")
vars <- c(vars100, "N_pt", "P_pt", "potassium_mehlich3_extractable", "cec_nh4_ph_7", "ph_h2o", "bulk_density_third_bar")


summary(gsn_df3[vars])

# gt0
test_gt0 <- sapply(gsn_df3[vars], gt0)
summary(test_gt0)
View(cbind(gsn_df3["pedon_key"], test_gt0))

# bt0and100
test_bt0and100 <- sapply(gsn_df3[vars100], bt0and100)
summary(test_bt0and100)
test_bt0and100 <- sapply(gsn_df3[vars100], bt0and100)
View(cbind(gsn_df3["pedon_key"], test_gt0))

#sum100
test_sum100 <- sum100(gsn_df3[vars100[1:3]])
test0 <- gsn_df3[!test_sum100 & complete.cases(gsn_df3[vars100[1:3]]), ]
summary(test_sum100)

test <- gsn_df3[c("pedon_key", var100)] |>
  within({
  ssc_total = rowSums(cbind(clay_total, sand_total, silt_total))
  check = as.integer(round(ssc_total)) == 100L
}) |>
  subset(complete.cases(cbind(clay_total, sand_total, silt_total)))



# fix units ----
# C, N, & P
vars <- c("total_carbon_ncs", "N_pt", "P_pt")
gsn_df3[vars] <- sapply(gsn_df3[vars], function(x) ifelse(x <= 0, 0, x))
summary(gsn_df3[vars])

gsn_df3$genhz <- aqp::generalize.hz(
  gsn_df3$hzn_desgn,
  new = c("O", "A", "BA|AB", "B", "C"),
  pattern = c("O", "^A", "AB|BA", "^B|^IIB|2B^|^IIIB|3B^", "^C|^IIIC|^4C"),
  ordered = TRUE
  )
aggregate(x = gsn_df3[vars], by = list(gsn_df3$genhz), FUN = function(x) format(summary(x), digits = 1))
as.data.table(gsn_df3)[, lapply(.SD, function(x) quantile(x, probs = seq(0, 1, 0.1), na.rm = TRUE)), .SDcols = vars, by = "genhz"]



# segment depths ----
gsn_seg_df <- segment(gsn_df3, hzdepcols = c("hzn_top", "hzn_bot"), intervals = c(0, 30, 60, 100))
gsn_seg_df <- gsn_seg_df |>
  subset(valid_dep_all == TRUE) |>
  transform(thk = hzn_bot - hzn_top)

vars <- c("clay_total", "silt_total", "sand_total", "total_carbon_ncs", "bulk_density_third_bar", "ph_h2o", "N_pt", "P_pt", "potassium_water_extractable", "cec_nh4_ph_7")
gsn_seg_df <- as.data.table(gsn_seg_df)[
  , lapply(.SD, function(x) weighted.mean(x, w = thk)), 
  .SDcols = vars, 
  by = c("pedon_key", "segment_id")
  ] |>
  as.data.frame()

# save(gsn_seg_df, file = file.path(fp, "gsn_seg_df.RData"))
load(file = file.path(fp, "gsn_seg_df.RData"))

test <- transform(gsn_seg_df, id = paste(pedon_key, segment_id))
seg_lo <- reshape(
  test, 
  direction = "long",
  idvar = "id",
  timevar = "variable", times = vars,
  v.names = "value", varying = vars
  )


ggplot(seg_lo, aes(x = value)) +
  geom_density(fill = "grey", alpha = 0.5) +
  facet_wrap(~ variable, scales = "free") +
  ggtitle("Distribution of GSN 0-30cm soil properties")

logit <- function(y) log(y / (100 - y))
hist(logit(gsn_seg_df[gsn_seg_df$segment_id == "000-030", ]$clay_total))



# transform skewed variables ----
vars <- c("total_carbon_ncs", "N_pt", "P_pt", "potassium_water_extractable")
nms <- paste0(vars, "_log")
gsn_seg_df[nms] <- lapply(gsn_seg_df[vars], function(x) log(x + 0.01))

# save(gsn_seg_df, file = file.path(fp, "gsn_seg_t_df.RData"))
load(file = file.path(fp, "gsn_seg_t_df.RData"))



