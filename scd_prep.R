
library(sf)
library(ggplot2)
library(dplyr)
library(data.table)


# load snapshots ----
fp <- "C:/Users/stephen.roecker/OneDrive - USDA/data/scd"
scd_l <- readRDS(file = paste0(fp, "/ncss-scd_sda_20230808.rds"))
# scd_l2 <- readRDS(file = file.path(fp, "ncss_labdata.rds"))
f <- readRDS(file = "C:/Users/stephen.roecker/Box/nasis-pedons/fetchNASIS_spc_20230926.rds")



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

s <- within(s, {
  peiid   = as.integer(peiid)
  obsdate = as.Date(obsdate, "%Y-%m-%d")
  obsyear = format(obsdate, "%Y") |> as.integer()
})

scd_l$combine_nasis_ncss <- merge(
  scd_l$combine_nasis_ncss, 
  s[c("peiid", "obsdate", "obsyear")],
  by.x = "pedoniid", by.y = "peiid",
  all.x = TRUE,
  sort = FALSE
)


scd_l$combine_nasis_ncss <- within(scd_l$combine_nasis_ncss, {
  samp_classdate2 = strptime(samp_classdate, "%e/%m/%Y %H:%M:%S %p")
  corr_classdate2 = strptime(corr_classdate, "%e/%m/%Y %H:%M:%S %p")
  SSL_classdate2  = strptime(SSL_classdate,  "%e/%m/%Y %H:%M:%S %p")
  site_obsdate2   = strptime(site_obsdate,   "%e/%m/%Y %H:%M:%S %p")
  
  samp_year = format(samp_classdate2, "%Y") |> as.integer()
  corr_year = format(corr_classdate2, "%Y") |> as.integer()
  SSL_year  = format(SSL_classdate2,  "%Y") |> as.integer()
  year      = apply(cbind(samp_year, corr_year, SSL_year), 1, min, na.rm = TRUE)
  year      = ifelse(is.na(obsyear), year, obsyear)
})



# coordinates ----
## rename columns ----
nm <- names(scd_l$site) 
idx <- grepl("longitude", nm)
names(scd_l$site)[idx] <- sub("longitude", "lon", nm[idx])
idx <- grepl("latitude", nm)
names(scd_l$site)[idx] <- sub("latitude",  "lat", nm[idx])
names(scd_l$site)[nm == "horizontal_datum_name"] <- "datum"

nm <- names(scd_l$site) 
var <- "std_decimal_degrees"
idx <- which(grepl(var, nm))
names(scd_l$site)[idx] <- gsub(var, "dd", nm[idx])

table(dms = complete.cases(scd_l$site[3:11]), dd = complete.cases(scd_l$site[12:13]))


## fix datum ----
scd_l$site <- within(scd_l$site, {
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
  y = (lat_degrees + lat_minutes / 60 + lat_seconds / 3600) * ifelse(lon_direction == "south", -1, 1)
})
table(datum = scd_l$site$datum, ellps =  scd_l$site$ellps, useNA = "always")


# estimate datum for missing
s2 <- aggregate(year ~ site_key, data = scd_l$combine_nasis_ncss, min, na.rm = TRUE)
scd_l$site <- merge(scd_l$site, s2, by = "site_key", all.x = TRUE, sort = FALSE)

test <- scd_l$site |>
  transform(xy = complete.cases(x, y),
            year = as.integer(year)
            ) |>
  aggregate(xy ~ year + datum, data = _, sum, na.rm = TRUE)

subset(test, year %in% 1920:2023) |>
  ggplot(aes(x = year, y = xy, col = datum)) +
  geom_line(size = 1, alpha = 0.7)


scd_l$site <- within(scd_l$site, {
  idx = (is.na(datum) | datum == "") & is.na(ellps)
  datum[idx  & year <= 2012]  = "NAD83"
  datum[idx & year  >   2012] = "WGS84"
  datum[is.na(datum) & is.na(ellps)] = "NAD83"
  idx = NULL
})
table(datum = scd_l$site$datum, ellps =  scd_l$site$ellps, useNA = "always")


## convert dms to dd ----
s2 <- {
  split(scd_l$site, ~ paste(datum, ellps), drop = TRUE) ->.;
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
scd_site <- merge(scd_l$site, s2[vars], by = "site_key", all.x = TRUE, sort = FALSE)
scd_site <- within(scd_site, {
  lon_dd = ifelse(is.na(lon_dd), X, lon_dd)
  lat_dd = ifelse(is.na(lat_dd), Y, lat_dd) 
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

# 5 digits = ~ 1m; 4 digits .0= ~ 11m; 3 digits = ~ 111m
test <- coord_prec2(scd_site$lon_dd, scd_site$lat_dd)
{apply(test, 1, max) >= 3} |> summary()

scd_site <- cbind(scd_site, test)



## project dd ----
scd_sf <- subset(scd_site, complete.cases(lon_dd, lat_dd))
# scd_sf <- subset(scd_site, complete.cases(longitude_std_decimal_degrees, latitude_std_decimal_degrees))
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

# saveRDS(scd_sf, file = "C:/Users/stephen.roecker/USDA/NSSC - SBS/projects/gsp-gsn/data/scd_sf.rds")



# remove duplicates ----
vars <- c("layer", "physical_properties", "chemical_properties", "")
sapply(scd_l[vars], function(x) sum(duplicated(x$labsampnum), na.rm = TRUE))

aggregate(siteiid ~ site_key, data = scd_l$combine_nasis_ncss, length) |> summary()
aggregate(pedoniid ~ pedon_key, data = scd_l$combine_nasis_ncss, length) |> summary()
aggregate(pedon_key ~ site_key, data = scd_l$combine_nasis_ncss, length) |> summary()
aggregate(pedlabsampnum ~ pedon_key, data = scd_l$combine_nasis_ncss, length) |> summary()





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
  inner_join(scd_l$combine_nasis_ncss, by = "pedon_key")
test2 <- test %>% 
  group_by(pedon_key) %>% 
  summarize(n_N = sum(!is.na(total_nitrogen_ncs)) > 0) %>%
  inner_join(scd_l$combine_nasis_ncss, by = "pedon_key") %>%
  select(pedon_key, n_N, samp_year, corr_year, SSL_year, site_year)


# tally the number of pedon per decade with 1 or more measurements for each soil property
test <- as.data.table(gsn_df)
vars <- c(vars_phys, vars_chem)
vars <- vars[!grepl("method", vars)]
test <- test[, lapply(.SD, function(x) any(!is.na(x))), .SDcols = vars, by = c("site_key", "pedon_key")]
vars <- c("site_key", "pedon_key", "year")
test <- merge(scd_l$combine_nasis_ncss[vars], test, by = vars[1:2], all.x = TRUE, sort = FALSE)
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
with(gsn_df, cor(total_nitrogen_ncs, total_nitrogen_pt, use = "complete.obs"))

# scatter plot
with(gsn_df, plot(total_nitrogen_ncs, total_nitrogen_pt))
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



