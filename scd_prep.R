
library(dplyr)
library(sf)


# load snapshot ----
fp <- "C:/Users/stephen.roecker/OneDrive - USDA/data/scd"
scd_l <- readRDS(file = paste0(fp, "/ncss-scd_sda_20230808.rds"))



# metadata
dm <- dm::dm_from_src(soilDB:::.openNASISchannel(), learn_keys = TRUE)
dm$system |> collect() |> subset(sysver == "Lab SDA Data Mart 1.0") |> as.data.frame()
md <- dm$attribute |> collect() |> 
  subset(sysiidref == 41045) |>
  merge(
    dm$uom |> collect() |> subset(select = c(uomiid, uomsym)), 
    by.x = "uomiidref", by.y = "uomiid", all.x = TRUE, sort = FALSE
    )
md <- md[c(ncol(md), 1:(ncol(md) - 1))]



# evaluate ----
test <- scd_l$chemical_properties %>%
  filter(!is.na(total_nitrogen_ncs)) %>%
  inner_join(scd_l$layer, by = "labsampnum") %>%
  inner_join(scd_l$combine_nasis_ncss, by = "pedon_key")
test2 <- test %>% 
  group_by(pedon_key) %>% 
  summarize(n_N = sum(!is.na(total_nitrogen_ncs)) > 0) %>%
  inner_join(scd_l$combine_nasis_ncss, by = "pedon_key") %>%
  select(pedon_key, n_N, samp_year, corr_year, SSL_year, site_year)



# convert dates ----
scd_l$combine_nasis_ncss <- within(scd_l$combine_nasis_ncss, {
  samp_classdate2 = strptime(samp_classdate, "%e/%m/%Y %H:%M:%S %p")
  corr_classdate2 = strptime(corr_classdate, "%e/%m/%Y %H:%M:%S %p")
  SSL_classdate2  = strptime(SSL_classdate,  "%e/%m/%Y %H:%M:%S %p")
  site_obsdate2   = strptime(site_obsdate,   "%e/%m/%Y %H:%M:%S %p")
  
  samp_year = format(samp_classdate2, "%Y")
  corr_year = format(corr_classdate2, "%Y")
  SSL_year  = format(SSL_classdate2,  "%Y")
  site_year = format(samp_classdate2, "%Y")
})



# extract variables ----
vars <- c("layer", "physical_properties", "chemical_properties", "")
sapply(scd_l[vars], function(x) sum(duplicated(x$labsampnum), na.rm = TRUE))


pat <- c("total_nitrogen|^phosphorus_|^new_zealand_phos|^potassium_|^cec|ph_h2o|ph_cacl2|ph_kcl|total_carbon|organic_carbon")
nms  <- names(scd_l$chemical_properties)
vars <- nms[grep(pat, nms)]
gsn_df <- scd_l$layer |>
  merge(scd_l$physical_properties, by = "labsampnum", all.x = TRUE, sort = FALSE) |>
  merge(scd_l$chemical_properties[c(1:3, match(vars, nms))], 
                                   by = "labsampnum", all.x = TRUE, sort = FALSE)


# uncode method codes ----
idx <- grep("method", names(gsn_df))
test <- lapply(idx, function(x) {
  factor(gsn_df[[x]], levels = scd_l$method_code$proced_code, labels = scd_l$method_code$proced_name) |> droplevels()
}) |> as.data.frame()
names(test) <- paste0(names(gsn_df)[idx], "_desc")
gsn_df <- cbind(gsn_df, test)



# remove duplicates ----


# coordinates ----

scd_sf <- subset(scd_l$site, complete.cases(longitude_std_decimal_degrees, latitude_std_decimal_degrees))
scd_sf <- st_as_sf(
  scd_sf,
  coords = c("longitude_std_decimal_degrees", "latitude_std_decimal_degrees")
  # crs = 4326
)
st_crs(scd_sf) <- 4326
st_bbox(scd_sf)  


# intersect with the GADM reference
world_bndy <- geodata::world(path = getwd(), resolution = 1) |> 
  st_as_sf() |>
  st_make_valid()

idx <- st_intersects(scd_sf, world_bndy) |> 
  lapply(function(x) x[1]) |> 
  unlist()
scd_sf <- scd_sf |> cbind(st_drop_geometry(world_bndy[idx, 1:2]))

table(USA = scd_sf$GID_0 == "USA", Count = !is.na(scd_sf$site_key), useNA = "always")



# nitrogen ----
nm  <- names(gsn_df)
idx <- which(grepl("nitrogen", nm) & !grepl("method", nm))
summary(gsn_df[idx])
aggregate(total_nitrogen_ncs ~ total_nitrogen_ncs_method_desc, data = gsn_df, function(x) {
  summary(x) |> round(2) |> t()
  })

table(gsn_df$total_nitrogen_ncs_method_desc, gsn_df$total_nitrogen_ncs_method)


## pedotransfer function from Tomer 2017 ----
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



# phosphorus ----
nm  <- names(gsn_df)
idx <- which(grepl("phosphorus", nm) & !grepl("method", nm))
summary(gsn_df[idx])


# subset non-missing data
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


## pedotransfer function for missing Bray 1 ----
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



