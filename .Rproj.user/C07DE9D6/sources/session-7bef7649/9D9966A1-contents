
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
## chemical properties
vars <- c("total_nitrogen|phosphorus_bray|phosphorus_mehlich|phosphorus_olwsen|potassium_mehlich|^cec|ph_h2o|ph_cacl2|ph_kcl|total_carbon|organic_carbon")
nms  <- names(scd_l$chemical_properties)
vars <- nms[grep(vars, nms)]
gsn_df <- scd_l$chemical_properties[c(1:3, match(vars, nms))]


# uncode method codes
idx <- grep("method", names(gsn_df))
test <- lapply(idx, function(x) {
  factor(gsn_df[[x]], levels = scd_l$method_code$proced_code, labels = scd_l$method_code$proced_name) |> droplevels()
}) |> as.data.frame()
names(test) <- paste0(names(gsn_df)[idx], "_uncoded")
gsn_df <- cbind(gsn_df, test)



# remove duplicates ----


## create bbox ----

scd_sf <- subset(scd_l$site, complete.cases(longitude_std_decimal_degrees, latitude_std_decimal_degrees))
scd_sf <- st_as_sf(
  scd_sf,
  coords = c("longitude_std_decimal_degrees", "latitude_std_decimal_degrees"),
  crs = "4326"
)                 
st_bbox(scd_sf)  



# Phosphorus ----
# subset non-missing data
nm  <- names(gsn_df)
idx <- which(grepl("phosphorus", nm) & !grepl("method", nm))
gsn_df[idx] > 0 ->.;
# .[rowSums(., na.rm = TRUE) > 0, ] |>
rowSums(., na.rm = TRUE) |> 
  table()
gsn_df_P <- gsn_df[rowSums(., na.rm = TRUE) > 0, c(1:3, idx)]

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

lm(phosphorus_bray1 ~ phosphorus_mehlich_3, data = gsn_df_P) |> summary()
lm(phosphorus_mehlich_3 ~ phosphorus_mehlich3_extractable, data = gsn_df_P) |> summary()

