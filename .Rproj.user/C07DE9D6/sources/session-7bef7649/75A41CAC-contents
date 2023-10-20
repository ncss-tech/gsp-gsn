# SDA ----

fp <- "C:/Users/stephen.roecker/OneDrive - USDA/data/scd"


## table names & keys ----
tbs <- c(
  "lab_analysis_procedure", 
  "lab_analyte", 
  "lab_area", 
  "lab_calculations_including_estimates_and_default_values", 
  "lab_chemical_properties",
  "lab_combine_nasis_ncss",
  "lab_layer",
  "lab_major_tr_elements_oxides",
  "lab_method_code",
  "lab_mineralogy_glass_count",
  "lab_mir",
  "lab_mir_wavelength",
  "lab_pedon",
  "lab_physical_properties",
  "lab_preparation",
  "lab_rosetta_key",
  "lab_site",
  "lab_webmap",
  "lab_xray_and_thermal"
)

td <- data.frame(tbs)
td$key[td$tbs == "lab_area"]        <- "area_key"
td$key[td$tbs == "lab_site"]        <- "site_key"
td$key[td$tbs == "lab_analyte"]     <- "analyte_key"
td$key[td$tbs == "lab_preparation"] <- "prep_key"
# td$key[td$tbs == "lab_rosetta_key"] <- "rosetta_key" # doesn't exist in SDA

idx <- c("lab_layer", "lab_mir", "lab_rosetta_key")
td$key[td$tbs %in% idx]             <- "layer_key"

idx <- c("lab_analysis_procedure", "lab_method_code")
td$key[td$tbs %in% idx]             <- "procedure_key"

idx <- c("lab_pedon", "lab_combine_nasis_ncss", "lab_webmap")
td$key[td$tbs %in% idx]             <- "pedon_key"

td$key[is.na(td$key)]               <- "labsampnum"

# missing tables
idx <- ! td$tbs %in% c("lab_mir", "lab_mir_wavelength", "lab_major_tr_elements_oxides")
td <- td[idx, ]



## table sizes ----

tbs_n <- by(td, td$tbs, function(x){
  cat("getting", x$tbs)
  tmp <- SDA_query(paste0(
    "SELECT 
    COUNT(*) n_records,
    COUNT(DISTINCT ", x$key, ") n_keys
    FROM ", x$tbs,
    ";"))
  tmp <- cbind(x, tmp)
})
td_n <- do.call("rbind", tbs_n)
row.names(td_n) <- NULL


## table columns ----
td_col <- by(td_n, td$tbs, function(x) {
  tmp <- SDA_query(paste0(
    "SELECT TOP(1) * FROM ", x$tbs, ";"
  ))[0, ]
  x$col <- list(tmp)
  return(x)
})
td_col <- do.call("rbind", td_col)
row.names(td_col) <- NULL


## query tables ----

scd <- by(td_col, td_col$tbs, function(x) {
  
  cat("getting", x$tbs, as.character(format(Sys.time(), "%Y/%m/%d %H:%M:%S")), "\n")
  
  n <- 2.5e4
  if (x$n_records > n) {
    x$chunk <- ceiling(x$n_records/n)
  } else x$chunk <- 1
  
  tmp <- lapply(seq(1, n * x$chunk, by = n), function(i) {
    cat("  getting offset ", i, as.character(format(Sys.time(), "%Y/%m/%d %H:%M:%S")), "\n")
    SDA_query(paste0(
      "SELECT * 
    FROM ", x$tbs, " 
    ORDER BY ", x$key, " 
    OFFSET ", i, " ROWS 
    FETCH NEXT 25000 ROWS ONLY;"
    ))
  })
  tmp <- do.call("rbind", tmp)
  
  write.csv(tmp, paste0(x$tbs, ".csv"), row.names = FALSE)
  
  return(NULL)
})



## save outputs ----
scd_l <- lapply(paste0(td$tbs, ".csv"), read.csv)
names(scd_l) <- gsub("lab_", "", td$tbs)
# saveRDS(scd_l, file = paste0(fp, "/ncss-scd_sda_20230808.rds"))
scd_l <- readRDS(file = paste0(fp, "/ncss-scd_sda_20230808.rds"))

lapply(names(scd_l), function(x) {
  cat("writing ", x, " to .sqlite \n")
  sf::write_sf(scd_l[[x]], dsn = paste0(fp, "/ncss-scd_sda_20230808.sqlite"), layer = x)
})


## check sqlite file
library(DBI)

con <- dbConnect(RSQLite::SQLite(), paste0("C:/Users/stephen.roecker/Downloads/NCSSLabDataMartSQLite/NCSSLabDataMartSQLite.sqlite3"))
(ldm_names <- dbListTables(con))


# test
scd_l <- readRDS(file = paste0(fp, "/ncss-scd_sda_20230808.rds"))

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


## create bbox ----

scd_sf <- subset(scd_l$site, complete.cases(longitude_std_decimal_degrees, latitude_std_decimal_degrees))
scd_sf <- st_as_sf(
  scd_sf,
  coords = c("longitude_std_decimal_degrees", "latitude_std_decimal_degrees"),
  crs = "4326"
)                 
st_bbox(scd_sf)  


library(dplyr)

test <- scd_l$chemical_properties %>%
  filter(!is.na(total_nitrogen_ncs)) %>%
  inner_join(scd_l$layer, by = "labsampnum") %>%
  inner_join(scd_l$combine_nasis_ncss, by = "pedon_key")
test2 <- test %>% 
  group_by(pedon_key) %>% 
  summarize(n_N = sum(!is.na(total_nitrogen_ncs)) > 0) %>%
  inner_join(scd_l$combine_nasis_ncss, by = "pedon_key") %>%
  select(pedon_key, n_N, samp_year, corr_year, SSL_year, site_year)


