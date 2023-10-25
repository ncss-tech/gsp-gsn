

# load packages
library(terra)
library(sf)
library(data.table)



# find files ----
# original predictions
lf <- list.files("./cov100/Predictionsv2", recursive = TRUE)
fn <- lf[grepl(".tif$", lf) & grepl("adj", lf) & !grepl("^trimmed", lf)]


# trimmed predictions
lf <- list.files("./cov100/Predictionsv2/trimmed", recursive = TRUE)
fn <- lf[grepl(".tif$", lf) & grepl("adj", lf)]



# parse filenames ----
fn_df <- lapply(fn, function(x) {
  x2 = {strsplit(x, "/")[[1]][2] ->.; strsplit(., "_")[[1]]}
  var = x2[1] |>
    # tolower() |> 
    unname()
  times = x2[which(grepl("[1:100]", x2) & grepl("x$", x2))] |> 
    unname() |> gsub("x", "", x = _) |> 
    as.integer()
  deps = x2[which(x2 == "cm") - 1] |> 
    unname() |> 
    as.integer()
  return(data.frame(var = var, times = times, dep = deps))
}) 
fn_df <- do.call("rbind", fn_df) |> cbind(fn = fn)
fn_df <- fn_df[order(fn_df$var, fn_df$dep), ] |>
  subset(var != "anylithicdpt")

lu <- list(
  solus = c("cec7", "ph", "claytotal", "silttotal", "sandtotal", "soc", "dbovendry"),
  gsn   = c("CEC",  "pH", "Clay",      "Silt",      "Sand",      "SOC", "BD")
)

fn_df <- subset(fn_df, var %in% lu$solus)
fn_df$var2 <- factor(fn_df$var, levels = lu$solus, labels = lu$gsn) |> as.character()

knitr::kable(fn_df)

# |    |var       | times| dep|fn                                                 |var2 |
#   |:---|:---------|-----:|---:|:--------------------------------------------------|:----|
#   |15  |cec7      |    10|   0|CEC7_gRPI_250k/cec7_r_10x_0_cm_2D_QRFadj.tif       |CEC  |
#   |20  |cec7      |    10|   5|CEC7_gRPI_250k/cec7_r_10x_5_cm_2D_QRFadj.tif       |CEC  |
#   |17  |cec7      |    10|  15|CEC7_gRPI_250k/cec7_r_10x_15_cm_2D_QRFadj.tif      |CEC  |
#   |19  |cec7      |    10|  30|CEC7_gRPI_250k/cec7_r_10x_30_cm_2D_QRFadj.tif      |CEC  |
#   |21  |cec7      |    10|  60|CEC7_gRPI_250k/cec7_r_10x_60_cm_2D_QRFadj.tif      |CEC  |
#   |16  |cec7      |    10| 100|CEC7_gRPI_250k/cec7_r_10x_100_cm_2D_QRFadj.tif     |CEC  |
#   |18  |cec7      |    10| 150|CEC7_gRPI_250k/cec7_r_10x_150_cm_2D_QRFadj.tif     |CEC  |
#   |22  |claytotal |     1|   0|Clay_gRPI_250k/claytotal_r_1x_0_cm_2D_QRFadj.tif   |Clay |
#   |27  |claytotal |     1|   5|Clay_gRPI_250k/claytotal_r_1x_5_cm_2D_QRFadj.tif   |Clay |
#   |24  |claytotal |     1|  15|Clay_gRPI_250k/claytotal_r_1x_15_cm_2D_QRFadj.tif  |Clay |
#   |26  |claytotal |     1|  30|Clay_gRPI_250k/claytotal_r_1x_30_cm_2D_QRFadj.tif  |Clay |
#   |28  |claytotal |     1|  60|Clay_gRPI_250k/claytotal_r_1x_60_cm_2D_QRFadj.tif  |Clay |
#   |23  |claytotal |     1| 100|Clay_gRPI_250k/claytotal_r_1x_100_cm_2D_QRFadj.tif |Clay |
#   |25  |claytotal |     1| 150|Clay_gRPI_250k/claytotal_r_1x_150_cm_2D_QRFadj.tif |Clay |
#   |1   |dbovendry |   100|   0|BD_gRPI_250k/dbovendry_r_100x_0_cm_2D_QRFadj.tif   |BD   |
#   |6   |dbovendry |   100|   5|BD_gRPI_250k/dbovendry_r_100x_5_cm_2D_QRFadj.tif   |BD   |
#   |3   |dbovendry |   100|  15|BD_gRPI_250k/dbovendry_r_100x_15_cm_2D_QRFadj.tif  |BD   |
#   |5   |dbovendry |   100|  30|BD_gRPI_250k/dbovendry_r_100x_30_cm_2D_QRFadj.tif  |BD   |
#   |7   |dbovendry |   100|  60|BD_gRPI_250k/dbovendry_r_100x_60_cm_2D_QRFadj.tif  |BD   |
#   |2   |dbovendry |   100| 100|BD_gRPI_250k/dbovendry_r_100x_100_cm_2D_QRFadj.tif |BD   |
#   |4   |dbovendry |   100| 150|BD_gRPI_250k/dbovendry_r_100x_150_cm_2D_QRFadj.tif |BD   |
#   |87  |sandtotal |     1|   0|Sand_gRPI_250k/sandtotal_r_1x_0_cm_2D_QRFadj.tif   |Sand |
#   |92  |sandtotal |     1|   5|Sand_gRPI_250k/sandtotal_r_1x_5_cm_2D_QRFadj.tif   |Sand |
#   |89  |sandtotal |     1|  15|Sand_gRPI_250k/sandtotal_r_1x_15_cm_2D_QRFadj.tif  |Sand |
#   |91  |sandtotal |     1|  30|Sand_gRPI_250k/sandtotal_r_1x_30_cm_2D_QRFadj.tif  |Sand |
#   |93  |sandtotal |     1|  60|Sand_gRPI_250k/sandtotal_r_1x_60_cm_2D_QRFadj.tif  |Sand |
#   |88  |sandtotal |     1| 100|Sand_gRPI_250k/sandtotal_r_1x_100_cm_2D_QRFadj.tif |Sand |
#   |90  |sandtotal |     1| 150|Sand_gRPI_250k/sandtotal_r_1x_150_cm_2D_QRFadj.tif |Sand |
#   |101 |silttotal |     1|   0|Silt_gRPI_250k/silttotal_r_1x_0_cm_2D_QRFadj.tif   |Silt |
#   |106 |silttotal |     1|   5|Silt_gRPI_250k/silttotal_r_1x_5_cm_2D_QRFadj.tif   |Silt |
#   |103 |silttotal |     1|  15|Silt_gRPI_250k/silttotal_r_1x_15_cm_2D_QRFadj.tif  |Silt |
#   |105 |silttotal |     1|  30|Silt_gRPI_250k/silttotal_r_1x_30_cm_2D_QRFadj.tif  |Silt |
#   |107 |silttotal |     1|  60|Silt_gRPI_250k/silttotal_r_1x_60_cm_2D_QRFadj.tif  |Silt |
#   |102 |silttotal |     1| 100|Silt_gRPI_250k/silttotal_r_1x_100_cm_2D_QRFadj.tif |Silt |
#   |104 |silttotal |     1| 150|Silt_gRPI_250k/silttotal_r_1x_150_cm_2D_QRFadj.tif |Silt |
#   |108 |soc       |  1000|   0|SOC_gRPI_250k/soc_r_1000x_0_cm_2D_QRFadj_bt.tif    |SOC  |
#   |113 |soc       |  1000|   5|SOC_gRPI_250k/soc_r_1000x_5_cm_2D_QRFadj_bt.tif    |SOC  |
#   |110 |soc       |  1000|  15|SOC_gRPI_250k/soc_r_1000x_15_cm_2D_QRFadj_bt.tif   |SOC  |
#   |112 |soc       |  1000|  30|SOC_gRPI_250k/soc_r_1000x_30_cm_2D_QRFadj_bt.tif   |SOC  |
#   |114 |soc       |  1000|  60|SOC_gRPI_250k/soc_r_1000x_60_cm_2D_QRFadj_bt.tif   |SOC  |
#   |109 |soc       |  1000| 100|SOC_gRPI_250k/soc_r_1000x_100_cm_2D_QRFadj_bt.tif  |SOC  |
#   |111 |soc       |  1000| 150|SOC_gRPI_250k/soc_r_1000x_150_cm_2D_QRFadj_bt.tif  |SOC  |
  
  

# compute weighted averages ----

fn_l <- split(fn_df, fn_df$var)
lapply(fn_l[2:6], function(x) {
  
  var  <- x$var[1]
  var2 <- x$var2[1]
  times <- x$times[1]
  dep  <- x$dep[1]
  
  cat("aggregating ", var, as.character(Sys.time()), "\n")
  
  tmp_r <- rast(file.path("./cov100/Predictionsv2/trimmed", x[, "fn"]))
  
  deps <- c(0, 5, 15, 30, 60, 100, 150, 200)
  wts  <- (deps[-1] - deps[-8])[1:3] / 30
  tmp_30cm <- (tmp_r[[1]] * wts[1] + tmp_r[[2]] * wts[2] + tmp_r[[3]] * wts[3]) / times
  
  writeRaster(
    tmp_30cm,
    filename = file.path("cov100/glosis", paste0("USA_GSNmap_", var2, "_Map30", ".tif"))
  )
})



# create tiles ----
tiles <- ext(clay_r) |> 
  vect() |> 
  st_as_sf() |> 
  st_make_grid() |> 
  st_as_sf()
st_crs(tiles) <-  st_crs(crs(clay_r, proj = TRUE))


# extract tiles
clay_fn <- subset(fn_df, var == "claytotal", select = fn)
clay_r  <- rast(file.path("./cov100/Predictionsv2/trimmed", clay_fn$fn))


clay_r2 <- crop(clay_r, vect(tiles[82, ]))
tmp_r <- c(clay_r2, clay_r2, clay_r2)


# linear interpolate ----
li_fun <- function(rs) {
  
  if (sum(!is.na(rs)) > 1) {
    mean(approx(
      x = c(0, 5, 15, 30),
      y = rs,
      xout = 0:30,
      na.rm = F
    )$y,
    na.rm = T
    ) 
  } else NA_real_
}

system.time(test <- app(clay_r[[1:4]], li_fun, cores = 15))

fn_l <- split(fn_df, fn_df$var)
lapply(fn_l, function(x) {
  
  var  <- x$var[1]
  var2 <- x$var2[1]
  times <- x$times[1]
  dep  <- x$dep[1]
  
  cat("aggregating ", var, as.character(Sys.time()), "\n")
  
  tmp_r <- rast(file.path("./cov100/Predictionsv2/trimmed", x[, "fn"]))
  
  tmp_30cm <- app(tmp_r[[1:4]], li_fun, cores = 15)
  
  writeRaster(
    tmp_30cm,
    filename = file.path("cov100/glosis", paste0("USA_GSNmap_", var2, "_Map30", "_linear.tif"))
  )
})



# mpspline ----

mps_fun <- function(rs) {
  suppressMessages(mpspline2::mpspline_one(
    site = data.frame(
      id = 1, 
      top   = c(0, 5, 15, 30, 60, 100, 150), 
      bot   = c(1, 6, 16, 31, 61, 101, 151),
      # value = c(10, 15, 25, 30, 10, 5, 2)
      value = rs
    ),
    d = c(0, 30)
  )$est_dcm)
}

system.time(test <- app(tmp_r, mps_fun, cores = 15))


tmp_df <- data.frame(
  id = 1:ncell(tmp_r), 
  values(tmp_r), 
  prop = "clay"
)
deps <- c(0, 5, 15, 30, 60, 100, 150)
vars <- paste0("z", deps)
names(tmp_df)[2:8] <- vars
# tmp_lo <- stats::reshape(
#   tmp_df, 
#   # idvar = c("id"), #"prop"),
#   direction = "long",
#   timevar = "variable", times = vars,
#   v.names = "value",    varying = vars
#   )
tmp_lo <- tmp_df |>
  as.data.table() |> 
  melt(id.vars = c("id", "prop"), measure.vars = vars) |>
  subset(!is.na(value))

tmp_lo$top <- tmp_lo$variable |>
  as.character() |>
  gsub("z", "", x = _) |> 
  as.integer()
tmp_lo$bot <- tmp_lo$top + 1
tmp_lo <- tmp_lo[, c("id", "top", "bot", "value", "prop")][order(id, top), ]



# data.table option 1
test2[1:100000, suppressMessages(mpspline(cbind(.BY, .SD), var_name = "value")), by = id] |> system.time()


# data.table option 2 with parallel apply
idx <- 1:1e6
system.time(tmp_s <- split(tmp_lo[idx], tmp_lo$id[idx]))

system.time(temp3 <- lapply(tmp_s[1:1000], function(x) {
  mpspline_one(x, var_name = "value", d = c(0, 30))
}))

future::plan(future::multisession, workers = 15)
system.time(temp4 <- future.apply::future_lapply(tmp_s, function(x) {
  suppressMessages(mpspline_one(x, var_name = "value", d = c(0, 30)))
}))
test3 <- do.call("rbind", lapply(mps30, function(x) x$est_dcm))
test3 <- data.frame(id = names(mps30), value = test3)



# Interps

system.time(test <- hsg_calc(tmp_r))









