
lf <- list.files("./cov100/Predictionsv2", recursive = TRUE)
fn <- lf[grepl(".tif$", lf) & grepl("adj", lf)]

vars <- sapply(fn, function(x) strsplit(x, "/|_")[[1]][2]) |> tolower() |> unname()
times <- sapply(fn, function(x) {
  x2 = strsplit(x, "/|_")[[1]]
  idx = which(grepl("[1:100]", x2) & grepl("x$", x2))
  x2[idx]
}) |> unname() |> gsub("x", "", x = _) |> as.integer()
deps <- sapply(fn, function(x) {
  x2 = strsplit(x, "/|_")[[1]]
  idx = which(x2 == "cm")
  x2 = x2[idx - 1]
}) |> unname() |> as.integer()
dat <- {
  data.frame(fn = fn, var = vars, multipler = times, depth = deps) ->.;
  .[order(.$var, .$dep), ]
}

idx <- which(dat$var == "clay")

clay <- rast(file.path("./cov100/Predictionsv2", dat[idx, "fn"]))

# tiles
tiles <- ext(clay) |> vect() |> st_as_sf() |> st_make_grid() |> st_as_sf()
st_crs(tiles) <-  st_crs(crs(clay, proj = TRUE))

tmp_r <- crop(clay, vect(tiles[82, ]))

test <- cbind(id = 1:ncell(tmp), values(tmp, dataframe = TRUE), prop = "clay")
deps <- c(0, 5, 15, 30, 60, 100, 150)
vars <- paste0("z", deps)
names(test)[2:8] <- vars
# test_lo <- stats::reshape(
#   test, 
#   # idvar = c("id"), #"prop"),
#   direction = "long",
#   timevar = "variable", times = vars,
#   v.names = "value",    varying = vars
#   )
test_lo <- data.table::melt(data.table::as.data.table(test), id.vars = c("id", "prop"), measure.vars = vars) |>
  as.data.frame()

test_lo$top <- as.character(test_lo$variable) |> gsub("z", "", x = _) |> as.integer()
test_lo$bot <- test_lo$top + 1
test2 <- test_lo[c("id", "top", "bot", "value", "prop")] |>
  poorman::arrange(id, top)


mps30 <- mpspline2::mpspline(test2, var = "value", d = c(0, 30))
test3 <- do.call("rbind", lapply(mps30, function(x) x$est_dcm))
test3 <- data.frame(id = names(mps30), value = test3)



