
# see Clay 2017 p. 112-113 DOI:10.2134/practicalmath2017.0022
# see Logsdon 2008 p. 199 DOI:10.2136/2008.soilsciencestepbystep
# see Halvin 2017 p. 31-32
# see Oyatokun 2023 DOI:10.9734/asrj/2023/v7i2127

# ppm = part per million
# mw = molecular weight (in grams or milligrams)
# eq = equilivalent charge

# ppm = mg/kg
# cmol/kg =  meq/100g


ppm_to_cmolkg <- function(ppm, mw, eq) {
  # ppm * (1 g / 1000 mg) * ( 1 mol / X mw in g) * (X cmol / 1 eq) = X cmol / kg
  ppm * (1 / 1000) * (1 / mw) * (eq / 1) * (100 / 1)
}

ppm_to_meq100g <- function(ppm, mw, eq) {
  ppm * (1/ mw) * (eq / 1) * (1 / (10 * 100)) * 100
}

cmolkg_to_ppm <- function(cmolkg, mw, eq) {
  # ppm * (1 g / 1000 mg) * ( 1 mol / X mw in g) * (X cmol / 1 eq) = X cmol / kg
  cmolkg * 1000 * mw * (1 / eq) * (1 / 100)
}


# Ca
ppm_to_cmolkg(2000, 40, 2)
ppm_to_meq100g(2000, 40, 2)

ppm_to_meqg <- 


ppm_to_cmolkg(1350, 40, 2)


ppm = 1350; mm = 40; eq = 2
