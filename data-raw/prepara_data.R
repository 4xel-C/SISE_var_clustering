# ============================================
# Dataset preparation for package vignette.
# ============================================

library(readxl)
library(dplyr)

# ---- autos ----
autos <- read_excel("data-raw/autos.xls")

# Convert to factor
autos <- autos %>%
  dplyr::mutate(across(where(is.character), as.factor))

# Convert to df
autos <- data.frame(autos)

# Save data file
usethis::use_data(autos, overwrite = TRUE)

# ---- players_22_cleaned ----
players <- read.csv("data-raw/players_22_cleaned.csv", stringsAsFactors = FALSE)

# Save data file
usethis::use_data(players, overwrite = TRUE)

# ---- Votes ----
vote <- read_excel("data-raw/vote.xls")

# Convert to factor
vote <- vote %>%
  dplyr::mutate(across(where(is.character), as.factor))

# Save data file
usethis::use_data(vote, overwrite = TRUE)

# ---- players_22_cleaned ----
canines <- read_excel("data-raw/races_canines.xls")

# Convert to factor
canines <- canines %>%
  dplyr::mutate(across(where(is.character), as.factor))

usethis::use_data(canines, overwrite = TRUE)

