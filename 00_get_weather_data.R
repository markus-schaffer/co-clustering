# Title: Fetch weather data from the Danish Meteorological Institute (DMI)
#
# Purpose: Use the API of DMI to fetch the outdoor temperature and radiation for
# the weather station located closet to the location of the SHM data
#
# To run this script one must obtain a free API key from DMI as described here:
# https://opendatadocs.dmi.govcloud.dk/en/Authentication
#
# Data files:
#   - None
#
# Function files:
#   - None
#
# Author: M. Schaffer
# Contact details: msch@build.aau.dk


# Import packages  -------------------------------------------------------
library(httr2)
library(data.table)
library(purrr)
library(lubridate)

# API parameters ---------------------------------------------------------
api_key <- "" # To be supplied by the user
station_id <- "06031"
time_range <- "2020-01-01T00:00:00+01:00/2021-12-31T23:59:59+01:00"
parameter_id <- c("temp_dry", "radia_glob")


# Define and run api ------------------------------------------------------

# Define and perform API calls
base_request <- request("https://dmigw.govcloud.dk/v2/metObs/collections/observation/items") |>
  req_headers("Accept" = "application/json", "X-Gravitee-Api-Key" = api_key) |>
  req_url_query(
    stationId = station_id,
    datetime = time_range,
    limit = "300000", # maximum allowed
  ) |>
  req_throttle(rate = 500 / 5)

full_requests <- map(parameter_id, ~ base_request |> req_url_query(parameterId = .x))
resp <- req_perform_sequential(full_requests, )


# Parse the response to long format data.tables
fn_parse <- function(res) {
  parsed <- resp_body_json(res, simplifyVector = T)
  if (resp_is_error(res)) {
    stop(paste(
      "DMI API request failed",
      resp_status(res),
      parsed$message
    ))
  }
  if (is.null(parsed$features$properties)) {
    result <- NULL
  } else {
    result <- parsed$features$properties |>
      setDT()
    result[, created := NULL]
    result <- dcast(result, observed + stationId ~ parameterId, value.var = "value")
    # Parse observation time
    result[, observed := fast_strptime(result$observed, "%Y-%m-%dT%H:%M:%S%z", tz = "UTC", lt = FALSE)]
  }
  return(result)
}
api_result <- map(resp, ~ fn_parse(.x))


# Combine different quantities into one data.table
api_result <- discard(api_result, is.null)
weather_data <- Reduce(function(...) merge(..., all = TRUE), api_result)
weather_data[, stationId := NULL]

# Sort and set time to local Danish time
setorder(weather_data, observed)
weather_data[, observed := with_tz(observed, tz = "Europe/Copenhagen")]

# Aggregate to daily values
weather_data[, day := floor_date(observed, unit = "day")]
weather_data <- weather_data[, .(day_mean_temp = mean(temp_dry), day_mean_rad = mean(radia_glob)), by = day]

# Add NA values for missing data
full_data <- data.table(day = floor_date(seq.POSIXt(weather_data[1, day], last(weather_data[, day]), by = "day"), unit = "day"))
weather_data <- merge(full_data, weather_data, by = "day", all = TRUE)

# Linearly interpolate nay missing values - i.e. if a whole day is missing
myna.approx <- function(x, y) {
  y[is.na(y)] <- approx(x = x[!is.na(y)], y = y[!is.na(y)], xout = x[is.na(y)])$y
  y
}
weather_data[, day_mean_temp := myna.approx(x = as.numeric(day), y = day_mean_temp)]
weather_data[, day_mean_rad := myna.approx(x = as.numeric(day), y = day_mean_rad)]

fwrite(weather_data, "data/00_weather_data.csv")
