# Demo of how to run the wipr tool thus far

# Load in functions
source("ForestedWetlands/wipr/wipr_no_fortran.R")

# Load in the DEM, training data, as well as other rasters and polygons
r <- terra::rast("PackForest/PF_DTM3.tif")
pts <- terra::vect("PackForest/PF_trainingdata.shp")

# Calculate the elevation derivatives
input_rast1 <- surface_met(r, len = 20, elev_dev = c("grad", "prof", "plan"))
input_rast2 <- surface_met(r, len = 100, elev_dev = c("grad", "prof", "plan"))
input_rast3 <- surface_met(r, len = 10, elev_dev = c("twi"))
input_rast <- c(input_rast1, input_rast2)

# Build the model
mod <- build_model(input_rast, train = pts)

# Find the estimated error of the model
CV_err(input_rast, train = pts)

# Calculate the probability raster using the models
prob <- run_model(mod, input_rast)

# Plot the probability rasters for each class
plot(prob["WET"])
