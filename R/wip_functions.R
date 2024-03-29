# Load packages
library(terra)
library(MultiscaleDTM)
library(randomForest)
library(caret)
library(nnet)

# -----------------------------------------------------------------------------
# Creates surface metrics from a set of external files

surface_met1 <- function(len, metrics = c("grad", "plan", "prof", "dev"),
                         dem_dir, exec_dir, out_dir=getwd(), re_sample = NA) {
  
  # Checking to see if directories exist
  if(!file.exists(dem_dir)) {
    stop("DEM directory does not exist!")
  }
  
  if(!dir.exists(exec_dir)) {
    stop("Executable Files directory does not exist!")
  }
  
  # Prepare inputs
  dem_dir <- normalizePath(dem_dir)
  out_dir <- normalizePath(out_dir)
  if(!endsWith(out_dir, "\\")) {
    out_dir <- paste0(out_dir, "\\")
  }
  exec_dir <- normalizePath(exec_dir)
  
  # Write input file
  file_name <- paste0(out_dir, "input_makeGrids.txt")
  file.create(file_name)
  
  writeLines(c("# Input file for makeGrids",
               "",
               paste0("DEM: ", dem_dir),
               paste0("SCRATCH DIRECTORY: ", out_dir),
               paste0("LENGTH SCALE: ", len)), con = file_name)
  
  if("grad" %in% metrics) {
    write(paste0("GRID: GRADIENT, OUTPUT FILE = ", out_dir, "grad", len, ".flt"),
          file = file_name, append = T) 
  }
  
  if("plan" %in% metrics) {
    write(paste0("GRID: PLAN CURVATURE, OUTPUT FILE = ", out_dir,
                 "plan", len), file = file_name, append = T)
  }
  
  if("prof" %in% metrics) {
    write(paste0("GRID: PROFILE CURVATURE, OUTPUT FILE = ", out_dir,
                 "prof", len), file = file_name, append = T)
  }
  
  # Run surface metrics sans DEV
  system(paste0(exec_dir, "\\makeGrids"), input = file_name)
  
  # Writing input file for DEV
  if ("dev" %in% metrics) {
    if(is.na(re_sample)) {
      stop("Set re_sample level")
    }
    
    # Prepare inputs
    file_name <- paste0(out_dir, "input_localRelief.txt")
    rad <- len / 2
    
    # Create and write input file
    file.create(file_name)
    writeLines(c("# Input file for LocalRelief",
                 "# Creating by surfaceMetrics.R",
                 paste0("# On ", Sys.time()),
                 paste0("DEM: ", dem_dir),
                 paste0("SCRATCH DIRECTORY: ", out_dir),
                 paste0("RADIUS: ", rad),
                 paste0("DOWN SAMPLE: ", re_sample),
                 paste0("SAMPLE INTERVAL: ", re_sample),
                 paste0("OUTPUT LOCAL RASTER: ", out_dir, "local", len)),
               con = file_name)
    
    # Run DEV in console
    system(paste0(exec_dir, "\\localRelief"), input = file_name)
  }
}

# -----------------------------------------------------------------------------
# Creates surface metrics in R

surface_met <- function(DEM, len, export = FALSE,
                        elev_dev = c("grad", "plan", "prof", "dev", "twi")) {
  # Checks if inputs are file names and loads them in
  if(is.character(DEM)) {
    if(!file.exists(DEM)) {
      stop("Cannot find DEM file")
    }
    DEM <- terra::rast(DEM)
  }
  # Sets up the resolution
  k <- round(len/terra::res(DEM)[1])
  if (k %% 2 == 0) {
    k <- k + 1
  }
  
  j <- k/2 + 0.5
  
  # Sets up the moving window which is used to calculate the metrics
  w_mat <- matrix(nrow = k, ncol = k)
  w_mat[c(rep(1,3), rep(j,3), rep(k,3)), rep(c(1,j,k), 3)] <- 1
  
  # Initialize the inputs for the model
  in_rast <- list()
  
  if("grad" %in% elev_dev) {
    j <- k/2 - 0.5
    
    xl.end <- matrix(c(1, rep(NA_real_, times=k-1)), ncol=k, nrow=1)
    xr.end <- matrix(c(rep(NA_real_, times=k-1), 1), ncol=k, nrow=1)
    
    x.mids <- matrix(NA_real_, ncol=k, nrow=j-1)
    
    xl.mid <- matrix(c(2, rep(NA_real_, times=k-1)), ncol=k, nrow=1)
    xr.mid <- matrix(c(rep(NA_real_, times=k-1), 2), ncol=k, nrow=1)
    
    xl.mat <- rbind(xl.end, x.mids, xl.mid, x.mids, xl.end)
    xr.mat <- rbind(xr.end, x.mids, xr.mid, x.mids, xr.end)
    
    yt.end <- matrix(c(1, rep(NA_real_, times=k-1)), ncol=1, nrow=k)
    yb.end <- matrix(c(rep(NA_real_, times=k-1), 1), ncol=1, nrow=k)
    
    y.mids <- matrix(NA_real_, ncol=j-1, nrow=k)
    
    yt.mid <- matrix(c(2, rep(NA_real_, times=k-1)), ncol=1, nrow=k)
    yb.mid <- matrix(c(rep(NA_real_, times=k-1), 2), ncol=1, nrow=k)
    
    yt.mat <- cbind(yt.end, y.mids, yt.mid, y.mids, yt.end)
    yb.mat <- cbind(yb.end, y.mids, yb.mid, y.mids, yb.end)
    
    dz.dx.l <- terra::focal(DEM, xl.mat, fun=sum, na.rm=T, na.policy = "omit")
    dz.dx.r <- terra::focal(DEM, xr.mat, fun=sum, na.rm=T, na.policy = "omit")
    dz.dy.t <- terra::focal(DEM, yt.mat, fun=sum, na.rm=T, na.policy = "omit")
    dz.dy.b <- terra::focal(DEM, yb.mat, fun=sum, na.rm=T, na.policy = "omit")
    
    wts.l <- terra::focal(!is.na(DEM), w=xl.mat, fun=sum, na.rm=TRUE,
                          na.policy = "omit")
    wts.r <- terra::focal(!is.na(DEM), w=xr.mat, fun=sum, na.rm=TRUE,
                          na.policy = "omit")
    wts.t <- terra::focal(!is.na(DEM), w=yt.mat, fun=sum, na.rm=TRUE,
                          na.policy = "omit")
    wts.b <- terra::focal(!is.na(DEM), w=yb.mat, fun=sum, na.rm=TRUE,
                          na.policy = "omit")
    dz.dx <- ((dz.dx.r/wts.r) - (dz.dx.l/wts.l))/(2*j*terra::xres(DEM))
    dz.dy <- ((dz.dy.t/wts.t) - (dz.dy.b/wts.b))/(2*j*terra::yres(DEM))
    
    grad <- sqrt(dz.dx^2 + dz.dy^2)
    in_rast <- c(in_rast, grad)
    
    names(in_rast)[length(in_rast)] <- paste0("grad", len)
  }
  
  if("plan" %in% elev_dev) {
    if ("prof" %in% elev_dev) {
      both <- MultiscaleDTM::Qfit(DEM, metrics = c("planc", "profc"),
                                  w = k, na.rm = T)
      in_rast <- c(in_rast, both[[1]], both[[2]])
      
      names(in_rast)[length(in_rast)-1] <- paste0("plan", len)
      names(in_rast)[length(in_rast)] <- paste0("prof", len)
    } else {
      plan <- MultiscaleDTM::Qfit(DEM, metrics = "planc", w = k, na.rm = T)
      in_rast <- c(in_rast, plan)
      
      names(in_rast)[length(in_rast)] <- paste0("plan", len)
    }
  } else if("prof" %in% elev_dev) {
    prof <- MultiscaleDTM::Qfit(DEM, metrics = "profc", w = k, na.rm = T)
    in_rast <- c(in_rast, prof)
    
    names(in_rast)[length(in_rast)] <- paste0("prof", len)
  }
  
  if("dev" %in% elev_dev) {
    dev <- (DEM - focal(DEM, w = w_mat, fun = "mean", na.rm = T, na.policy = "omit")) / focal(DEM, w = w_mat, fun = "sd", na.rm = T, na.policy = "omit") 
    in_rast <- c(in_rast, dev)
    
    names(in_rast)[length(in_rast)] <- paste0("dev", len)
  }
  
  if("twi" %in% elev_dev) {
    topidx <- topmodel::topidx(terra::as.matrix(DEM), res = terra::res(DEM)[1])
    twi <- terra::setValues(DEM, topidx$atb)
    terra::values(twi) <- ifelse(terra::values(twi) < 0, 0, terra::values(twi))
    twi <- terra::focal(twi, w = w_mat, mean, na.rm = T, na.policy = "all")
    
    in_rast <- c(in_rast, twi)
    
    names(in_rast)[length(in_rast)] <- paste0("twi", len)
  }
  
  # Exports the surface metrics
  if(export) {
    for(i in 1:length(in_rast)) {
      writeRaster(in_rast[[i]],
                  filename = paste0(names(in_rast[i]), len, ".tif"))
    }
  }
  return(in_rast)
}

# -----------------------------------------------------------------------------
# Build training points

build_train_pts <- function(region_poly, wet_poly, multi_class = FALSE,
                            wet_types = c("Freshwater Forested/Shrub Wetland",
                                          "Freshwater Emergent Wetland",
                                          "Freshwater Pond",
                                          "Estuarine and Marine Wetland",
                                          "Riverine", "Lake",
                                          "Estuarine and Marine Deepwater",
                                          "Other"), wet_field = "WETLAND_TY",
                            sample_points = c(50, 150), export = FALSE) {
  
  # Loads in polygons if input is a file name
  if(is.character(wet_poly[[1]])) {
    temp_poly <- list()
    for(i in 1:length(wet_poly)) {
      temp_poly[[i]] <- terra::vect(wet_poly) 
    }
    wet_poly <- temp_poly
  }
  
  if(is.character(region_poly)) {
    region_poly <- terra::vect(region_poly)
  }
  
  # Filters the wetland polygons to only include wanted types
  wet_poly <- wet_poly[unlist(wet_poly[[wet_field]]) %in% wet_types]
  if(length(wet_poly) == 0) {
    stop("No wetlands to sample!")
  }
  
  # Cropping the wetland polygon(s) to the overall region
  wet_poly <- terra::project(wet_poly, region_poly)
  wet_poly <- terra::crop(wet_poly, region_poly)
  
  # Checks if output is supposed to be more than two classes before proceeding
  if(multi_class) {
    # Initialize parameters
    train_crds <- NULL
    train_atts <- c()
    wet_samp <- sample_points[1]
    up_samp <- sample_points[2]
    
    # Sample points for each wetland class
    for(i in 1:length(wet_types)) {
      temp_poly <- wet_poly[unlist(wet_poly[[wet_field]]) == wet_types[i]]
      
      # Checking if polygons of that type of wetland exist
      if(length(temp_poly) == 0) {
        warning(paste0(wet_types[i], " not found!"))
      }
      else {
        wet_crds <- NULL
        samp_wet_pts <- terra::spatSample(temp_poly, wet_samp)
        coords <- terra::crds(samp_wet_pts)
        wet_crds <- rbind(wet_crds, coords)
        
        num_coords <- nrow(coords)
        while(num_coords < wet_samp) {
          new_points <- terra::spatSample(temp_poly,
                                          wet_samp-(num_coords))
          new_crds <- terra::crds(new_points)
          wet_crds <- rbind(wet_crds, new_crds)
          num_coords <- num_coords + nrow(new_crds)
        }
        
        train_crds <- rbind(train_crds, wet_crds)
        train_atts <- c(train_atts, rep(wet_types[i], wet_samp))
      }
    }
    
    # Sample points from non-wetland areas
    up_poly <- terra::erase(region_poly, wet_poly)
    samp_up_pts <- terra::spatSample(up_poly, up_samp)
    up_crds <- terra::crds(samp_up_pts)
    
    # Create the points
    train_crds <- rbind(train_crds, up_crds)
    train_atts <- c(train_atts, rep("UPL", up_samp))
    train_atts <- data.frame(class = factor(train_atts))
    pts <- terra::vect(train_crds, atts = train_atts,
                       crs = terra::crs(region_poly))
  } else {
    
    # Sample the wetland points
    for(i in 1:length(wet_types)) {
      temp_poly <- wet_poly[unlist(wet_poly[[wet_field]]) == wet_types[i]]
      
      # Checking if polygons of that type of wetland exist
      if(length(temp_poly) == 0) {
        warning(paste0(wet_types[i], " not found!"))
      }
    }
    
    samp_wet_points <- terra::spatSample(wet_poly, sample_points[1])
    wet_crds <- terra::crds(samp_wet_points)
    num_wet <- nrow(wet_crds)
    while(num_wet < sample_points[1]) {
      new_points <- terra::spatSample(temp_poly,
                                      num_points[i]-(num_wet))
      new_crds <- terra::crds(new_points)
      wet_crds <- rbind(wet_crds, new_crds)
      num_wet <- num_wet + nrow(new_crds)
    }
    
    # Sample points from non-wetland areas
    up_poly <- terra::erase(region_poly, wet_poly)
    samp_up_pts <- terra::spatSample(up_poly, sample_points[2])
    up_crds <- terra::crds(samp_up_pts)
    num_up <- nrow(up_crds)
    
    # Create the points
    train_crds <- rbind(wet_crds, up_crds)
    train_atts <- data.frame(class = factor(c(rep("WET", num_wet),
                                              rep("UPL", num_up))))
    pts <- terra::vect(train_crds, atts = train_atts,
                       crs = terra::crs(region_poly))
  }
  
  # Return the points and exports them, if desired
  if(export) {
    terra:writeVector(pts, filename = "trainingdata.shp")
  }
  
  return(pts)
}
# -----------------------------------------------------------------------------
# Builds the model

build_model <- function(in_rasts, poly_inputs = list(), train, ref_raster,
                        model_type = "forest", model_params = list(ntree = 200),
                        class_field_name = "class") {
  
  # Checking if input rasters are file names, then load them in
  if(is.character(in_rasts[1])) {
    temp_rast <- rep(list(), length(in_rasts))
    for(i in 1:length(in_rasts)) {
      temp_rast[[i]] <- terra::rast(in_rasts[i])
    }
    names(temp_rast) <- in_rasts
    in_rasts <- temp_rast
  }
  
  # Checks if there are any polygon inputs
  if(length(poly_inputs) > 0) {
    
    # Checking to see the polygon inputs are filenames
    if(is.character(poly_inputs[[1]])) {
      temp_poly <- rep(list(), length(poly_inputs))
      for(i in 1:length(poly_inputs)) {
        if(!file.exists(poly_inputs[[i]])) {
          stop(paste0("Cannot find poly input file:", poly_inputs[i]))
        }
        temp_poly[[i]] <- terra::vect(poly_inputs[i]) 
      }
      names(temp_poly) <- poly_inputs
      poly_inputs <- temp_poly
    }
    
    # Rasterize polygon inputs
    for(i in 1:length(poly_inputs)) {
      vr_name <- names(poly_inputs)[i]
      temp_rast <- terra::rasterize(poly_inputs[i], ref_raster, field = vr_name)
      in_rasts <- c(in_rasts, temp_rast)
    }
  }
  
  # Ensure that all inputs are covering the same area
  print("Formatting inputs")
  for(i in 1:length(in_rasts)) {
    in_rasts[[i]] <- terra::project(in_rasts[[i]], ref_raster)
    in_rasts[[i]] <- terra::crop(in_rasts[[i]], ref_raster)
  }
  
  
  # Set up training data
  print("Setting up training data")
  train <- terra::project(train, ref_raster)
  df_train <- data.frame(class = factor(as.vector(unlist(train[[class_field_name]]))))
  for(i in 1:length(in_rasts)) {
    vals <- terra::extract(in_rasts[[i]], train, ID = F)
    df_train <- cbind(df_train, vals)
  }
  df_train <- na.omit(df_train)
  colnames(df_train) <- c("class", names(in_rasts))
  
  # Build the model
  print("Building model")
  if(model_type == "forest"){
    mod <- randomForest::randomForest(class ~ ., data = df_train, 
                                      ntree = model_params$ntree)
  } else if (model_type == "tree") {
    mod <- randomForest::randomForest(class ~ ., data = df_train, ntree = 1)
  } else if(model_type == "glm") {
    if(length(levels(df_train$class)) > 2) {
      mod <- nnet::multinom(class ~ ., data = df_train)
    } else {
      mod <- glm(class ~ ., data = df_train, family = "binomial")
    }
  } else if(model_type == "knn") {
    mod <- caret::knn3(formula = class ~ ., data = df_train, k = model_params$k)
  } else {
    stop("Incorrect model type")
  }
  
  print("Done!")
  return(mod)
}


# -----------------------------------------------------------------------------
# Runs the model

run_model <- function(mod, in_rasts = list(), poly_inputs = list(), ref_raster,
                      model_type = "forest", class_rast = FALSE,
                      export = FALSE) {
  
  # Checking if inputs are file names, then load them in
  if(is.character(in_rasts[1])) {
    temp_rast <- rep(list(), length(in_rasts))
    for(i in 1:length(in_rasts)) {
      temp_rast[[i]] <- terra::rast(in_rasts[i])
    }
    names(temp_rast) <- in_rasts
    in_rasts <- temp_rast
  }
  
  if(length(poly_inputs) > 0) {
    if(is.character(poly_inputs[[1]])) {
      temp_poly <- rep(list(), length(poly_inputs))
      for(i in 1:length(poly_inputs)) {
        if(!file.exists(poly_inputs[[i]])) {
          stop(paste0("Cannot find poly input file:", poly_inputs[i]))
        }
        temp_poly[[i]] <- terra::vect(poly_inputs[i]) 
      }
      names(temp_poly) <- poly_inputs
      poly_inputs <- temp_poly
    }
    for(i in 1:length(poly_inputs)) {
      vr_name <- names(poly_inputs)[i]
      temp_rast <- terra::rasterize(poly_inputs[i], ref_raster, field = vr_name)
      in_rasts <- c(in_rasts, temp_rast)
    }
  }
  
  # Ensure that all inputs are covering the same area
  print("Formatting inputs")
  for(i in 1:length(in_rasts)) {
    in_rasts[[i]] <- terra::project(in_rasts[[i]], ref_raster)
    in_rasts[[i]] <- terra::crop(in_rasts[[i]], ref_raster)
  }
  
  # Stacks the rasters on top of each other to create one raster
  print("Stacking rasters")
  input_raster <- in_rasts[[1]]
  if(length(in_rasts) > 1) {
    for(i in 2:length(in_rasts)) {
      input_raster <- c(input_raster, in_rasts[[i]])
    }
  }
  names(input_raster) <- names(in_rasts)
  
  # Run the model
  print("Running model")
  
  if(class_rast) {
    if(isTRUE(mod$call[[1]] == "glm")) {
      output <- terra::predict(input_raster, mod, na.rm = T,
                               type = "response")
      vals <- terra::values(output)
      vals <- ifelse(vals > 0.5, "WET", "UPL")
      terra::values(output) <- vals
      
    } else {
      output <- terra::predict(input_raster, mod, na.rm = T)
    }
    
  } else {
    if(isTRUE(mod$call[[1]] == "glm")) {
      output <- terra::predict(input_raster, mod, na.rm = T,
                               type = "response")
    } else {
      output <- terra::predict(input_raster, mod, na.rm = T, type = "prob")
    }
  }  
  
  if(export) {
    for(i in 1:length(output)) {
      file_name <- paste0(names(input_raster)[i], "prob.tif")
      terra::writeRaster(output[[i]], filename = file_name)
    }
  }
  
  print("Done!")
  return(output)
}

# -----------------------------------------------------------------------------
# Tests model

CV_err <- function(in_rasts, poly_inputs = list(), ref_raster,
                   model_type = "forest", model_params = list(ntree = 200), 
                   train, kfold= 5, class_field_name = "class") {
  
  # Checking if inputs are file names, then load them in
  if(is.character(in_rasts[1])) {
    temp_rast <- rep(list(), length(in_rasts))
    for(i in 1:length(in_rasts)) {
      temp_rast[[i]] <- terra::rast(in_rasts[i])
    }
    names(temp_rast) <- in_rasts
    in_rasts <- temp_rast
  }
  
  if(length(poly_inputs) > 0) {
    if(is.character(poly_inputs[[1]])) {
      temp_poly <- rep(list(), length(poly_inputs))
      for(i in 1:length(poly_inputs)) {
        if(!file.exists(poly_inputs[[i]])) {
          stop(paste0("Cannot find poly input file:", poly_inputs[i]))
        }
        temp_poly[[i]] <- terra::vect(poly_inputs[i]) 
      }
      names(temp_poly) <- poly_inputs
      poly_inputs <- temp_poly
    }
  }
  
  # Convert the polygons into rasters
  if(length(poly_inputs) > 0) {
    for(i in 1:length(poly_inputs)) {
      vr_name <- names(poly_inputs)[i]
      temp_rast <- terra::rasterize(poly_inputs[i], ref_raster, field = vr_name)
      in_rasts <- c(in_rasts, temp_rast)
    }
  }
  
  # Ensure that all inputs are covering the same area
  for(i in 1:length(in_rasts)) {
    in_rasts[[i]] <- terra::project(in_rasts[[i]], ref_raster)
    in_rasts[[i]] <- terra::crop(in_rasts[[i]], ref_raster)
  }
  
  # Set up training data
  train <- terra::project(train, ref_raster)
  df_train <- data.frame(class = factor(as.vector(unlist(train[[class_field_name]]))))
  for(i in 1:length(in_rasts)) {
    vals <- terra::extract(in_rasts[[i]], train, ID = F)
    df_train <- cbind(df_train, vals)
  }
  df_train <- na.omit(df_train)
  colnames(df_train) <- c("class", names(in_rasts))
  
  k <- kfold
  test_err <- c()
  index <- sample(k, nrow(df_train), replace = T)
  
  for(i in 1:k) {
    train_df <- df_train[index != i,]
    test_df <- df_train[index == i,]
    y_test <- test_df$class
    
    if(model_type == "forest"){
      mod <- randomForest::randomForest(class ~ ., data = train_df, 
                                        ntree = model_params$ntree)
    } else if (model_type == "tree") {
      mod <- randomForest::randomForest(class ~ ., data = train_df, ntree = 1)
    } else if(model_type == "glm") {
      if(length(levels(df_train$class)) > 2) {
        mod <- nnet::multinom(class ~ ., data = train_df)
      } else {
        mod <- glm(class ~ ., data = train_df, family = "binomial")
      }
    } else if(model_type == "knn") {
      mod <- caret::knn3(formula = class ~ ., data = train_df,
                         k = model_params$k)
    } else {
      stop("Incorrect model type")
    }
    
    if(model_type == "glm") {
      pred <- predict(mod, newdata = test_df, type = "response")
    } else {
      pred <- predict(mod, newdata = test_df)
    }
    test_err[i] <- mean(pred != y_test)
  }
  mean_err <- mean(test_err)
  ci_err <- round(100 * (mean_err + c(-1, 1)*qnorm(0.975)*sd(test_err)/k), 1)
  print(paste0("Test Error Estimate: ", round(mean_err * 100, 1), "%"))
  print(paste0("95% Confidence Interval: [", ci_err[1], ", ", ci_err[2], "]"))
}
