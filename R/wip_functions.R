# Load packages
library(terra)
library(MultiscaleDTM)
library(randomForest)
library(caret)
library(nnet)

# -----------------------------------------------------------------------------
# Creates surface metrics

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
  k <- round(len/res(DEM)[1])
  if (k %% 2 == 0) {
    k <- k + 1
  }
  
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
    dev <- (DEM - focal(DEM, w = k, fun = "mean", na.rm = T, na.policy = "omit")) / focal(DEM, w = k, fun = "sd", na.rm = T, na.policy = "omit") 
    in_rast <- c(in_rast, rast_dev)
    
    names(in_rast)[length(in_rast)] <- paste0("dev", len)
  }
  
  if("twi" %in% elev_dev) {
    topidx <- topmodel::topidx(terra::as.matrix(DEM), res = res(DEM)[1])
    a <- terra::setValues(DEM, topidx$area)
    twi <- a / tan(terra::terrain(DEM, unit = "radians"))
    values(twi) <- ifelse(values(twi) < 0, 0, values(twi))
    twi <- terra::focal(twi, w = k, mean, na.rm = T, na.policy = "omit")
    
    in_rast <- c(in_rast, twi)
    
    names(in_rast)[length(in_rast)] <- paste0("twi", len)
  }
  
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

build_train_pts <- function() {
  
}

# -----------------------------------------------------------------------------
# Builds the model

build_model <- function(in_rasts, poly_inputs = list(), train, 
                        model_type = "forest", model_params = list(ntree = 200),
                        class_field_name = "Class") {
  
  # Checking if inputs are file names, then load them in
  for(i in 1:length(in_rasts)) {
    if(is.character(in_rasts[[i]])) {
      in_rasts[i] <- terra::rast(in_rasts[[i]])
    }
  }
   
  if(length(poly_inputs) > 0) {
    for(i in 1:length(poly_inputs)) {
      if(is.character(poly_inputs[[i]])) {
        if(!file.exists(poly_inputs[[i]])) {
          stop(paste0("Cannot find poly input file:", poly_inputs[i]))
        }
        poly_inputs[i] <- terra::vect(poly_inputs[[i]]) 
      }
    }
  }
  
  # Changing polygons to rasters
  if(length(poly_inputs) > 0) {
    for(i in 1:length(poly_inputs)) {
      vr_name <- names(poly_inputs)[i]
      temp_rast <- terra::rasterize(poly_inputs[i], DEM, field = vr_name)
      in_rast <- c(in_rast, temp_rast)
    }
  }
  
  # Set up training data
  df_train <- data.frame(class = factor(as.vector(unlist(train[[class_field_name]]))))
  for(i in 1:length(in_rasts)) {
    vals <- terra::extract(in_rasts[[i]], train, ID = F)
    df_train <- cbind(df_train, vals)
  }
  df_train <- na.omit(df_train)
  colnames(df_train) <- c("class", names(in_rasts))
  
  # Build the model
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

  return(mod)
}


# -----------------------------------------------------------------------------
# Runs the model

run_model <- function(mod, in_rasts = list(), poly_inputs = list(),
                      model_type = "forest") {
  for(i in 1:length(in_rasts)) {
    if(is.character(in_rasts[[i]])) {
      if(!file.exists(in_rasts[[i]])) {
        stop(paste0("Cannot find poly input file:", in_rasts[i]))
      }
      in_rasts[i] <- terra::rast(in_rasts[[i]])
    }
  }
  
  if(length(poly_inputs) > 0) {
    for(i in 1:length(poly_inputs)) {
      if(is.character(poly_inputs[[i]])) {
        if(!file.exists(poly_inputs[[i]])) {
          stop(paste0("Cannot find poly input file:", poly_inputs[i]))
        }
        poly_inputs[i] <- terra::vect(poly_inputs[[i]]) 
      }
    }
  }
  
  if(length(poly_inputs) > 0) {
    for(i in 1:length(poly_inputs)) {
      vr_name <- names(poly_inputs)[i]
      temp_rast <- terra::rasterize(poly_inputs[i], DEM, field = vr_name)
      in_rasts <- c(in_rasts, temp_rast)
    }
  }
  
  input_raster <- in_rasts[[1]]
  if(length(in_rasts) > 1) {
    for(i in 2:length(in_rasts)) {
      ext(in_rasts[[i]]) <- ext(in_rasts[[1]])
      input_raster <- c(input_raster, in_rasts[[i]])
    }
  }
  names(input_raster) <- names(in_rasts)
  
  # Run the model
  if(model_type == "glm") {
    prob_rast <- terra::predict(input_raster, mod, na.rm = T, type = "response")
  } else {
    prob_rast <- terra::predict(input_raster, mod, na.rm = T, type = "prob")
  }
  return(prob_rast)
}

# -----------------------------------------------------------------------------
# Tests model

CV_err <- function(in_rasts, poly_inputs = list(), 
                   model_type = "forest", model_params = list(ntree = 200), 
                   train, kfold= 5, class_field_name = "Class") {
  for(i in 1:length(in_rasts)) {
    if(is.character(in_rasts[[i]])) {
      in_rasts[i] <- terra::rast(in_rasts[[i]])
    }
  }
  
  if(length(poly_inputs) > 0) {
    for(i in 1:length(poly_inputs)) {
      if(is.character(poly_inputs[[i]])) {
        if(!file.exists(poly_inputs[[i]])) {
          stop(paste0("Cannot find poly input file:", poly_inputs[i]))
        }
        poly_inputs[i] <- terra::vect(poly_inputs[[i]]) 
      }
    }
  }
  
  if(length(poly_inputs) > 0) {
    for(i in 1:length(poly_inputs)) {
      vr_name <- names(poly_inputs)[i]
      temp_rast <- terra::rasterize(poly_inputs[i], DEM, field = vr_name)
      in_rast <- c(in_rast, temp_rast)
    }
  }
  
  # Set up training data
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

# -----------------------------------------------------------------------------
# Attempts to do whole WIP tool without fortran

dem_to_prob <- function(DEM, len, poly_inputs = list(),
                 elev_dev = c("grad", "plan", "prof", "dev", "twi"),
                 train, class_field_name = "Class", wet_name = "WET", 
                 export_prob = FALSE, prob_rast_name = NA, cv_err = FALSE) {
  
  # Checks if inputs are file names and loads them in
  if(is.character(DEM)) {
    if(!file.exists(DEM)) {
      stop("Cannot find DEM file")
    }
    DEM <- terra::rast(DEM)
  }
  
  if(is.character(train)) {
    if(!file.exists(train)) {
      stop("Cannot find training points file")
    }
    train <- terra::vect(train)
  }
  
  if(length(poly_inputs) > 0) {
    for(i in 1:length(poly_inputs)) {
      if(is.character(poly_inputs[[i]])) {
        if(!file.exists(poly_inputs[[i]])) {
          stop(paste0("Cannot find poly input file:", poly_inputs[i]))
        }
        poly_inputs[[i]] <- terra::vect(poly_inputs[[i]]) 
      }
    }
  }
  
  # Sets up the resolution
  k <- round(len/res(DEM)[1])
  if (k %% 2 == 0) {
    k <- k + 1
  }
  
  # Initialize the inputs for the model
  print("Creating input metrics")
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
    in_rast <- c(in_rast, grad = grad)
  }
  
  if("plan" %in% elev_dev) {
    if ("prof" %in% elev_dev) {
      both <- MultiscaleDTM::Qfit(DEM, metrics = c("planc", "profc"),
                                  w = k, na.rm = T)
      in_rast <- c(in_rast, plan = both[[1]], prof = both[[2]])
    } else {
      plan <- MultiscaleDTM::Qfit(DEM, metrics = "planc", w = k, na.rm = T)
      in_rast <- c(in_rast, plan = plan)
    }
  } else if("prof" %in% elev_dev) {
    prof <- MultiscaleDTM::Qfit(DEM, metrics = "profc", w = k, na.rm = T)
    in_rast <- c(in_rast, prof = prof)
  }
  
  if("dev" %in% elev_dev) {
    dev <- (DEM - focal(DEM, w = k, fun = "mean", na.rm = T, na.policy = "omit")) / focal(DEM, w = k, fun = "sd", na.rm = T, na.policy = "omit") 
    in_rast <- c(in_rast, dev = dev)
  }
  
  if("twi" %in% elev_dev) {
    topidx <- topmodel::topidx(terra::as.matrix(DEM), res = res(DEM)[1])
    a <- terra::setValues(DEM, topidx$area)
    twi <- a / tan(terra::terrain(DEM, unit = "radians"))
    values(twi) <- ifelse(values(twi) < 0, 0, values(twi))
    twi <- terra::focal(twi, w = k, mean, na.rm = T, na.policy = "omit")
    
    in_rast <- c(in_rast, twi = twi)
  }
  
  if(length(poly_inputs) > 0) {
    for(i in 1:length(poly_inputs)) {
      vr_name <- names(poly_inputs)[i]
      temp_rast <- terra::rasterize(poly_inputs[i], DEM, field = vr_name)
      in_rast <- c(in_rast, temp_rast)
    }
  }
  
  # Set up training data
  df_train <- data.frame(class = factor(as.vector(unlist(train[[class_field_name]]))))
  for(i in 1:length(in_rast)) {
    vals <- terra::extract(in_rast[[i]], train, ID = F)
    df_train <- cbind(df_train, vals)
  }
  df_train <- na.omit(df_train)
  colnames(df_train) <- c("class", names(in_rast))
  
  # Build the model
  print("Building and running model")
  
  mod <- randomForest::randomForest(class ~ ., data = df_train, ntree = 200)
  
  # Find estimated test error using cross-validation
  if(cv_err) {
    k <- 5
    index <- sample(k, nrow(df_train), replace = T)
    
    for(i in 1:k) {
      train_df <- df_train[index != i,]
      test_df <- df_train[index == i,]
      y_test <- test_df$class
      
      mod1 <- randomForest::randomForest(mod$formula, data = train_df,
                                         ntree = 200)
      pred <- predict(mod1, newdata = test_df)
      test_err[i] <- mean(pred != y_test)
    }
    print(paste0("Test Error Estimate: ", mean(test_err)))
  }
  
  # Move the rasters into one multilayered raster
  input_raster <- in_rast[[1]]
  if(length(in_rast) > 1) {
    for(i in 2:length(in_rast)) {
      input_raster <- c(input_raster, in_rast[[i]])
    }
  }
  names(input_raster) <- names(in_rast)
  
  # Run the model
  prob_rast <- terra::predict(input_raster, mod, na.rm = T, type = "prob")
  wet_prob <- prob_rast[wet_name]
  
  # Return the wetlands probability raster as the output
  if(export_prob) {
    if(!is.na(prob_rast_name)) {
      file_name <- paste0(prob_rast_name, ".tif")
    } else {
      file_name = "wetProb.tif"
    }
    terra::writeRaster(wet_prob, filename = file_name)
    print(paste0("Created Wetlands Probability Raster: ", file_name))
  }
  return(wet_prob)
}

# -----------------------------------------------------------------------------
# Surface metrics function using makeGrids (fortran)
surface_met1 <- function(len, metrics = c("grad", "plan", "prof", "dev"),
                         dem_dir, exec_dir, out_dir=getwd(), re_sample = NA) {
  
  # Checking to see if directories exist
  if(!file.exists(dem_dir)) {
    stop("DEM directory does not exist!")
  }
  
  if(!dir.exists(exec_dir)) {
    stop("Executable Files directory does not exist!")
  }
  
  # if(!endsWith(dem_dir, ".flt")) {
  #   r <- terra::rast(dem_dir)
  #   dem_dir <- paste0(out_dir, "/temp.flt")
  #   writeRaster(r, dem_dir, overwrite = T)
  #   print("Created .flt file")
  # }
  
  # Prepare inputs
  dem_dir <- normalizePath(dem_dir)
  # dem_dir <- substr(dem_dir, 1, nchar(dem_dir)-4)
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
