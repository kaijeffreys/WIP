# Development surface metrics function
#
# Contains gradient, planar curvature, profile curvature, local relief, and twi
#
# For gradient, it is calculated in two directions, north-south and east-west
#
# For types of curvature, the function iterates through each cell, and using
# nine points surrounding the point (including the point itself), it fits a
# quadratic model y = ax^2 + by^2 + cxy + dx + ey, and uses the coefficients # from that model to calculate the types of curvature

surface_met_dev <- function(DEM, len, export = FALSE,
                            met = c("grad", "plan", "prof", "dev", "twi")) {
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
  
  # Calculation of gradient
  if("grad" %in% met) {
    print("Running gradient")
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
  
  # Calculation of curvature
  if("prof" %in% met | "plan" %in% met) {
    print("Running curvature")
    
    # Extracting values of elevation as well as the x,y coordinates
    foc_vals <- terra::focalValues(dem, w = w_mat)
    foc_x <- terra::focalValues(init(dem, "x"), w = w_mat)
    foc_y <- terra::focalValues(init(dem, "y"), w = w_mat)
    foc_x2 <- foc_x ^ 2
    foc_y2 <- foc_y ^ 2
    foc_xy <- foc_x * foc_y
    
    # Set number of iterations and initialize vectors of values
    n_iter <- nrow(foc_vals)
    
    if("plan" %in% met) {
      plan_vec <- rep(NA, n_iter)
    }
    
    if("prof" %in% met) {
      prof_vec <- c(NA, n_iter)
    }
    
    # Condition to check if there are enough focal values to run a model,
    # speeds up process
    cond <- rowSums(!is.na(foc_vals)) > 5
    
    # Set up progress bar
    pb <- txtProgressBar(max = n_iter, width = 50, style = 3)
    
    # Iterate through each cell
    for(i in (1:n_iter)[cond]) {
      
      # Fits a quadratic model for the cell using the points specified in
      # moving window
      mod <- mgcv::bam(elev ~ x2 + y2 + xy + x + y,
                         data = list(elev = foc_vals[i,],x = foc_x[i,],
                                     y = foc_y[i,],x2 = foc_x2[i,],
                                     y2 = foc_y2[i,],xy = foc_xy[i,]),
                         drop.intercept = T, method = "GCV.Cp")
      
      # Extract the coefficients  
      coeff <- mod$coefficients
      
      # Calculate curvature values using the coefficients
      if("plan" %in% met) {
          plan_vec[i] <- -2*(coeff[[1]]*(coeff[[5]]^2) - coeff[[3]]*coeff[[4]]*coeff[[5]] + coeff[[2]]*coeff[[4]]^2)/((coeff[[4]]^2+coeff[[5]]^2) * sqrt(1+coeff[[4]]^2+coeff[[5]]^2))
      }
        
      if("prof" %in% met) {
          prof_vec[i] <- (-2 * (coeff[[1]]*coeff[[4]]^2 + coeff[[3]]*coeff[[4]]*coeff[[5]] + coeff[[2]]*coeff[[5]]^2)) / ((coeff[[4]]^2 + coeff[[5]]^2)*(1 + coeff[[4]]^2 + coeff[[5]]^2)^1.5)
      }
      
      # Update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Set the values of the output rasters for curvature
    if("plan" %in% met) {
      plan <- terra::setValues(DEM, plan_vec)
      in_rast <- c(in_rast, plan)
      names(in_rast)[length(in_rast)] <- paste0("plan", len)
    }
    
    if("prof" %in% met) {
      prof <- terra::setValues(DEM, prof_vec)
      in_rast <- c(in_rast, prof)
      names(in_rast)[length(in_rast)] <- paste0("prof", len)
    }
  }
  
  # Calculation of local relief
  if("dev" %in% met) {
    print("Running Local Relief")
    dev <- (DEM - focal(DEM, w = w_mat, fun = "mean", na.rm = T, na.policy = "omit")) / focal(DEM, w = w_mat, fun = "sd", na.rm = T, na.policy = "omit") 
    in_rast <- c(in_rast, dev)
    
    names(in_rast)[length(in_rast)] <- paste0("dev", len)
  }
  
  # Calculation of TWI
  if("twi" %in% met) {
    print("Running TWI")
    topidx <- topmodel::topidx(terra::as.matrix(DEM), res = terra::res(DEM)[1])
    twi <- terra::setValues(DEM, topidx$atb)
    terra::values(twi) <- ifelse(terra::values(twi) < 0, 0, terra::values(twi))
    twi <- terra::focal(twi, w = w_mat, mean, na.rm = T, na.policy = "all")
    
    in_rast <- c(in_rast, twi)
    names(in_rast)[length(in_rast)] <- paste0("twi", len)
  }
  return(in_rast)
}
