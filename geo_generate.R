library(tidyverse)
library(sn)
library(png)
library(colorspace)
library(reshape2)
library(cetcolor)
library(caret)

MAP_DIR = 'img_maps'

PEAK_AMPL_DEFAULT_FN = function(n) rgamma(n,3,1/20)
# relative to amplitude
PEAK_RADIUS_DEFAULT_FN = function(n,h=1)  h*(1.25+rgamma(n, 2, 1/2))

# determines relative scale of x-y dimensions (pre-rotation)
PEAK_STRETCH_DEFAULT_FN = function(n) exp((rbeta(n, 0.8, 0.8)-0.5)*2)

#skewness (positive = right-tailed, negative = left-tailed)
PEAK_SKEWNESS_DEFAULT_FN = function(n) rnorm(n) * 2

#distance between features

# INCOMPLETE
generate_geo_grid <- function(
  lat,
  lon,
  n_peaks=0,
  n_valleys=0,
  n_saddles=0,
  peak_ampl_distro = PEAK_AMPL_DEFAULT_FN,
  peak_radius_distro =PEAK_RADIUS_DEFAULT_FN,
  peak_stretch_distro = PEAK_STRETCH_DEFAULT_FN
){

  grid = matrix(0,ncol=lon,nrow=lat)
  # generate peaks, valleys, and saddles (saddles will be mix of 4 peaks/
  # valleys near each other)
  n_feature_centers = n_peaks + n_valleys + n_saddles


  rotations = 2*pi(runif(n_feature_centers))

  peak_ampls = peak_ampl_distro(n_peaks)
  valley_ampls = -peak_ampl_distro(n_valleys)
  saddle_ampls = peak_ampl_distro(n_saddles * 4) * (1-2*rbinom(4 * n_saddles,1,0.5))

}

## THIS NEEDS TESTING
build_matrix_from_posterior <- function(x,name,column=ncol(x),debug=F){
  rn = rownames(x)
  rrn_pattern = sprintf('^%s\\[(\\d+),(\\d+)\\]',name)
  rrn = rn[grepl(rrn_pattern,rn)]
  if (column == -1){
    column = ncol(x)
  }
  if (debug){
    print(rrn_pattern)
    print(length(rrn))
    print(unique(gsub('[0-9]','',rrn)))
  }
  indices = (str_match_all(rrn,rrn_pattern) %>% do.call(what=rbind))[,-1] %>%
    apply(MARGIN=2,as.numeric)

  i1_range = range(indices[,1])
  i2_range = range(indices[,2])

  mat = matrix(x[rrn,column],nrow=i1_range[2],ncol=i2_range[2],byrow=TRUE)
  mat
}

set_0_to_na <- function(x){
  x = ifelse(x==0,NA,x)
}


# COMPLETE
load_map <- function(idx=1,plot=FALSE, scale=100){
  if (is.numeric(idx)){
    X = readPNG(file.path(MAP_DIR,sprintf('map%s.png',idx))) * scale
  } else {
    X = readPNG(idx) * scale
  }
  avg = (X[,,1] + X[,,2] + X[,,3])/3
  if (plot){
    df = reshape2::melt(avg) %>% transmute(lat=rev(Var1),lon=Var2,z=value)
    ggplot(df)+geom_raster(aes(y=lat,x=lon,fill=z)) +
      scale_fill_gradientn(colors = cet_pal(7,'l15'))
  }
  X
}

# diagnostic function
plot_matrix <- function(mat, pal=cet_pal(7,'inferno')){
  require(cetcolor)
  mm = reshape2::melt(mat)
  ggplot(mm) +
    geom_raster(aes(x=Var2,y=-Var1,fill=value)) +
    scale_fill_gradientn(colors=pal) +
    theme_bw()
}




# COMPLETE
row_grid_walk <- function(
    map,
    y_breaks,
    x_breaks,
    y_first=TRUE,
    start_sides = c(-1, -1),
    margins = list(x=c(0,0),y=c(0,0)),
    time_drift_coefs=c(0,0),
    debug=FALSE,
    return_instructions=FALSE
){
  ydim = nrow(map)
  xdim = ncol(map)

  nyb = length(y_breaks)
  nxb = length(x_breaks)

  # calculate # of x and y positions
  # -1 = start at 1; 1 = start at end

  instruction_list = list()
  idx = 1

  printif <- function(x) if(debug) print(x)
  sprintif <- function(x,...) if (debug) (print(sprintf(x,...)))

  if (y_first){
    printif('Starting on y')
    side = start_sides[1]
    x = y_breaks[1]
    y = ifelse(side==-1,1,ydim)
    yi = 1
    while (yi < nyb){
      printif(sprintf('start of yi=%s < nyb loop',yi))
      # vertical sweep
      if (side == -1){
        # starting on top
        printif('start on top')
        x2 = x
        y2 = ydim-margins$y[2]
        instruction_list[[idx]] = data.frame(x=x,y=y:y2)
        idx = idx + 1
        y = y2
        side = 1
        printif(instruction_list[[length(instruction_list)]])

      } else {
        # starting on bottom
        printif('start on bottom')
        x2 = x
        y2 = 1+margins$y[2]
        instruction_list[[idx]] = data.frame(x=x,y=y:y2)
        idx = idx + 1
        y = y2
        side = -1
        printif(instruction_list[[length(instruction_list)]])
      }

      # move along margin if not last break

      x2 = y_breaks[yi + 1]

      direction = sign(x2-x)
      adiff = abs(x2-x)

      # horizontal shift
      last_shift_offset = ifelse(yi == nyb -1, 1, 0)
      if (direction == 1 & adiff > 1){
        printif('move to right')
        # move to right
        instruction_list[[idx]] = data.frame(x=(x+1):(x2-1 + last_shift_offset),y=y)
        idx = idx + 1
        printif(instruction_list[[length(instruction_list)]])

      } else if (direction == -1 & adiff > 1){
        printif('move to left')
        # move to left
        instruction_list[[idx]] = data.frame(x=(x-1):(x2+1 - last_shift_offset),y=y)
        idx = idx + 1
        printif(instruction_list[[length(instruction_list)]])
      } else {
        printif('no extra side movement')
        # moving back-and-forth or between adjacent columns
      }
      x = x2
      yi = yi + 1
    }
    # last move before swapping to X movement
    printif('last column')
    if (start_sides[2] == -1){
      printif('need to move to left')
      x2 = 1 + margins$x[1]
      direction = -1
      side = -1
    } else {
      printif('need to move to right')
      x2 = xdim - margins$x[2]
      direction = 1
      side = 1
      #printif(instruction_list[[length(instruction_list)]])
    }
    if (x != x2){
      # draw to X margin
      printif(sprintf('drawing to margin on idx=%s, from %s to %s',idx,x+direction,x2))
      instruction_list[[idx]] = data.frame(x=(x+direction):x2,y=y)
      idx = idx + 1
      printif(instruction_list[[length(instruction_list)]])
      x = x2
    } else {
      printif('Already at X margin, skipping step')
    }

    # need to shift y value to correct starting spot
    if (abs(y-x_breaks[1]) > 1){
      printif('Need to shift y value for first x loop')
      y2 = x_breaks[1]
      if (y2 < y){
        printif('shift y up')
        direction = -1
      } else {
        printif('shift y down')
        direction = 1
      }
      instruction_list[[idx]] = data.frame(x=x2, y=(y+direction):(y2-direction))
      idx = idx + 1
      printif(instruction_list[[length(instruction_list)]])
    }

    # x loop
    xi = 1
    #direction = direction * -1

    printif('starting x loops')
    while (xi < nxb){
      printif(sprintf('starting xi=%s < nxb loop',xi))
      y = x_breaks[xi]#y2
      #y2 = x_breaks[xi]

      # horizontal sweep
      if (side == -1){
        sprintif('on left side, idx=%s',idx)
        # starting on left
        y2 = y
        x2 = xdim-margins$x[2]
        instruction_list[[idx]] = data.frame(y=y,x=x:x2)
        idx = idx + 1
        x = x2
        side = 1
        printif(instruction_list[[length(instruction_list)]])

      } else {
        sprintif('on right side, idx=%s',idx)
        # starting on right
        y2 = y
        x2 = 1+margins$x[2]
        instruction_list[[idx]] = data.frame(y=y,x=x:x2)
        idx = idx + 1
        x = x2
        side = -1
        printif(instruction_list[[length(instruction_list)]])
      }

      y2 = x_breaks[xi + 1]

      direction = sign(y2-y)
      adiff = abs(y2-y)
      # vertical shift
      if (direction == 1 & adiff > 1){
        # move downwards
        sprintif('moving down, idx=%s',idx)
        instruction_list[[idx]] = data.frame(y=(y+1):(y2-1),x=x)
        idx = idx + 1
        printif(instruction_list[[length(instruction_list)]])

      } else if (direction == -1 & adiff > 1){
        # move upwards
        sprintif('moving up, idx=%s', idx)
        instruction_list[[idx]] = data.frame(y=(y-1):(y2+1),x=x)
        idx = idx + 1
        printif(instruction_list[[length(instruction_list)]])
      } else {
        # moving back-and-forth or between adjacent rows
        printif('no vertical movement')
      }
      y = y2
      xi = xi + 1
    }

    # last one to exit
    sprintif('last row, idx=%s',idx)
    y = y2
    y2 = x_breaks[xi]

    if (nxb %% 2 == 0){
      end_side = start_sides[2]
      printif('ends on same side')
    } else {
      end_side = -start_sides[2]
      printif('ends on opposite side')
    }

    if (end_side == 1){
      # to right
      printif('ends on right')
      y2 = ydim
      direction = 1
      side = 1
    } else {
      printif('ends on left')
      y2 = 1
      direction = -1
      side = -1
    }
    printif(sprintf('drawing to edge on idx=%s',idx))
    instruction_list[[idx]] = data.frame(y=y:y2,x=x)
    printif(instruction_list[[length(instruction_list)]])


  } else {
    # let's just call the function recursively and then rotate the result
    # it's not the bottleneck of this exercise
    printif('calling transposed function')
    #stop('Need to test y-first mode')
    transposed = row_grid_walk(
      t(map),
      y_breaks=x_breaks,
      x_breaks=y_breaks,
      y_first=TRUE,
      start_sides = rev(start_sides),
      margins = list(x=margins$y,y=margins$x),
      time_drift_coefs=time_drift_coefs,
      debug=debug,
      return_instructions=TRUE
    )

    instruction_list = transposed

  }

  # todo: build matrices and return those instead
  # todo: build symmetric x-first mode

  if (return_instructions){
    printif('returning basic instructions')
   return(instruction_list)
  }

  # build full data frame
  printif('building full instructions list')
  full_instructions = do.call(rbind,instruction_list) %>%
    mutate(
      step = 1:n()
    )

  if (!y_first){
    printif('transposing full instructions')
    full_instructions = full_instructions %>% mutate(x=full_instructions$y,y=full_instructions$x)
  }

  N = nrow(full_instructions)
  if (prod(time_drift_coefs) != 0){
    printif('building time drift')
    time_drift = time_drift_coefs[1] * arima.sim(list(ar=time_drift_coefs[2]),n=N)
  } else {
    printif('skipping time drift')
    time_drift = 0
  }

  printif('assigning map values to full instructions')
  full_instructions = full_instructions %>%
    rowwise() %>%
    mutate(
      raw_value = map[y,x]
    ) %>%
    ungroup() %>%
    mutate(
      value = raw_value + time_drift
    )

  mask = 0 * map
  printif('building mask')
  for (i in 1:nrow(full_instructions)){
    x = full_instructions$x[i]
    y = full_instructions$y[i]
    mask[y,x] = 1
  }

  printif('returning values')
  return(
    list(
      map=map,
      instruction_list = instruction_list,
      full_instructions = full_instructions,
      time_drift = time_drift,
      mask=mask,
      masked_map = mask * map,
      inputs = list(
        margins=margins,
        y_first=y_first,
        start_sides=start_sides,
        time_drift_coefs=time_drift_coefs,
        x_breaks=x_breaks,
        y_breaks=y_breaks
      )
    )
  )
}


if (FALSE){
  test_map = matrix(1:300, nrow=20,ncol=15)
  rgw = row_grid_walk(
    test_map,
    c(4,8,7,14),
    c(3,20,18,1,14),
    start_sides = c(-1,1),
    margins=list(x=c(1,1),y=c(1,1)),
    debug=TRUE,
    y_first=FALSE,
    time_drift_coefs = c(1,0.995)
  )
}

# returns list of lists
# uses character indices, starting at '0'
# 'y' (vertical) is considered x_1, 'x' (horizontal) is considered x_2
# first element is the object itself (0th derivative)
build_derivatives <- function(mat,max_degree=1){
  # note: y gets 1st index (row #), x gets 2nd (col #)
  nr = nrow(mat)
  nc = ncol(mat)
  derivatives_list = list()
  derivatives_list[['0']] = list('0'=mat)
  # define bands
  ym = 2:(nr-1)
  xm = 2:(nc-1)

  # loop
  for (D in 1:max_degree){
    cD = as.character(D)
    cpD = as.character(D-1)
    current_derivatives = list()
    for (i in 0:D){
      ci = as.character(i)
      cpi = as.character(i-1)
      n_down = i
      n_up = D-i
      # first: y-shift first of previous
      # center: y-shift from equivalent of previous index
      # last: x-shift last of previous

      # middle band: average up/down (y) or left/right (x)
      # edges: 1-sided slope
      if (i != D){
        dmat = matrix(NA, nrow=nr,ncol=nc)
        # first: y-shift first of previous
        pmat = derivatives_list[[cpD]][[ci]]
        # top
        dmat[1,] = pmat[2,]-pmat[1,]
        # bot
        dmat[nr,] = pmat[nr,]-pmat[nr-1,]
        # middle bands
        dmat[ym,] = (pmat[ym+1,] - pmat[ym-1,])/2
      } else {
        dmat = matrix(NA, nrow=nr,ncol=nc)
        # first: x-shift last of previous
        pmat = derivatives_list[[cpD]][[cpi]]
        # top
        dmat[,1] = pmat[,2]-pmat[,1]
        # bot
        dmat[,nc] = pmat[,nc]-pmat[,nc-1]
        # middle bands
        dmat[,xm] = (pmat[,xm+1] - pmat[,xm-1])/2
      }
      current_derivatives[[ci]] = dmat
    }
    derivatives_list[[cD]] = current_derivatives
  }
  return(derivatives_list)
}


create_knn_init <- function(data, n_chains=1, noise_func=function(n) 0, noise_sigma=1/100){
  if (n_chains > 1){
    return(lapply(
      1:n_chains,
      create_knn_init,
      data=data,
      n_chains=1,
      noise_func = function(n) rnorm(n) * noise_sigma
    ))
  }
  data = path_data
  ymax = nrow(data$map)
  xmax = ncol(data$map)
  df = data$full_instructions
  # x, y, value
  knn_model = knnreg(value ~ x + y, data=df)

  full_df = expand.grid(x=1:xmax,y=1:ymax)

  full_df$value = predict(knn_model, full_df)

  new_mat = reshape2::acast(full_df, y ~ x, fun.aggregate=mean)
  data_mask = data$mask
  data_masked_map = data$masked_map

  new_mat = (new_mat+noise_func(ymax*xmax)) * (1-data_mask) + data_masked_map
  # estimate derivative matrix
  ym1 = 2:(ymax-1)
  xm1 = 2:(xmax-1)

  map_xd_init = matrix(0,nrow=ymax,ncol=xmax)
  # left edge
  map_xd_init[,1] = new_mat[,2] - new_mat[,1]
  # right edge
  map_xd_init[,xmax] = new_mat[,xmax] - new_mat[,xmax-1]
  # everything else
  map_xd_init[,xm1] = 1/2 * (new_mat[,xm1+1] - new_mat[,xm1-1])

  map_yd_init = matrix(0,nrow=ymax,ncol=xmax)
  # top edge
  map_yd_init[2,] = new_mat[2,] - new_mat[1,]
  # bot edge
  map_yd_init[ymax,] = new_mat[ymax,] - new_mat[ymax-1,]
  # everything else
  map_yd_init[ym1,] = 1/2 * (new_mat[ym1+1,] - new_mat[ym1-1,])

  # use y derivative in x direction
  map_xyd_init = 0 * map_yd_init
  # left edge
  map_xyd_init[,1] = map_yd_init[,2] - map_yd_init[,1]
  # right edge
  map_xyd_init[,xmax] = map_yd_init[,xmax] - map_yd_init[,xmax-1]
  # everything else
  map_xyd_init[,xm1] = 1/2 * (map_yd_init[,xm1+1] - map_yd_init[,xm1-1])


  list(
    map = new_mat,
    map_1d_y = map_yd_init,
    map_1d_x = map_xd_init,
    map_2d_xy = map_xyd_init
  )
}



# INCOMPLETE
# (x,y) => (-y, x) (when comparing to normal cartesian coordiantes)
generate_walked_path_on_grid <- function(map, n_steps=3000, time_drift_coef=0, start=c(1,1),
                                         step_controller = DEFAULT_STEP_CONTROLLER, step_scale=1){
  if (time_drift_coef > 0){
    time_error = arima.sim(list(ar=time_drift_coef),n=n_steps)
  } else{
    time_error = rep(0,n_steps)
  }

  ymax = nrow(map)
  xmax = ncol(map)

  #
  records = matrix(ncol=4,nrow=n_steps)

  records[1,] = c(1,start,time_error[1] + map[start[1],start[2]])
  for (i in 2:n_steps){

  }
}
