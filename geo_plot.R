library(tidyverse)
library(cetcolor)
library(rayshader)
library(gganimate)
library(gridExtra)
library(RColorBrewer)
#require(reshape2)
source('geo_generate.R')

# testing
if (FALSE){
  STAN_GEO_DIR = 'stan_geo_files'
  STAN_GEO_PLOT_DIR = 'geo_plots'
  data1 = readRDS(file.path(STAN_GEO_DIR,'geo_estimate_01_v3.RDS'))
  data2 = readRDS(file.path(STAN_GEO_DIR,'geo_estimate_02_02_v3.RDS'))

  name='TEST'
  name_nice = 'test, but nice'
  name2=''
  directory=file.path(STAN_GEO_PLOT_DIR,'test')

  width_default=840
  height_default=840
  width_3d=4
  height_3d=4
  res=200
  scale_3d=310
  frames_3d=720
  fps_3d=30
  zoom_3d=0.6
  fov_3d=30
  topo_pal = colorRampPalette(rev(brewer.pal(11, "Spectral")))
}

set_0_to_na <- function(x){
  x = ifelse(x==0,NA,x)
}

plot_matrix <- function(
    mat,
    pal=cet_pal(7,'inferno'),
    zrange=range(mat),
    zname='Height'
){
  mm = reshape2::melt(mat)
  ggplot(mm) +
    geom_raster(aes(x=Var2,y=-Var1,fill=value)) +
    scale_fill_gradientn(zname, colors=pal,limits=zrange) +
    theme_bw()
}

# 2-color gradient for easier reading
plot_matrix_centered <- function(
    mat,
    zrange=range(mat),
    zname='Height',
    ...
){
  mm = reshape2::melt(mat)
  ggplot(mm) +
    geom_raster(aes(x=Var2,y=-Var1,fill=value)) +
    scale_fill_gradient2(zname, midpoint=0, limits=zrange, ...) +
    theme_bw()
}

# removes more labels
plot_matrix2 <- function(...,keep_legend=1){
  if (keep_legend==1){
    p <- plot_matrix(...) +
      xlab('') + ylab('') +
      theme(legend.title=element_blank())
  } else if (keep_legend==2) {
    p <- plot_matrix(...) +
      xlab('') + ylab('')
  } else {
    p <- plot_matrix(...) +
      xlab('') + ylab('') +
      theme(legend.position='none')
  }
  return(p)
}

# 2-color graadient for easier comparing
plot_matrix2c <- function(...,keep_legend=1){
  if (keep_legend==1){
    p <- plot_matrix_centered(...) +
      xlab('') + ylab('') +
      theme(legend.title=element_blank())
  } else if (keep_legend == 2) {
    p <- plot_matrix_centered(...) +
      xlab('') + ylab('')
  } else {
    p <- plot_matrix_centered(...) +
      xlab('') + ylab('') +
      theme(legend.position='none')
  }
  return(p)
}

plot_compare <- function(
  data1,
  data2,
  path_data,
  name,
  directory,
  width_default=840,
  height_default=840,
  grid_factor=1.5,
  grid_res_factor=0.5,
  ...
){
  if (is.character(data1)){
    data1 = readRDS(data1)
  }
  if (is.character(data2)){
    data2 = readRDS(data2)
  }

  dir.create(directory,recursive=T,showWarnings=FALSE)

  # consts
  topo_pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  inferno_pal = cet_pal(7,'inferno')

  # unpack data
  fit1=data1$fit
  params1=data1$stan_params
  knn_init1=data1$knn_init
  control1=data1$stan_ctrl
  time1 = data1$time

  fit2=data2$fit
  params2=data2$stan_params
  knn_init2=data2$knn_init
  control2=data2$stan_ctrl
  time2 = data2$time

  map_derivatives = build_derivatives(path_data$map,2)

  estimate1 = get_posterior_mean(fit1)
  smooth_grid1 = build_matrix_from_posterior(estimate1,'granular_map')
  smooth_grid_x1 = build_matrix_from_posterior(estimate1,'granular_map_dx')
  smooth_grid_y1 = build_matrix_from_posterior(estimate1,'granular_map_dy')
  smooth_grid_xy1 = build_matrix_from_posterior(estimate1,'granular_map_dxy')

  col = 1
  estimate_mat1 = build_matrix_from_posterior(estimate1,'map')
  estimate_mat_y1 = build_matrix_from_posterior(estimate1,'map_1d_y')
  estimate_mat_x1 = build_matrix_from_posterior(estimate1,'map_1d_x')
  estimate_mat_xy1 = build_matrix_from_posterior(estimate1,'map_2d_xy')

  estimate2 = get_posterior_mean(fit2)
  smooth_grid2 = build_matrix_from_posterior(estimate2,'granular_map')
  smooth_grid_x2 = build_matrix_from_posterior(estimate2,'granular_map_dx')
  smooth_grid_y2 = build_matrix_from_posterior(estimate2,'granular_map_dy')
  smooth_grid_xy2 = build_matrix_from_posterior(estimate2,'granular_map_dxy')

  col = 1
  estimate_mat2 = build_matrix_from_posterior(estimate2,'map')
  estimate_mat_y2 = build_matrix_from_posterior(estimate2,'map_1d_y')
  estimate_mat_x2 = build_matrix_from_posterior(estimate2,'map_1d_x')
  estimate_mat_xy2 = build_matrix_from_posterior(estimate2,'map_2d_xy')


  #
  granularity = (dim(smooth_grid1)[1]-1)/(dim(estimate_mat1)[1]-1)
  granularity_trans <- function(x) (x-1)/granularity + 1
  granularity_trans2 <- function(x) -((x+1)/granularity + 1)

  ppng <- function(x,...){
    png(file.path(directory, sprintf('%s_%s.png',x,name)), width=width_default, height=height_default, res=res,...)
  }

  gppng <- function(x,...){
    png(file.path(directory, sprintf('%s_%s.png',x,name)), width=grid_factor*width_default, height=grid_factor * height_default, res=res*grid_res_factor,...)
  }

  # smooth comparisons
  ppng('plot_comparison')
  plt_cmp <- plot_matrix2c(smooth_grid2-smooth_grid1,keep_legend=2,zname='ΔHeight') +
    ggtitle('Difference between models', subtitle='2nd & 3rd-derivative smoothed model and \n1st-derivative-only model') +
    scale_y_continuous(label=granularity_trans2) +
    scale_x_continuous(label=granularity_trans)
  print(
    plt_cmp
  )
  dev.off()

  ppng('plot_comparison_dx')
  plt_cmp_x <- plot_matrix2c(smooth_grid_x2-smooth_grid_x1,keep_legend=2,zname=bquote(paste(delta,f[x]))) +
    ggtitle(bquote(paste(f[x],' difference between models')), subtitle='2nd & 3rd-derivative smoothed model and \n1st-derivative-only model') +
    scale_y_continuous(label=granularity_trans2) +
    scale_x_continuous(label=granularity_trans)
  print(
    plt_cmp_x
  )
  dev.off()

  ppng('plot_comparison_dy')
  plt_cmp_y <-   plot_matrix2c(smooth_grid_y2-smooth_grid_y1,keep_legend=2,zname=bquote(paste(delta,f[y]))) +
    ggtitle(bquote(paste(f[y],' difference between models')), subtitle='2nd & 3rd-derivative smoothed model and \n1st-derivative-only model') +
    scale_y_continuous(label=granularity_trans2) +
    scale_x_continuous(label=granularity_trans)

  print(
    plt_cmp_y
  )
  dev.off()

  ppng('plot_comparison_dxy')
  plt_cmp_xy <-   plot_matrix2c(smooth_grid_xy2-smooth_grid_xy1,keep_legend=2,zname=bquote(paste(delta,f[xy]))) +
    ggtitle(bquote(paste(f[xy],' difference between models')), subtitle='2nd & 3rd-derivative smoothed model and \n1st-derivative-only model') +
    scale_y_continuous(label=granularity_trans2) +
    scale_x_continuous(label=granularity_trans)
  print(
    plt_cmp_xy
  )
  dev.off()

  # "improvements"
  ppng('improvement')
  plt_imp <-   plot_matrix2c(
    abs(path_data$map-estimate_mat1) - abs(path_data$map-estimate_mat2),
    keep_legend=2,
    zname=bquote(paste(Delta,'|',epsilon[f],'|'))
  ) +
    ggtitle('Improvement in estimates of models',subtitle='2nd and 3rd-derivative-smoothed model vs. \n1st-derivative-only model') +
    scale_y_continuous(label=function(x)-x)
  print(
    plt_imp
  )
  dev.off()

  ppng('improvement_dx')
  plt_imp_x <-   plot_matrix2c(
    abs(map_derivatives[['1']][['1']]-estimate_mat_x1) - abs(map_derivatives[['1']][['1']]-estimate_mat_x2),
    keep_legend=2,
    zname=bquote(paste(Delta,'|',epsilon[f[x]],'|'))
  ) +
    ggtitle(bquote(paste('Improvement in ',f[x],' estimates of models')),
            subtitle='2nd and 3rd-derivative-smoothed model vs. \n1st-derivative-only model') +
    scale_y_continuous(label=function(x)-x)
  print(
    plt_imp_x
  )
  dev.off()

  ppng('improvement_dy')
  plt_imp_y <-   plot_matrix2c(
    abs(map_derivatives[['1']][['0']]-estimate_mat_y1) - abs(map_derivatives[['1']][['0']]-estimate_mat_y2),
    keep_legend=2,
    zname=bquote(paste(Delta,'|',epsilon[f[y]],'|'))
  ) +
    ggtitle(bquote(paste('Improvement in ',f[y],' estimates of models')),
            subtitle='2nd and 3rd-derivative-smoothed model vs. \n1st-derivative-only model') +
    scale_y_continuous(label=function(x)-x)
  print(
    plt_imp_y
  )
  dev.off()

  ppng('improvement_dxy')
  plt_imp_xy <-   plot_matrix2c(
    abs(map_derivatives[['2']][['1']]-estimate_mat_xy1) - abs(map_derivatives[['2']][['1']]-estimate_mat_xy2),
    keep_legend=2,
    zname=bquote(paste(Delta,'|',epsilon[f[xy]],'|'))
  ) +
    ggtitle(bquote(paste('Improvement in ',f[xy],' estimates of models')),
            subtitle='2nd and 3rd-derivative-smoothed model vs. \n1st-derivative-only model') +
    scale_y_continuous(label=function(x)-x)
  print(
    plt_imp_xy
  )
  dev.off()

  # grid arrange

  gppng('comparison_grid')
  print(
    grid.arrange(plt_cmp,plt_cmp_x,plt_cmp_y,plt_cmp_xy)
  )
  dev.off()

  gppng('improvement_grid')
  print(
    grid.arrange(plt_imp,plt_imp_x,plt_imp_y,plt_imp_xy)
  )
  dev.off()


}

plot_geo <- function(
    data,
    name, # describes algorithm
    name_nice, # for human reading
    name2='', # describes map
    directory='.',
    width_default=840,
    height_default=840,
    width_3d=4,
    height_3d=4,
    res=200,
    scale_3d=310,
    frames_3d=720,
    fps_3d=30,
    zoom_3d=0.6,
    fov_3d=30,
    topo_pal = colorRampPalette(rev(brewer.pal(11, "Spectral"))),
    return_plots=TRUE,
    skip3d=FALSE,
    ...
){
  name2_pad = ifelse(name2=='','',sprintf('(%s)',name2))

  if (is.character(data)){
    data = readRDS(data)
  }

  ppng <- function(x,...){
    png(file.path(directory, sprintf('%s_%s.png',x,name)), width=width_default, height=height_default, res=res,...)
  }

  pnpng <- function(x,...){
    png(file.path(directory, sprintf('%s_%s_%s.png',x,name,name2)), width=width_default, height=height_default, res=res,...)
  }

  # consts
  #topo_pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  inferno_pal = cet_pal(7,'inferno')

  # unpack data
  fit=data$fit
  params=data$stan_params
  knn_init=data$knn_init
  control=data$stan_ctrl
  time = data$time

  # expand values

  estimate = get_posterior_mean(fit)
  smooth_grid = build_matrix_from_posterior(estimate,'granular_map')
  smooth_grid_x = build_matrix_from_posterior(estimate,'granular_map_dx')
  smooth_grid_y = build_matrix_from_posterior(estimate,'granular_map_dy')
  smooth_grid_xy = build_matrix_from_posterior(estimate,'granular_map_dxy')

  col = 1
  estimate_mat = build_matrix_from_posterior(estimate,'map')
  estimate_mat_y = build_matrix_from_posterior(estimate,'map_1d_y')
  estimate_mat_x = build_matrix_from_posterior(estimate,'map_1d_x')
  estimate_mat_xy = build_matrix_from_posterior(estimate,'map_2d_xy')

  map = path_data$map
  map_derivatives = build_derivatives(map,2)
  map_dx = map_derivatives[['1']][['1']]
  map_dy = map_derivatives[['1']][['0']]
  map_dxy = map_derivatives[['2']][['0']]

  eprox = abs(estimate_mat-path_data$map)
  kprox = abs(knn_init[[col]]$map - path_data$map)
  proxdiff = kprox-eprox

  granularity = (dim(smooth_grid)[1]-1)/(dim(estimate_mat)[1]-1)
  granularity_trans <- function(x) (x-1)/granularity + 1
  granularity_trans2 <- function(x) -((x+1)/granularity + 1)

  # global z scale
  get_range <- function(...){
    do.call(range,list(na.rm=T,...))
  }
  widen_range <- function(x,f=1.04){
    d = x[2] - x[1]
    m = mean(x)
    c(m-d * f/2,m+d*f/2)
  }
  cat_sublists <- function(x,subname){
    do.call(rbind, lapply(x,function(x)x[[subname]]))
  }
  z_range = get_range(smooth_grid,estimate_mat,path_data$map,cat_sublists(knn_init,'map')) %>% widen_range()
  dx_range = get_range(smooth_grid_x, estimate_mat_x,cat_sublists(knn_init,'map_1d_x'),map_dx) %>% widen_range()
  dy_range = get_range(smooth_grid_y,estimate_mat_y,cat_sublists(knn_init,'map_1d_y'),map_dy) %>% widen_range()
  dxy_range = get_range(smooth_grid_xy,estimate_mat_xy,cat_sublists(knn_init,'map_2d_xy'),map_dxy) %>% widen_range()

  # expand plot functions
  plot_matrix_z <- function(...) plot_matrix2(zrange=z_range,pal=topo_pal(4),...)
  plot_matrix_dx <- function(...) plot_matrix2c(zrange=dx_range,zname=bquote(f[x]),...)
  plot_matrix_dy <- function(...) plot_matrix2c(zrange=dy_range,zname=bquote(f[y]),...)
  plot_matrix_dxy <- function(...) plot_matrix2c(zrange=dxy_range,zname=bquote(f[xy]),high=scales::muted('green'),low=scales::muted('orange'),...)

  # create dir
  dir.create(directory,showWarnings=FALSE,recursive=TRUE)

  sample_plots <- list(
    plot_matrix(estimate_mat) + ggtitle("ESTIMATE"),
    plot_matrix(path_data$map) + ggtitle("ACTUAL"),
    plot_matrix(estimate_mat-path_data$map) + ggtitle("DIFFERENCE"),
    plot_matrix(set_0_to_na(path_data$masked_map)) + ggtitle("MASKED"),
    plot_matrix(knn_init[[col]]$map) + ggtitle("KNN IMPUTE"),
    plot_matrix(estimate_mat - set_0_to_na(path_data$masked_map)) + ggtitle('Estimate Minus Masked'),
    plot_matrix(estimate_mat - knn_init[[col]]$map) + ggtitle('ESTIMATE vs. KNN INIT DIFF'),
    plot_matrix(proxdiff) + ggtitle('Improvement of Spline vs. KNN')
  )

  pnpng('sample_plots')
  do.call(grid.arrange,sample_plots)
  dev.off()

  # unused?
  nm_app <- function(x){
    if (name != ''){
      paste(x,' - ',name)
    } else {
      x
    }
  }

  # plot actual
  ppng('actual')
  (plt_actual <- plot_matrix_z(path_data$map,keep_legend=2) +
    ggtitle('Full Map') + scale_y_continuous(label=function(x)-x)
  )
  print(plt_actual)
  dev.off()

  ppng('actual_dx')
  (plt_actual_dx <- plot_matrix_dx(map_dx,keep_legend=2) +
      ggtitle(bquote(paste('Map, ',f[x]))) + scale_y_continuous(label=function(x)-x)
  )
  print(plt_actual_dx)
  dev.off()

  ppng('actual_dy')
  (plt_actual_dy <- plot_matrix_dy(map_dy,keep_legend=2) +
      ggtitle(bquote(paste('Map, ',f[y]))) + scale_y_continuous(label=function(x)-x)
  )
  print(plt_actual_dy)
  dev.off()

  ppng('actual_dxy')
  (plt_actual_dxy <- plot_matrix_dxy(map_dxy,keep_legend=2) +
      ggtitle(bquote(paste('Map, ',f[xy]))) + scale_y_continuous(label=function(x)-x)
  )
  print(plt_actual_dxy)
  dev.off()


  # plot path
  ppng('path')
  (plt_path<-plot_matrix_z(set_0_to_na(path_data$masked_map),keep_legend=2) +
    ggtitle('Sampled Map')
  )
  print(plt_path)
  dev.off()

  # plot knn estimate
  ppng('knn')
  (plt_knn <- plot_matrix_z(knn_init[[col]]$map,keep_legend=2) +
    ggtitle(paste('KNN Estimate'),subtitle=name_nice) + scale_y_continuous(label=function(x)-x)
  )
  print(plt_knn)
  dev.off()

  ppng('knn_noname')
  (plt_knn_nn <- plot_matrix_z(knn_init[[col]]$map,keep_legend=2) +
      ggtitle('KNN Estimate') + scale_y_continuous(label=function(x)-x)
  )
  print(plt_knn_nn)
  dev.off()

  ppng('knn_dx')
  (plt_knn_x <- plot_matrix_dx(knn_init[[col]]$map_1d_x,keep_legend=2) +
      ggtitle(bquote(paste('KNN Estimate, ',f[x])),subtitle=name_nice) + scale_y_continuous(label=function(x)-x)
  )
  print(plt_knn_x)
  dev.off()

  ppng('knn_dy')
  (plt_knn_y <- plot_matrix_dy(knn_init[[col]]$map_1d_y,keep_legend=2) +
      ggtitle(bquote(paste('KNN Estimate, ',f[x])),subtitle=name_nice) + scale_y_continuous(label=function(x)-x)
  )
  print(plt_knn_y)
  dev.off()

  ppng('knn_dxy')
  (plt_knn_xy <- plot_matrix_dxy(knn_init[[col]]$map_2d_xy,keep_legend=2) +
      ggtitle(bquote(paste('KNN Estimate, ',f[x])),subtitle=name_nice) + scale_y_continuous(label=function(x)-x)
  )
  print(plt_knn_xy)
  dev.off()

  # plot estimate
  ppng('estimate')
  (plt_est <- plot_matrix2(estimate_mat,zrange=z_range,keep_legend=2) +
    ggtitle(paste("Spline Estimate"),subtitle=name_nice) + scale_y_continuous(label=function(x)-x)
  )
  print(plt_est)
  dev.off()

  ppng('estimate_dx')
  (plt_est_x <- plot_matrix_dx(estimate_mat_x,keep_legend=2) +
      ggtitle(bquote(paste('Estimate, ',f[x])),subtitle=name_nice) + scale_y_continuous(label=function(x)-x)
  )
  print(plt_est_x)
  dev.off()

  ppng('estimate_dy')
  (plt_est_y <- plot_matrix_dy(estimate_mat_y,keep_legend=2) +
      ggtitle(bquote(paste('Estimate, ',f[y])),subtitle=name_nice) + scale_y_continuous(label=function(x)-x)
  )
  print(plt_est_y)
  dev.off()

  ppng('estimate_dxy')
  (plt_est_xy <- plot_matrix_dxy(estimate_mat_xy,keep_legend=2) +
      ggtitle(bquote(paste('Estimate, ',f[xy])),subtitle=name_nice) + scale_y_continuous(label=function(x)-x)
  )
  print(plt_est_xy)
  dev.off()

  # plot smooth
  ppng('estimate_smooth')
  (plt_sest <- plot_matrix_z(smooth_grid,keep_legend=2) +
      ggtitle(paste('Spline Interpolation'),subtitle=name_nice) +
    scale_x_continuous(label=granularity_trans) +
    scale_y_continuous(label=granularity_trans2)
  )
  print(plt_sest)
  dev.off()

  ppng('estimate_smooth_dx')
  (plt_sest_x <- plot_matrix_dx(smooth_grid_x,keep_legend=2) +
    ggtitle(bquote(paste('Spline Interpolation, ',f[x])),subtitle=name_nice)+
      scale_x_continuous(label=granularity_trans) +
      scale_y_continuous(label=granularity_trans2)
  )
  print(plt_sest_x)
  dev.off()

  ppng('estimate_smooth_dy')
  (plt_sest_y <- plot_matrix_dy(smooth_grid_y,keep_legend=2) +
    ggtitle(bquote(paste('Spline Interpolation, ',f[y])),subtitle=name_nice)+
      scale_x_continuous(label=granularity_trans) +
      scale_y_continuous(label=granularity_trans2)
  )
  print(plt_sest_y)
  dev.off()

  ppng('estimate_smooth_dxy')
  (plt_sest_xy <- plot_matrix_dxy(smooth_grid_xy,keep_legend=2) +
    ggtitle(bquote(paste('Spline Interpolation, ',f[xy])),subtitle=name_nice)+
      scale_x_continuous(label=granularity_trans) +
      scale_y_continuous(label=granularity_trans2)
  )
  print(plt_sest_xy)
  dev.off()


  plot_gg2 <- function(x){
    plot_gg(x,width=width_3d,height=height_3d,scale=scale_3d)
  }

  render_snap <- function(fn){
    render_snapshot(file.path(directory,sprintf('%s_%s',fn,name)))
  }
  render_mov <- function(fn){
    render_movie(
      file.path(directory,sprintf('%s_%s.mp4',fn,name)),
      frames=frames_3d,
      fps=fps_3d,
      zoom=zoom_3d,
      fov=fov_3d
    )
  }

  multi_3d <- function(x,fn,mode='',debug=F){
    if (mode == 'z'){
      add = add_topo_color(mode)
    } else {
      add = NULL
    }
    x = x + add + theme_void() + theme(legend.position='none')
    if (debug){
      print(x)
      return(T)
    }

    plot_gg2(x)
    render_snap(fn)
    render_mov(fn)
  }

  add_topo_color <- function(x){
    r = get(sprintf('%s_range',x),envir=parent.env(environment()))
    scale_fill_gradientn(colors=topo_pal(4),limits=r)
  }
  if (!skip3d){
    # 3d actual
    multi_3d(plt_actual,'actual3d',mode='z')
    multi_3d(plt_actual_dx,'actual_dx3d',mode='dx')
    multi_3d(plt_actual_dy,'actual_dy3d',mode='dy')
    multi_3d(plt_actual_dxy,'actual_dxy3d',mode='dxy')

    # 3d estimate
    multi_3d(plt_est,'estimate3d',mode='z')
    multi_3d(plt_est_x,'estimate3d_x', mode='dx')
    multi_3d(plt_est_y,'estimate3d_y',mode='dy')
    multi_3d(plt_est_xy,'estimate3d_xy',mode='dxy')


    # 3d path
    multi_3d(plt_path,'pathsample3d', mode='z')

    # 3d knn estimate
    multi_3d(plt_knn,'knn3d',mode='z')
    multi_3d(plt_knn_x,'knn3d_x', mode='dx')
    multi_3d(plt_knn_y,'knn3d_y',mode='dy')
    multi_3d(plt_knn_xy,'knn3d_xy',mode='dxy')

    # 3d smooth
    multi_3d(plt_sest,'estimate3d',mode='z')
    multi_3d(plt_sest_x,'estimate3d_x', mode='dx')
    multi_3d(plt_sest_y,'estimate3d_y',mode='dy')
    multi_3d(plt_sest_xy,'estimate3d_xy',mode='dxy')
  }

  ## comparisons...
  #eprox = abs(estimate_mat-path_data$map)
  #kprox = abs(knn_init[[col]]$map - path_data$map)
  kediff = abs(knn_init[[col]]$map - estimate_mat)
  #proxdiff = kprox-eprox

  # knn vs. actual
  ppng('diffs_knn_actual')
  (plt_diff_knn_actual<- plot_matrix2(kprox,zname='|Δz|',keep_legend=2) +
    ggtitle('Magnitude of Difference between KNN and Actual Map')
  )
  print(plt_diff_knn_actual)
  dev.off()

  # knn vs. estimate
  ppng('diffs_knn_est')
  (plt_diff_knn_est <- plot_matrix2(kediff,zname='|Δz|',keep_legend=2) +
    ggtitle('Magnitude of Difference between KNN and Spline Fit',subtitle=name_nice)  + scale_y_continuous(label=function(x)-x)
  )
  print(plt_diff_knn_est)
  dev.off()

  # actual vs. estimate
  ppng('diffs_actual_est')
  (plt_diff_actual_est<- plot_matrix2(eprox,zname='|Δz|',keep_legend=2) +
    ggtitle('Magnitude of Difference between Spline Fit and Actual Map',subtitle=name_nice)  + scale_y_continuous(label=function(x)-x)
  )
  print(plt_diff_actual_est)
  dev.off()

  # improvements in estimate
  ppng('diffs_diff_knn_est')
  (plt_diffs_diff_knn_est <- plot_matrix2c(proxdiff,zname='Δ|Δz|',keep_legend=2) +
    ggtitle(paste0('Improvements in Estimate vs. KNN\n',name_nice),subtitle=bquote(paste('|',epsilon[KNN],'|','-','|',epsilon[spline],'|'))) +
       scale_y_continuous(label=function(x)-x)
  )
  print(plt_diffs_diff_knn_est)
  dev.off()

  # steepness
  map_steepness = sqrt(map_dx^2 + map_dy^2)
  knn_steepness = sqrt(knn_init[[col]]$map_1d_x^2 + knn_init[[col]]$map_1d_y^2)
  est_steepness = sqrt(estimate_mat_x^2 + estimate_mat_y^2)
  smooth_steepness = sqrt(smooth_grid_x^2 + smooth_grid_y^2)

  steepness_zrange = c(0, max(map_steepness,knn_steepness,est_steepness,smooth_steepness))
  # map

  ppng('steepness_map')
  (plt_steepness_actual <- plot_matrix2(map_steepness,zname=bquote(sqrt(f[x]^2+f[y]^2)),keep_legend=2,zrange=steepness_zrange) +
    ggtitle('"Steepness" of map')  + scale_y_continuous(label=function(x)-x)
  )
  print(plt_steepness_actual)
  dev.off()

  # knn

  ppng('steepness_knn')
  (plt_steepness_knn <- plot_matrix2(knn_steepness,zname=bquote(sqrt(f[x]^2+f[y]^2)),keep_legend=2,zrange=steepness_zrange) +
      ggtitle('"Steepness" of KNN estimate')  + scale_y_continuous(label=function(x)-x)
  )
  print(plt_steepness_knn)
  dev.off()

  # estimate

  ppng('steepness_est')
  (plt_steepness_est <- plot_matrix2(est_steepness,zname=bquote(sqrt(f[x]^2+f[y]^2)),keep_legend=2,zrange=steepness_zrange) +
      ggtitle('"Steepness" of spline fit',subtitle=name_nice)  + scale_y_continuous(label=function(x)-x)
  )
  print(plt_steepness_est)
  dev.off()

  # smooth

  ppng('steepness_smooth')
  (plt_steepness_smooth <- plot_matrix2(smooth_steepness,zname=bquote(sqrt(f[x]^2+f[y]^2)),keep_legend=2,zrange=steepness_zrange) +
      ggtitle('"Steepness" of interpolated spline fit',subtitle=name_nice) +
      scale_x_continuous(label=granularity_trans) +
      scale_y_continuous(label=granularity_trans)
  )
  print(plt_steepness_smooth)
  dev.off()

  # 2nd derivative steepness
  # map only for now until I add granular stats for the stan code
  d2_steepness = sqrt(
    Reduce('+',lapply(map_derivatives[['2']],function(x)x^2))
  )
  ppng('steepness_map_d2')
  (plt_steepness_actual_d2 <- plot_matrix2(d2_steepness,zname=bquote(sqrt(f[xx]^2+2*f[xy]^2+f[yy]^2)),keep_legend=2,zrange=steepness_zrange) +
      ggtitle('"Steepness" of map')  + scale_y_continuous(label=function(x)-x)
  )
  print(plt_steepness_actual_d2)
  dev.off()

  # performance statistics
  # R^2 vs. KNN
  ss = sum((map-mean(map))^2)
  knn_r2 = 1-sum(knn_init[[col]]$map - map)^2/ss
  est_r2 = 1-sum(estimate_mat - map)^2/ss

  # 1st derivatives

  # mixed 2nd

  # return plots
  # I am too lazy to a
  plots = list()
  for (object in ls()){
    if (identical(class(get(object)), class(plt_est_y))){
      plots[[object]] = get(object)
    }
  }
  if (return_plots){
    return(plots)
  }
}


plot_animated_walk <- function(path_data,name,directory,...){
  p <- ggplot(path_data$full_instructions) +
    geom_tile(aes(x=x,y=-y,fill=raw_value)) +
    transition_manual(step,cumulative=TRUE) +
    scale_fill_gradientn('Height', colors=topo_pal(4)) +
    xlab('') + ylab('') +
    ggtitle('Path Sampling Walk') +
    theme_bw() +
    scale_y_continuous(label=function(x)-x)

  anim_save(
    file.path(directory,sprintf("walk_%s.gif",name)),
    animation=animate(p,...),
    ...
  )
}
