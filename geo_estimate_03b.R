library(rstan)
library(gridExtra)
source('geo_generate.R')


SLUG = 'v3'

PLOT_DIR='geo_plots'
dir.create(PLOT_DIR,showWarnings=FALSE)

STAN_GEO_DIR = 'stan_geo_files'
dir.create(STAN_GEO_DIR,showWarnings=FALSE)

# sample data

map = load_map(5)[,,1]

path_data = row_grid_walk(
  map,
  x_breaks=c(3,5,15,18,20,25,30,40,50,60,61,65),
  y_breaks=c(4,8,16,24,32,40,47,48,55,2),
  y_first=TRUE,
  margins=list(y=c(1,1),x=c(1,1)),
  debug=F,
  return_instructions=F
)

# diagnostic function
plot_matrix <- function(mat, pal=cet_pal(7,'inferno')){
  require(cetcolor)
  mm = reshape2::melt(mat)
  ggplot(mm) +
    geom_raster(aes(x=Var2,y=-Var1,fill=value)) +
    scale_fill_gradientn(colors=pal) +
    theme_bw() +
    scale_y_continuous(label=function(x) -x)
}

create_matrix_diffs <- function(mat){
  ydim = nrow(mat)
  xdim = ncol(mat)

}

## THIS NEEDS TESTING
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

knn_init = create_knn_init(path_data,n_chains=4)

if (FALSE){
  plot_matrix(path_data$masked_map)

  plot_matrix(knn_init[[1]]$map)
}

stan_params = list(
  known_map=path_data$masked_map,
  mask=path_data$mask,
  xdim = ncol(map),
  ydim = nrow(map),
  W_known = 20,
  W_derivative = 1,
  sigma_dscale=30,
  granularity=20
)

stan_ctrl = list(max_treedepth=6,adapt_delta=0.8)

# do not save pars or warmup, since dimensionality makes the size ridiculous
t1 = Sys.time()
print(t1)
options(mc.cores=4)
fit = stan(
  file='geo_estimate_01_alt.stan',
  data = stan_params,
  model_name='geo_estimate_v1',
  verbose=TRUE,
  #diagnostic_file='stan_geo_diagnostic',
  chains=4,
  cores=4,
  save_warmup=FALSE,
  pars = c('sum_slope_squared','penalty','penalty_err'),
  init=knn_init,
  control=stan_ctrl
)
t2 = Sys.time()
print(t2-t1)
estimate = get_posterior_mean(fit)
if (TRUE){
  saveRDS(
    list(fit=fit,params=stan_params,knn_init=knn_init,control=stan_ctrl,time=t2-t1,path_data=path_data),
    file.path(STAN_GEO_DIR,sprintf('geo_estimate_01_%s.RDS',SLUG))
  )
}






col = 1
estimate_mat = build_matrix_from_posterior(estimate,'map',-1)
eprox = abs(estimate_mat-path_data$map)
kprox = abs(knn_init[[col]]$map - path_data$map)
proxdiff = kprox-eprox
plots <- list(
  plot_matrix(estimate_mat) + ggtitle("ESTIMATE"),
  plot_matrix(path_data$map) + ggtitle("ACTUAL"),
  plot_matrix(estimate_mat-path_data$map) + ggtitle("DIFFERENCE"),
  plot_matrix(set_0_to_na(path_data$masked_map)) + ggtitle("MASKED"),
  plot_matrix(knn_init[[col]]$map) + ggtitle("KNN IMPUTE"),
  plot_matrix(estimate_mat - set_0_to_na(path_data$masked_map)) + ggtitle('Estimate Minus Masked'),
  plot_matrix(estimate_mat - knn_init[[col]]$map) + ggtitle('ESTIMATE vs. KNN INIT DIFF'),
  plot_matrix(proxdiff) + ggtitle('Improvement of Spline vs. KNN')
)

do.call(grid.arrange,plots)

smooth_grid = build_matrix_from_posterior(estimate,'granular_map')
plot_matrix(smooth_grid)

library(rayshader)
library(RColorBrewer)
topo_pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
(topo_plot <- plot_matrix(smooth_grid,pal=topo_pal(4)))
(topo_3d <- plot_gg(topo_plot,width=4,height=4,scale=310))

(topo_3d <- plot_gg(plot_matrix(map,pal=topo_pal(4)),width=4,height=4,scale=310))

render_snapshot(file.path(PLOT_DIR,sprintf('smooth_01_snapshot_topo,_%s',SLUG)))
render_movie(file.path(PLOT_DIR,sprintf('smooth_01_animated_topo_%s.mp4',SLUG)),
             frames=720,fps=30,zoom=0.6,fov=30)
