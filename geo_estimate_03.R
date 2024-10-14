library(rstan)
library(gridExtra)
source('geo_generate.R')

#

SLUG = 'v3_mtd6'

PLOT_DIR='geo_plots'
dir.create(PLOT_DIR,showWarnings=FALSE)

STAN_GEO_DIR = 'stan_geo_files'
dir.create(STAN_GEO_DIR,showWarnings=FALSE)


XY_mat = expand.grid(
  x=seq(-1,5,0.1),
  y=seq(-1,5,0.1)
) %>%
  mutate(
    Val=x*y
  ) %>%
  acast(
    x~y
  )


# sample data

map = load_map(5)[,,1]

path_data2 = row_grid_walk(
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
    theme_bw()
}

# INCOMPLETE
create_matrix_diffs <- function(mat){
  ydim = nrow(mat)
  xdim = ncol(mat)

}

knn_init2 = create_knn_init(path_data2,n_chains=4)

if (FALSE){
  plot_matrix(path_data$masked_map)

  plot_matrix(knn_init[[1]]$map)
}

stan_params2 = list(
  known_map=path_data$masked_map,
  mask=path_data$mask,
  xdim = ncol(map),
  ydim = nrow(map),
  W_known = 20,
  W_derivative = 1,
  sigma_dscale=80,
  granularity=20,
  W_derivative2 = 1,
  W_derivative3 = 0.5
)

stan_ctrl2 = list(max_treedepth=6,adapt_delta=0.8)

# do not save pars or warmup, since dimensionality makes the size ridiculous
t1 = Sys.time()
print(t1)
options(mc.cores=4)
fit2 = stan(
  file='geo_estimate_02.stan',
  data = stan_params2,
  model_name='geo_estimate_v2',
  verbose=TRUE,
  #diagnostic_file='stan_geo_diagnostic',
  chains=4,
  cores=4,
  save_warmup=FALSE,
  pars = c('sum_slope_squared','penalty','penalty_err'),
  init=knn_init,
  control=stan_ctrl2
)
t2 = Sys.time()
print(t2-t1)
estimate2 = get_posterior_mean(fit2)
if (TRUE){
  saveRDS(
    list(fit=fit2,params=stan_params2,knn_init=knn_init2,control=stan_ctrl2,time=t2-t1,path_data=path_data2),
    file.path(STAN_GEO_DIR,sprintf('geo_estimate_02_02_%s.RDS',SLUG))
  )
}



col = 2
estimate_mat2 = build_matrix_from_posterior(estimate2,'map',-1)
eprox2 = abs(estimate_mat2-path_data2$map)
kprox2 = abs(knn_init2[[col]]$map - path_data$map)
proxdiff2 = kprox2-eprox2
plots <- list(
  plot_matrix(estimate_mat2) + ggtitle("ESTIMATE"),
  plot_matrix(path_data2$map) + ggtitle("ACTUAL"),
  plot_matrix(estimate_mat2-path_data2$map) + ggtitle("DIFFERENCE"),
  plot_matrix(set_0_to_na(path_data2$masked_map)) + ggtitle("MASKED"),
  plot_matrix(knn_init2[[col]]$map) + ggtitle("KNN IMPUTE"),
  plot_matrix(estimate_mat2 - set_0_to_na(path_data2$masked_map)) + ggtitle('Estimate Minus Masked'),
  plot_matrix(estimate_mat2 - knn_init2[[col]]$map) + ggtitle('ESTIMATE vs. KNN INIT DIFF'),
  plot_matrix(proxdiff2) + ggtitle('Improvement of Spline vs. KNN')
)

do.call(grid.arrange,plots)

smooth_grid2 = build_matrix_from_posterior(estimate2,'granular_map')
plot_matrix(smooth_grid2)

library(rayshader)
library(RColorBrewer)
topo_pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
(topo_plot2 <- plot_matrix(smooth_grid2,pal=topo_pal(4)))
(topo_3d2 <- plot_gg(topo_plot2,width=4,height=4,scale=310))
render_snapshot(file.path(PLOT_DIR,sprintf('smooth_02_snapshot_topo_%s',SLUG)))
render_movie(file.path(PLOT_DIR,sprintf('smooth_02_animated_topo_%s.mp4',SLUG)),
             frames=720,fps=30,zoom=0.6,fov=30)
#(topo_3d <- plot_gg(plot_matrix(map,pal=topo_pal(4)),width=4,height=4,scale=310))

#plot_gg(plot_matrix(set_0_to_na(path_data2$masked_map)),width=4,height=4,scale=310)

