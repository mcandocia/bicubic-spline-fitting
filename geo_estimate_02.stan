// this version takes in the known map values and a (1,0) mask
// this also controls 2nd and 3rd derivative magnitudes to reduce sharp bumps
data {
  int<lower=1> ydim;
  int<lower=1> xdim;
  matrix[ydim,xdim] mask;
  matrix[ydim,xdim] known_map;
  // adjust strength of penalties for mis-estimates
  real W_known;
  // different criteria for known vs. unkown points
  // note: a squared 1D incline for a given slope would give
  //       a penalty of W * slope^2, so the below W
  // should be less than the one above
  // this will avoid the optimizer being unable to decide whether to minimize
  // slope or to get a proper fit
  real W_derivative;
  real W_derivative2;
  real W_derivative3;
  // this should be roughly the maximum slope you see
  real<lower=0> sigma_dscale;
  // controls generated output
  int<lower=1> granularity;
}
// not particularly necessary, but just keeping it here in case needed later
transformed data {
  real<lower=0> W_dscale;
  real<lower=0> step_size;
  int n_y_tiles;
  int n_x_tiles;
  n_y_tiles = 1 + granularity * (ydim-1);
  n_x_tiles = 1 + granularity * (xdim-1);
  step_size = 1./granularity;
  W_dscale = 0.1;
}
// the map itself is being modeled
parameters {
  matrix[ydim,xdim] map;
  matrix[ydim,xdim] map_1d_y;
  matrix[ydim,xdim] map_1d_x;
  matrix[ydim,xdim] map_2d_xy;
}
// this is where bulk of math happens
transformed parameters {
  matrix[ydim-1,xdim-1] f_00;
  matrix[ydim-1,xdim-1] fdy_00;
  matrix[ydim-1,xdim-1] fdx_00;
  matrix[ydim-1,xdim-1] fdxy_00;
  matrix[ydim-1,xdim-1] f_01;
  matrix[ydim-1,xdim-1] fdy_01;
  matrix[ydim-1,xdim-1] fdx_01;
  matrix[ydim-1,xdim-1] fdxy_01;
  matrix[ydim-1,xdim-1] f_10;
  matrix[ydim-1,xdim-1] fdy_10;
  matrix[ydim-1,xdim-1] fdx_10;
  matrix[ydim-1,xdim-1] fdxy_10;
  matrix[ydim-1,xdim-1] f_11;
  matrix[ydim-1,xdim-1] fdy_11;
  matrix[ydim-1,xdim-1] fdx_11;
  matrix[ydim-1,xdim-1] fdxy_11;
  matrix[ydim-1,xdim-1] a00;
  matrix[ydim-1,xdim-1] a01;
  matrix[ydim-1,xdim-1] a02;
  matrix[ydim-1,xdim-1] a03;
  matrix[ydim-1,xdim-1] a10;
  matrix[ydim-1,xdim-1] a11;
  matrix[ydim-1,xdim-1] a12;
  matrix[ydim-1,xdim-1] a13;
  matrix[ydim-1,xdim-1] a20;
  matrix[ydim-1,xdim-1] a21;
  matrix[ydim-1,xdim-1] a22;
  matrix[ydim-1,xdim-1] a23;
  matrix[ydim-1,xdim-1] a30;
  matrix[ydim-1,xdim-1] a31;
  matrix[ydim-1,xdim-1] a32;
  matrix[ydim-1,xdim-1] a33;
  matrix[ydim,xdim] estimate;
  real sum_slope_squared;
  real penalty;
  real penalty_err;
  // define blocks
  f_00 = block(map,1,1,ydim-1,xdim-1);
  fdy_00 = block(map_1d_y,1,1,ydim-1,xdim-1);
  fdx_00 = block(map_1d_x,1,1,ydim-1,xdim-1);
  fdxy_00 = block(map_2d_xy,1,1,ydim-1,xdim-1);
  f_01 = block(map,1,2,ydim-1,xdim-1);
  fdy_01 = block(map_1d_y,1,2,ydim-1,xdim-1);
  fdx_01 = block(map_1d_x,1,2,ydim-1,xdim-1);
  fdxy_01 = block(map_2d_xy,1,2,ydim-1,xdim-1);
  f_10 = block(map,2,1,ydim-1,xdim-1);
  fdy_10 = block(map_1d_y,2,1,ydim-1,xdim-1);
  fdx_10 = block(map_1d_x,2,1,ydim-1,xdim-1);
  fdxy_10 = block(map_2d_xy,2,1,ydim-1,xdim-1);
  f_11 = block(map,2,2,ydim-1,xdim-1);
  fdy_11 = block(map_1d_y,2,2,ydim-1,xdim-1);
  fdx_11 = block(map_1d_x,2,2,ydim-1,xdim-1);
  fdxy_11 = block(map_2d_xy,2,2,ydim-1,xdim-1);
  // calculate (meta-programmed :D)
  a00=f_00;
  a01=fdx_00;
  a02=-3*f_00 + 3*f_01 - 2*fdx_00 - fdx_01;
  a03=2*f_00 - 2*f_01 + fdx_00 + fdx_01;
  a10=fdy_00;
  a11=fdxy_00;
  a12=-2*fdxy_00 - fdxy_01 - 3*fdy_00 + 3*fdy_01;
  a13=fdxy_00 + fdxy_01 + 2*fdy_00 - 2*fdy_01;
  a20=-3*f_00 + 3*f_10 - 2*fdy_00 - fdy_10;
  a21=-3*fdx_00 + 3*fdx_10 - 2*fdxy_00 - fdxy_10;
  a22=9*f_00 - 9*f_01 - 9*f_10 + 9*f_11 + 6*fdx_00 + 3*fdx_01 - 6*fdx_10 - 3*fdx_11 + 4*fdxy_00 + 2*fdxy_01 + 2*fdxy_10 + fdxy_11 + 6*fdy_00 - 6*fdy_01 + 3*fdy_10 - 3*fdy_11;
  a23=-6*f_00 + 6*f_01 + 6*f_10 - 6*f_11 - 3*fdx_00 - 3*fdx_01 + 3*fdx_10 + 3*fdx_11 - 2*fdxy_00 - 2*fdxy_01 - fdxy_10 - fdxy_11 - 4*fdy_00 + 4*fdy_01 - 2*fdy_10 + 2*fdy_11;
  a30=2*f_00 - 2*f_10 + fdy_00 + fdy_10;
  a31=2*fdx_00 - 2*fdx_10 + fdxy_00 + fdxy_10;
  a32=-6*f_00 + 6*f_01 + 6*f_10 - 6*f_11 - 4*fdx_00 - 2*fdx_01 + 4*fdx_10 + 2*fdx_11 - 2*fdxy_00 - fdxy_01 - 2*fdxy_10 - fdxy_11 - 3*fdy_00 + 3*fdy_01 - 3*fdy_10 + 3*fdy_11;
  a33=4*f_00 - 4*f_01 - 4*f_10 + 4*f_11 + 2*fdx_00 + 2*fdx_01 - 2*fdx_10 - 2*fdx_11 + fdxy_00 + fdxy_01 + fdxy_10 + fdxy_11 + 2*fdy_00 - 2*fdy_01 + 2*fdy_10 - 2*fdy_11;
  // assign (ydim-1)x(xdim-1) block the a00 value
  estimate[1:(ydim-1),1:(xdim-1)] = a00;
  // asign sum(aY0) coefficients to bottom
  estimate[ydim,1:(xdim-1)] = to_row_vector(block(a00 + a10 + a20 + a30,ydim-1,1,1,xdim-1));
  // assign sum(a0X) coefficients to right
  estimate[1:(ydim-1),xdim] = to_vector(block(a00 + a01 + a02 + a03,1,xdim-1,ydim-1,1));
  // assign sum(aYX) coefficients to bottom right corner
  estimate[ydim,xdim] = a00[ydim-1,xdim-1] + a01[ydim-1,xdim-1] + a02[ydim-1,xdim-1] + a03[ydim-1,xdim-1] +
    a10[ydim-1,xdim-1] + a11[ydim-1,xdim-1] + a12[ydim-1,xdim-1] + a13[ydim-1,xdim-1] +
    a20[ydim-1,xdim-1] + a21[ydim-1,xdim-1] + a22[ydim-1,xdim-1] + a23[ydim-1,xdim-1] +
    a30[ydim-1,xdim-1] + a31[ydim-1,xdim-1] + a32[ydim-1,xdim-1] + a33[ydim-1,xdim-1];
  // derivatives are already included in parameterization and are not part of any
  // estimation function, so we can just write in penalty term  here

  sum_slope_squared = sum(a01.*a01 + 2*a01.*a02 + 2*a01.*a03 + a01.*a11 + a01.*a12 + a01.*a13 + 2*a01.*a21/3 + 2*a01.*a22/3 + 2*a01.*a23/3 + a01.*a31/2 + a01.*a32/2 + a01.*a33/2 + 4*a02.*a02/3 + 3*a02.*a03 + a02.*a11 + 4*a02.*a12/3 + 3*a02.*a13/2 + 2*a02.*a21/3 + 8*a02.*a22/9 + a02.*a23 + a02.*a31/2 + 2*a02.*a32/3 + 3*a02.*a33/4 + 9*a03.*a03/5 + a03.*a11 + 3*a03.*a12/2 + 9*a03.*a13/5 + 2*a03.*a21/3 + a03.*a22 + 6*a03.*a23/5 + a03.*a31/2 + 3*a03.*a32/4 + 9*a03.*a33/10 + a10.*a10 + a10.*a11 + 2*a10.*a12/3 + a10.*a13/2 + 2*a10.*a20 + a10.*a21 + 2*a10.*a22/3 + a10.*a23/2 + 2*a10.*a30 + a10.*a31 + 2*a10.*a32/3 + a10.*a33/2 + 2*a11.*a11/3 + 7*a11.*a12/6 + 16*a11.*a13/15 + a11.*a20 + 7*a11.*a21/6 + a11.*a22 + 9*a11.*a23/10 + a11.*a30 + 16*a11.*a31/15 + 9*a11.*a32/10 + 4*a11.*a33/5 + 29*a12.*a12/45 + 4*a12.*a13/3 + 2*a12.*a20/3 + a12.*a21 + 16*a12.*a22/15 + 13*a12.*a23/12 + 2*a12.*a30/3 + 9*a12.*a31/10 + 14*a12.*a32/15 + 14*a12.*a33/15 + 26*a13.*a13/35 + a13.*a20/2 + 9*a13.*a21/10 + 13*a13.*a22/12 + 83*a13.*a23/70 + a13.*a30/2 + 4*a13.*a31/5 + 14*a13.*a32/15 + 176*a13.*a33/175 + 4*a20.*a20/3 + 4*a20.*a21/3 + 8*a20.*a22/9 + 2*a20.*a23/3 + 3*a20.*a30 + 3*a20.*a31/2 + a20.*a32 + 3*a20.*a33/4 + 29*a21.*a21/45 + 16*a21.*a22/15 + 14*a21.*a23/15 + 3*a21.*a30/2 + 4*a21.*a31/3 + 13*a21.*a32/12 + 14*a21.*a33/15 + 8*a22.*a22/15 + 47*a22.*a23/45 + a22.*a30 + 13*a22.*a31/12 + 47*a22.*a32/45 + a22.*a33 + 289*a23.*a23/525 + 3*a23.*a30/4 + 14*a23.*a31/15 + a23.*a32 + 36*a23.*a33/35 + 9*a30.*a30/5 + 9*a30.*a31/5 + 6*a30.*a32/5 + 9*a30.*a33/10 + 26*a31.*a31/35 + 83*a31.*a32/70 + 176*a31.*a33/175 + 289*a32.*a32/525 + 36*a32.*a33/35 + 18*a33.*a33/35);  // now some sensible controls to make sure the derivatives don't blow up when being sampled

  // penalty
  penalty = 0;
  penalty_err = 0;
  penalty_err += W_known * sum((known_map - estimate)^2 .* mask);
  penalty += penalty_err;

  // calculate the integral of the sum of absolute value of all derivatives squared
  // same as sqrt(dx^2 +dy^2)^2 which is nice mathematically
  // for second derivatives 2 * (dxy)^2 is added
  // for third derivatives 3 * (dxy^2)^2 + 3 * (dx^2y) is added
  // there is a combinatorics derivation I don't completely understand
  // but the example rotated equations checked out
  penalty += W_derivative * sum_slope_squared;
  penalty += W_derivative2 * sum(4*a02.*a02 + 12*a02.*a03 + 4*a02.*a12 + 6*a02.*a13 + 8*a02.*a22/3 + 4*a02.*a23 + 2*a02.*a32 + 3*a02.*a33 + 12*a03.*a03 + 6*a03.*a12 + 12*a03.*a13 + 4*a03.*a22 + 8*a03.*a23 + 3*a03.*a32 + 6*a03.*a33 + 2*a11.*a11 + 4*a11.*a12 + 4*a11.*a13 + 4*a11.*a21 + 4*a11.*a22 + 4*a11.*a23 + 4*a11.*a31 + 4*a11.*a32 + 4*a11.*a33 + 4*a12.*a12 + 10*a12.*a13 + 4*a12.*a21 + 22*a12.*a22/3 + 9*a12.*a23 + 4*a12.*a31 + 104*a12.*a32/15 + 42*a12.*a33/5 + 38*a13.*a13/5 + 4*a13.*a21 + 9*a13.*a22 + 66*a13.*a23/5 + 4*a13.*a31 + 42*a13.*a32/5 + 12*a13.*a33 + 4*a20.*a20 + 4*a20.*a21 + 8*a20.*a22/3 + 2*a20.*a23 + 12*a20.*a30 + 6*a20.*a31 + 4*a20.*a32 + 3*a20.*a33 + 4*a21.*a21 + 22*a21.*a22/3 + 104*a21.*a23/15 + 6*a21.*a30 + 10*a21.*a31 + 9*a21.*a32 + 42*a21.*a33/5 + 232*a22.*a22/45 + 176*a22.*a23/15 + 4*a22.*a30 + 9*a22.*a31 + 176*a22.*a32/15 + 13*a22.*a33 + 272*a23.*a23/35 + 3*a23.*a30 + 42*a23.*a31/5 + 13*a23.*a32 + 578*a23.*a33/35 + 12*a30.*a30 + 12*a30.*a31 + 8*a30.*a32 + 6*a30.*a33 + 38*a31.*a31/5 + 66*a31.*a32/5 + 12*a31.*a33 + 272*a32.*a32/35 + 578*a32.*a33/35 + 1734*a33.*a33/175);
  penalty += W_derivative3 * sum(36*a03.*a03 + 36*a03.*a13 + 24*a03.*a23 + 18*a03.*a33 + 12*a12.*a12 + 36*a12.*a13 + 24*a12.*a22 + 36*a12.*a23 + 24*a12.*a32 + 36*a12.*a33 + 48*a13.*a13 + 36*a13.*a22 + 90*a13.*a23 + 36*a13.*a32 + 432*a13.*a33/5 + 12*a21.*a21 + 24*a21.*a22 + 24*a21.*a23 + 36*a21.*a31 + 36*a21.*a32 + 36*a21.*a33 + 32*a22.*a22 + 84*a22.*a23 + 36*a22.*a31 + 84*a22.*a32 + 108*a22.*a33 + 384*a23.*a23/5 + 36*a23.*a31 + 108*a23.*a32 + 924*a23.*a33/5 + 36*a30.*a30 + 36*a30.*a31 + 24*a30.*a32 + 18*a30.*a33 + 48*a31.*a31 + 90*a31.*a32 + 432*a31.*a33/5 + 384*a32.*a32/5 + 924*a32.*a33/5 + 4896*a33.*a33/35);
  // for some reason you can't use a matrix with the *_lpdf functions, so manually calculating instead
  penalty += W_dscale * sum(map_1d_y^2)/sigma_dscale^2;
  penalty += W_dscale * sum(map_1d_x^2)/sigma_dscale^2;
  penalty += W_dscale * sum(map_2d_xy^2)/sigma_dscale^2;
}

model {
  target += -penalty;
}

generated quantities {


  // will generate 4 separate interpolated maps
  matrix[n_y_tiles,n_x_tiles] granular_map;
  matrix[n_y_tiles,n_x_tiles] granular_map_dy;
  matrix[n_y_tiles,n_x_tiles] granular_map_dx;
  matrix[n_y_tiles,n_x_tiles] granular_map_dxy;

  real y;
  real x;

  //
  for (i in 1:ydim){
    for (j in 1:(xdim)){
      if (j != xdim){
        if (i != ydim){
          // no edges
          for (k in 1:granularity){
            for (l in 1:granularity){
              y = step_size * (k-1);
              x = step_size * (l-1);
              granular_map[granularity * (i-1) + k][granularity * (j-1) + l] = a00[i,j] + a01[i,j]*x + a02[i,j]*x^2 + a03[i,j]*x^3 + a10[i,j]*y + a11[i,j]*x*y + a12[i,j]*x^2*y + a13[i,j]*x^3*y + a20[i,j]*y^2 + a21[i,j]*x*y^2 + a22[i,j]*x^2*y^2 + a23[i,j]*x^3*y^2 + a30[i,j]*y^3 + a31[i,j]*x*y^3 + a32[i,j]*x^2*y^3 + a33[i,j]*x^3*y^3;
              granular_map_dy[granularity * (i-1) + k][granularity * (j-1) + l] = a10[i,j] + a11[i,j]*x + a12[i,j]*x^2 + a13[i,j]*x^3 + 2*a20[i,j]*y + 2*a21[i,j]*x*y + 2*a22[i,j]*x^2*y + 2*a23[i,j]*x^3*y + 3*a30[i,j]*y^2 + 3*a31[i,j]*x*y^2 + 3*a32[i,j]*x^2*y^2 + 3*a33[i,j]*x^3*y^2;
              granular_map_dx[granularity * (i-1) + k][granularity * (j-1) + l] = a01[i,j] + 2*a02[i,j]*x + 3*a03[i,j]*x^2 + a11[i,j]*y + 2*a12[i,j]*x*y + 3*a13[i,j]*x^2*y + a21[i,j]*y^2 + 2*a22[i,j]*x*y^2 + 3*a23[i,j]*x^2*y^2 + a31[i,j]*y^3 + 2*a32[i,j]*x*y^3 + 3*a33[i,j]*x^2*y^3;
              granular_map_dxy[granularity * (i-1) + k][granularity * (j-1) + l] = a11[i,j] + 2*a12[i,j]*x + 3*a13[i,j]*x^2 + 2*a21[i,j]*y + 4*a22[i,j]*x*y + 6*a23[i,j]*x^2*y + 3*a31[i,j]*y^2 + 6*a32[i,j]*x*y^2 + 9*a33[i,j]*x^2*y^2;
            }
          }
        } else {
          // y edge
            y = 1;
            for (l in 1:granularity){
              // drop i by 1
                x = step_size * (l-1);
                granular_map[n_y_tiles][granularity * (j-1) + l] = a00[i-1,j] + a01[i-1,j]*x + a02[i-1,j]*x^2 + a03[i-1,j]*x^3 + a10[i-1,j]*y + a11[i-1,j]*x*y + a12[i-1,j]*x^2*y + a13[i-1,j]*x^3*y + a20[i-1,j]*y^2 + a21[i-1,j]*x*y^2 + a22[i-1,j]*x^2*y^2 + a23[i-1,j]*x^3*y^2 + a30[i-1,j]*y^3 + a31[i-1,j]*x*y^3 + a32[i-1,j]*x^2*y^3 + a33[i-1,j]*x^3*y^3;
                granular_map_dy[n_y_tiles][granularity * (j-1) + l] = a10[i-1,j] + a11[i-1,j]*x + a12[i-1,j]*x^2 + a13[i-1,j]*x^3 + 2*a20[i-1,j]*y + 2*a21[i-1,j]*x*y + 2*a22[i-1,j]*x^2*y + 2*a23[i-1,j]*x^3*y + 3*a30[i-1,j]*y^2 + 3*a31[i-1,j]*x*y^2 + 3*a32[i-1,j]*x^2*y^2 + 3*a33[i-1,j]*x^3*y^2;
                granular_map_dx[n_y_tiles][granularity * (j-1) + l] = a01[i-1,j] + 2*a02[i-1,j]*x + 3*a03[i-1,j]*x^2 + a11[i-1,j]*y + 2*a12[i-1,j]*x*y + 3*a13[i-1,j]*x^2*y + a21[i-1,j]*y^2 + 2*a22[i-1,j]*x*y^2 + 3*a23[i-1,j]*x^2*y^2 + a31[i-1,j]*y^3 + 2*a32[i-1,j]*x*y^3 + 3*a33[i-1,j]*x^2*y^3;
                granular_map_dxy[n_y_tiles][granularity * (j-1) + l] = a11[i-1,j] + 2*a12[i-1,j]*x + 3*a13[i-1,j]*x^2 + 2*a21[i-1,j]*y + 4*a22[i-1,j]*x*y + 6*a23[i-1,j]*x^2*y + 3*a31[i-1,j]*y^2 + 6*a32[i-1,j]*x*y^2 + 9*a33[i-1,j]*x^2*y^2;
            }
        }
      } else {
         // evaluate x = x_dim and y != ydim
        if (i != ydim){
            x = 1;
            for (k in 1:granularity){
              // drop j by 1
                y = step_size * (k-1);
                granular_map[granularity * (i-1) + k][n_x_tiles] = a00[i,j-1] + a01[i,j-1]*x + a02[i,j-1]*x^2 + a03[i,j-1]*x^3 + a10[i,j-1]*y + a11[i,j-1]*x*y + a12[i,j-1]*x^2*y + a13[i,j-1]*x^3*y + a20[i,j-1]*y^2 + a21[i,j-1]*x*y^2 + a22[i,j-1]*x^2*y^2 + a23[i,j-1]*x^3*y^2 + a30[i,j-1]*y^3 + a31[i,j-1]*x*y^3 + a32[i,j-1]*x^2*y^3 + a33[i,j-1]*x^3*y^3;
                granular_map_dy[granularity * (i-1) + k][n_x_tiles] = a10[i,j-1] + a11[i,j-1]*x + a12[i,j-1]*x^2 + a13[i,j-1]*x^3 + 2*a20[i,j-1]*y + 2*a21[i,j-1]*x*y + 2*a22[i,j-1]*x^2*y + 2*a23[i,j-1]*x^3*y + 3*a30[i,j-1]*y^2 + 3*a31[i,j-1]*x*y^2 + 3*a32[i,j-1]*x^2*y^2 + 3*a33[i,j-1]*x^3*y^2;
                granular_map_dx[granularity * (i-1) + k][n_x_tiles] = a01[i,j-1] + 2*a02[i,j-1]*x + 3*a03[i,j-1]*x^2 + a11[i,j-1]*y + 2*a12[i,j-1]*x*y + 3*a13[i,j-1]*x^2*y + a21[i,j-1]*y^2 + 2*a22[i,j-1]*x*y^2 + 3*a23[i,j-1]*x^2*y^2 + a31[i,j-1]*y^3 + 2*a32[i,j-1]*x*y^3 + 3*a33[i,j-1]*x^2*y^3;
                granular_map_dxy[granularity * (i-1) + k][n_x_tiles] = a11[i,j-1] + 2*a12[i,j-1]*x + 3*a13[i,j-1]*x^2 + 2*a21[i,j-1]*y + 4*a22[i,j-1]*x*y + 6*a23[i,j-1]*x^2*y + 3*a31[i,j-1]*y^2 + 6*a32[i,j-1]*x*y^2 + 9*a33[i,j-1]*x^2*y^2;
            }
         } else {
            // evaluate x = x_dim and y == ydim
            // drop i and j by 1
            x = 1;
            y = 1;
            granular_map[n_y_tiles][n_x_tiles] = a00[i-1,j-1] + a01[i-1,j-1]*x + a02[i-1,j-1]*x^2 + a03[i-1,j-1]*x^3 + a10[i-1,j-1]*y + a11[i-1,j-1]*x*y + a12[i-1,j-1]*x^2*y + a13[i-1,j-1]*x^3*y + a20[i-1,j-1]*y^2 + a21[i-1,j-1]*x*y^2 + a22[i-1,j-1]*x^2*y^2 + a23[i-1,j-1]*x^3*y^2 + a30[i-1,j-1]*y^3 + a31[i-1,j-1]*x*y^3 + a32[i-1,j-1]*x^2*y^3 + a33[i-1,j-1]*x^3*y^3;
            granular_map_dy[n_y_tiles][n_x_tiles] = a10[i-1,j-1] + a11[i-1,j-1]*x + a12[i-1,j-1]*x^2 + a13[i-1,j-1]*x^3 + 2*a20[i-1,j-1]*y + 2*a21[i-1,j-1]*x*y + 2*a22[i-1,j-1]*x^2*y + 2*a23[i-1,j-1]*x^3*y + 3*a30[i-1,j-1]*y^2 + 3*a31[i-1,j-1]*x*y^2 + 3*a32[i-1,j-1]*x^2*y^2 + 3*a33[i-1,j-1]*x^3*y^2;
            granular_map_dx[n_y_tiles][n_x_tiles] = a01[i-1,j-1] + 2*a02[i-1,j-1]*x + 3*a03[i-1,j-1]*x^2 + a11[i-1,j-1]*y + 2*a12[i-1,j-1]*x*y + 3*a13[i-1,j-1]*x^2*y + a21[i-1,j-1]*y^2 + 2*a22[i-1,j-1]*x*y^2 + 3*a23[i-1,j-1]*x^2*y^2 + a31[i-1,j-1]*y^3 + 2*a32[i-1,j-1]*x*y^3 + 3*a33[i-1,j-1]*x^2*y^3;
            granular_map_dxy[n_y_tiles][n_x_tiles] = a11[i-1,j-1] + 2*a12[i-1,j-1]*x + 3*a13[i-1,j-1]*x^2 + 2*a21[i-1,j-1]*y + 4*a22[i-1,j-1]*x*y + 6*a23[i-1,j-1]*x^2*y + 3*a31[i-1,j-1]*y^2 + 6*a32[i-1,j-1]*x*y^2 + 9*a33[i-1,j-1]*x^2*y^2;
          }
          // end j = xdim
      }

      // end x(j) loop
    }

    // end y(i) loop
  }



}
