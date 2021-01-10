# distance and neighbour matrix setup
# library(plyr) 

# returns a distance matrix from set of x,y coordinates 
distance_matrix = function(x_points, y_points){
  
  dist_mat = matrix(NA, ncol = length(x_points), nrow = length(y_points)) 

  for(row in seq_along(x_points)){
  
    dist_mat[row, ] = sqrt((x_points - x_points[row]) ^ 2 + (y_points - y_points[row]) ^ 2)
      
  }

  return(dist_mat)

}

# returns a distance matrix from two sets of x,y coordinates, return a length(x1) by length(x2) matrix 
distance_matrix2 = function(x1, y1, x2, y2){
  
  dist_mat = matrix(NA, nrow = length(x1), ncol = length(x2)) 

  for(row in seq_along(x1)){
  
    dist_mat[row, ] = sqrt((x2 - x1[row]) ^ 2 + (y2 - y1[row]) ^ 2)
      
  }

  return(dist_mat)

}

# I want a slightly different output, so that the matrix has length(data$sur) rows 
# and a coloumn for each other individuals alive that year. fill in to make square with NAs
# return 1 matrix with distances and another with heights and vector giving row
# lengths

dist_height = function(xLoc, yLoc, height, year){

  max_num_inyear = max(tapply(year, INDEX = year, FUN = length)) - 1
  dist_mat = matrix(1, nrow = length(xLoc), ncol = max_num_inyear)
  height_mat = matrix(0, nrow = length(xLoc), ncol = max_num_inyear) 
  num_contemp = numeric(length(xLoc))
  
  for(i in seq_along(xLoc)){
      
    ind_vec = year == year[i] & !is.na(height)
    ind_vec[i] = FALSE
    
    num_contemp[i] = sum(ind_vec)
    
    dist_mat[i, 1:num_contemp[i]] = sqrt((xLoc[ind_vec] - xLoc[i]) ^ 2 + (yLoc[ind_vec] - yLoc[i]) ^ 2)
      
    height_mat[i, 1:num_contemp[i]] = height[ind_vec] 
    
  }
  
  return(list(dist = dist_mat, height_mat = height_mat, num_contemp = num_contemp))

}

# set up the distance from every indivudal to every other indivdual within neoghbourhood of the target
# also record the response variable for each of those individuals and the number of neighbours
dist_res_neigh = function(xLoc, yLoc, year, win_r){

  dist_all = distance_matrix(xLoc, yLoc)
  dist_mat = matrix(1000, nrow = length(xLoc), ncol = length(xLoc) - 1)
  ind_mat = matrix(0, nrow = length(xLoc), ncol = length(xLoc) - 1) 
  num_contemp = numeric(length(xLoc))
  
  for(i in seq_along(xLoc)){
  
    ind_vec = dist_all[i, ] < win_r
    ind_vec[i] = FALSE
    
    num_contemp[i] = sum(ind_vec)
    
    dist_mat[i, 1:num_contemp[i]] = dist_all[i, ind_vec]
      
    ind_mat[i, 1:num_contemp[i]] = which(ind_vec) 
    
  }
  
  return(list(dist = dist_mat, ind_mat = ind_mat, num_contemp = num_contemp))

}

# setup a set of knot points, remove those that are more than iso_dist from nearest data point  
knot_coords = function(x_min, x_max, y_min, y_max, res, dat_x, dat_y, iso_dist){
  
  cords = expand.grid(seq(x_min - res, x_max + res, res), seq(y_min - res, y_max + res, res))
  
  colnames(cords) = c('x_cord', 'y_cord')
  
  # take out a bunch of points that are too isolated to be useful, to stop 
  # extrapolationg to far from the data
  # dist from each knot point to each data point
  is_close = logical(length = length(cords[,1])) 
  
  for(row in seq_along(cords[,1])){
  
    is_close[row] = ifelse(min(sqrt((dat_x - cords[row, 'x_cord']) ^ 2 + (dat_y - cords[row, 'y_cord']) ^ 2)) < iso_dist,
      TRUE, FALSE)
      
  }
  
  cords = cords[is_close, ]
  
  return(cords)

}

# setup a set of knot points, use a different stratergy to try and minimise the 
# the number of knot points used for a given set of 
knot_coords2 = function(dat_x, dat_y, min_dist){
  
  # set every data location to be a knot 
  knot_locs = data.frame(dat_x, dat_y)
  # remove duplicated locations
  knot_locs = knot_locs[!duplicated(knot_locs),]
  
  knot_locs_cl = matrix(NA, ncol = 2, nrow = dim(knot_locs)[1])
  count = 1
  search = TRUE
  while(search){
  
    #make a distance matrix between knots 
    knot_dist = sqrt((knot_locs[,1] - knot_locs[1, 1]) ^ 2 + (knot_locs[,2] - knot_locs[1, 2]) ^ 2)
    # take the mean locations a cluster
    cl_ind = which(knot_dist <= min_dist)
    knot_locs_cl[count,] = c(mean(knot_locs[cl_ind, 1]), mean(knot_locs[cl_ind, 2]))
    
    # remove the used ones from the set
    knot_locs = knot_locs[-cl_ind, ]
    
    count = count + 1
    if(dim(knot_locs)[1] < 1) search = FALSE
  
  }
  
  return(knot_locs_cl[!is.na(knot_locs_cl[,1]), ])

}

# take a distance matrix and return only the entries less than win_r in the form of a sparse matrix
distmat_2_sparse = function(full_distmat, win_r){

  mat_rows = dim(full_distmat)[1]
  mat_cols = dim(full_distmat)[2]
  
  reduced_mat = matrix(NA, nrow = mat_rows, ncol = mat_cols)
  
  # find the entires that are less than win_r
  for(i in 1:mat_rows){
  
    in_win = full_distmat[i, ] <= win_r
    reduced_mat[i, in_win] = full_distmat[i, in_win]
  
  }
  
  # data structures to hold the sparse matrix
  tot_zeros_lower_tri = sum(is.na(reduced_mat[lower.tri(reduced_mat)]))
  tot_ent_lower_tri = sum(!is.na(reduced_mat[lower.tri(reduced_mat)]))
  zero_ent_ind = matrix(NA, nrow = 2, ncol = tot_zeros_lower_tri)
  ent_ind = matrix(NA, nrow = 2, ncol = tot_ent_lower_tri)
  dists = numeric(length = dim(ent_ind)[2])
  
  # only store the lower trianle in the sparse mat
  count_z = 1
  count_e = 1
  for(i in 2:mat_rows){
  
    #zero entries
    inds = which(is.na(reduced_mat[i, 1:(i - 1)]))
    num_zeros = length(inds)
    if(num_zeros > 0){
    
      zero_ent_ind[1, count_z:(count_z + num_zeros - 1)] = i
      zero_ent_ind[2, count_z:(count_z + num_zeros - 1)] = inds
      count_z = count_z + num_zeros
    
    }
    
    # non-zero entries
    inds = which(!is.na(reduced_mat[i, 1:(i - 1)]))
    num_ents = length(inds)
    if(num_ents > 0){
    
      ent_ind[1, count_e:(count_e + num_ents - 1)] = i
      ent_ind[2, count_e:(count_e + num_ents - 1)] = inds
      dists[count_e:(count_e + num_ents - 1)] = reduced_mat[i, inds] 
      count_e = count_e + num_ents
    
    }
  }
  
  return(list(zero_ind = zero_ent_ind, ent_ind = ent_ind, dist = dists))
  
}

# slightly different version that does not do the reflection, to be used with
# output from distance_matrix2(), which may not be symetrical 
distmat_2_sparse2 = function(full_distmat, win_r){

  mat_rows = dim(full_distmat)[1]
  mat_cols = dim(full_distmat)[2]
  
  reduced_mat = matrix(NA, nrow = mat_rows, ncol = mat_cols)
  
  # find the entires that are less than win_r
  for(i in 1:mat_rows){
  
    in_win = full_distmat[i, ] <= win_r
    reduced_mat[i, in_win] = full_distmat[i, in_win]
  
  }
  
  # data structures to hold the sparse matrix
  tot_zeros = sum(is.na(reduced_mat))
  tot_ent = sum(!is.na(reduced_mat))
  zero_ent_ind = matrix(NA, nrow = 2, ncol = tot_zeros)
  ent_ind = matrix(NA, nrow = 2, ncol = tot_ent)
  dists = numeric(length = dim(ent_ind)[2])
  
  # only store the lower trianle in the sparse mat
  count_z = 1
  count_e = 1
  for(i in 1:mat_rows){
  
    #zero entries
    inds = which(is.na(reduced_mat[i, ]))
    num_zeros = length(inds)
    if(num_zeros > 0){
    
      zero_ent_ind[1, count_z:(count_z + num_zeros - 1)] = i
      zero_ent_ind[2, count_z:(count_z + num_zeros - 1)] = inds
      count_z = count_z + num_zeros
    
    }
    
    # non-zero entries
    inds = which(!is.na(reduced_mat[i, ]))
    num_ents = length(inds)
    if(num_ents > 0){
    
      ent_ind[1, count_e:(count_e + num_ents - 1)] = i
      ent_ind[2, count_e:(count_e + num_ents - 1)] = inds
      dists[count_e:(count_e + num_ents - 1)] = reduced_mat[i, inds] 
      count_e = count_e + num_ents
    
    }
  }
  
  return(list(zero_ind = zero_ent_ind, ent_ind = ent_ind, dist = dists))
  
}

# find any locations in test_loc within dist of target location targ_lco
find_near = function(targ_X, targ_Y, test_X, test_Y, dist){

  dist_targ_2_test = sqrt((targ_X - test_X) ^ 2 + (targ_Y - test_Y) ^ 2)
  
  return(which(dist_targ_2_test <= dist))

}





