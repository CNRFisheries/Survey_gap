#BIMAC interpolation algorithms - created by Gianpaolo Coro - gianpaolo.coro@cnr.it
#rm(list=ls(all=TRUE))
start_time <- Sys.time()
library(raster)
library ( R2jags )
library ( coda )
library(plyr)
library(dplyr)
library(digest)


#DEFAULT EXAMPLE: Global scale temperature
#ASC files definitions
#punctual_data_file = "../data/solemon/BIMAC_2022_Solea_solea.csv"
#solea<-read.csv("../data/other data/haul_location.csv", sep=";")
punctual_data_file = paste0("../BIMAC_input/",specie_b,"/",year_b,"/",specie_b,"_hauls_year_", year_b,".csv")
#solea<-read.csv("../data/other data/haul_location.csv", sep=";")

min_x_in_raster<-12#min(solea$lon)#-(2*resolution)
max_x_in_raster<-16#max(solea$lon)#+(2*resolution)
min_y_in_raster<-41.6#min(solea$lat)#-(2*resolution)
max_y_in_raster<-46#max(solea$lat)#+(2*resolution)

#min_x_in_raster<-13.1
#max_x_in_raster<-18.8
#min_y_in_raster<-51.3
#max_y_in_raster<-56.6

depth_file = "../data/gebco_30sec_8.asc"
analysis_depth_min<--12000 #discard all points over this limit
analysis_depth_max<-50

#OUTPUT files will be: 
#prior_data_output = "./BIMAC_output/BIMAC_IDW_prior_no_advection_2022.asc"
#posterior_data_output = "./BIMAC_output/BIMAC_interpolation_no_advection_2022.asc"
#posterior_data_output_sd = "./BIMAC_output/BIMAC_interpolation_sd_no_advection_2022.asc"
prior_data_output = paste0("../BIMAC_output/",specie_b,"/",year_b,"/","/BIMAC_IDW_prior_no_advection_",year_b,"_",specie_b,".asc")
posterior_data_output = paste0("../BIMAC_output/",specie_b,"/",year_b,"/","/BIMAC_interpolation_no_advection_",year_b,"_",specie_b,".asc")
posterior_data_output_sd = paste0("../BIMAC_output/",specie_b,"/",year_b,"/","/BIMAC_interpolation_sd_no_advection_",year_b,"_",specie_b,".asc")


automatic_extent<-F
fast_solving<-T
resolution<-0.1

moving_average_points<-0
smooth = F
if (moving_average_points==0)
  smooth = F

#Step 1 - Create matrices out of the input files
cat("Reading data file\n")
punctual_data<-read.csv(file = punctual_data_file)
names(punctual_data)<-c("x","y","value")


if (automatic_extent){
  min_x_in_raster<-min(punctual_data$x,na.rm = T)
  max_x_in_raster<-max(punctual_data$x,na.rm = T)
  min_y_in_raster<-min(punctual_data$y,na.rm = T)
  max_y_in_raster<-max(punctual_data$y,na.rm = T)
}

adapt_depth<-function(depth_file, resolution, min_x_in_raster, min_y_in_raster, max_x_in_raster, max_y_in_raster){
  caching_string<-paste0(resolution,";",min_x_in_raster,";",min_y_in_raster,";",max_x_in_raster,";",max_y_in_raster)
  sha<-sha1(caching_string)
  cachedfile<-paste0(depth_file,"_",sha)
if (!file.exists(cachedfile)){
    
  asc_file<-raster(depth_file)
  xseq<-seq(from=min_x_in_raster+(resolution/2),to=max_x_in_raster-(resolution/2),by=resolution)
  yseq<-seq(from=min_y_in_raster+(resolution/2),to=max_y_in_raster-(resolution/2),by=resolution)
  grid_of_points<-expand.grid(x = xseq, y = yseq) #combine the x and y coordinate to generate pairs
  grid_values<-extract(x=asc_file,y=grid_of_points,method='simple') #extract raster values for the observations and the grid
  grid_val<-grid_of_points
  grid_val$value<-grid_values
  ypoints<-unique(grid_of_points$y)
  xpoints<-unique(grid_of_points$x)
  ncol_r<-length(xpoints)
  nrow_r<-length(ypoints)
  values<-matrix(nrow = nrow_r,ncol = ncol_r,data = -9999)
  row_counter<-1
  for (y_c in 1:(nrow_r)){
    yp<-ypoints[y_c]
    row_rast<-grid_val[which(grid_val$y == yp),]
    row_rast<-row_rast[order(row_rast$x),]
    values[(nrow_r-row_counter+1),]<-row_rast$value[1:(ncol_r)]
    row_counter<-row_counter+1
  }
  save(values,file=cachedfile)
}else{
  cat("...Loading cached file",cachedfile,"\n")
  load(cachedfile)
}
  return(values)
}

depth_matrix<-adapt_depth(depth_file, resolution, min_x_in_raster, min_y_in_raster, max_x_in_raster, max_y_in_raster)

#Step 2 - Retrieve points on land as positive bathymetry points
cat("Retrieving land\n")
land<-which( (as.vector(t(depth_matrix))>analysis_depth_max) | (as.vector(t(depth_matrix))<analysis_depth_min))

if (length(land)==0) {
  land<-c()
  cat("NOTE: NO LAND LOCATIONS ARE PRESENT\n")
}else if (length(land)==length(as.vector(t(depth_matrix)))) {
  stop("ERROR: No point at the given depth is present in the area, please change the analysis depth input") 
}

#Step 3 - building data grid
cat("Building observation grid\n")
#create the grid
xseq<-seq(from=min_x_in_raster+(resolution/2),to=max_x_in_raster-(resolution/2),by=resolution)
yseq<-seq(from=min_y_in_raster+(resolution/2),to=max_y_in_raster-(resolution/2),by=resolution)
grid_of_points<-expand.grid(x = xseq, y = yseq) #combine the x and y coordinate to generate pairs
grid_val<-grid_of_points
grid_val$value<-0
#read the observation points
ypoints<-unique(grid_of_points$y)
xpoints<-unique(grid_of_points$x)
ncol_r<-length(xpoints)
nrow_r<-length(ypoints)

#delete points falling out of the velocity bounding box
cat("Filtering observation points\n")
punctual_data_within_area<-punctual_data[which(punctual_data$x>min_x_in_raster & punctual_data$x<max_x_in_raster),]
punctual_data_within_area<-punctual_data_within_area[which(punctual_data_within_area$y>min_y_in_raster & punctual_data_within_area$y<max_y_in_raster),]
if (dim(punctual_data_within_area)[1]==0){
  stop("ERROR: NO POINT PRESENT IN THE ANALYSIS AREA\n")
}
#calculate the matrix indices of the observations
punctual_data_within_area$x_i<-round((punctual_data_within_area$x-min_x_in_raster)/resolution)+1
punctual_data_within_area$y_i<-round((max_y_in_raster-punctual_data_within_area$y)/resolution)+1
n_punctual_data<-dim(punctual_data_within_area)[1]
data_matrix<-matrix(nrow = nrow_r,ncol = ncol_r,data = NA)
#assign the observations to matrix elements
distinct_observation_pairs<-distinct(data.frame(x_i=punctual_data_within_area$x_i,y_i=punctual_data_within_area$y_i))
distinct_observation_pairs_xyv<-list()

cat("Transforming observations into a matrix through IDW\n")
inserted_observations<-1
for (pair in 1:nrow(distinct_observation_pairs)){
      xi = distinct_observation_pairs[pair,]$x_i
      yi = distinct_observation_pairs[pair,]$y_i
      pk<-punctual_data_within_area[which(punctual_data_within_area$x_i==xi & punctual_data_within_area$y_i == yi),]
      dist<-sqrt( ( (xi-pk$x)*(xi-pk$x) ) + ( (yi-pk$y)*(yi-pk$y) ) ) 
      if (length(which(dist==0))>0)
        dist[which(dist==0)] = resolution/0.1
      
      inv_dist<-1/(dist)
      avg<-sum((inv_dist*pk$value)/sum(inv_dist))
      if (yi<=nrow(data_matrix) && xi<=ncol(data_matrix)){
        data_matrix[yi,xi]<-avg
        #index to coordinates
        x_m<-((xi-1)*resolution)+min_x_in_raster
        y_m<-max_y_in_raster-((yi-1)*resolution)
        distinct_observation_pairs_xyv[[inserted_observations]]<-data.frame(x_m=x_m,y_m=y_m,value=avg)
        inserted_observations<-inserted_observations+1
      }
}
min_real_value<-min(data_matrix,na.rm = T)
max_real_value<-max(data_matrix,na.rm = T)
distinct_observation_pairs_xyv <- ldply(distinct_observation_pairs_xyv, data.frame)

#Step 4 - calculate the average proximity between the points to assess a radius for prior value averaging
cat("Calculating average proximity between the points\n")
min_distances<-sapply(1:nrow(distinct_observation_pairs_xyv), function(i){
  x<-distinct_observation_pairs_xyv$x_m[i]
  y<-distinct_observation_pairs_xyv$y_m[i]
  ddx<-(x-distinct_observation_pairs_xyv$x_m[-i])
  ddy<-(y-distinct_observation_pairs_xyv$y_m[-i])
  dist<-sqrt( (ddx*ddx) + (ddy*ddy))
  return(min(dist))
},simplify = T)

#use the geometric mean
log_mean_min_dist<-mean(log(min_distances[which(min_distances>0)]))
log_sd_min_dist<-sd(log(min_distances[which(min_distances>0)]))
upper_limit_min_distances<-exp(log_mean_min_dist+1.96*log_sd_min_dist)
upper_limit_min_distances_index<-round(upper_limit_min_distances/resolution)
cat("Proximity range is",upper_limit_min_distances,"deg","(=",upper_limit_min_distances_index,"indices)","\n")

#Step 5 - calculate prior values of the quantity through inverse weighting
cat("Inverse weighting\n")
#data_matrix_prefilled is the matrix of prior values
data_matrix_prefilled<-data_matrix

on_land<-function(i,j,resolution,min_x_in_raster,data_matrix,land){
  x<-((j-1)*resolution)+min_x_in_raster
  y<-max_y_in_raster-((i-1)*resolution)
  index_in_vector<-((i-1)*ncol(data_matrix))+j
  if ( (length(land)>0) && (index_in_vector%in%land) )
    return(T)
  else
    return(F)
}

cat("   Filling gaps around the observations\n")
#for each NA element (=without observation) retrieve observation points within the proximity range; calculate the inverse-weighted values and assign the weighted average
data_matrix_prefilled<-sapply(1:nrow_r, function(i){
  data_matrix_prefilled_row<-vector()
  for (j in 1:ncol_r){
    element<-data_matrix[i,j]
    #if the element is not an observation then process it
    if (is.na(element)){
      #retrieve the x,y coordinates of the cell
      x<-((j-1)*resolution)+min_x_in_raster
      y<-max_y_in_raster-((i-1)*resolution)
      #assign -9999 to land points
      if (on_land(i,j,resolution,min_x_in_raster,data_matrix,land)){
        avg<- -9999
      }
      else{
        #calculate the inverse distance with respect to all points
        ddx<-(x-distinct_observation_pairs_xyv$x_m)
        ddy<-(y-distinct_observation_pairs_xyv$y_m)
        dist<-sqrt( (ddx*ddx) + (ddy*ddy))
        inv_dist<-1/(dist)
        #select the distances within the proximity range
        good_points<-which(dist<=upper_limit_min_distances)
        if (length(good_points)>0){
          #if there are points within the proximity range, calculate the inverse-distance weighted average and assign it to the current cell
          avg<-sum((inv_dist[good_points]*distinct_observation_pairs_xyv$value[good_points]))/sum(inv_dist[good_points])
        }else{ 
          #if no proximity point was found force a hole in the matrix
          avg<- NA
        }
      }
      data_matrix_prefilled_row[j]<-avg
    }else{
      data_matrix_prefilled_row[j]<-element
    }
  }
  return(data_matrix_prefilled_row)
},simplify = T)
#adjust the matrix after the sapply, which returns the rows as columns
data_matrix_prefilled<-t(data_matrix_prefilled)

#Step 6 - Fill the holes' matrix holes through an iterative gap-filling process
cat("   Averaging gaps\n")
cat("Iterating.. ")
still_elements_to_fill<-T
iterations<-1
#iterate until no holes are present
while (still_elements_to_fill){
  still_elements_to_fill<-F
  elements_to_fill<-0
  elements_to_fill_prev<-0
  for (i in 1:nrow_r){
    for (j in 1:ncol_r){
      element<-data_matrix_prefilled[i,j]
      #if the current element is a hole then process it
      if (is.na(element) && !(on_land(i,j,resolution,min_x_in_raster,data_matrix,land))){
        avg<-NA
        #retrieve all cells within a proximity distance (this time taken as a matrix index)
        di<-upper_limit_min_distances_index#1 #distance as an index
        i_0 = max(1,(i-di))
        i_1 = min(nrow_r,(i+di))
        j_0 = max(1,(j-di))
        j_1 = min(ncol_r,(j+di))
        i_c = ((i_1-i_0)/2)+1
        j_c = ((j_1-j_0)/2)+1
        
        submatrix<-as.vector(data_matrix_prefilled[i_0:i_1,j_0:j_1])
        submatrix_dist<-as.vector(data_matrix_prefilled[i_0:i_1,j_0:j_1])
        #calculate the inverse-distance values depending on how far the points are from the current point
        submatrix_vector_index<-1
        for (gy in j_0:j_1){
          dy = abs(gy-j)*resolution
          for (gx in i_0:i_1){
            dx = abs(gx-i)*resolution
            id = 1/((dx*dx)+(dy*dy))
            if (is.na(data_matrix_prefilled[gx,gy]) || (data_matrix_prefilled[gx,gy]==-9999) ){
                id = 1
            }
            submatrix_dist[submatrix_vector_index]<-id
            submatrix_vector_index=submatrix_vector_index+1
          }
        }
        
        #retrieve only non-hole and non-land points
        valid_indices<-which(!is.na(submatrix) & submatrix!=-9999)
        if (length(valid_indices)>0){
          #if there are points around fill the hole with the inverse-distance weighted average
          sm_valid<-submatrix[valid_indices]
          sm_valid_d<-submatrix_dist[valid_indices]
          avg<-sum(sm_valid_d*sm_valid)/sum(sm_valid_d)
        }else{
          #if there hole remains then alert the loop that holes are still present
          still_elements_to_fill<-T
          #increase the number of holes
          elements_to_fill<-elements_to_fill+1
        }
        data_matrix_prefilled[i,j]<-avg
      }
    }
  }
  cat(iterations,"(",elements_to_fill,"); ")
  if ((elements_to_fill_prev>0) && (elements_to_fill<elements_to_fill_prev)){
    elements_to_fill_prev<-elements_to_fill
  }else{
    upper_limit_min_distances_index<-upper_limit_min_distances_index+1
  }
  #increase the number of iterations
  iterations<-iterations+1
}
cat("\n")

if (smooth){
#smooth
  data_matrix_prefilled_smoothed<-sapply(1:nrow_r, function(i){
    #take moving_average_points-elements around each point in a row
    min_i<-max(1,(i-moving_average_points))
    max_i<-min(nrow_r,(i+moving_average_points))
    smoothed<-integer(ncol_r) #inizialise with zeros
    smoothed<-sapply(1:ncol_r, function(j){
      #take moving_average_points-elements around each point along the column
      min_j<-max(1,(j-moving_average_points))
      max_j<-min(ncol_r,(j+moving_average_points))
      #extract the submatrix
      sub_matrix<-data_matrix_prefilled[min_i:max_i,min_j:max_j]
      #transform the submatrix into a vector and delete NA values (land points)
      sub_vector<-as.vector(sub_matrix)
      valid_indices_subvector<-which(!is.na(sub_vector) & sub_vector!=-9999)
      #if points exist report the submatrix mean for the current matrix element; otherwise report -9999
      if (length(valid_indices_subvector)>0){
        avg<-mean(sub_matrix[valid_indices_subvector])
      }
      else
        avg<--9999
      return(avg)
    },simplify = T)
    
    return(smoothed)
  },simplify = T)
data_matrix_prefilled<-t(data_matrix_prefilled_smoothed)
}
cat("Saving the inverse weighted map\n")
ro <- raster(ncol=ncol_r, nrow=nrow_r)
length(values(ro))
data_matrix_prefilled_vector<-as.vector(t(data_matrix_prefilled))
data_matrix_prefilled_vector_to_save<-data_matrix_prefilled_vector

if (length(land)>0){
  data_matrix_prefilled_vector_to_save[land]<--9999
}
values(ro)<-data_matrix_prefilled_vector_to_save
extent(ro)<-extent( min_x_in_raster+(resolution/2),
                    max_x_in_raster-(resolution/2),
                    min_y_in_raster+(resolution/2),
                    max_y_in_raster-(resolution/2))
NAvalue(ro)<- -9999
writeRaster(ro, filename=prior_data_output, format="ascii",overwrite=TRUE)


#Step 8 - solve the advection equation and recalculate the grid-matrix values
cat("Solving advection equation...\n")
#data_matrix_postfilled is the posterior values' matrix
data_matrix_postfilled<-matrix(nrow = nrow_r,ncol = ncol_r,data = 1)
data_matrix_postfilled<-as.vector(data_matrix_postfilled)
N<-length(data_matrix_postfilled)
#calculate the standard deviation of the observation values
prior_sd<-sd(data_matrix,na.rm=T)
data_vector<-as.vector(t(data_matrix))
if (length(which(!is.na(data_vector)))==0){
  stop("ERROR IN THE DATA: NO REAL OBSERVATION PRESENT IN THE AREA!")
}
if (length(which(is.na(data_vector)))==0){
  stop("ERROR IN THE DATA: THERE ARE NO POINTS TO RECONSTRUCT!")
}
#record the indexes of the observation values in the data matrix
is_true_value<-!is.na(data_vector)
real_values_idx<-which(!is.na(data_vector))
#record the indexes of the non-observation values to recalculate
values_to_estimate_idx<-which(is.na(data_vector))

#exclude land points from the data
land_no_values<-c()
if (length(land)>0){
  #exclude land points from the non-observation values indexes
  values_to_estimate_idx<-values_to_estimate_idx[-which(values_to_estimate_idx%in%land)]
  #exclude possible observation values falling on land from the land indexes - these should be used as likelihood data
  true_in_land<-which(land%in%real_values_idx)
  if (length(true_in_land)>0){
    land_no_values<-land[-true_in_land]
  }else
    land_no_values<-land
  #set all land values to zero
  data_matrix_prefilled_vector[land_no_values]<-0  
}

#Set the input parameters of the model
#prior_sd = standard deviation of the observation values
#data_matrix_prefilled_vector = priors' vector, zero for land values
#min_real_value = minimum observation value
#max_real_value = maximum observation value
#real_values_idx = vector indices of the real observations
#values_to_estimate_idx = vector indices of the values to estimate
#land_no_values = indices of land locations excluding those containing real observations
jags.data <- list ("data_matrix_prefilled_vector",
                   "min_real_value","max_real_value",
                   "prior_sd",
                   "real_values_idx",
                   "values_to_estimate_idx",
                   "land_no_values")
#P = posterior distribution values
jags.params <- c("P")

#Note: we use differential equations like https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/advection.pdf
Model = "
model {
#priors

#force posterior distribution values to 0
for (g in land_no_values){
  P[g]<-0
}

#initialise the individual-location posterior distribution values in the real observations
for (g in real_values_idx){
  P[g] ~ dunif (
      (min_real_value),
      (max_real_value))
}

#initialise the individual-location posterior distribution values in the IDW-estimated observations
for (g in values_to_estimate_idx){
  P[g] ~ dunif (
      (min_real_value),
      (max_real_value))
}

#fit the individual distribution values to real observations in the real-value locations - this accounts for error in the observations
for (g in real_values_idx){
  invsigma[g] <- pow(0.1,-2)
  data_matrix_prefilled_vector[g] ~ dnorm(P[g],invsigma[g])
}

#fit the other individual-location values to inverse-weighted estimated observations using a larger standard deviation
for (g in values_to_estimate_idx){
  invsigma_na[g] <- pow(prior_sd,-2)
  data_matrix_prefilled_vector[g] ~ dnorm (P[g],invsigma_na[g])
}

}"

end_time_idw <- Sys.time()
cat("Elapsed Data Preparation + IDW\n")
print(end_time_idw-start_time)

#write the BUGS model to a file
JAGSFILE =" r2ssb.bug "
cat (Model , file = JAGSFILE )

if (fast_solving){
  Nchains = 1
  Nburnin = 10
  Niter = 100
  Nthin = 5 
}else{
  Nchains = 1 # number of Markov chains - to account for non - ergodic convergence
  Nburnin = 100 # burn -in iterations - n. of initial iterations to discard
  Niter = 1000 # total n. of iterations
  Nthin = 10 # thinning - take every 10 samples to lower the dependency among the samples
}

#Run the Gibbs sampling
start_time_jags <- Sys.time()
jagsfit <- jags ( data = jags.data , working.directory=NULL , inits =NULL , jags.params,
                  model.file = JAGSFILE , n.chains=Nchains, n.thin=Nthin , n.iter=Niter , n.burnin=Nburnin )
end_time_jags <- Sys.time()
cat("Elapsed Jags\n")
print(end_time_jags-start_time_jags)

#retrieve the posterior vector
dmp<-jagsfit$BUGSoutput$sims.list$P
#the posterior vector has one column for each P element, samples from the gibbs sampling are contained in each column
#the mean extract the optimal value according to the Monte Carlo integration theorem
dmpC<-colMeans(dmp)
dmpV<-apply(dmp, 2, sd)
#set land points to NA
if(length(land)>0){
  dmpC[land]<-NA
  dmpV[land]<-NA
}

#Step 9 - distribution smoothing
cat("Smoothing\n")
data_matrix_postfilled<-matrix(dmpC, nrow = nrow_r,ncol = ncol_r,byrow = T)
data_matrix_postfilled_sd<-matrix(dmpV, nrow = nrow_r,ncol = ncol_r,byrow = T)

if (smooth){
  data_matrix_postfilled_smoothed<-sapply(1:nrow_r, function(i){
    #take moving_average_points-elements around each point in a row
    min_i<-max(0,(i-moving_average_points))
    max_i<-min(nrow_r,(i+moving_average_points))
    smoothed<-integer(ncol_r) #inizialise with zeros
    smoothed<-sapply(1:ncol_r, function(j){
      #take moving_average_points-elements around each point along the column
      min_j<-max(1,(j-moving_average_points))
      max_j<-min(ncol_r,(j+moving_average_points))
      #extract the submatrix
      sub_matrix<-data_matrix_postfilled[min_i:max_i,min_j:max_j]
      #transform the submatrix into a vector and delete NA values (land points)
      sub_vector<-as.vector(sub_matrix)
      valid_indices_subvector<-which(!is.na(sub_vector))
      #if points exist report the submatrix mean for the current matrix element; otherwise report -9999
      if (length(valid_indices_subvector)>0){
        avg<-mean(sub_matrix[valid_indices_subvector])
      }
      else
        avg<--9999
      return(avg)
    },simplify = T)
    return(smoothed)
  },simplify = T)
  data_matrix_postfilled_smoothed_sd<-sapply(1:nrow_r, function(i){
    #take moving_average_points-elements around each point in a row
    min_i<-max(0,(i-moving_average_points))
    max_i<-min(nrow_r,(i+moving_average_points))
    smoothed<-integer(ncol_r) #inizialise with zeros
    smoothed<-sapply(1:ncol_r, function(j){
      #take moving_average_points-elements around each point along the column
      min_j<-max(1,(j-moving_average_points))
      max_j<-min(ncol_r,(j+moving_average_points))
      #extract the submatrix
      sub_matrix<-data_matrix_postfilled_sd[min_i:max_i,min_j:max_j]
      #transform the submatrix into a vector and delete NA values (land points)
      sub_vector<-as.vector(sub_matrix)
      valid_indices_subvector<-which(!is.na(sub_vector))
      #if points exist report the submatrix mean for the current matrix element; otherwise report -9999
      if (length(valid_indices_subvector)>0){
        avg<-mean(sub_matrix[valid_indices_subvector])
      }
      else
        avg<--9999
      return(avg)
    },simplify = T)
    return(smoothed)
  },simplify = T)
  
  #transpose since sapply returns rows as columns
  data_matrix_postfilled_smoothed<-t(data_matrix_postfilled_smoothed)
  data_matrix_postfilled_smoothed_sd<-t(data_matrix_postfilled_smoothed_sd)
  #substitute the smoothed matrix to the posterior matrix
  data_matrix_postfilled<-data_matrix_postfilled_smoothed
  data_matrix_postfilled_sd<-data_matrix_postfilled_smoothed_sd
  cat("Saving smoothed map...\n")
}

#Step 10 - save the posterior probability matrix
ro_p <- raster(ncol=ncol_r, nrow=nrow_r)
ro_p_sd <- raster(ncol=ncol_r, nrow=nrow_r)
values_vec_post<-as.vector(t(data_matrix_postfilled))
values_vec_post_sd<-as.vector(t(data_matrix_postfilled_sd))
if (length(land)>0){
  values_vec_post[land]<--9999
  values_vec_post_sd[land]<--9999
}

values(ro_p)<-values_vec_post
extent(ro_p)<-extent( min_x_in_raster+(resolution/2),
                      max_x_in_raster-(resolution/2),
                      min_y_in_raster+(resolution/2),
                      max_y_in_raster-(resolution/2))
NAvalue(ro_p)<- -9999

writeRaster(ro_p, filename=posterior_data_output, format="ascii",overwrite=TRUE)


values(ro_p_sd)<-values_vec_post_sd
extent(ro_p_sd)<-extent( min_x_in_raster+(resolution/2),
                         max_x_in_raster-(resolution/2),
                         min_y_in_raster+(resolution/2),
                         max_y_in_raster-(resolution/2))
NAvalue(ro_p_sd)<- -9999
writeRaster(ro_p_sd, filename=posterior_data_output_sd, format="ascii",overwrite=TRUE)

prior<-as.vector(t(data_matrix_prefilled))
posterior<-as.vector(t(data_matrix_postfilled))

#evaluate the difference between posterior and prior probabilities
if (length((land)>0)){
  prior[land]<-NA
  posterior[land]<-NA
}
diff<-posterior-prior
if (length((land)>0)){
  diff<-abs(diff[-which(is.na(posterior-prior))])
  diff_rel<-diff/abs(posterior[-land])
}else{
  diff<-abs(diff)
  diff_rel<-diff/abs(posterior)
}
cat("Maximum absolute discrepancy between prior and posterior distributions",max(diff),"\n")
cat("Mean relative discrepancy",100*mean(diff_rel),"%\n")
cat("Sd of the discrepancy",sd(diff),"\n")
end_time <- Sys.time()
cat("Done.\n")
cat("Elapsed.\n")
print(end_time-start_time)
