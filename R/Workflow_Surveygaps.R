library(raster)
library(sqldf)
feature_selection=T
species<-c("Genus_species","Genus_species") 
years<-c(2019,2020)

for (specie in species){
 cat(paste0("species in progress ",specie,"\n"))
 data_path=paste0("../data/",specie,".csv")
 missing_haul_path=("../data/missing_hauls.csv")
 data<-read.csv(data_path)
 haul_m<-read.csv(missing_haul_path)
 
 for (year in years){
   cat(paste0("year in progress ",year,"\n"))
   haul_y<-haul_m[which(haul_m$year==year),"station"]

  ### Data preparation SSA###
  for (i in 1:length(haul_y)){
    cat(paste0("generating ssa input for haul ", haul_y[i],"\n"))
    data_2<-data[which(data$station==haul_y[i] & data$year<year),]
    data_3<-data_2[,c("year","biomindex")]
    path=paste0("../ssa_input")
    path_1=paste0("../ssa_input/",specie)
    ssa_input_folder=paste0("./ssa_input/",specie,"/",year)
    if(!dir.exists(path)) dir.create(path=path)
    if(!dir.exists(path_1)) dir.create(path=path_1)
    if(!dir.exists(ssa_input_folder)) dir.create(path=ssa_input_folder)
    write.csv(data_3,paste0(ssa_input_folder,"/",specie,"_haul_", haul_y[i],".csv"),row.names=FALSE)
    
    path_3=paste0("../ssa_output")
    path_4=paste0("../ssa_output/",specie)
    ssa_output_folder=paste0("../ssa_output/",specie,"/",year)
    if(!dir.exists(path_3)) dir.create(path=path_3)
    if(!dir.exists(path_4)) dir.create(path=path_4)
    if(!dir.exists(ssa_output_folder)) dir.create(path=ssa_output_folder)
    
   ### Run SSA #####
    cat(paste0("generating ssa output for haul ", haul_y[i],"\n"))
    command<-paste0("java -cp ../java/ssa.jar it.cnr.timeseries.analysis.test.Main ",
                    paste0(ssa_input_folder,"/",specie,"_haul_", haul_y[i],".csv")," year biomindex 1 0.01 6 true")
    SSA_execution<-system(command, intern = T,
                          ignore.stdout = FALSE, ignore.stderr = FALSE,
                          wait = TRUE, input = NULL, show.output.on.console = TRUE,
                          minimized = FALSE, invisible = TRUE)
    execution_success<-(length(which(grepl(pattern="OK SSA",x=SSA_execution)))>0)
    from_file<-c("./forecast.csv")
    to_file<-(paste0(ssa_output_folder,"/","forecast_","haul_", haul_y[i],".csv"))
    file.rename(from_file,to_file)
    from_file<-c("./forecast.png")
    to_file<-(paste0(ssa_output_folder,"/","forecast_","haul_", haul_y[i],".png"))
    file.rename(from_file,to_file)
    from_file<-c("./eigenvalues.png")
    to_file<-(paste0(ssa_output_folder,"/","eigenvalues_","haul_", haul_y[i],".png"))
    file.rename(from_file,to_file)
    from_file<-c("./signal_samples_u.png")
    to_file<-(paste0(ssa_output_folder,"/","signal_samples_u_","haul_", haul_y[i],".png"))
    file.rename(from_file,to_file)
    file.remove("eigenvalues.csv","signal_u.csv","signal_u.png","summary.txt")
    unlink("./*.tmp")
    if (!file.exists(paste0(ssa_output_folder,"/","forecast_","haul_", haul_y[i],".csv"))){
      df<-data.frame(time=year,value=0)
      write.csv(df,paste0(ssa_output_folder,"/","forecast_","haul_", haul_y[i],".csv"),row.names=FALSE)
    }# End if
  }# End haul_y and Run SSA
  
   ### Data preparation BIMAC ###
   cat(paste0("generating BIMAC for year ", year,", ", species,"\n"))
   data_4<-data[which(data$year==year),c("lon","lat","biomindex")]
   data_4_b<-data[which(data$year==year),c("station","lat", "lon","biomindex")]
   path_6=paste0("../BIMAC_input")
   path_7=paste0("../BIMAC_input/",specie)
   BIMAC_input_folder=paste0("../BIMAC_input/",specie,"/",year)
   if(!dir.exists(path_6)) dir.create(path=path_6)
   if(!dir.exists(path_7)) dir.create(path=path_7)
   if(!dir.exists(BIMAC_input_folder)) dir.create(path=BIMAC_input_folder)
   write.csv(data_4,paste0(BIMAC_input_folder,"/",specie,"_hauls_year_", year,".csv"), row.names=F)  
  
   path_9=paste0("../BIMAC_output")
   path_10=paste0("../BIMAC_output/",specie)
   BIMAC_output_folder=paste0("../BIMAC_output/",specie,"/",year)
   if(!dir.exists(path_9)) dir.create(path=path_9)
   if(!dir.exists(path_10)) dir.create(path=path_10)
   if(!dir.exists(BIMAC_output_folder)) dir.create(path=BIMAC_output_folder)
  
   ### Run BIMAC #####
   cat(paste0("run BIMAC for year ", year,"\n")) 
   specie_b<<-specie[1]
   year_b<<-year[1]
   source("./BIMAC_no_advection.R")
  
   #### Data preparation MaxEnt ##
   cat(paste0("generating MaxEnt for year ", year,", ", species,"\n")) 
   treshold<-max(data_4$biomindex)*0.05
   data_5_b<-data[which(data$year==year & data$biomindex > treshold),c("station","species","lat", "lon")]  
   data_6<-setdiff(data_4_b$station, data_5_b$station)
   cat(paste0(" cutted_hauls ",data_6, "\n"))
   data_5<-data[which(data$year==year & data$biomindex > treshold),c("species","lat", "lon")]
   data_5$species<-specie
   path_12=paste0("../MaxEnt_input")
   path_13=paste0("../MaxEnt_input/",specie)
   MaxEnt_input_folder=paste0("../MaxEnt_input/",specie,"/",year)
   if(!dir.exists(path_12)) dir.create(path=path_12)
   if(!dir.exists(path_13)) dir.create(path=path_13)
   if(!dir.exists(MaxEnt_input_folder)) dir.create(path=MaxEnt_input_folder)
   write.csv(data_5,paste0(MaxEnt_input_folder,"/",specie,"_presence_year_", year,".csv"), row.names=F)
  
   path_15=paste0("../MaxEnt_output")
   path_16=paste0("../MaxEnt_output/",specie)
   MaxEnt_output_folder=paste0("../MaxEnt_output/",specie,"/",year)
   if(!dir.exists(path_15)) dir.create(path=path_15)
   if(!dir.exists(path_16)) dir.create(path=path_16)
   if(!dir.exists(MaxEnt_output_folder)) dir.create(path=MaxEnt_output_folder)
  
   ### Run MAxEnt #####
   cat(paste0("run MaxEnt all variables for year ", year,", ", species,"\n")) 
   presence_data<-paste0(MaxEnt_input_folder,"/",specie,"_presence_year_", year,".csv")
   env_data_folder<-paste0("../data/Environmental_inputs/MaxEnt_",year)
   maxent_out<-paste0(MaxEnt_output_folder,"/")
   prevalence<-0.5
   command<-paste0("java -jar ../java/max_ent_cyb.jar \"",presence_data,"\" \"",
                   env_data_folder, "/\" \"",maxent_out,"\" ",prevalence)
   maxent_execution<-system(command, intern = T,
                           ignore.stdout = FALSE, ignore.stderr = FALSE,
                           wait = TRUE, input = NULL, show.output.on.console = TRUE,
                           minimized = FALSE, invisible = TRUE)
   execution_success<-(length(which(grepl(pattern="OK MaxEnt",x=maxent_execution)))>0)
  
   if (feature_selection){
    ## select important variables ##
    treshold_maxent<-read.csv(paste0(MaxEnt_output_folder,"/","maxentResults.csv"))
    treshold_1<-t(treshold_maxent)
    treshold_2<-as.data.frame(treshold_1)
    treshold_3<- cbind(variable = rownames(treshold_2), treshold_2)
    rownames(treshold_3) <- 1:nrow(treshold_2)
    colnames(treshold_3)<-c("variable", "value")
    treshold_data<-data.frame(treshold_3)
    treshold_4<-treshold_data[grepl('contribution', treshold_data$variable),]
    value<-as.numeric(treshold_4$value)
    variable<-as.character(treshold_4$variable)
    treshold_5<-data.frame(variable,value)
    treshold_6<-max(treshold_5$value)*0.05 ##treshold for all values over the 95% of the maximum value
    treshold_7<-treshold_5[which(treshold_5$value > treshold_6),] 
    treshold_7$variable<-gsub(".contribution",".asc",treshold_7$variable)
    temp_environment_folder=paste0("Environmental_inputs/MaxEnt_",year,"/","refine_variable","_",specie)
    if(!dir.exists(temp_environment_folder)) dir.create(path=temp_environment_folder)
     for (file in treshold_7$variable){
      file_origin<-paste0("../data/Environmental_inputs/MaxEnt_",year,"/",file)
      file_destination<-paste0(temp_environment_folder,"/")
      file.copy(file_origin,file_destination)
      } ## End select important variables

   ### second run maxent refine variables####
   cat(paste0("run MaxEnt refine variables for year ", year,", ", species,"\n")) 
   MaxEnt_output_folder_refine=paste0("./MaxEnt_output/",specie,"/",year,"_refine_variables","/")
   if(!dir.exists(MaxEnt_output_folder_refine)) dir.create(path=MaxEnt_output_folder_refine)
   presence_data<-paste0(MaxEnt_input_folder,"/",specie,"_presence_year_", year,".csv")
   env_data_folder_refine<-paste0(temp_environment_folder,"/")
   maxent_out_refine<-paste0(MaxEnt_output_folder_refine,"/")
   prevalence<-0.5
   command<-paste0("java -jar ../java/max_ent_cyb.jar \"",presence_data,"\" \"",
                  env_data_folder_refine, "/\" \"",maxent_out_refine,"\" " ,prevalence)
   maxent_execution<-system(command, intern = T,
                           ignore.stdout = FALSE, ignore.stderr = FALSE,
                           wait = TRUE, input = NULL, show.output.on.console = TRUE,
                           minimized = FALSE, invisible = TRUE)
   execution_success<-(length(which(grepl(pattern="OK MaxEnt",x=maxent_execution)))>0)
   }## End if second run MaxEnt
  
   ######### MAKE HBIE INPUT ################
   cat(paste0("generating HBIE input for year ", year,", ", species,"\n")) 
   
   ### exctract values from asc MaxEnt#####
   haul_m_2<-haul_m[which(haul_m$year==year),]
   if (feature_selection){
     asc_file<-raster(paste0(MaxEnt_output_folder_refine,"/",specie,".asc"))
   }else{
     asc_file<-raster(paste0(MaxEnt_output_folder,"/",specie,".asc"))
   }#End if
   min_x_in_raster<-asc_file@extent[1]
   max_x_in_raster<-asc_file@extent[2]
   min_y_in_raster<-asc_file@extent[3]
   max_y_in_raster<-asc_file@extent[4]
   resolution <- res(asc_file)[1]
   xseq<-seq(from=min_x_in_raster+(resolution/2),to=max_x_in_raster-(resolution/2),by=resolution)
   yseq<-seq(from=min_y_in_raster+(resolution/2),to=max_y_in_raster-(resolution/2),by=resolution)
   grid_of_points<-expand.grid(x = xseq, y = yseq)
   #extract values filling NA with the nearest neighbor#
   crs(asc_file)<- "+proj=utm + zone=33 +datum=WGS84 +no_defs"
   dist<-distance(asc_file)
   direct <- direction(asc_file, from=FALSE)
   rna <- is.na(asc_file)
   na.x <- init(rna, 'x')
   na.y <- init(rna, 'y')
   co.x <- na.x + dist * sin(direct)
   co.y <- na.y + dist * cos(direct)
   co <- cbind(co.x[], co.y[])
   NAVals <- raster::extract(asc_file, co, method='simple')
   r.NAVals <- rna
   r.NAVals[] <- NAVals
   r.filled <- cover(x=asc_file, y= r.NAVals)
   grid_of_points_extracted_fill_NA <- raster::extract(x = r.filled, y = grid_of_points,method='simple')
   grid_of_points$value<-grid_of_points_extracted_fill_NA
   crs(asc_file)<- "+proj=longlat +datum=WGS84 +no_defs"
   coordinate_at_res<-function(origin, coordinate, resolution){
     times<-round((coordinate-origin)/resolution)
     coordinate<-(times*resolution)+origin
     return (coordinate)
   }#End coordinate_at_res MaxEnt
   haul_m_2$lon_res<-coordinate_at_res(origin = min_x_in_raster+(resolution/2),
                                       coordinate = haul_m_2$lon, resolution = resolution)
   haul_m_2$lat_res<-coordinate_at_res(origin = min_y_in_raster+(resolution/2),
                                       coordinate = haul_m_2$lat, resolution =  resolution)
   haul_m_2$ecological<-0
   #i=1
   for (i in 1:nrow(haul_m_2)){
     extracted_value<-grid_of_points[which(grid_of_points$x==haul_m_2$lon_res[i] 
                                           & grid_of_points$y==haul_m_2$lat_res[i]),"value"]
     haul_m_2$ecological[i]<-extracted_value[1]
   }#End for haul_m_2 MaxEnt
  
  ### exctract values from asc BIMAC IDW #####
   asc_file<-raster(paste0(BIMAC_output_folder,"/","BIMAC_IDW_prior_no_advection_",year,"_",specie,".asc"))
   max_x_in_raster<-asc_file@extent[2]
   min_y_in_raster<-asc_file@extent[3]
   max_y_in_raster<-asc_file@extent[4]
   resolution <- res(asc_file)[1]
   xseq<-seq(from=min_x_in_raster+(resolution/2),to=max_x_in_raster-(resolution/2),by=resolution)
   yseq<-seq(from=min_y_in_raster+(resolution/2),to=max_y_in_raster-(resolution/2),by=resolution)
   grid_of_points<-expand.grid(x = xseq, y = yseq)
   grid_of_points_extracted<-extract(x=asc_file,y=grid_of_points,method='simple')
   grid_of_points$value<-grid_of_points_extracted   
   coordinate_at_res<-function(origin, coordinate, resolution){
    times<-round((coordinate-origin)/resolution)
    coordinate<-(times*resolution)+origin
    return (coordinate)
   }#End coordinate_at_res BIMAC IDW
   haul_m_2$lon_res<-coordinate_at_res(origin = min_x_in_raster+(resolution/2),
                                       coordinate = haul_m_2$lon, resolution = resolution)
   haul_m_2$lat_res<-coordinate_at_res(origin = min_y_in_raster+(resolution/2),
                                       coordinate = haul_m_2$lat, resolution = resolution)
   haul_m_2$spatial<-0
   #i=1
    for (i in 1:nrow(haul_m_2)){
       extracted_value<-grid_of_points[which(grid_of_points$x==haul_m_2$lon_res[i] 
                                             & grid_of_points$y==haul_m_2$lat_res[i]),"value"]
       haul_m_2$spatial[i]<-extracted_value[1]
    }#End for haul_m_2 BIMAD IDW
  
   ### exctract values from SSA #####
   cat(paste0("Estract values from SSA for year ", year,", ", species,"\n"))
   haul_m_2$temporal<-0
   #i<-1
   for (i in 1:length(haul_y)){
     temporal_1<-read.csv(paste0(ssa_output_folder,"/","forecast_","haul_", haul_y[i],".csv"))
     temporal_1[is.na(temporal_1)]<-0
     temporal_value<-tail(temporal_1$value, n=1)
     if (temporal_value[1]<0) temporal_value[1]=haul_m_2$spatial[i]
     if (is.nan(temporal_value[1])) temporal_value[1]=0
     haul_m_2$temporal[i]<-temporal_value[1]
   }#End for haul_y

   ### exctract values from asc BIMAC sd #####
   asc_file<-raster(paste0(BIMAC_output_folder,"/","BIMAC_interpolation_sd_no_advection_",year,"_",specie,".asc"))
   min_x_in_raster<-asc_file@extent[1]
   max_x_in_raster<-asc_file@extent[2]
   min_y_in_raster<-asc_file@extent[3]
   max_y_in_raster<-asc_file@extent[4]
   resolution <- res(asc_file)[1]
   xseq<-seq(from=min_x_in_raster+(resolution/2),to=max_x_in_raster-(resolution/2),by=resolution)
   yseq<-seq(from=min_y_in_raster+(resolution/2),to=max_y_in_raster-(resolution/2),by=resolution)
   grid_of_points<-expand.grid(x = xseq, y = yseq)
   grid_of_points_extracted<-extract(x=asc_file,y=grid_of_points,method='simple')
   grid_of_points$value<-grid_of_points_extracted   
   coordinate_at_res<-function(origin, coordinate, resolution){
     times<-round((coordinate-origin)/resolution)
     coordinate<-(times*resolution)+origin
     return (coordinate)
   }#End coordinate_at_res BIMAC sd
   haul_m_2$lon_res<-coordinate_at_res(origin = min_x_in_raster+(resolution/2),coordinate = haul_m_2$lon, resolution = resolution)
   haul_m_2$lat_res<-coordinate_at_res(origin = min_y_in_raster+(resolution/2),coordinate = haul_m_2$lat, resolution = resolution)
   haul_m_2$spatialerror<-0
   #i=1
    for (i in 1:nrow(haul_m_2)){
       extracted_value<-grid_of_points[which(grid_of_points$x==haul_m_2$lon_res[i] 
                                             & grid_of_points$y==haul_m_2$lat_res[i]),"value"]
       haul_m_2$spatialerror[i]<-extracted_value[1]
    }#End for haul_m_2 BIMAC sd
  
   hbie_input<-haul_m_2[,-c(2,3,4)]
   colnames(hbie_input)<-c("station","lon", "lat", "ecological","spatial","temporal","spatialerror")
   hbie_input[is.na(hbie_input)] <- 0
   hbie_input_folder=paste0("./hbie_input")
   if(!dir.exists(hbie_input_folder)) dir.create(path=hbie_input_folder)
   write.csv(hbie_input, paste0(hbie_input_folder,"/",specie,"_hbie_input_", year,".csv"), row.names=F)
   hbie_output_folder=paste0("./hbie_output")
   if(!dir.exists(hbie_output_folder)) dir.create(path=hbie_output_folder)
  
   ##### RUN HBIE ##########
   #########################
   cat(paste0("Run HBIE for year ", year,", ", species,"\n"))
   fileIdx = paste0(hbie_input_folder,"/",specie,"_hbie_input_", year,".csv")
   fileIdxOutput = paste0(hbie_output_folder,"/",specie,"_hbie_output_", year,".csv")
   #threshold MaxEnt percent omission rate########
    if (feature_selection){
      MaxEnt_htlm<- paste0(MaxEnt_output_folder_refine,"/",specie,".html")
    }else{
      MaxEnt_htlm<- paste0(MaxEnt_output_folder,"/",specie,".html")
    }#End if
   read_MaxEnt_htlm<-readChar(MaxEnt_htlm, file.info(MaxEnt_htlm)$size)
   omission_rate<-gsub(".*>Training omission rate</th><tr align=center><td>","",read_MaxEnt_htlm)
   omission_rate<- gsub("</td><td>Fixed cumulative value 1.*","",omission_rate)    
   omission_rate<- as.double(gsub(".*</td><td>","",omission_rate))
   maxEnt_Percent_Omission_Rate<-omission_rate
   cat(paste0("Omission rate MaxEnt for year ", year,", ", species," ="),maxEnt_Percent_Omission_Rate, "\n")
   
   hbie_output<-read.csv(file = fileIdx)
   alpha = 1
   beta = 1
   penalty = 0.4
   hbie_output$biomass_est_mean<-0
   hbie_output$biomass_est_low<-0
   hbie_output$biomass_est_high<-0
  
    for (i in 1:nrow(hbie_output)){
     haul<-hbie_output[i,]
     w<-1
     if (haul$ecological<maxEnt_Percent_Omission_Rate)
      w<-penalty
      b_mean<-w*( ((haul$spatial*alpha)+(haul$temporal*beta))/(alpha+beta))
      b_low<-w*( ( ( (haul$spatial-3*haul$spatialerror)*alpha)+(haul$temporal*beta))/(alpha+beta))
      b_low<-max(b_low,0)
      b_high<-w*( (((haul$spatial+3*haul$spatialerror)*alpha)+(haul$temporal*beta))/(alpha+beta))
      hbie_output$biomass_est_mean[i]<-b_mean
      hbie_output$biomass_est_low[i]<-b_low
      hbie_output$biomass_est_high[i]<-b_high
    } # End RUN HBIE
   write.csv(hbie_output,file=fileIdxOutput,row.names = F)
  
   ##### RUN biomass calculator ##########
   #######################################
   # Data preparation
   cat(paste0("Run biomass calculator for year ", year,", ", species,"\n"))
   biomass_calc<-data[which(data$year==year),]
   K_hauls<-biomass_calc[,c("station","lat","lon","biomindex")]
   colnames(K_hauls)<-c("station","lat","lon","biom_index")
   U_hauls<-read.csv(paste0(hbie_output_folder,"/",specie,"_hbie_output_", year,".csv"))
   U_hauls_mean<-U_hauls[,c("station","lat","lon","biomass_est_mean")]
   colnames(U_hauls_mean)<-c("station","lat","lon","biom_index")
   U_K_hauls_mean<-rbind(K_hauls,U_hauls_mean)
   U_hauls_low<-U_hauls[,c("station","lat","lon","biomass_est_low")]
   colnames(U_hauls_low)<-c("station","lat","lon","biom_index")
   U_K_hauls_low<-rbind(K_hauls,U_hauls_low)
   U_hauls_high<-U_hauls[,c("station","lat","lon","biomass_est_high")]
   colnames(U_hauls_high)<-c("station","lat","lon","biom_index")
   U_K_hauls_high<-rbind(K_hauls,U_hauls_high)
  
   path_18=paste0("../biom_index_input/")
   path_19=paste0("../biom_index_input/",specie)
   biom_index_input_folder=paste0("../biom_index_input/",specie,"/",year)
   if(!dir.exists(path_18)) dir.create(path=path_18)
   if(!dir.exists(path_19)) dir.create(path=path_19)
   if(!dir.exists(biom_index_input_folder)) dir.create(path=biom_index_input_folder)
   write.csv(K_hauls,paste0(biom_index_input_folder,"/",specie,"_known_", year,".csv"), row.names=F)
   write.csv(U_K_hauls_mean,paste0(biom_index_input_folder,"/",specie,"_mean_", year,".csv"), row.names=F)
   write.csv(U_K_hauls_low,paste0(biom_index_input_folder,"/",specie,"_low_", year,".csv"), row.names=F)
   write.csv(U_K_hauls_high,paste0(biom_index_input_folder,"/",specie,"_high_", year,".csv"), row.names=F)
  
   path_21=paste0("../biom_index_output/")
   path_22=paste0("../biom_index_output/",specie)
   biom_index_output_folder=paste0("../biom_index_output/",specie,"/",year)
   if(!dir.exists(path_21)) dir.create(path=path_21)
   if(!dir.exists(path_22)) dir.create(path=path_22)
   if(!dir.exists(biom_index_output_folder)) dir.create(path=biom_index_output_folder)
  
   indexes=data.frame(ncol=2,nrow=4)
   colnames(indexes)<-c("method","biom_index")
   tag_index<-1
  
   level_tag<-c("known","mean","low","high")
  
    for (tag in level_tag){
    #HaulBiomass <- read.csv("HaulBiomass.csv") # Biomass Index by Haul
      fileIdx = paste0(biom_index_input_folder,"/",specie,"_",tag,"_", year,".csv")
      data_calculator<-read.csv(file = fileIdx)
  
      calcBiomass<-function(HaulBiomass){
      StrataWeight <- read.csv("../data/StrataWeight.csv") # Area of strata
      HaulData <- read.csv("../data/HaulData.csv") # Swept Areas
    
      # Compute Index
      Index=merge(HaulBiomass, HaulData, by=c("year", "Station"))
      Index$BiomRaw=Index$BiomIndex*Index$SweptArea # return from Biomass Index to raw weight of species in the haul
      Index=Index[!is.na(Index$BiomRaw),]
      Index=aggregate(list(Biomass=Index$BiomRaw, SweptArea=Index$SweptArea), 
                    by = list(year=Index$year, Stratum=Index$Stratum), 
                    FUN=sum) # summarize weight and swept area by year and stratum
      Index$BiomassStratum=Index$Biomass/Index$SweptArea # Stratum Index 
      Index=merge(Index, StrataWeight, by="Stratum") 
      Index$BiomassStratumWeighted=Index$BiomassStratum*Index$StratumWeight 
      # Weight Stratum Index by the relative area of the stratum
      Index=aggregate(list(Index=Index$BiomassStratumWeighted), by=list(year=Index$year), FUN=sum) # Get Index
      return (Index)
      }#End calcBiomass
      
      indexes[tag_index,"method"]<-paste0("biomass_",tag)
      HaulBiomass<-data.frame(year=year,Station=data_calculator$station,BiomIndex=data_calculator$biom_index)
      TotalIndex<-calcBiomass(HaulBiomass)
      output_biom_calc<-(paste0("biomass_",tag," = ",TotalIndex[2]))
      print(paste0("biomass_",tag," = ",TotalIndex[2]))
      indexes[tag_index,"biom_index"]<-TotalIndex[2]
      tag_index=tag_index+1

    } # End Tag
   write.csv(indexes,paste0(biom_index_output_folder,"/indexes.csv"),row.names=F)

}# end species
}# end year


