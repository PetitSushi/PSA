library(sf)
library(tmap)
library(dplyr)
library(sp)
library(RColorBrewer)
library(spdep)
library(janitor)
library(ggplot2)
library(stringr)
library(spatialreg)

manchester_lsoa <-  st_read('raw/BoundaryData/england_lsoa_2011.shp')

st_crs(manchester_lsoa)

lsoa_WGS84 <- st_transform(manchester_lsoa, 4326)
st_crs(lsoa_WGS84)

plot(st_geometry(lsoa_WGS84))

rm(manchester_lsoa)


#  Read Attribute data into R
# List all 12 files on monthly crime data for Greater Manchester
listcrimedata <- fs::dir_info(here::here("raw", "greater_manchester_crime_data_2018"))%>%
  dplyr::filter(str_detect(path, ".csv")) %>%
  dplyr::select(path)%>%
  pull()%>%
  as.character()

# upload files and create dataframe 
data_csv = plyr::ldply(listcrimedata, read_csv)
data_csv

data_csv <- clean_names(data_csv)

#Filter out to select burglary
burglary_year <- filter(data_csv, crime_type == "Burglary")

#Transform into spatial object
burglary_spatial = st_as_sf(burglary_year, coords = c("longitude", "latitude"), 
                            crs = 4326, agr = "constant")

#Remove redundant non spatial burglary object
rm(burglary_year)
rm(data_csv)

#Select burglaries that intersect with the Manchester city LSOA map.
bur_mc <- st_intersects(lsoa_WGS84, burglary_spatial)

bur_mc <- burglary_spatial[unlist(bur_mc),]
#Check results
plot(st_geometry(bur_mc))

#Remove redundant objects
rm(burglary_spatial)

#Point in polygon spatial operation (be patient this can take time)
burglaries_per_lsoa <- bur_mc %>% 
  st_join(lsoa_WGS84, ., left = FALSE) %>% 
  count(lsoa_code)

#  rename the column with the count of burglaries (n) into something more meaningful
burglaries_per_lsoa <- rename(burglaries_per_lsoa, burglary = n)


# set CRS to BNG 27700 to get maps
st_crs(burglaries_per_lsoa)

lsoa_WGS84 <- st_transform(burglaries_per_lsoa, 27700)

# Plot with tmap 
# chloropleth map with quantile breaks

# tmap_mode("view")
tmap_mode("plot")

# quantile 
# could set breaks 1, 13, 20, 28, 39 (close to jenks)
tm_shape(burglaries_per_lsoa) + 
  tm_fill("burglary", style = "quantile", palette = "Reds", title = "Burglary counts") +
  tm_borders(alpha = 0.1) +
  tm_scale_bar(position=c(0.6,0.05), text.size=0.4)+
  tm_compass(north=0, position=c(0.65,0.15), text.size = 0.5)+
  tm_layout(main.title = "Burglary counts across Manchester's lsoas", main.title.size = 0.7 , 
            legend.outside =  TRUE, 
            legend.position = c("right", "bottom"), legend.title.size = 0.8)+
  tm_credits("(c) UK Data Service and Data Police UK", position=c(0.0,0.0))

# jenks breaks 
tm_shape(burglaries_per_lsoa) + 
  tm_fill("burglary", style = "jenks", palette = "Reds",  title = "Burglary counts") +
  tm_borders(alpha = 0.1) +
  tm_scale_bar(position=c(0.6,0.05), text.size=0.4)+
  tm_compass(north=0, position=c(0.65,0.15), text.size = 0.5)+
  tm_layout(main.title = "Burglary counts across Manchester's wards", main.title.size = 0.7 , legend.outside =  TRUE, 
            legend.position = c("right", "bottom"), legend.title.size = 0.8)+
  tm_credits("(c) UK Data Service and Data Police UK", position=c(0.0,0.0))


### Onto Spatial Autocorrelation analysis 
## Define neighbours

# calculate the centroids of all Wards in Manchester 
coordsW <- burglaries_per_lsoa %>%
  st_centroid()%>%
  st_geometry()

plot(coordsW,axes=TRUE)

#create a neighbours list
# Queen Case 
MWard_nb <- burglaries_per_lsoa %>%
  poly2nb(., queen=T)

#plot them
plot(MWard_nb, st_geometry(coordsW), col="red")
#add a map underneath
plot(burglaries_per_lsoa$geometry, add=T)

#create a spatial weights object from these weights
Mward.lw <- MWard_nb %>%
  nb2listw(., style="C")

head(Mward.lw$neighbours)

# Rook case 
MWard_nb_Rook <- burglaries_per_lsoa %>%
  poly2nb(., queen=F)

#plot them
plot(MWard_nb_Rook, st_geometry(coordsW), col="red")
#add a map underneath
plot(burglaries_per_lsoa$geometry, add=T)

#create a spatial weights object from these weights
Mward.lw_Rook <- MWard_nb_Rook %>%
  nb2listw(., style="C")


## Global Moran's I queen case
I_MWard_Global_Density <- burglaries_per_lsoa %>%
  pull(burglary) %>%
  as.vector()%>%
  moran.test(., Mward.lw)

I_MWard_Global_Density

# Global moran with queen case = 0.428 rather clustered, nothing dramatic 

## Global Moran's I Rook case
I_MWard_Global_Density_Rook <- burglaries_per_lsoa %>%
  pull(burglary) %>%
  as.vector()%>%
  moran.test(., Mward.lw_Rook)

I_MWard_Global_Density_Rook
# pretty much same result with Rook case:  0.429

rm(Mward.lw_Rook)
rm(MWard_nb_Rook)
rm(I_MWard_Global_Density_Rook)

# Monte Carlo for more robust model
mc_model <- moran.mc(burglaries_per_lsoa$burglary, Mward.lw, nsim=599)

# what do we get?
mc_model


# Getis-Ord global G statistic
G_MWard_Global_Density <- 
  burglaries_per_lsoa %>%
  pull(burglary) %>%
  as.vector()%>%
  globalG.test(., Mward.lw)

G_MWard_Global_Density
# Global G stat = 5.347451e-03    & expectation =  3.558719e-03
# The General G statistic = G > expected, so high values are tending to cluster


# Local Moran

I_MWard_Local_Density <- burglaries_per_lsoa %>%
  pull(burglary) %>%
  as.vector()%>%
  localmoran(., Mward.lw)%>%
  as_tibble()

# have a look at output (the localMoran object)
slice_head(I_MWard_Local_Density, n=5)

# copy the I score (column 1) and the z-score standard deviation (col 4)) back into the burglary df
burglaries_per_lsoa <- burglaries_per_lsoa %>%
  mutate(bur_localI =as.numeric(I_MWard_Local_Density$Ii))%>%
  mutate(bur_localIz =as.numeric(I_MWard_Local_Density$Z.Ii))

# plot a map of the local Moran’s I outputs

breaks1<-c(-1000,-2.58,-1.96,-1.65,1.65,1.96,2.58,1000)
MoranColours<- rev(brewer.pal(8, "RdGy"))

tm_shape(burglaries_per_lsoa) +
  tm_polygons("bur_localIz",
              style="fixed",
              breaks=breaks1,
              palette=MoranColours,
              midpoint=NA, title="Local Moran's I")+
  tm_scale_bar(position=c(0.6,0.05), text.size=0.4)+
  tm_compass(north=0, position=c(0.65,0.15), text.size = 0.5)+
  tm_layout(main.title = "Local Moran's I, Burglaries in Manchester", main.title.size = 0.9, legend.outside =  TRUE, legend.outside.position =  "right" )


##  Getis Ord G*i statisic for hot and cold spots 
Gi_MWard_Local_Density <- burglaries_per_lsoa %>%
  pull(burglary) %>%
  as.vector()%>%
  localG(., Mward.lw)

head(Gi_MWard_Local_Density)

burglaries_per_lsoa <- burglaries_per_lsoa %>%
  mutate(density_G = as.numeric(Gi_MWard_Local_Density))

GIColours<- rev(brewer.pal(8, "RdBu"))

tmap_mode("view")
#now plot G*i statisic on a  map
tm_shape(burglaries_per_lsoa) +
  tm_polygons("density_G",
              style="fixed",
              breaks=breaks1,
              palette=GIColours,
              midpoint=NA, title="Gi*, Burglaries ") + 
  tm_scale_bar(position=c(0.6,0.05), text.size=0.4)+
  tm_compass(north=0, position=c(0.65,0.15), text.size = 0.5)+
tm_layout(main.title = "Getis Ord G*i statistic, Burglaries in Manchester", 
          main.title.size = 0.8, legend.outside =  TRUE, legend.outside.position =  "right" )


# LISA 

# use the localmoran() function to store values
local_moran_manc_burglaries <- localmoran(burglaries_per_lsoa$burglary, Mward.lw)

# rescale that variable
burglaries_per_lsoa$scale_n_burglaries <- scale(burglaries_per_lsoa$burglary)

# create a spatial lag variable and save it to a new column
burglaries_per_lsoa$lag_scale_n_burglaries <- lag.listw(Mward.lw, burglaries_per_lsoa$scale_n_burglaries)

# run all of these through to assign variables
# for the statistical significance version assign a level 
# of statistical significance for the p value, column 5 of the local moran model
sig_level <- 0.1

# version with significance value
burglaries_per_lsoa$quad_sig <- ifelse(burglaries_per_lsoa$scale_n_burglaries > 0 & 
                                         burglaries_per_lsoa$lag_scale_n_burglaries > 0 & 
                                         local_moran_manc_burglaries[,5] <= sig_level, 
                                       "high-high", 
                                       ifelse(burglaries_per_lsoa$scale_n_burglaries <= 0 & 
                                                burglaries_per_lsoa$lag_scale_n_burglaries <= 0 & 
                                                local_moran_manc_burglaries[,5] <= sig_level, 
                                              "low-low", 
                                              ifelse(burglaries_per_lsoa$scale_n_burglaries > 0 & 
                                                       burglaries_per_lsoa$lag_scale_n_burglaries <= 0 & 
                                                       local_moran_manc_burglaries[,5] <= sig_level, 
                                                     "high-low", 
                                                     ifelse(burglaries_per_lsoa$scale_n_burglaries <= 0 & 
                                                              burglaries_per_lsoa$lag_scale_n_burglaries > 0 & 
                                                              local_moran_manc_burglaries[,5] <= sig_level, 
                                                            "low-high",
                                                            ifelse(local_moran_manc_burglaries[,5] > sig_level, 
                                                                   "not-significant", 
                                                                   "not-significant")))))

# version without significance check
burglaries_per_lsoa$quad_non_sig <- ifelse(burglaries_per_lsoa$scale_n_burglaries > 0 & 
                                             burglaries_per_lsoa$lag_scale_n_burglaries > 0, 
                                           "high-high", 
                                           ifelse(burglaries_per_lsoa$scale_n_burglaries <= 0 & 
                                                    burglaries_per_lsoa$lag_scale_n_burglaries <= 0, 
                                                  "low-low", 
                                                  ifelse(burglaries_per_lsoa$scale_n_burglaries > 0 & 
                                                           burglaries_per_lsoa$lag_scale_n_burglaries <= 0, 
                                                         "high-low", 
                                                         ifelse(burglaries_per_lsoa$scale_n_burglaries <= 0 & 
                                                                  burglaries_per_lsoa$lag_scale_n_burglaries > 0,
                                                                "low-high",NA))))


# plot the results without the satistical significance
ggplot(burglaries_per_lsoa, aes(x = scale_n_burglaries, 
                                y = lag_scale_n_burglaries, 
                                color = quad_non_sig)) +
  geom_vline(xintercept = 0) + # plot vertical line
  geom_hline(yintercept = 0) + # plot horizontal line
  xlab("Scaled Burglaries (n)") +
  ylab("Lagged Scaled Burglaries (n)") +
  labs(colour="Relative to neighbours") +
  geom_point()

# plot the results with the satistical significance
ggplot(burglaries_per_lsoa, aes(x = scale_n_burglaries, 
                                y = lag_scale_n_burglaries, 
                                color = quad_sig)) +
  geom_vline(xintercept = 0) + # plot vertical line
  geom_hline(yintercept = 0) + # plot horizontal line
  xlab("Scaled Burglaries (n)") +
  ylab("Lagged Scaled Burglaries (n)") +
  labs(colour="Relative to neighbours") +
  geom_point()

# map all of the results 
tmap_mode("plot")

tm_shape(burglaries_per_lsoa) +
  tm_fill(col = "quad_non_sig", title = "Similarity to neighbours (non sig)") + 
  tm_scale_bar(position=c(0.6,0.05), text.size=0.4)+
  tm_compass(north=0, position=c(0.65,0.15), text.size = 0.5)+
  tm_layout(main.title = "Local Moran, Relative similarity to neighbours, non stat. sig", 
            main.title.size = 0.8, legend.outside =  TRUE, legend.outside.position =  "right" )


# map only the statistically significant results here
tm_shape(burglaries_per_lsoa) +
  tm_fill(col = "quad_sig", title = "Similarity to neighbours (sig)") + 
  tm_scale_bar(position=c(0.6,0.05), text.size=0.4)+
  tm_compass(north=0, position=c(0.65,0.15), text.size = 0.5)+
  tm_layout(main.title = "Local Moran, Relative similarity to neighbours, stat. sig only", 
            main.title.size = 0.8, legend.outside =  TRUE, legend.outside.position =  "right" )


# Onto II - regression to infer about demographics and  burglaries in Manchester

# upload imd

imd_file <- read.csv('raw/imd/imd_2019_allscores.csv')

###cleaning

##  remove columns
# inspect the dataframe: column names
names(imd_file)

imd_sub <- select(imd_file, LSOA.code..2011., LSOA.name..2011., Index.of.Multiple.Deprivation..IMD..Score, Index.of.Multiple.Deprivation..IMD..Rank..where.1.is.most.deprived. ,
                  Index.of.Multiple.Deprivation..IMD..Decile..where.1.is.most.deprived.10..of.LSOAs.)

# select only manchester

#Filter out to select Manchester ward only
imd_manchester <- imd_sub %>%
  dplyr::filter(str_detect(LSOA.name..2011.,  "^Manchester"))

# remove redundant files
rm(imd_file)
rm(imd_sub)

imd_manchester <- clean_names(imd_manchester)

# join imd to sf lsoa file
st_crs(lsoa_WGS84)

imd_manchester <- rename(imd_manchester, lsoa11cd = lsoa_code_2011)
#burglaries_per_lsoa <- rename(burglaries_per_lsoa, lsoa11cd = lsoa_code)

burglaries_imd <- dplyr::left_join(burglaries_per_lsoa,imd_manchester,by=c("lsoa11cd"))

burglaries_imd <- rename(burglaries_imd, imd_dec = index_of_multiple_deprivation_imd_decile_where_1_is_most_deprived_10_of_lso_as)
burglaries_imd <- rename(burglaries_imd, imd_score = index_of_multiple_deprivation_imd_score)
burglaries_imd <- rename(burglaries_imd, imd_rank = index_of_multiple_deprivation_imd_rank_where_1_is_most_deprived)

rm(imd_manchester)

### EDA

## stat IMD
# summary imd deciles
summary(burglaries_imd$imd_dec)

# summary imd scores
summary(burglaries_imd$imd_score)

# boxplot imd deciles
# variable imd_decile = where_1_is_most_deprived_10_of_lsoas
ggplot(burglaries_imd, aes(x=imd_dec))+
  geom_boxplot()

# boxplot imd scores
ggplot(burglaries_imd, aes(x=imd_score))+
  geom_boxplot()

ggplot(burglaries_imd, aes(x=imd_score)) +
  geom_histogram(binwidth = 5)


# map values of IMD decile scores, to see which wards are most deprived (1=most deprived)
breaksdeciles<-c(1,2,3,4,5,6,7,8,9,10)

tmap_mode("view")
tmap_mode("plot")

tm_shape(burglaries_imd) + 
  tm_polygons('imd_dec', breaks=breaksdeciles, palette='-Reds',  title="IDM Deciles") +
  tm_scale_bar(position=c(0.6,0.05), text.size=0.4)+
  tm_compass(north=0, position=c(0.65,0.15), text.size = 0.5)+
  tm_layout(main.title = "Index of Multiple Deprivation deciles across Manchester", main.title.size = 0.9, legend.outside =  TRUE, legend.outside.position =  "right" )


## stat Burglaries
# boxplot of burglaries
ggplot(burglaries_imd, aes(x=burglary))+
  geom_boxplot()

# summary
summary(burglaries_imd$burglary)

# get output ?
ggplot(burglaries_imd, aes(x=burglary)) +
  geom_histogram(binwidth = 5)

# burglary data looks very  skewed, cube transform the burglary data 
# and inspect the transformed data again
burglaries_imd$burglary_cube <- abs(burglaries_imd$burglary)^(1/3)

# summary and box plot cubed burglary data
summary(burglaries_imd$burglary_cube)

# get output? 
ggplot(burglaries_imd, aes(x=burglary_cube))+
  geom_boxplot()

ggplot(burglaries_imd, aes(x=burglary_cube)) +
  geom_histogram(binwidth = 0.5)

# visual relationship between input (independent) variables and our output (dependent) variable 
# imd decile vs burglary number cubed
ggplot(burglaries_imd,aes(x=imd_dec, y=burglary_cube)) +
  geom_point()
# The IMD deciles, which are in fact ordinal data, do not show any obvious pattern

# imd score vs burglary count cubed
# get output
ggplot(burglaries_imd,aes(x=imd_score, y=burglary_cube)) +
  geom_point()



# imd score vs burglary count  
ggplot(burglaries_imd,aes(x=imd_score, y=burglary)) +
  geom_point()

# Correlation : imd decile vs burglary cubed
cor(burglaries_imd$imd_dec,burglaries_imd$burglary_cube)
# = 0.06502851

# Correlation : imd score vs burglary cubed 
cor(burglaries_imd$imd_score,burglaries_imd$burglary_cube)
# -0.09068189

# Correlation : imd score vs burglary
cor(burglaries_imd$imd_score,burglaries_imd$burglary)
# -0.1279016

# test  hypothesis using an Ordinary Least Squares (OLS) regression model.
# IMD deciles
# linear model
burglary_model <- lm(burglary_cube ~ imd_dec, data=burglaries_imd)

# get the results
summary(burglary_model)

#IMD scores 
burglary_model2 <- lm(burglary_cube ~ imd_score, data=burglaries_imd)

# get the results
summary(burglary_model2)
tidy(burglary_model2)

#IMD scores  & burglaries numbers
burglary_model3 <- lm(burglary ~ imd_score, data=burglaries_imd)
# get the results
summary(burglary_model3)
tidy(burglary_model3)

# Accounting for spatial heterogeneity
# extract the residuals from the model and assign to our spatial dataset
burglaries_imd$residuals <- residuals(burglary_model)
burglaries_imd$residuals2 <- residuals(burglary_model2)


# extract the predicted values from the model and assign to our spatial dataset
burglaries_imd$predicted <- fitted(burglary_model)
burglaries_imd$predicted2 <- fitted(burglary_model2)

# example observed, residual, predicted of first LSOA in our data
burglaries_imd[1,c(1,13,15,16)]

# standardise
burglaries_imd$sd_mean <- (burglaries_imd$predicted - mean(burglaries_imd$predicted)) / sd(burglaries_imd$predicted)

burglaries_imd$sd_mean2 <- (burglaries_imd$predicted2 - mean(burglaries_imd$predicted2)) / sd(burglaries_imd$predicted2)

summary(burglaries_imd$sd_mean)
summary(burglaries_imd$sd_mean2)

# inspect the result

tm_shape(burglaries_imd) + 
  tm_fill('sd_mean', style='jenks', palette='-RdBu') +
  tm_borders(alpha = 0.1) 


st_crs(burglaries_imd)
st_crs(burglaries_imd) = 4326
# get output 

summary(burglaries_imd$sd_mean2)

tm_shape(burglaries_imd) + 
  tm_fill('sd_mean2', style='jenks', palette='-RdBu') +
  tm_borders(alpha = 0.1) 

breakssd  <-c(-2.3,-1.4,-0.6,0.2,0.9,2)

tm_shape(burglaries_imd) + 
  tm_polygons("sd_mean2", style="fixed", breaks=breakssd,
              palette='-RdBu',
              midpoint=NA, title="Stand. Dev. Mean")+
  tm_scale_bar(position=c(0.6,0.05), text.size=0.4)+
  tm_compass(north=0, position=c(0.65,0.15), text.size = 0.5)+
  tm_layout(main.title = "Residuals Standard Deviation from the mean", 
            main.title.size = 0.8, legend.outside =  TRUE, legend.outside.position =  "right" )


# execute a Moran’s test for the regression residuals
lm.morantest(burglary_model, Mward.lw, alternative='two.sided')

lm.morantest(burglary_model2, Mward.lw, alternative='two.sided')


# row standardised weight matrix 
Mward_row_wm <- nb2listw(MWard_nb, style='W')

# lagrange mulitplier
lm.LMtests(burglary_model, Mward_row_wm, test = c('LMerr','LMlag','RLMerr','RLMlag'))


# execute the spatially lagged model
# this can take a little bit of time, be patient!
burglary_model_lag <- lagsarlm(burglary_cube ~ imd_dec, data=burglaries_imd, listw=Mward_row_wm)

#  inspect the results
summary(burglary_model_lag)

# execute the spatial error model
# this can take a little bit of time, be patient!
burglary_model_err <- errorsarlm(burglary_cube ~ imd_dec, data=burglaries_imd, listw=Mward_row_wm)

burglary_model_err2 <- errorsarlm(burglary_cube ~ imd_score, data=burglaries_imd, listw=Mward_row_wm)

# let's inspect the results
summary(burglary_model_err)
tidy(burglary_model_err)

summary(burglary_model)

# model 2
summary(burglary_model2)
tidy(burglary_model2)

summary(burglary_model_err2)
tidy(burglary_model_err2)
