library(rgdal)
library(raster)
library(rNOMADS)

library(RSQLite)
library(dplyr)
library(dbplyr)

library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
library(plot3D)
library(ecmwfr)
library(keyring)
library(maptools)
library(ggplot2)

#https://rpubs.com/boyerag/297592

# 2. Import  temperature 


# wf_set_key(user="juan.vesga-gaviria@lshtm.ac.uk", service="cds")
wf_set_key(user="100002", key="2342ea70-4952-4578-ba44-3bb678c32da5", service="cds")


request <- list(
  variable = "soil_temperature_level_1",
  year = c("2016","2017","2018","2019","2020"),
  month = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"),
  day = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31"),
  time = "10:00",
  # area = c(35, 62.1, 34, 62.5),
  area = c(-27.78, 24.06, -29.20, 26.15),
  format = "netcdf",
  dataset = "reanalysis-era5-land",
  target = "download.nc"
)

request <- list("dataset_short_name" = "reanalysis-era5-land",
                "product_type"   = "reanalysis",
                "variable"       = "soil_temperature_level_1",
                "year"           = c("2016","2017","2018","2019","2020"),
                "month"          = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"),
                "day"            = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31"),
                "time"           = "10:00",
                "area"           = "-27.78/24.06/-29.20/26.15",
                "format"         = "netcdf",
                "target"         = "download.nc")



nc_stl1 <- wf_request(user = "100002",
                      request = request,   
                      transfer = TRUE,  
                      path = "~",
                      verbose = TRUE)


nstl1 <- nc_open(nc_stl1)
r_stl1 <- brick(nc_stl1) # brick overlay all the temporal maps we downloaded in one object easy to manipulate

# mask  dataset with France shapefile
# r_stl1_fra <- mask(r_stl1,FRA) 
# plot(r_stl1_fra[[3]]) # check if this works
# plot(FRA, add=TRUE) # perfect overlay of raster and shp

id_cell <- raster(r_stl1) 
id_cell[] <- 1:ncell(r_stl1[[1]])
summary(id_cell)
cat("\n", "Number of cells in raster: ", ncell(id_cell), "\n")  
mw <- as(r_stl1[[1]], "SpatialPixelsDataFrame")
class(mw)

grid_center <- coordinates(mw)
centroids <- SpatialPoints(grid_center, proj4string = crs(r_stl1))

############################################################################

par(mfrow=c(1,1))
plot(r_stl1[[1]])
plot(centroids,add=T, pch=6, cex=.2)


# Soil Temperature
T_data <- extract(r_stl1, centroids)
class(T_data)
str(T_data)
T_data[1:10, 1:5] # row: id cell number, col: variable


#convert into a dataframe
T_datat <- t(T_data) # transpose the matrix
Temp_df <- data.frame(T_datat)
# prepare labels for columns
year = c("2016","2017","2018","2019","2020")
month = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
# day = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31")
# time = c("00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00", "23:00")



lab <- seq(from=as.Date("2016-01-01"), to=as.Date("2020-12-31"), by="day")
names_T <- c(paste("T", sep="_", lab))
Temp_df$median <- apply(Temp_df, 1, median, na.rm = T)
Temp_df$IQR <- apply(Temp_df, 1, FUN = IQR, na.rm = T)
#Temp_df$max <- mapply(Temp_df,(Temp_df$X1:Temp_df$X121), FUN = max)
#Temp_df$min <- apply(Temp_df,1, FUN = min)
Temp_df$time <- c(1:nrow(Temp_df))
Temp_df$vector_days <- lab
#save dataset
temp_2017_2020 <- Temp_df[,c("median", "vector_days")]
temp_2017_2020$Temperature <- as.numeric(temp_2017_20200$median) - 273.15




write.csv(temp_2017_2020,here("data","soil_temp_2016_2020_FreState_SA.csv"), row.names=TRUE)


# Compare temperature range with Afghanistan
temp_day2= read.csv(here("data",'soil_temp_05_14_v2.csv'), 
                    header=TRUE, sep=',')

plot(seq(1,length(temp_2017_2020$Temperature)), 
     temp_2017_2020$Temperature,"l",col="blue",
     ylab=("Celsius"),xlab=("days"),ylim=(c(0,40)),main="Soil temperaature range")
lines(seq(1,length(temp_day2$Temperature)),
      temp_day2$Temperature,"l",col="red")
legend("topright", inset=c(-0.8,0), legend=c("Free State, SA", "Herat, Afghanistan"),
       col=c("blue", "red"), lty=1:2, cex=0.8)



#to export all the brick object

if (require(ncdf4)) {	
  rnc <- writeRaster(r_stl1, filename=here("data","r_stl1SA.nc"), format="CDF", overwrite=TRUE)   
}


# To bring it back  
FRA_stl1 <- brick(here("data","r_stl1SA.nc"))




