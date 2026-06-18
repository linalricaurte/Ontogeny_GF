############################
#### figure 1: tracking data
############################

m(list = ls())
library(tidyverse)
library(lubridate)
library(amt)  
library(sf)  
library(rnaturalearth)
library(plotly)

#### figure 1: tracking data

##### read all GPS TARCK DATA SET
data<-read.csv("Data/AB/GF_resampled1h_fv_v20251021.csv", header = TRUE, dec = ".", sep=",") #I create this txt using the 02_data_process within the refinemnt folder


data$date_time <- as.POSIXct(strptime(data$date_time, format="%Y-%m-%d %H:%M:%S"), tz='UTC')
data$date <- as.Date(data$date_time)
## order dataframe chronologically per device
data<-data[order(data$dev,data$date_time),]

unique_individuals <- data%>%
  group_by(dev) %>%
  dplyr::summarise(unique_count = n_distinct(dev))

unique_individuals <- data%>%
  group_by(type2) %>%
  dplyr::summarise(unique_count = n_distinct(dev))


meta <- read.csv('Metadata/Metadata_PR_forR_v20241212_V3.csv', header = TRUE, dec = ".", sep=",")
library(plyr)
meta$dev <- revalue(meta$dev, c("8444"="008444"))
meta$dev <- revalue(meta$dev, c("5727"="005727"))
meta$dev <- revalue(meta$dev, c("6118"="006118"))
names(meta)[names(meta)=="country"] <- "origin"
meta$colony<-as.factor(meta$colony)

### remove dead birds and adult bird

meta$dev <- sprintf("%9s",meta$dev)
meta$dev<- gsub('^\\s+', '0', meta$dev)
meta$dev <- paste("B",meta$dev,sep="")

### remove dead birds =6 but for the staging paper 7 because I removed 233156 who died in ALgerie,  tag malfunction/not useful data=6 

### remove two tags deployed in birds that died, stop working and resident birds: in total 18 tags in the list and adult bird because we dont know the exact age
meta <- meta[ !(meta$dev %in% c("B0005727", "B0006118", "B0221225", "B0221228", "B0221235","B0233156",
                                "B0233157", "B0234766", "B0234768", "B0234772", 
                                "B0234772_b", "B0234780", "B0234782", "B0234782_b", 
                                "B0MALA30", "B0ISPR01", "B0ISPR02", "B0ISPR03", "B0MALA03")), ]

data <- data[ !(data$dev %in% c("B0233156")), ]



meta <- meta[order(meta$dev),]
data <- data[order(data$dev),]

unique_individuals <- data %>%
  group_by(type2) %>%
  dplyr::summarise(unique_count = n_distinct(dev))


unique_individuals <- meta %>%
  group_by(Tagging.year) %>%
  dplyr::summarise(unique_count = n_distinct(dev))

### Add tagging year column 
### Subset to tracks males and females

filtered_data <- meta %>%
  filter(Tagging.year == 2015)
data15 <- subset(data,data$dev %in% c("B0MALA01","B0MALA02","B0MALA04"))
table1<-as.data.frame(table(data15$dev))

filtered_data <- meta %>%
  filter(Tagging.year == 2016)
data16 <- subset(data,data$dev %in% c("B0FLAI01","B0FLAI02","B0FLAI03a","B0FLAI04" ,"B0FLAI05" ,
                                      "B0ISPR04","B0ISPR05" ,"B0ISPR06" ,"B0ISPR07" ,"B0ISPR08" ,"B0ISPR09" ,"B0ISPR10" ,"B0ISPR11a",
                                      "B0ISPR12" ,"B0ISPR13" ,"B0ISPR14" ,"B0ISPR15" ,"B0ISPR16"  ,"B0ISPR17" ,"B0ISPR18" ,"B0MALA05", 
                                      "B0MALA06" ,"B0MALA11" ,"B0MALA13" ,"B0MALA14" ,"B0MALA31" ,"B0MALA34" ,"B0SPOO01" ,"B0SPOO02", 
                                      "B0SPOO03" ,"B0SPOO04"))

table2<-as.data.frame(table(data16$dev))

filtered_data <- meta %>%
  filter(Tagging.year == 2017)

data17 <- subset(data,data$dev %in% c("B0FLAI03b","B0ISPR11b","B0MALA35"))

table3<-as.data.frame(table(data17$dev))


filtered_data <- meta %>%
  filter(Tagging.year == 2021)

data21 <- subset(data,data$dev %in% c("B0008444"))
table5<-as.data.frame(table(data21$dev))

filtered_data <- meta %>%
  filter(Tagging.year == 2022)

data22 <- subset(data,data$dev %in% c("B0221218","B0221220","B0221222","B0221224" ,"B0221226" ,"B0221227",
                                      "B0221229" ,"B0221230" , "B0221231" ,  "B0221236" , "B0221237" , "B04159A610", "B076C9AC10","B078C9AC10" ,"B07DC9AC10" ,"B084C9AC10", "B0WSIK67",  "B0WSIK82" ))
table6<-as.data.frame(table(data22$dev))

filtered_data <- meta %>%
  filter(Tagging.year == 2023)

data23 <- subset(data,data$dev %in% c("B0234764","B0234765" ,"B0234767" ,"B0234769", "B0234770", "B0234771" ,"B0234773" ,"B0234775","B0234776" ,"B0234778", "B0234779",  "B0234781", "B0234783", "B0234784", "B0234785" ,"B0234786","B0234788"))
table7<-as.data.frame(table(data23$dev))

filtered_data <- meta %>%
  filter(Tagging.year == 2024)

data24 <- subset(data,data$dev %in% c("B0233140", "B0233141" ,"B0233146", "B0233147","B0233148" ,"B0233150","B0233154" ,"B0234774","B0234777" ,"B0234787"))
table8<-as.data.frame(table(data24$dev))

### Add tag year column
data15$Tagging.year<-NA
data15$Tagging.year<-"2015"
data16$Tagging.year<-NA
data16$Tagging.year<-"2016"
data17$Tagging.year<-NA
data17$Tagging.year<-"2017"


data21$Tagging.year<-NA
data21$Tagging.year<-"2021"

data22$Tagging.year<-NA
data22$Tagging.year<-"2022"

data23$Tagging.year<-NA
data23$Tagging.year<-"2023"

data24$Tagging.year<-NA
data24$Tagging.year<-"2024"

### to bind two lists
data<-rbind(data15,data16,data17,data21,data22,data23,data24)

rm(table2,table3,table5,table6,table7,table8)

# Calculate distance of bird to center point of colony
source('Rfunctions/sidescript_pt2pt_fxns_v20230911.R')

# The pt2pt.range function returns distance in km 

# Calculate distance of bird to center point of colony
#source('sidescript_pt2pt_fxns.R')
xx <- pt2pt.range(data$LAT,data$LON,data$colony.lat,data$colony.long,threshold=100)
data$dist.to.colony <- xx$range.distance
data$near.colony <- xx$is.in.range


### Format date_time
data$Tagging.year2 <- as.Date(paste(data$Tagging.year, "-07-01", sep = ""), format = "%Y-%m-%d")
data$Year.tag <- format(data$Tagging.year2, "%Y")

### 
cycle.starts <- data %>%
  group_by(Tagging.year) %>%
  dplyr::summarize(cycle.start = as.Date(paste(min(year),"-09-01 00:00:01",sep=''))) 
cycle.starts <- as.data.frame(cycle.starts)

data <- merge(data,cycle.starts,all.x=TRUE)
data$cycle.start <- as.POSIXct(strptime(paste(data$cycle.start,"00:00:01",sep=" "), format="%Y-%m-%d %H:%M:%S"), tz='UTC')
data$cycle.time <- as.numeric(difftime(data$date_time,data$cycle.start,units="days"))
rm(cycle.starts)

data<-data[order(data$dev,data$date_time),]

###read data colonies
colonies<-read.csv("Data/AB/GF_colonies_fig1.csv", header = TRUE, dec = ".", sep=",") ### I create this txt using script 02_data_process and it includes movement metrics
colonies <- unique(colonies[,c("colony","colony.long","colony.lat")])
colonies <- colonies[order(colonies$colony.long),]
colonies$lab <- c("1","2","3","4")

####map study area and flamingo tracks
longlimits <- c(min(data$LON)-0.5,max(data$LON)+0.5)
latlimits <- c(min(data$LAT)-0.5,max(data$LAT)+0.5)

#grouping variable


world <- ne_countries(scale = "medium", returnclass = "sf")

###################################################
#####     PRODUCE FIG 2 Tracking data        ####     
###################################################
##### Fig 2a. Mapping movements

data <- data[order(data$dev,data$date_time),]

library(ggspatial)

p <- ggplot() +
  # Add countries
  geom_sf(data = world, fill = "grey50", col = "grey5", size = .9) +
  # Zoom map
  coord_sf(xlim = longlimits, ylim = latlimits, expand = FALSE) +
  theme_bw()
p

p1 <- p +
  annotation_scale(location = "bl",width_hint = 0.3,unit_category = "metric",text_col = "black",line_col = "white",
  bar_cols = c("black", "white")) +
  geom_text(data = data.frame(x = min(longlimits) + 2, y = max(latlimits) - 1.5), aes(x = x, y = y), label = "(a)", size = 4, fontface = "plain", hjust = 0) +
  xlab("Longitude") + ylab("Latitude") +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(-0.2, "cm"),
    panel.background  = element_rect(fill = "#1a5276"),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.title.x.top  = element_text(size = 13, face = "plain"),
    axis.text.x.top   = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.x       = element_text(size = 12, face = "plain"),
    axis.text.y       = element_text(size = 12, face = "plain"),
    axis.title        = element_text(size = 12, face = "plain"),
    axis.line         = element_line(linewidth = .4),
    axis.text         = element_text(size = 12),
    legend.title      = element_text(size = 12, face = "plain"),
    legend.text       = element_text(size = 10, face = "plain"),
    legend.position   = 'none',
    plot.margin       = unit(c(0,0,0,0), "cm")
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 3, title.position = 'left',
                          override.aes = list(size = 2)),
    size   = guide_legend(order = 2, nrow = 2, title.position = 'left'),
    fill   = guide_legend(order = 3, nrow = 3, title.position = 'left',
                          override.aes = list(alpha = 1)),
    alpha  = "none"
  )
p1

#### 
p2<-p1+
  geom_path(data=data,aes(x=LON,y=LAT,col=colony,group=dev),size=.4,alpha=.4)+
  # add colony coordinates
  geom_point(data=colonies,aes(x=colony.long,y=colony.lat,col=colony),size=3,shape=21,fill='transparent',stroke=1.8,alpha=.6)+
  geom_text(data = colonies, aes(x = colony.long + 1, y = colony.lat + 0.4, label = lab, color = colony), size = 5.5,
            fontface = "bold") +
  scale_color_manual(values = c("Aigues-Mortes" = "#FFC900", "Molentargius" = "orangered", "Comacchio" = "#56B4E9", "Margherita di Savoia" = "#AA4499"))+
  scale_fill_manual(values = c("Aigues-Mortes" = "#FFC900", "Molentargius" = "orangered", "Comacchio" = "#56B4E9", "Margherita di Savoia" = "#AA4499"))
p2

###############################
### Fig 2b. Time series
## annotate cycle since post-fledging date
data <- data[order(data$dev, data$date_time),]
data$md <- round(data$cycle.time)

x_pos <- 0
y_pos <- 4060

p3 <- ggplot() +
  
  geom_path(data = data, aes(x = md, y = dist.to.colony, col = colony, group = paste(dev, year)),
            linewidth = 0.7, alpha = 0.6) +
  scale_color_manual(values = c("Aigues-Mortes"= "#FFC900","Molentargius"= "orangered","Comacchio" = "#56B4E9",
    "Margherita di Savoia" = "#AA4499")) +
  scale_x_continuous(
    breaks = c(0, 365, 730, 1095, 1460, 1825, 2190, 2555, 2920),
    labels = c("1yo","2yo","3yo","4yo","5yo","6yo","7yo","8yo","9yo"),
    limits = c(0, 2900),
    position = 'bottom'
  ) +
  scale_y_continuous(limits = c(0, 4100)) +
  annotate("text", x = x_pos, y = y_pos,label = "(b)", size = 4, fontface = "plain", hjust = 0) +
  theme_classic() +
  xlab("Calendar year") + ylab("Distance to the natal colony (km)") +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(-0.2, "cm"),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.text.x       = element_text(size = 12, face = "plain"),
    axis.text.y       = element_text(size = 12, face = "plain"),
    axis.title.x      = element_text(size = 13, face = "plain"),
    axis.title.y      = element_text(size = 13, face = "plain"),
    legend.title      = element_text(size = 12, face = "plain"),
    legend.text       = element_text(size = 10, face = "plain"),
    legend.position   = 'none'
  )
p3

##################
#### Fig 2c. Traking data
unique_individuals <- data %>%
  group_by(type2) %>%
  dplyr::summarise(unique_count = n_distinct(dev))

unique_individuals <- data %>%
  group_by(type2  , colony) %>%
  dplyr::summarise(unique_count = n_distinct(dev))

# Use aggregate function to calculate the sum of 'dev' values for each 'year'
sum_dev <- aggregate(unique_count ~ type2, data =unique_individuals , FUN = sum)

# Now, create labels with year and sum of dev values
labels <- paste(sum_dev$type2, "\n(", sum_dev$unique_count, ")", sep = "")

library(paletteer)
# Define the color palette
my_palette <- paletteer_c("grDevices::YlGnBu", 13)

# Load required libraries
library(grid)
library(jpeg)

# Function to load a JPEG image
get_jpeg <- function(filename) {
  grid::rasterGrob(jpeg::readJPEG(filename), interpolate = TRUE)
}

# Load the image
gf <- get_jpeg("Figures_2025/Picture2.jpg")

x_pos_p4 <- 0.6
y_pos_p4  <- 86

p4 <- ggplot(data = unique_individuals, aes(x = factor(type2), y = unique_count, fill = factor(colony))) +
  geom_bar(stat = "identity", alpha = 0.8) +
  labs(x = "Age", y = "Number of individuals", fill = "Colony") +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(limits = c(0, 86), breaks = c(0, 20, 40, 60, 80)) +
  scale_fill_manual(values = c("Aigues-Mortes"="#FFC900", "Molentargius"= "orangered","Comacchio"= "#56B4E9", "Margherita di Savoia" = "#AA4499")) +
  annotate("text", x = x_pos_p4, y = y_pos_p4,
           label = "(c)", size = 4, fontface = "plain", hjust = 0) +
  theme_classic() +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(-0.2, "cm"),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.text.x       = element_text(size = 12, face = "plain"),
    axis.text.y       = element_text(size = 12, face = "plain"),
    axis.title.x      = element_text(size = 13, face = "plain"),
    axis.title.y      = element_text(size = 13, face = "plain"),
    legend.title      = element_text(size = 12, face = "plain"),
    legend.text       = element_text(size = 10, face = "plain"),
    legend.position   = "none",
    strip.text        = element_text(size = 12, face = "plain"),
    strip.background  = element_rect(fill = "white")
  ) +
  guides(
    alpha = "none",
    size  = guide_legend(nrow = 2)
  )

p4

gg <- p4 + annotation_custom(gf, xmin = 12, xmax = 3, ymin = 85, ymax = 16)
gg

######
#### Fig. 2d. Tracking duration
# Load the image
#install.packages("jpeg")
library(jpeg)
#gf2 <- readJPEG("Figures_2025/flamingos3.jpg")

gf2 <- get_jpeg("Figures_2025/flamingos3.jpg")

x_pos_p5 <- 4
y_pos_p5  <- nlevels(factor(data$dev)) * 0.97

p5 <- ggplot() +
  geom_point(data = data, aes(x = md, y = reorder(dev, cycle.time), color = colony),
             size = 0.4, stroke = 0.5, alpha = 0.8) + xlim(0, 3000) +
  scale_color_manual(values = c("Aigues-Mortes"= "#FFC900","Molentargius"= "orangered","Comacchio"= "#56B4E9","Margherita di Savoia" = "#AA4499"))+
  annotate("text", x = -Inf, y = Inf,label = "(d)", size = 4, fontface = "plain",hjust = -0.5, vjust = 1.5) +
  labs(x = "Tracking duration (days)", y = "Individual devices", fill = "Colony") +
  theme_classic() +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(-0.2, "cm"),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.text.x       = element_text(size = 12, face = "plain"),
    axis.text.y       = element_blank(),
    axis.title.x      = element_text(size = 13, face = "plain"),
    axis.title.y      = element_text(size = 13, face = "plain"),
    legend.title      = element_text(size = 12, face = "plain"),
    legend.text       = element_text(size = 10, face = "plain"),
    legend.position   = "none",
    strip.text        = element_blank(),
    strip.background  = element_rect(fill = "white")
  ) +
  guides(
    alpha = "none",
    size  = guide_legend(nrow = 2)
  )
p5

gg2 <- p5 + 
  annotation_custom(
    gf2, 
    xmin = 1000, xmax = 3000, 
    ymin = 0, ymax = 50 
  )

gg2
library(patchwork)
pfin <- (p2 | p3) / (gg | gg2)
ggsave(plot= pfin,filename = './Figures_2025/fig_2.tiff',dpi=300,width= 14,height= 11)


