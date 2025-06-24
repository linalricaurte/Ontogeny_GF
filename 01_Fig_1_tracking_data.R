############################
#### figure 1: tracking data
############################

#### read all 1-h resampled GPS tracking data set


data<-read.csv("Data/GF-fULL_resampled_v20250212.csv", header = TRUE, dec = ".", sep=",") #


data$date_time <- as.POSIXct(strptime(data$date_time, format="%Y-%m-%d %H:%M:%S"), tz='UTC')
data$date <- as.Date(data$date_time)

## order dataframe chronologically per device
data<-data[order(data$dev,data$date_time),]

unique_individuals <- data%>%
  group_by(dev) %>%
  dplyr::summarise(unique_count = n_distinct(dev))

### read metadata

meta <- read.csv('Metadata/Metadata_PR_forR_v20241212_V4.csv', header = TRUE, dec = ".", sep=",")

#### format metadata dev
meta$dev <- revalue(meta$dev, c("8444"="008444"))
meta$dev <- revalue(meta$dev, c("5727"="005727"))
meta$dev <- revalue(meta$dev, c("6118"="006118"))
names(meta)[names(meta)=="country"] <- "origin"
meta$colony<-as.factor(meta$colony)



meta$dev <- sprintf("%9s",meta$dev)
meta$dev<- gsub('^\\s+', '0', meta$dev)
meta$dev <- paste("B",meta$dev,sep="")

####Filter out birds with no useful data
### remove dead birds =7 including 233156 who died in Algeria,  6 tags malfunction/not useful, 5 resident birds and 1 adult bird

meta <- meta[ !(meta$dev %in% c("B0005727", "B0006118", "B0221225", "B0221228", "B0221235","B0233156",
                                "B0233157", "B0234766", "B0234768", "B0234772", 
                                "B0234772_b", "B0234780", "B0234782", "B0234782_b", 
                                "B0MALA30", "B0ISPR01", "B0ISPR02", "B0ISPR03", "B0MALA03")), ]

data <- data[ !(data$dev %in% c("B0233156")), ] #### remove 233156 who died in Algeria from the df



meta <- meta[order(meta$dev),]
data <- data[order(data$dev),]


### Add tagging year column to df

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
                                      "B0221229" ,"B0221230" , "B0221231" ,  "B0221236" , "B0221237" , "B4159A610", "B76C9AC10","B78C9AC10" ,"B7DC9AC10" ,"B84C9AC10", "B0WSIK67",  "B0WSIK82" ))
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

xx <- pt2pt.range(data$LAT,data$LON,data$colony.lat,data$colony.long,threshold=100)
data$dist.to.colony <- xx$range.distance
data$near.colony <- xx$is.in.range

### Format date_time for plotting
data$Tagging.year2 <- as.Date(paste(data$Tagging.year, "-07-01", sep = ""), format = "%Y-%m-%d")
data$Year.tag <- format(data$Tagging.year2, "%Y")

### compute cycle start for plotting
cycle.starts <- data %>%
  group_by(Tagging.year) %>%
  dplyr::summarize(cycle.start = as.Date(paste(min(year),"-09-01 00:00:00",sep=''))) 
cycle.starts <- as.data.frame(cycle.starts)

data <- merge(data,cycle.starts,all.x=TRUE)
data$cycle.start <- as.POSIXct(strptime(paste(data$cycle.start,"00:00:00",sep=" "), format="%Y-%m-%d %H:%M:%S"), tz='UTC')
data$cycle.time <- as.numeric(difftime(data$date_time,data$cycle.start,units="days"))
rm(cycle.starts)


## order df chronologically per device
data<-data[order(data$dev,data$date_time),]


#### define map limits based on flamingo tracks
longlimits <- c(min(data$LON)-0.5,max(data$LON)+0.5)
latlimits <- c(min(data$LAT)-0.5,max(data$LAT)+0.5)


### Use the rnaturalearth package to download the world map from Natural Earth data
world <- ne_countries(scale = "medium", returnclass = "sf")

###read data colonies
colonies<-read.csv("Data/GF_colonies_fig1.csv", header = TRUE, dec = ".", sep=",") ### I create this txt using script 02_data_process and it includes movement metrics
colonies <- unique(colonies[,c("colony","colony.long","colony.lat")])
colonies <- colonies[order(colonies$colony.long),]
colonies$lab <- c("1","2","3","4")


###################################################
#####     Produce figure 1                     ####     
###################################################

### Map Figure 1a
# Create the base plot

p <- ggplot() +
  # Add countries
  geom_sf(data = world, fill = "grey50", col = "grey5", size = .9) +

  coord_sf(xlim = longlimits, ylim = latlimits, expand = FALSE) +
  # Add theme
  theme_bw()

# Plot the map
p


p1 <- p+ theme_bw() + 

  xlab("") + ylab ("") +
  theme_dark()+
  theme(
    plot.background = element_rect(fill = 'transparent', colour = NA),
    panel.background = element_rect(fill = 'black', colour = NA),
    legend.position   = 'none',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x.top = element_text(size = 13, face = 'bold'),  # Add x-axis title on top
    axis.text.x.top = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1),  # Adjust x-axis text on top
    axis.line         = element_line(size=.4),
    axis.text		      = element_text(size=12),
    
    plot.margin=unit(c(0,0,0,0), "cm"), 
    axis.title		      = element_text(size=12,face='bold'))+
  guides(colour= guide_legend(order =1,nrow =3,title.position='left', override.aes = list(size = 2)),
         size = guide_legend(order =2,nrow =2,title.position='left'),
         fill=guide_legend(order=3,nrow =3,title.position='left',override.aes= list(alpha = 1 )),
         alpha=FALSE)



p1 

p2<-p1+
  
  ## Add tracks
  geom_path(data=data,aes(x=LON,y=LAT,col=colony,group=dev),size=.4,alpha=.4)+
  
  # add colony coordinates
  geom_point(data=colonies,aes(x=colony.long,y=colony.lat,col=colony),size=3,shape=21,fill='transparent',stroke=1.8,alpha=.6)+
  geom_text(data = colonies, aes(x = colony.long + 1, y = colony.lat + 0.4, label = lab, color = colony), size = 5.5,
            fontface = "bold") +
  
  scale_color_manual(values = c("Aigues-Mortes" = "#FFC900", "Molentargius" = "orangered", "Comacchio" = "#56B4E9", "Margherita di Savoia" = "#AA4499"))+
  
scale_fill_manual(values = c("Aigues-Mortes" = "#FFC900", "Molentargius" = "orangered", "Comacchio" = "#56B4E9", "Margherita di Savoia" = "#AA4499"))

p2


###Figure 1b
## annotate cycle since post-fledging date
data <- data[order(data$dev, data$date_time),]
data$md <- round(data$cycle.time)

#######################################
#plot for timing
p3 <- ggplot()+
  
  geom_path(data=data,aes(x=md,y=dist.to.colony,col=colony,group=paste(dev,year)),size=0.7,alpha=0.6)+
  
  # colours
  
  scale_color_manual(values = c("Aigues-Mortes" = "#FFC900", "Molentargius" = "orangered", "Comacchio" = "#56B4E9", "Margherita di Savoia" = "#AA4499"))+
  #scale_x_continuous(breaks=c(0,366,732,1098,1464,1830,2196,2562),labels=c("1yo","2yo", "3yo", "4yo","5yo","6yo", "7yo","8yo"),limits=c(0,2980),position='bottom')+
  scale_x_continuous(breaks=c(0,365,730,1095,1460,1825,2190,2555,2920),labels=c("1yo","2yo", "3yo", "4yo","5yo","6yo", "7yo","8yo", "9yo"),limits=c(0,2900),position='bottom')+
  
  scale_y_continuous(limits=c(0,4000))+
  theme_classic()+
  xlab("Calendar year") + ylab ("Distance to the natal colony (km)") + 
  theme(plot.background	= element_rect(fill='transparent',colour='NA'),
        legend.position   = 'none',
        panel.grid.minor.x   = element_blank(),
        panel.grid.major.x   = element_blank(),
        panel.grid.minor.y   = element_line(size=.2,linetype='dashed',colour='grey'),
        panel.grid.major.y   = element_line(size=.2,linetype='dashed',colour='grey'),
        axis.text		      = element_text(size=12),
        axis.title		      = element_text(size=13))

p3

### Figure 1c

unique_individuals <- data %>%
  group_by(type2  , colony) %>%
  dplyr::summarise(unique_count = n_distinct(dev))


# Use aggregate function to calculate the sum of 'dev' values for each 'year'
sum_dev <- aggregate(unique_count ~ type2, data =unique_individuals , FUN = sum)

# Now, create labels with year and sum of dev values
labels <- paste(sum_dev$type2, "\n(", sum_dev$unique_count, ")", sep = "")


# Define the color palette
my_palette <- paletteer_c("grDevices::YlGnBu", 13)


# Function to load a JPEG image
get_jpeg <- function(filename) {
  grid::rasterGrob(jpeg::readJPEG(filename), interpolate = TRUE)
}

# Load the image
gf <- get_jpeg("Figures_2024/photos/Picture2.jpg")


p4<-ggplot(data=unique_individuals,aes(x = factor(type2), y= unique_count, fill = factor(colony))) +
  geom_bar(stat = "identity", alpha=0.8)+ 
  labs( x = "Age", y = "Number of individuals", fill = "type2")+
  scale_x_discrete(labels = labels) +  
  scale_fill_manual(values = c("Aigues-Mortes" = "#FFC900", "Molentargius" = "orangered", "Comacchio" = "#56B4E9", "Margherita di Savoia" = "#AA4499"))+
  theme(panel.background = element_rect(fill = alpha("#FFFFFF",0.2))) + 
  theme(panel.grid.major	= element_blank(),
        panel.grid.minor	= element_blank(),
        axis.text	    = element_text(size=12),
        axis.title		    = element_text(size=13),
        axis.line = element_line(colour = "black"),
        legend.position  = "none",
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=10),
        strip.text		    = element_text(size=12,face='bold'),
        strip.background = element_rect(fill='white'))+
  guides(alpha=FALSE,
         size=guide_legend(nrow =2))


p4


# Define the plot limits if not defined already
latlimits <- range(unique_individuals$unique_count)


# Add falcon image to plot
gg <- p4 + annotation_custom(gf, xmin = 5+3, xmax = 1+3, ymin = 57+40, ymax = latlimits[1])
gg

gf2 <- readJPEG("Figures_2024/photos/flamingos3.jpg")

gf2 <- get_jpeg("Figures_2024/photos/flamingos3.jpg")

### Figure 1c

p5 <- ggplot()+
  ## Add tracks
  geom_point(data=data,aes(x=md,y=reorder(dev,cycle.time), color=colony),size=0.4,stroke = 0.5, alpha=.8)+
  xlim(0,3000)+
  scale_color_manual(values = c("Aigues-Mortes" = "#FFC900", "Molentargius" = "orangered", "Comacchio" = "#56B4E9", "Margherita di Savoia" = "#AA4499"))+
  labs( x = "Tracking duration (days)", y = "Individual devices", fill = "type2")+
  theme(panel.background = element_rect(fill = alpha("#FFFFFF",0.2))) + 
  theme(panel.grid.major	= element_blank(),
        panel.grid.minor	= element_blank(),
        axis.text		    = element_text(size=12),
        axis.title		    = element_text(size=13),
        axis.text.y		    = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position  = "",
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=10),
        strip.text		    = element_blank(),
        strip.background = element_rect(fill='white'))+
  guides(alpha=FALSE,
         size=guide_legend(nrow =2))

p5

gg2 <- p5 + 
  annotation_custom(
    gf2, 
    xmin = 1000, xmax = 3000,  
    ymin = -10, ymax = 50      
  )

gg2


# Create the top row with four panels
top_row <- cowplot::plot_grid(p2, p3, labels = c("A", "B"),
                              label_size = 16, ncol = 2,rel_heights = c(2, 1))

top_row


bottom_row <- cowplot::plot_grid(gg, gg2, labels = c("C", "D"),
                              label_size = 16, ncol = 2)

bottom_row

# Combine the top row and the bottom panel
pfin <- cowplot::plot_grid(top_row, bottom_row, label_size = 15, 
                           nrow = 2) # Adjust heights as needed

# Display the final plot
pfin


ggsave(plot=pfin,filename='./Figures_2025/fig_1_v20250318.tiff',dpi=300,width=14,height=11)
