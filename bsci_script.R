#-----------------------------
# Install and load libraries
#-----------------------------
pkg_list = c('terra', 'randomForest', 'ggplot2', 'dtwclust', 'sf', 'caret')

installed_packages <- pkg_list %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(pkg_list[!installed_packages])
}

#-----------------------------
# Load Packages
#-----------------------------
lapply(pkg_list, function(p) {require(p,
                                      character.only = TRUE,
                                      quietly=TRUE)})
#
# install.packages('raster')
# library(raster)
# install.packages('cluster')
# library("cluster")
# install.packages('randomForest')
# library("randomForest")
# install.packages('ggplot2')
# library(ggplot2)
# install.packages('data.table')
# library(data.table)
# library(cluster)
# install.packages('clusterCrit')
# library(clusterCrit)
# install.packages('dtwclust')
# library(dtwclust)
# library(rgdal)
# library(dplyr)
# install.packages('ppclust')
# library(ppclust)
# library(dtw)
# install.packages("terra")
# library("terra")
# install.packages("sf")
# library("sf")
# install.packages("stars")
# library("stars")
# install.packages("rgdal")
# library("rgdal")
# install.packages("caTools")
# library("caTools")
# install.packages("caret")
# library("caTools")
# install.packages("RStoolbox")
# library("caTools")
# install.packages("reshape")
# library("reshape")
# install.packages("sClust")
# library("sClust")
# install.packages("eclust")
# library("eclust")
# install.packages("factoextra")
# library("factoextra")

#-----------------------------
# Load rasters into stack
#-----------------------------
Data_dir <-  "../" # Edit to reflect the path on your computer
tif_list <-  list.files(Data_dir, pattern=".tif", full.names = TRUE)
stk <-  rast(tif_list)
print("Dimensions of stack: "); print(dim(stk))
print(paste("Extent of raster:", ext(stk)))
print(paste("Number of cells:", ncell(stk)))

# Have a look
plot(stk$`MOD13A3.A2021305.h20v05.006.2021338162541.psrpgs_000501811712.1_km_monthly_NDVI-1_km_monthly_NDVI`)

#-----------------------------
# Take a random sample of points 
#-----------------------------
#Use a random sampling of points instead of the full raster
sample_pts <- spatSample(stk, 1000, na.rm = TRUE, as.points=TRUE)
length(sample_pts[,1])
# Get only the values (for clustering:
sample_df = as.data.frame(sample_pts)

#-----------------------------
# Time series clustering
#-----------------------------
# Even with the very small sample of points, this takes a long time!
clust <-  tsclust(sample_df, k=3)
# and join back to the points

sample_result = cbind(sample_pts, data.frame(clust@cluster))
plot(sample_result, col=sample_result$clust.cluster)

#-----------------------------
# Done
#-----------------------------



# bsci_test_list<-list.files(path = "E:/NORTHERN_NEGEV/PROCESSING/NDVI/testpam/ndvi_pam", pattern= ".tif", all.files=TRUE, full.names=FALSE)
# dir<-"E:/NORTHERN_NEGEV/PROCESSING/NDVI/testpam/ndvi_pam/"
# bsci_raster1 <- lapply(paste0(dir,bsci_test_list), raster)
# bsci_rasterstack<- stack(bsci_raster1)
# bsci_rasterstack<- rast(bsci_rasterstack)
#
# shp_nn<-shapefile("E:/NORTHERN_NEGEV/ARCMAP_docs/shapefile/Export_Output_forest.shp")
# full_bsci<-mask(bsci_rasterstack, shp_nn)


years<-c(1984:2011, 2013:2020)
names(full_bsci)<-as.character(years)

#get values for each pixel
v <- getValues(full_bsci)
i <- which(!is.na(v))
v <- na.omit(v)


#pam clustering
set.seed(123)
vx<-v[sample(nrow(v), 10000)]
distdtw<-dist(vx, method = "dtw_basic")
clust_pam_bsci<-pam(distdtw, 3, diss=T)
clust_pam_bsci$silinfo$clus.avg.widths
clust_pam_bsci$silinfo$avg.width


#RF model
library(randomForest)
library(caTools)
library(caret)

set.seed(122)
vx<-as.data.frame(vx)
vx$sample<-sample.split(vx$X1984, SplitRatio = 0.7)
vx$cluster<-clust_pam_ci$clustering
pam_train<-subset(vx, sample==TRUE)
pam_test<-subset(vx, sample==FALSE)

cluster_train<-pam_train$cluster
cluster_test<-pam_test$cluster

pam_train1<-pam_train[,-c(31,33)]
pam_test1<-pam_test[,-c(31,33)]

rf_train <- randomForest(pam_train1,as.factor(cluster_train),ntree = 1000)

rf_raster_ci<- predict(object=bsci_val, model=rf_train)
print(rf_train)

pam_test$pred<-predict(rf_train,pam_test)
pam_test$pred<-as.factor(pam_test$pred)
pam_test$cluster<-as.factor(pam_test$cluster)
confusionMatrix(pam_test$pred, pam_test$cluster)

dev.off()
#visualisation

library(RStoolbox)
fullbsci<-as.data.frame(fullbsci)
clustering_pam3<-ggplot(fullbsci, geom_raster = T)+
  scale_fill_discrete(labels=paste("Cluster", 1:3), guide_legend(reverse=T), na.translate = FALSE,)+
  scale_color_manual(labels=paste("Cluster", 1:3),values=c("palegreen3", "dodgerblue3", "orchid3"),
                      scale_y_continuous(limits = c(-1, 1)),
                      aesthetics = c("colour", "colour"),
                      na.translate = FALSE)+
  theme_bw()+
  labs(y="", x = "")+
  theme(legend.title = element_blank())

fullbsci

###
#plot(rf_raster_ci, col=brewer.pal(5,"Set1"))

####################################################################################################
library(data.table)
library(ggplot2)

vx1<-vx[,-c(32,33)]
data_plotci <- data.table(melt(data.table(cluster = as.factor(clust_pam_ci$clustering),
                                          vx1)))
data_plotci[, Time := rep(c(1987:2010,2013:2018), each = nrow(vx))]
data_plotci[, ID := rep(1:nrow(vx), ncol(vx))]

#centroids
med_val<-vx[c(clust_pam_ci$medoids),]
med_val$cluster<-c("Cluster 1", "Cluster 2", "Cluster 3")
centers_ci_ci <- data.table(melt(med_val))
labs<-c("Cluster 1", "Cluster 2", "Cluster 3")
levels(centers_ci_ci$cluster)<-labs
centers_ci_ci$Time<- rep(c(1984:2011, 2013:2020), each=3)

#plot
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
par(mar=c(10,6,10,4), mai=c(8, 1, 2, 1))
ci_over_time<-ggplot(centers_ci_ci, aes(Time, value, group = cluster)) +
  geom_line(data = centers_ci_ci, aes(Time, value),
            color="darkgreen", alpha = 0.70, size = 1) +
  labs(x = "\n Time \n", y = "\nCI\n") + ylim(c(0,1))+
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        strip.text = element_text(size=18))+
  facet_wrap(~cluster, ncol = 1, scales = "free_y")





ci_over_time1<-ggplot(centers_ci_ci, aes( x = Time, y = value, group = cluster ,colour=variable)) +
  geom_line(data = centers_ci_ci$cluster=="Cluster 1", aes(Time, value),
            color = "palegreen3", alpha = 0.70, size = 1) +
  geom_line(data = centers_ci_ci$cluster=="Cluster 2", aes(Time, value),
            color = "dodgerblue3", alpha = 0.70, size = 1) +
  geom_line(data = centers_ci_ci$cluster=="Cluster 3", aes(Time, value),
            color = "orchid3", alpha = 0.70, size = 1) +
  labs(x = "\n Time \n", y = "\nCI\n") + ylim(c(0,1))+
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        strip.text = element_text(size=18))

plot(centers_ci_ci$Time, centers_ci_ci$variable)

####
#summary statistics
as.numeric(med_val)
mean(as.numeric(med_val))
mean(as.numeric(med_val[2,-32]))
mean(as.numeric(med_val[3,-32]))

sd(as.numeric(med_val[1,-32]))
sd(as.numeric(med_val[2,-32]))
sd(as.numeric(med_val[3,-32]))

sd(as.numeric(med_val[1,-32]))/mean((as.numeric(med_val[1,-32])))
sd(as.numeric(med_val[2,-32]))/mean((as.numeric(med_val[2,-32])))
sd(as.numeric(med_val[3,-32]))/mean((as.numeric(med_val[3,-32])))

pacf(as.numeric(med_val[1,-32]), lag.max = 15)
pacf(as.numeric(med_val[2,-32]), lag.max = 15)
pacf(as.numeric(med_val[3,-32]), lag.max = 15)

cor(t(med_val$cluster[,-32]))

t.test(x=as.numeric(med_val[2,-32]), y=as.numeric(med_val[3,-32]), paired = T)
t.test(x=as.numeric(med_val[2,-32]), y=as.numeric(med_val[1,-32]), paired = T)
t.test(x=as.numeric(med_val[1,-32]), y=as.numeric(med_val[3,-32]), paired = T)

ks.test(x=as.numeric(med_val[2,-32]), y=as.numeric(med_val[1,-32]))
ks.test(x=as.numeric(med_val[3,-32]), y=as.numeric(med_val[1,-32]))
ks.test(x=as.numeric(med_val[2,-32]), y=as.numeric(med_val[3,-32]))

##################################################################################################################
library(corrplot)
#change to medoids object, not csv file
cluster_meds<-read.csv("G:/klil1/artical/cluster_medoids.csv")
med_cor<-cor(cluster_meds[,c(2:4)])

tiff("G:/klil1/artical/correlation_plot.tif", width = 30, height=28, res=100, units="cm")
corrplot(med_cor, method="color", addCoef.col = "black", type="upper",
         tl.col="black", diag=T, tl.cex=2, number.cex = 1.5, cl.cex=1.5)
dev.off()
