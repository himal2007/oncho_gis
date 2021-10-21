setwd("C:/Users/User/OneDrive - LA TROBE UNIVERSITY/oncho_git/onchocerciasis/espen_ethiopia_data/analysis/210723_analysis_writing/")


# Load libraries ----------------------------------------------------------
suppressMessages({
  library(tidyverse);
  library(sf);
  library(tmap);
  library(raster);
  library(RColorBrewer);
  library(INLAutils);
  ## for bivariate maps
  library(classInt);
  library(rgdal);
  library(dismo);
  library(XML);
  library(maps);
  library(sp);
  
})

extrafont::loadfonts(device="win")

# Load data ---------------------------------------------------------------
data <- read.csv("data/eth_espen_nodule_covs_1km.csv")
data$prevalence <- data$CASES/data$N *100

m_0 <- getData(name = "GADM", country = "ETH", level = 0)
m_0 <- m_0 %>% st_as_sf()
m_1 <- getData(name = "GADM", country = "ETH", level = 1)
m_1 <- m_1 %>% st_as_sf()
mypalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))

load("docs/INLA_model_run_5km.RData")


# Bivariate map functions-----------------------------------------------------------

colmat<-function(nquantiles=10, upperleft=rgb(0,150,235, maxColorValue=255), upperright=rgb(130,0,80, maxColorValue=255), bottomleft="grey", bottomright=rgb(255,230,15, maxColorValue=255), xlab="x label", ylab="y label"){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3)
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]}

col.matrix<-colmat(nquantiles=10, upperleft="blue", upperright="yellow", bottomleft="green", bottomright="red", xlab="My x label", ylab="My y label")

col.matrix<-colmat(nquantiles=4)

## Bivariate map function
### START COPYING HERE ###
bivariate.map<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=10){
  quanmean<-getValues(rasterx)
  temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
  brks<-with(temp, quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r1<-within(temp, quantile <- cut(quanmean, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr<-data.frame(r1[,2]) 
  quanvar<-getValues(rastery)
  temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
  brks<-with(temp, quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r2<-within(temp, quantile <- cut(quanvar, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr2<-data.frame(r2[,2])
  as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
  col.matrix2<-colormatrix
  cn<-unique(colormatrix)
  for(i in 1:length(col.matrix2)){
    ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
  cols<-numeric(length(quantr[,1]))
  for(i in 1:length(quantr[,1])){
    a<-as.numeric.factor(quantr[i,1])
    b<-as.numeric.factor(quantr2[i,1])
    cols[i]<-as.numeric(col.matrix2[b,a])}
  r<-rasterx
  r[1:length(r)]<-cols
  return(r)}


# Plot bivariate maps -----------------------------------------------------

my.colors = colorRampPalette(c("white","lightblue", "yellow","orangered", "red"))
plot(r_prev_mean,frame.plot=F,axes=F,box=F,add=F,legend.width=.5,legend.shrink=.5,col=my.colors(255)) 
plot(st_geometry(m_0), add= T)

plot(r_prev_sd,frame.plot=F,axes=F,box=c("#D62F48CA", "#FFFFFF", "#FFFFFF")F,add=F,legend.width=.5,legend.shrink=.5,col=my.colors(255)) 
plot(st_geometry(m_0), add= T)

bivmap<-bivariate.map(r_prev_mean,r_prev_sd, colormatrix=col.matrix, nquantiles=4)

png(filename = "docs/bivariate_index.png", width = 4, height = 4, units = "in", res= 1000, bg = "transparent")
col.matrix<-colmat(nquantiles=4, upperleft="lightgreen", upperright="#F0842B5E", 
                   bottomleft="#6EA2F0", bottomright="#B31212D4", xlab="Prevalence", ylab="Uncertainty" )
dev.off()

#B31212D4 - red
#2F69D6CA - blue
#F0842B5E - orange

png(filename = "docs/bivariate_map.png", width = 10, height = 10, units = "in", res= 1000, bg = "transparent")
bivmap<-bivariate.map(r_prev_mean,r_prev_sd, colormatrix=col.matrix, nquantiles=4)
plot(bivmap,frame.plot=F,axes=F,box=F,add=F,legend=F,col=as.vector(col.matrix))
plot(st_geometry(m_0), add= T, lwd = 4)
dev.off()

c("#FF000033", "#FFFFFF", "#FFFFFF")
c("#EBB9B9D6", "#E8A190", "#D12A2A")
c("#F5C4C4", "#F29191", "#F58282", "#FA4141")
c("#2A6EEB", "#6EA2F0", "#B8CCF2", "#ABC5F5") - Blues
