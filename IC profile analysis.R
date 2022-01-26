# Title: Connectivity Index Longitudinal analysis

# Author details: Author: Loris Torresani & Guillaume Piton,

# Script and data info: This script performs a profile analysis on connectivity index (Cavalli et al.,2013)  and its components

# Copyright statement: This script is the product of Loris Torresani & Guillaume Piton

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


# Packages used
# install packages("ggplot2")
# install.packages("gridExtra")
# install.packages("tidyverse")
# install.packages("lmtest")
# install.packages("psych")
# install.packages("GGally")
# install.packages("ggpubr")
# install.packages("nls.multstart")
# install.packages("cowplot")
# install.packages("egg")
# install.packages("wesanderson")
# install.packages("devtools")
# install.packages("ellipsis")

# Load package
library("GGally")
library("ggplot2")
library("dplyr")
library("gridExtra")
library("foreign")
library("lmtest")
library("psych")
library("ggpubr")
library("stats")
library("compareGroups")
library("tidyverse")
library("rstatix")
library("Raster")
library("nls.multstart")
library("cowplot") # 1.0.0
library("egg")# 0.4.5
library("wesanderson")
library("devtools")

# INPUTS
# Set working directory to your location
# !!! WARNING !!!   do not write "/" as last characte

setwd("user defined working directory")

# Input dataset -> csv file with IC components extracted along a point sampling profile - user defined 
csvname<-"insert csv filename.csv"

#Read the data from the csv table
DATA <- data.frame(read.csv(csvname, header = TRUE, sep = ";",fileEncoding="UTF-8-BOM"))   #foreign::read.dbf (way to open a dbf file

#
###Sorting the data----  
DATA<-DATA[order(DATA$sorted),] #sort data from top to bottom or viceversa according to a user defined scale


###To compute distance----
DATA$Absc<-rep(NA,length(Tab$sorted))
DATA$Absc[1]<-0
for (i in 2:length(DATA$sorted))
{
  DATA$Absc[i]<-DATA$Absc[i-1]+((DATA$POINT_X[i]-DATA$POINT_X[i-1])^2+(DATA$POINT_Y[i]-DATA$POINT_Y[i-1])^2)^0.5
}

DATA$Absc<-max(DATA$Absc)-DATA$Absc
DATA$Domain <- NA

#---------------------------------------------------------------------

###Split profile into user defined domains----

######## defines interval of distance for create homogeneous domain. Distance refers to the "Absc" column in which the user 
#can use its personal value (expressed in "m")

# DATA$Domain<-Reach_mean$Domain<-NA
DATA$Domain[DATA$Absc<269]<-"A"
DATA$Domain[DATA$Absc>=269 & DATA$Absc<990]<-"B"
DATA$Domain[DATA$Absc>=990 & DATA$Absc<1250]<-"C"
DATA$Domain[DATA$Absc>=1250 & DATA$Absc<1960]<-"D"

DATA$Domain<-as.factor(DATA$Domain)

#---------------------------------------------------------------------

# detrend IC---

#to Check linear trend
LinearTrend=lm(DATA$CI~DATA$Absc)

#to Check a nonlinear trend
# define the equation
ICequation <- function(IC_0, Aa, B, Absc) { # equatio for fitting the IC and distance data
  
  eq <- IC_0 + Aa* (Absc **(B))
  
  return(eq)
}

nonLinearTrend<-nls_multstart(CI ~ ICequation(IC_0,Aa,B,Absc = Absc),
              data = DATA,
              iter = 500000,
              start_lower = c(IC_0=-10,Aa=-5, B=-1),
              start_upper = c(IC_0=10,Aa=5, B=4),
              lower = c(IC_0 = -10,Aa=-4, B = 0),
              trace = T,
              supp_errors = 'N')

#Check nonlinear model
print(nonLinearTrend)
formula(nonLinearTrend)

#Check linear model
LinearTrend=lm(DATA$CI~DATA$Absc)
print(LinearTrend)

#Control plot for trend verification
plot(DATA$Absc,DATA$CI,type="l");lines(DATA$Absc,predict(nonLinearTrend),col="red")

ggplot(data=DATA)+geom_line(aes(Absc,CI,color=Domain))+scale_colour_viridis_d("Domain Code",option="D")+theme_bw()+
  facet_grid(rows = vars(Domain))+
  # geom_boxplot(data=Reach,aes(Abscm,CIm,color=Domainm))
  geom_abline(aes(intercept=coef(LinearTrend)[1], slope=coef(LinearTrend)[2]))+
  geom_line(aes(x=DATA$Absc,y=predict(nonLinearTrend)),col="black")

#---------------------------------------------------------------------

###Plotting boxplots of different components---

#Boxplot weighting factor

ggplot(DATA,aes(Domain,W, fill= Domain))+
  geom_boxplot(varwidth=T,col="grey",lwd=0.2)+scale_fill_viridis_d("Domain Code",option="D")+theme_bw()+
  labs( x = "Domain code",y = "W [-]" 
  ) +guides(fill = "none")+
  theme(text = element_text(size=15))+ # ,caption=paste("(a)")
  # ,title = paste0("Connectivity index between domains")
  ylim(0,+1)
  ggsave("Boxplot_W.png", width = 3.14, height = 3.14)

  #Boxplot Roughness Index
  ggplot(DATA,aes(Domain,RI, fill= Domain))+
  geom_boxplot(varwidth=T,col="grey",lwd=0.2)+scale_fill_viridis_d("Domain Code",option="D")+theme_bw()+
  labs( x = "Domain code",y = "RI [m]" 
  ) +guides(fill = "none")+
  theme(text = element_text(size=15))+ # ,caption=paste("(a)")
  # ,title = paste0("Connectivity index between domains")
  ylim(0,+1)
  ggsave("Boxplot_RI.png", width = 3.14, height = 3.14)

  #Boxplot Slope
  ggplot(DATA,aes(Domain,S, fill= Domain))+
  geom_boxplot(varwidth=T,col="grey",lwd=0.2)+scale_fill_viridis_d("Domain Code",option="D")+theme_bw()+
  labs( x = "Domain code",y = "Slope [m/m]" 
  ) +guides(fill = "none")+
  theme(text = element_text(size=15))+ # ,caption=paste("(a)")
  # ,title = paste0("Connectivity index between domains")
  ylim(0,+1)
  ggsave("Boxplot_S.png", width = 3.14, height = 3.14)

  #Boxplot Upslope Component
  ggplot(DATA,aes(Domain,D_up, fill= Domain))+
  geom_boxplot(varwidth=T,col="grey",lwd=0.2)+scale_fill_viridis_d("Domain Code",option="D")+theme_bw()+
  labs( x = "Domain code",y = "D_up [-]" 
  ) +guides(fill = "none")+
  theme(text = element_text(size=15)) # ,caption=paste("(a)")
  # ,title = paste0("Connectivity index between domains")
   ggsave("Boxplot_D_Up.png", width = 3.14, height = 3.14)


#__________________________________________________________________ 

###Plotting IC and geomorphic components----

Profile<-ggplot()+ #Define data to plot
  geom_line(data=DATA,aes(x=Absc,y=Z,colour=Domain)) +scale_colour_viridis_d("Domain Code",option="D")+theme_bw()+
  theme(legend.position = "top")+
  theme(legend.direction = "horizontal")+
  guides(colour = guide_legend(nrow = 1))+#Where is the legend
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()
        ,plot.tag.position=c(0.97,0.4))+
  # annotate("text", x= 900 , y= 1800, label= "black dots indicate check dams") +#Remove axis
  labs(y = "Z [m.a.s.l]",tag=paste("(a)")
       # ,title = "Longitudinal profile"
       ,colour="Domain")+
  geom_linerange(data=DATA,aes(x=Absc,ymin=Z-Structures*50,ymax=Z-Structures*55),lwd=1.1)+
  geom_linerange(data=DATA,aes(x=Absc,ymin=min(Z)-Structures*20,ymax=min(Z)+Structures*10),lwd=1.1)

Profile

Slope<-ggplot(data = DATA,aes(x=Absc,y=S, color=Domain))+ #Define data to plot
  # geom_hline(aes(yintercept = 0,lwd=1))+ #Plot horizontal line at y =0
  scale_colour_viridis_d("Domain Code",option="D")+theme_bw()+
  geom_line()+
  #Plot a line of your data
  # scale_y_log10()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()
        ,legend.position = "none",plot.tag.position=c(0.97,0.5))+
  labs( x = "Distance [m]",y = "Slope [m/m]",tag=paste("(g)")
        # ,title = paste0("Slope")
  )
Slope


CI<-ggplot(data=DATA)+
  geom_line(aes(x=Absc,y=CI, color=Domain))+ #plot area
  scale_colour_viridis_d("Domain Code",option="D")+theme_bw()+
  #theme(legend.position = "top")+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()
        ,legend.position = "none",plot.tag.position=c(0.97,0.5))+
  labs( x = "Distance [m]",y = "IC [-]",tag=paste("(c)")
        # ,title = "Connectivity index"
  )+
  geom_line(aes(x=Absc,y=predict(nonLinearTrend)),col="grey")
  
CI


Area<-ggplot(data=DATA)+
  geom_line(aes(x=Absc,y=A/10^6, color=Domain))+ #plot area
  scale_colour_viridis_d("Domain Code",option="D")+theme_bw()+
  #theme(legend.position = "top")+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()
        ,legend.position = "none",plot.tag.position=c(0.97,0.5))+
  labs( x = "Distance [m]",y = expression("A [km"^2~"]") ,tag=paste("(b)")
        # ,title = "Contributing Area"
  )
Area

W_factor<-ggplot(data=DATA)+
  geom_line(aes(x=Absc,y=W, color=Domain))+ #plot area
  scale_colour_viridis_d("Domain Code",option="D")+theme_bw()+
  #theme(legend.position = "top")+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()
        ,legend.position = "none",plot.tag.position=c(0.97,0.5))+
  labs( x = "Distance [m]",y = bquote(W[x]~"[-]"),tag=paste("(d)")
        # ,title = "Weighting factor"
  )
W_factor
# 
RI<-ggplot(data=DATA)+
  geom_line(aes(x=Absc,y=RI, color=Domain))+ #plot area
  scale_colour_viridis_d("Domain Code",option="D")+theme_bw()+
  #theme(legend.position = "top")+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()
        ,legend.position = "none",plot.tag.position=c(0.97,0.5))+
  labs( x = "Distance [m]",y = bquote(RI[x]~"[-]"),tag=paste("(d)")
        # ,title = "Weighting factor"
  )
RI

D_up<-ggplot(data=DATA)+
  geom_point(aes(x=Absc,y=D_up/A^0.5, color=Domain))+ #plot area
  scale_colour_viridis_d("Domain Code",option="D")+theme_bw()+
  #theme(legend.position = "top")+
  
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()
        ,legend.position = "none",plot.tag.position=c(0.97,0.5))+
  labs( x = "Distance [m]",y = bquote("Dup [-]"),tag=paste("(e)")
        #,title = "Upslope component")
  )+geom_pointrange(aes(x=Absc,ymax=Structures*D_up/A^0.5,y=0.05*Structures,ymin=0.05*Structures),lwd=0.2,col="grey")+
  geom_point(aes(x=Absc,y=Structures*0.05),size=0.5)+
  coord_cartesian(ylim = c(0.05,1))
  # annotate("text", x= 900 , y= 0.08, label= "black dots indicate check dams")
  
D_up


D_dn<-ggplot(data=DATA)+
  geom_line(aes(x=Absc,y=D_down, color=Domain))+ #plot area
  scale_colour_viridis_d("Domain Code",option="D")+theme_bw()+
  #theme(legend.position = "top")+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()
        ,legend.position = "none",plot.tag.position=c(0.97,0.5))+
  labs( x = "Distance [m]",y = bquote("Ddn [-]"),tag=paste("(f)")
        #colour="Domain",
        #,title = "Downslope component")
  )
D_dn


png(paste0("Plot",csvname,".png"), width = 19, height = 30,units="cm",res=350)
{cowplot::plot_grid(Profile + theme(plot.margin = unit(c(-0.1,0.3,0.05,-0.5), "cm")
                                    ,legend.box.spacing=unit(c(0,0,0,0),"cm")
                                    ,legend.box.margin=unit(c(0,0,0,0),"cm")),
                   Area+theme(plot.margin=unit(c(0,0,0,-0.5), "cm")), 
                   CI+theme(plot.margin=unit(c(0,0,0,0.5), "cm")), 
                   W_factor+theme(plot.margin=unit(c(0,0,0,-0.5), "cm")),
                   D_up+theme(plot.margin=unit(c(0,0,0,-0.5), "cm")),
                   D_dn+theme(plot.margin=unit(c(0,0,0,-0.5), "cm")),
                   Slope+theme(plot.margin=unit(c(0,0,0,-0.5), "cm")),
                   nrow = 7,
                   align = c("hv"),
                   rel_heights = c(1,1,1)
                   )}

dev.off()

