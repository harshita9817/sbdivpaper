#A
library(tidyverse)
library(patchwork)
library(data.table)
library(ggsci)
library(scales)
library(ggrepel)
library(ggformula)
library(maps)
library(geosphere)
library(scatterpie)
library(dplyr)
library(tidyverse)
library(readxl)
library(rMVP)
library(RColorBrewer)
library(randomcoloR)
#setwd("~/Documents/Sorghumseedimage/351_biallelic/GWAS_biallelic/SuppWorldfigure/")
world <- map_data("world")
data=readxl::read_excel("ithaseverthingyoucanthinkof.xlsx")
color=data[,c(1,27)]
b=as.data.frame(table(data$Country))
colnames(color)[1]="Pinumber"
cn <- data.frame(region=unique(data$Country))
cn <- world[world$region %in% cn$region,]
cn <- cn %>% 
  group_by(region) %>% 
  group_modify(~ data.frame(centroid(cbind(.x$long, .x$lat))))
colnames(b)[1] <- "region"
b <- merge(b, cn, by="region")
colnames(b)[2] <- "n"
world <- world[world$lat>-55,]
genotypeID=data$Genotype
color=color %>%
  filter(Pinumber %in% genotypeID)
color=color[order(color$Pinumber, decreasing = TRUE),]
color=color[!duplicated(color$Pinumber),]
colnames(data)[1]="Pinumber"
data2=merge(data, color, by="Pinumber")
data2$Race.y=gsub("Mixed", "Unspecified", data2$Race.y)
data2$Race.y=gsub("Others", "Unspecified", data2$Race.y)

data2=data2[!is.na(data2$Country),]
data2=data2 %>% 
  group_by(Country, Race.y) %>%
  summarise(count=n())
data2=data2%>%
  pivot_wider(names_from = Race.y, values_from = count)
data2[is.na(data2)] <- 0
summary(data2)
data3=data2[,-1]
data3= (data3/rowSums(data3))*100
data2=cbind(data2[,1], data3)
country=b$region
data2=data2%>% 
  filter(Country %in% country)
country=data2$Country
b=b%>% 
  filter(region %in% country)
colnames(data2)[1]="region"
b=merge(b,data2, by="region")

b$radius=5*abs(rnorm(b$n))/2
b$n[b$n>=10]=10

theme(default)

b$n=(b$n)/2

#runthistogetthemap
p1 <- ggplot(world, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), fill = "#69b3a2", colour = "lightgray", size = 0.1) + 
  theme(panel.background = element_rect(fill = 'lightblue'), legend.position = "top") +
  geom_scatterpie(data = b, aes(x = lon, y = lat, group = region, r = n), 
                  cols = c( "Breeding-lines","Durra","Guinea","Bicolor","Caudatum","Kafir","Unspecified"), 
                  sorted_by_radius = TRUE, color = NA, size = 0.1) + 
  coord_cartesian(xlim = c(-10, 150), ylim = c(-30, 65)) +  # Limiting longitude and latitude ranges
  xlab("Longitude") + 
  ylab("Latitude") +
  geom_scatterpie_legend(b$n, x = -140, y = -70, n = 5) +
  scale_fill_manual(values = c("darkblue", "purple", "darkred", "pink", "darkorange", "darkgreen", "yellow"), 
                    name = "Subpopulation")
  theme(
    title = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12), 
    legend.text = element_text(size = 12, colour = "black"),
    legend.position = "right",
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.title = element_text(angle = 90),
    legend.title.align = 0.5
  ) +
  guides(fill = guide_legend(title.position = "left", title.hjust = 0.5))+
  theme(legend.position = 'None')+
  coord_fixed(ratio = 1)  # Set the aspect ratio to 1 for equal axes

###########################################################################################################################
#B
#this will be pca for the different race

dfrace <- data[, c("Race", "PCA.1", "PCA.2", "PCA.3")]

# Handle missing values if any
dfrace <- na.omit(dfrace)

# Create a color palette for sub-populations
groups <- c("Durra","Guinea","Bicolor","Caudatum","Kafir","Unspecified")
groups.col <- c("darkblue", "purple", "darkred", "pink", "darkorange", "darkgreen")
color_palette <- setNames(groups.col, groups)

# Plot PC1 vs PC2 with different colors for each sub-population
p2<- ggplot(dfrace, aes(x = PCA.1, y = PCA.2, color = Race)) +
  geom_point(size = 2) +
  scale_color_manual(values = color_palette) +
  labs(#title = "Populations based on Race", 
       x = "PC 1", 
       y = "PC 2",
       color = "Race") +
  theme_minimal()+
  theme(
    legend.position = 'None', 
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),  # Add black lines for x and y axes
    axis.text.x = element_text(color = "black"),  # Change color of x-axis labels to black
    axis.text.y = element_text(color = "black")   # Change color of y-axis labels to black
  ) +
  coord_fixed(ratio = 1)  # Set the aspect ratio to 1 for equal axes  # Set the aspect ratio to 1 for equal axes
p3<- ggplot(dfrace, aes(x = PCA.2, y = PCA.3, color = Race)) +
  geom_point(size = 2) +
  scale_color_manual(values = color_palette) +
  labs(#title = "Populations based on Race", 
       x = "PC 2", 
       y = "PC 3",
       color = "Race") +
  theme_minimal()+
  theme(
    legend.position = 'None', 
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),  # Add black lines for x and y axes
    axis.text.x = element_text(color = "black"),  # Change color of x-axis labels to black
    axis.text.y = element_text(color = "black")   # Change color of y-axis labels to black
  ) +
  coord_fixed(ratio = 1)  # Set the aspect ratio to 1 for equal axes  # Set the aspect ratio to 1 for equal axes
p4<- ggplot(dfrace, aes(x = PCA.1, y = PCA.3, color = Race)) +
  geom_point(size = 2) +
  scale_color_manual(values = color_palette) +
  labs(#title = "Populations based on Race", 
       x = "PC 1", 
       y = "PC 3",
       color = "") +
  theme_minimal()+
  theme(
    legend.position = 'None', 
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),  # Add black lines for x and y axes
    axis.text.x = element_text(color = "black"),  # Change color of x-axis labels to black
    axis.text.y = element_text(color = "black")   # Change color of y-axis labels to black
  ) +
  coord_fixed(ratio = 1)  # Set the aspect ratio to 1 for equal axes # Set the aspect ratio to 1 for equal axes
################################################################################################

#this will be pca for the different race

dfregions <- data[, c("Regions", "PCA.1", "PCA.2", "PCA.3")]

# Handle missing values if any
dfregions <- na.omit(dfregions)

# Create a color palette for sub-populations
groupsregions <- c("Asia","Africa","Americas","Arabia/MiddleEast","Unspecified")
groups.colreg <- c("yellow", "blue", "green", "orange", "black", "gray")
color_palette <- setNames(groups.colreg, groupsregions)

# Plot PC1 vs PC2 with different colors for each sub-population
p5<- ggplot(dfregions, aes(x = PCA.1, y = PCA.2, color = Regions)) +
  geom_point(size = 2) +
  scale_color_manual(values = color_palette) +
  labs(#title = "Populations based on Regions", 
       x = "PC 1", 
       y = "PC 2",
       color = "Regions") +
  theme_minimal()+
  theme(
    legend.position = 'None', 
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),  # Add black lines for x and y axes
    axis.text.x = element_text(color = "black"),  # Change color of x-axis labels to black
    axis.text.y = element_text(color = "black")   # Change color of y-axis labels to black
  ) +
  coord_fixed(ratio = 1)  # Set the aspect ratio to 1 for equal axes  # Set the aspect ratio to 1 for equal axes

p6<- ggplot(dfregions, aes(x = PCA.2, y = PCA.3, color = Regions)) +
  geom_point(size = 2) +
  scale_color_manual(values = color_palette) +
  labs(#title = "Populations based on Race", 
    x = "PC 2", 
    y = "PC 3",
    color = "Regions") +
  theme_minimal()+
  theme(
    legend.position = 'None', 
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),  # Add black lines for x and y axes
    axis.text.x = element_text(color = "black"),  # Change color of x-axis labels to black
    axis.text.y = element_text(color = "black")   # Change color of y-axis labels to black
  ) +
  coord_fixed(ratio = 1)  # Set the aspect ratio to 1 for equal axes  # Set the aspect ratio to 1 for equal axes
p7<- ggplot(dfregions, aes(x = PCA.1, y = PCA.3, color = Regions)) +
  geom_point(size = 2) +
  scale_color_manual(values = color_palette) +
  labs(#title = "Populations based on Race", 
    x = "PC 1", 
    y = "PC 3",
    color = "") +
  theme_minimal()+
theme(
    legend.position = 'None', 
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),  # Add black lines for x and y axes
    axis.text.x = element_text(color = "black"),  # Change color of x-axis labels to black
    axis.text.y = element_text(color = "black")   # Change color of y-axis labels to black
  ) +
  coord_fixed(ratio = 1)  # Set the aspect ratio to 1 for equal axes # Set the aspect ratio to 1 for equal axes


###########################################################################################
#fixing plots:
  library(gridExtra)

# Define the layout
fig2top <- plot_grid(p2, p3, p4, nrow = 1, labels = c('B', 'C', 'D'), rel_widths = rep(1/3, 3), rel_heights = rep(1/3, 3))
fig2bottom <- plot_grid(p5, p6, p7, nrow = 1, labels = c('E', 'F', 'G'), rel_widths = rep(1/3, 3), rel_heights = rep(1/3, 3))

fig2 <- plot_grid(
  plot_grid(p1, ncol = 1, labels = 'A', rel_heights = c(0.5, 0.5)), 
  plot_grid(fig2top, fig2bottom, nrow = 2), 
  nrow = 2
)

# Plot the final layout
fig3 <- plot_grid(fig2, nrow = 1, rel_heights = c(0.9, 0.1))

# Display the plot
print(fig3)

ggsave('fig5.png', plot = fig3, width = 9.28, height = 9.08, units = 'in', dpi=1000)









