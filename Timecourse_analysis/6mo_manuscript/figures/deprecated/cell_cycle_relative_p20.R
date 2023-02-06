library(monocle3)
library(tidyverse)
library(here)
library(tricycle)
library(circular)
library(RColorBrewer)
library(ggnewscale)

source(here("scripts/accessory_functions/monocle_mods.R"))

cds <- readRDS(here("tc_neuro_p20.rds"))

theta.v <- cds$tricyclePosition

#do density on each cell type separately 
color_var.v <- cds$cell_type_aggregate %>% as.factor()
palette.v <- brewer.pal(nlevels(color_var.v), "Set1")

density_list <- map(levels(color_var.v), function(cv){
  d <- density.circular(circular(theta.v[which(color_var.v == cv)]), bw = 10)
  
  return(data.frame(x = as.numeric(d$x), y = d$y, color = cv))
})

density_df <- do.call(rbind, density_list)

#do density on all cells to get mean
d <- density.circular(theta.v, bw = 10)
all_df <- data.frame(x = as.numeric(d$x), y = d$y)
all_df[,"color"] <- "all"

density_df_all <- rbind(density_df, all_df)

ggplot(density_df_all, aes(x = x, y = y, color = color)) + geom_line()

#Densities sum to one
map(levels(color_var.v), function(cv){
  
  sub_den <- density_df[density_df$color == cv,]
  pracma::trapz(sub_den$x, sub_den$y)
  
})

relative_density_df <-map(levels(color_var.v), function(cv){
  
  sub_den <- density_df[density_df$color == cv,]
  data.frame(x = sub_den$x, y = sub_den$y - all_df$y, color = cv)
  
}) %>% do.call(rbind, .)

ggplot(relative_density_df, aes(x = x, y = y, color = color)) + geom_line() +
  scale_color_manual(values = palette.v,name = "cell_type", labels = levels(color_var.v)) +
  scale_x_continuous(limits = c(0,2*pi), breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
                     labels = paste0(c(0, 0.5, 1, 1.5, 2), "π"), name = "θ") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ylab("Relative density to overall mean")

#Densities sum to zero
map(levels(color_var.v), function(cv){
  
  sub_den <- relative_density_df[density_df$color == cv,]
  pracma::trapz(sub_den$x, sub_den$y)
  
})

png(here("plots/6mo_manuscript/tricycle_relativedensity_p20.png"), width = 700, height = 450)
ggplot(relative_density_df, aes(x = x, y = y, color = color)) + geom_line() +
  scale_color_manual(values = palette.v,name = "cell_type", labels = levels(color_var.v)) +
  scale_x_continuous(limits = c(0,2*pi), breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
                     labels = paste0(c(0, 0.5, 1, 1.5, 2), "π"), name = "θ") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ylab("Relative density to overall mean")
dev.off()




