library(monocle3)
library(dplyr)
library(here)
library(ggplot2)
library(ggdist)
library(ggpubr)

lmmp <- readRDS(here("6month_LMMP.rds"))

cell_types_select <- c("NENS","MENS")
dat <- pData(lmmp) %>% as.data.frame() %>% mutate(cell_type = stringr::str_replace(cell_type, "NENS", "NC-neurons"))

pl1 <- ggplot(filter(dat, cell_type %in% cell_types_select), aes(x = cell_type, y = Total_UMIs)) + 
  geom_boxplot() + 
  theme(
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        plot.title = element_text(hjust= 0.5)) + 
  ggtitle("UMI distribution in Kulkarni et al.")  

ulrika_filename <- "/data/users/jared/ulrika_ens/juvenile/shared_processing/data/juvenile_ENS.rds"
ulrika_ens <- readRDS(ulrika_filename)

ulrika_dat <- pData(ulrika_ens) %>% as.data.frame()
ulrika_dat["NC-neurons"] <- "NC-neurons"

pl2 <- ggplot(ulrika_dat, aes(x = "NC-neurons", y = Total_UMIs)) + 
  geom_boxplot() + 
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 13),
    plot.title = element_text(hjust= 0.5)) + 
  ggtitle("UMI distribution in Morarch et al.") + 
  ylab("") + 
  xlab("")

fig <- ggpubr::ggarrange(pl1, pl2)

ggsave(here("./plots/supp_figures/morarch_umi.png"))
