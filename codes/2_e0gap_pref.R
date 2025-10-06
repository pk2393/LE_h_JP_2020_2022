rm(list=ls())

library(tidyverse);library(magrittr);library(ggplot2);library(grid);library(gridExtra)
palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

data<-readr::read_csv(file="./data_prefectures/LE_summ_pref_GH.csv")


data%<>%
  mutate(gap2221.total=total2022-total2021, 
         gap2120.total=total2021-total2020, 
         gap2019.total=total2020-total2019
         )%>%
  mutate(Prefecture=factor(stringr::str_to_title(pref)), .after=pref)%>%
  dplyr::arrange(gap2221.total)

x_order=rev(c(data$Prefecture))

#############
############# Total 
#############

data%>%
  ggplot(aes(x=Prefecture, y=gap2221.total, fill=gap2221.total>0) ) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("TRUE" = palette[11], "FALSE" = palette[2])) +
  scale_x_discrete(limits=x_order)+
  coord_flip() +ylim(-1.25,1.25)+
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(
      fill = "white"),
    legend.position = "none",
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.title.x = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  theme(axis.text = element_text(size = 12))+
  geom_hline(yintercept=c(0))+
  labs(title = "(C) 2021-22")+
  ylab("change (years)") +
  xlab("Prefecture")->gap2221


data%>%
  ggplot(aes(x=Prefecture, y=gap2120.total, fill=gap2120.total>0) ) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("TRUE" = palette[11], "FALSE" = palette[2])) +
  scale_x_discrete(limits=x_order)+
  coord_flip() +ylim(-1.25,1.25)+
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(
      fill = "white"),
    legend.position = "none",
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.title.x = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  theme(axis.text = element_text(size = 12))+
  geom_hline(yintercept=c(0))+
  labs(title = "(B) 2020-21")+
  ylab("change (years)") +
  xlab("Prefecture")->gap2120

data%>%
  ggplot(aes(x=Prefecture, y=gap2019.total, fill=gap2019.total>0) ) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("TRUE" = palette[11], "FALSE" = palette[2])) +
  scale_x_discrete(limits=x_order)+
  coord_flip() +ylim(-1.25,1.25)+
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(
      fill = "white"),
    legend.position = "none",
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.title.x = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  theme(axis.text = element_text(size = 12))+
  geom_hline(yintercept=c(0))+
  labs(title = "(A) 2019-20")+
  ylab("change (years)") +
  xlab("Prefecture")->gap2019


combined_plot <- arrangeGrob(gap2019, gap2120, gap2221, ncol=3, 
                                        top=textGrob("Life expectancy changes (years) by prefecture, Total Population", 
                                                     gp = gpar(fontface = "bold", fontsize = 20)))

ggsave("./figures/fig1_pref_e0_gap_GH.tiff", plot = combined_plot, dpi = 300, width = 15, height = 10, 
       compression = "lzw")
