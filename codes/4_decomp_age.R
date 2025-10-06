rm(list=ls())
library(tidyverse);library(magrittr)

load(file="life_table_2000_2022.RData")
palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")



decomp_age_func=function(year,  sex){
  
  Jp_next <- lt_national_for_arriaga[[paste0(year, ".", sex)]]
  Jp_prev <- lt_national_for_arriaga[[paste0(year-1, ".", sex)]]
  
  
  Jp_next$Age[22]<-"100+"
  Jp_prev$Age[22]<-"100+"
  
  ###################################################################################
  JData=data.frame(Age = Jp_prev$Age,
                   Jp_prev.lx = Jp_prev$lx,
                   Jp_prev.Lx = Jp_prev$Lx,
                   Jp_prev.Tx=Jp_prev$Tx,
                   
                   Jp_next.lx = Jp_next$lx,
                   Jp_next.Lx = Jp_next$Lx,
                   Jp_next.Tx = Jp_next$Tx
  )
  
  JData %<>% mutate(A=rep(NA, nrow(.)))
  for(i in 1:nrow(JData)){
    JData$A[i]=(JData$Jp_prev.lx[i]/100000)*(((JData$Jp_next.Lx[i])/(JData$Jp_next.lx[i]))-((JData$Jp_prev.Lx[i])/(JData$Jp_prev.lx[i])))
  }
  
  JData %<>% mutate(B=rep(NA,nrow(JData)))
  for(i in 1: nrow(JData)){
    JData$B[i] = (JData$Jp_next.Tx[i+1]/100000)*(((JData$Jp_prev.lx[i])/(JData$Jp_next.lx[i]))-((JData$Jp_prev.lx[i+1])/(JData$Jp_next.lx[i+1])))
  }
  
  JData %<>%
    mutate(
      nCx=case_when(
        is.na(B)==T ~ A, 
        T ~ A+B))
  
  
  ### save csv
  readr::write_csv(JData,paste0("./data_cause/Cx_GH", year-2000, ".", sex, ".csv"))
  
  ### create figure
  JData %>%
    select(Age, nCx)%>% 
    mutate(Age = factor(Age,Age))->JData.plot
  
  
  plot_age=ggplot(JData.plot,aes(x=Age, y=`nCx`, fill=nCx>0) ) +
    geom_bar(stat="identity") + 
    scale_fill_manual(values = c("TRUE" = palette[11], "FALSE" = palette[2])) +
    coord_flip() + ylim(-.1,.1)+
    
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(
        fill = "white"),
      
      legend.position = "none",
      axis.line = element_line(colour = "black")
    ) +
    theme(title=element_text(size = 24, color="black"))+
    theme(axis.title.x=element_text(size = 24, color="black"))+
    theme(axis.title.y=element_text(size = 24, color="black"))+
    theme(axis.ticks.x=element_line(color="black"))+
    theme(axis.text=element_text(size=24))+
    xlab("Age Group (years)") +
    ylab(paste0("Contribution by Age Group (years)"))
  
  
  
  res_list=list(
    data = JData.plot, 
    plot = plot_age)
  
  return(res_list)
  
}


decomp_age_func(2022, "Total")$plot+labs(title="(C) 2021-22")->decomp_age_total_2021_2022
decomp_age_func(2021, "Total")$plot+labs(title="(B) 2020-21")->decomp_age_total_2020_2021
decomp_age_func(2020, "Total")$plot+labs(title="(A) 2019-20")->decomp_age_total_2019_2020
decomp_age_func(2019, "Total")

library(grid);library(gridExtra)
combined_plot <- arrangeGrob(decomp_age_total_2019_2020, 
                             decomp_age_total_2020_2021, 
                             decomp_age_total_2021_2022,
                             ncol=3, 
                             top=textGrob("Contribution to Life Expectancy Changes (years) by Age Group: Total Population", 
                                          gp = gpar(fontface = "bold", fontsize = 32)))

plot(combined_plot)
ggsave("./figures/fig3_decomp_age_GH.tiff", plot = combined_plot, dpi = 300, width = 24, height = 12, 
       compression = "lzw")


#### male ####
decomp_age_func(2022, "Male")$plot+labs(title="(C) 2021-22")->decomp_age_male_2021_2022
decomp_age_func(2021, "Male")$plot+labs(title="(B) 2020-21")->decomp_age_male_2020_2021
decomp_age_func(2020, "Male")$plot+labs(title="(A) 2019-20")->decomp_age_male_2019_2020
decomp_age_func(2019, "Male")

combined_plot_male <- arrangeGrob(decomp_age_male_2019_2020, 
                                  decomp_age_male_2020_2021, 
                                  decomp_age_male_2021_2022,
                                  ncol=3, 
                                  top=textGrob("Contribution to Life Expectancy Changes (years) by Age Group: Male Population", 
                                               gp = gpar(fontface = "bold", fontsize = 32)))
plot(combined_plot_male)

ggsave("./figures/Suppl_fig3_decomp_age_male_GH.tiff", plot = combined_plot_male, dpi = 300, width = 24, height = 12, 
       compression = "lzw")


#### female ####
decomp_age_func(2022, "Female")$plot+labs(title="(C) 2021-22")->decomp_age_female_2021_2022
decomp_age_func(2021, "Female")$plot+labs(title="(B) 2020-21")->decomp_age_female_2020_2021
decomp_age_func(2020, "Female")$plot+labs(title="(A) 2019-20")->decomp_age_female_2019_2020
decomp_age_func(2019, "Female")

combined_plot_female <- arrangeGrob(decomp_age_female_2019_2020, 
                                    decomp_age_female_2020_2021, 
                                    decomp_age_female_2021_2022,
                                    ncol=3, 
                                    top=textGrob("Contribution to Life Expectancy Changes (years) by Age Group: Female Population", 
                                                 gp = gpar(fontface = "bold", fontsize = 32)))
plot(combined_plot_female)

ggsave("./figures/Suppl_fig3_decomp_age_female_GH.tiff", plot = combined_plot_male, dpi = 300, width = 24, height = 12, 
       compression = "lzw")


library(writexl)

decomp_age_func(2020, "Total")$data -> data_decomp_age_total_2020
decomp_age_func(2021, "Total")$data -> data_decomp_age_total_2021
decomp_age_func(2022, "Total")$data -> data_decomp_age_total_2022
data_decomp_age_total_2020 %>%
  cbind(data_decomp_age_total_2021$nCx) %>% 
  cbind(data_decomp_age_total_2022$nCx) -> data_decomp_age_total

colnames(data_decomp_age_total) <- c("Age", paste0("nCx_", 2020:2022))
write_xlsx(data_decomp_age_total, "./figures/age_decomp_total_GH.xlsx")

decomp_age_func(2020, "Male")$data -> data_decomp_age_male_2020
decomp_age_func(2021, "Male")$data -> data_decomp_age_male_2021
decomp_age_func(2022, "Male")$data -> data_decomp_age_male_2022

data_decomp_age_male_2020 %>%
  cbind(data_decomp_age_male_2021$nCx) %>% 
  cbind(data_decomp_age_male_2022$nCx) -> data_decomp_age_male

colnames(data_decomp_age_male) <- c("Age", paste0("nCx_", 2020:2022))
write_xlsx(data_decomp_age_male, "./figures/age_decomp_male_GH.xlsx")


decomp_age_func(2020, "Female")$data->data_decomp_age_female_2020
decomp_age_func(2021, "Female")$data->data_decomp_age_female_2021
decomp_age_func(2022, "Female")$data->data_decomp_age_female_2022
data_decomp_age_female_2020%>%cbind(data_decomp_age_female_2021$nCx)%>%cbind(data_decomp_age_female_2022$nCx)->data_decomp_age_female
colnames(data_decomp_age_female)<-c("Age", paste0("nCx_", 2020:2022))


write_xlsx(data_decomp_age_female, "./figures/age_decomp_female_GH.xlsx")
