rm(list=ls())
library(tidyverse);library(magrittr);library(lubridate);library(broom)


df.le<-readr::read_csv(file="./data_prefectures/LE_summ_pref_GH.csv")
pref.info<-readr::read_csv("pref_info.csv")
#https://www.stat.go.jp/data/nihon/02.htm

df.le<-df.le%>%
  cbind(pref.info[,2])%>%
  mutate(gap21=total2021-total2020, 
         gap22=total2022-total2021)

for(year in 2021:2022){
  int2021=lubridate::interval(start=paste0(year, "/1/1"), end=paste0(year, "/12/31"))
  
  readr::read_csv(file="./data_covid/newly_confirmed_cases_daily.csv")%>%
    mutate(Date=as.Date(Date))%>%
    dplyr::filter(Date %within% int2021)->df.covid.cases
  
  colSums(df.covid.cases[, -c(1:2)])%>%as.data.frame()%>%
    mutate(pref=rownames(.))->covid.cases2021
  colnames(covid.cases2021)=paste0("cases", year)
  
  df.le%<>%
    cbind(covid.cases2021)
  
  readr::read_csv(file="./data_covid/severe_cases_daily.csv")%>%
    mutate(Date=as.Date(Date))%>%
    dplyr::filter(Date %within% int2021)->df.severe.cases
  
  df.severe.cases[is.na(df.severe.cases)]<-0
  df.severe.cases[,-1] <- data.frame(lapply(df.severe.cases[,-1], function(x) as.numeric(as.character(x))))
  
  colSums(df.severe.cases[, -c(1:2)])%>%as.data.frame()%>%
    mutate(pref=rownames(.))->severe.cases2021
  colnames(severe.cases2021)=paste0("severe.personday", year)
  
  df.le%<>%
    cbind(severe.cases2021)
  
  
  
  readr::read_csv(file="./data_covid/number_of_deaths_daily.csv")%>%
    mutate(Date=as.Date(Date))%>%
    dplyr::filter(Date %within% int2021)->df.death
  
  colSums(df.death[, -c(1:2)])%>%as.data.frame()%>%
    mutate(pref=rownames(.))->death2021
  colnames(death2021)=paste0("death", year)
  
  df.le%<>%
    cbind(death2021)
  
}


df.le[, 9:14]<-df.le[, 9:14]/df.le[,6]

readr::write_csv(df.le, file="le_covid_compare_GH.csv")


palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
##################################
### covid cases per 100,000
##################################
dot_size = 1.5
p1 <- ggplot(df.le, aes(x=log(cases2021), y=gap21)) +
  geom_point(aes(color='2020-21'), size=dot_size, shape=17) +
  theme_minimal()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.ticks=element_line(color="black"),
        axis.text.x = element_text(hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))+
  labs(color='Year', x="Cases per 100,000 (log-scale)", y="Life expectancy change (years)")

p1

p2 <- p1 + geom_point(aes(x=log(cases2022), y=gap22, color='2021-22'), size=dot_size) +
  theme_minimal()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.ticks=element_line(color="black"),
        axis.text.x = element_text(hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))

p2+geom_hline(yintercept=0, linetype="dashed")+
  scale_color_manual(values = c('2020-21' = "black", '2021-22' = palette[2]))+
  labs(title="(A)")->plot_cases


lm.cases2021=lm(gap21 ~ log(cases2021), data=df.le)
lm.cases2021 %>% broom::tidy(conf.int = T, conf.level = .95) %>%
  mutate(  across(c(estimate, std.error, conf.low, conf.high),
           ~ sprintf(paste0("%.", 3, "f"), .x)))
# # A tibble: 2 × 7
# term           estimate std.error statistic p.value conf.low conf.high
# <chr>          <chr>    <chr>         <dbl>   <dbl> <chr>    <chr>    
#   1 (Intercept)    -0.661   0.220         -3.00 0.00438 -1.105   -0.218   
# 2 log(cases2021) -0.104   0.044         -2.38 0.0215  -0.192   -0.016   

lm.cases2022=lm(gap22 ~ log(cases2022), data=df.le)
lm.cases2022 %>% broom::tidy(conf.int = T, conf.level = .95)%>%
  mutate(  across(c(estimate, std.error, conf.low, conf.high),
                  ~ sprintf(paste0("%.", 3, "f"), .x)))
# # A tibble: 2 × 7
# term           estimate std.error statistic p.value conf.low conf.high
# <chr>          <chr>    <chr>         <dbl>   <dbl> <chr>    <chr>    
#   1 (Intercept)    -0.649   0.314        -2.07   0.0447 -1.281   -0.016   
# 2 log(cases2022) -0.103   0.200        -0.513  0.611  -0.506   0.301   

##################################
### severe person-day per 100,000
##################################

p3 <- ggplot(df.le, aes(x=log(severe.personday2021), y=gap21)) +
  geom_point(aes(color='2020-21'), size=dot_size, shape=17) +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.ticks=element_line(color="black"),
        axis.text.x = element_text(hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))+
  labs(color='Year', x="Person-day under Intensive Care per 100,000 (log-scale)", y="Life expectancy change (years)")

p4 <- p3 + geom_point(aes(x=log(severe.personday2022), y=gap22, color='2021-22'), size=dot_size) +
  theme_minimal()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.ticks=element_line(color="black"),
        axis.text.x = element_text(hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))

p4+geom_hline(yintercept=0, linetype="dashed")+
  scale_color_manual(values = c('2020-21' = "black", '2021-22' = palette[2]))+
  labs(title="(B)")->plot_icu

lm.severe2021=lm(gap21 ~ log(severe.personday2021), data=df.le)
lm.severe2021 %>% broom::tidy(conf.int = T, conf.level = .95)%>%
  mutate(  across(c(estimate, std.error, conf.low, conf.high),
                  ~ sprintf(paste0("%.", 3, "f"), .x)))
# # A tibble: 2 × 7
# term                      estimate std.error statistic p.value conf.low conf.high
# <chr>                     <chr>    <chr>         <dbl>   <dbl> <chr>    <chr>    
#   1 (Intercept)               -0.710   0.219         -3.25 0.00219 -1.150   -0.270   
# 2 log(severe.personday2021) -0.082   0.031         -2.63 0.0118  -0.146   -0.019  

lm.severe2022=lm(gap22 ~ log(severe.personday2022), data=df.le)
lm.severe2022 %>% broom::tidy(conf.int = T, conf.level = .95)%>%
  mutate(  across(c(estimate, std.error, conf.low, conf.high),
                  ~ sprintf(paste0("%.", 3, "f"), .x)))
# # A tibble: 2 × 7
# term                      estimate std.error statistic  p.value conf.low conf.high
# <chr>                     <chr>    <chr>         <dbl>    <dbl> <chr>    <chr>    
#   1 (Intercept)               -0.874   0.243         -3.60 0.000781 -1.363   -0.386   
# 2 log(severe.personday2022) -0.052   0.033         -1.60 0.116    -0.118   0.014  

##################################
### death per 100,000
##################################

p5 <- ggplot(df.le, aes(x=log(death2021), y=gap21)) +
  geom_point(aes(color='2020-21'), size=dot_size, shape=17) +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.ticks=element_line(color="black"),
        axis.text.x = element_text(hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))+
  labs(color='Year', x="COVID-19 Deaths per 100,000 (log-scale)", y="Life expectancy change (years)")

p6 <- p5 + geom_point(aes(x=log(death2022), y=gap22, color='2021-22'), size=dot_size) +
  theme_minimal()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.ticks=element_line(color="black"),
        axis.text.x = element_text(hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))

p6+geom_hline(yintercept=0, linetype="dashed")+
  scale_color_manual(values = c('2020-21' = "black", '2021-22' = palette[2]))+
  labs(title="(C)")->plot_death

lm.death2021=lm(gap21 ~ log(death2021), data=df.le)
lm.death2021%>% broom::tidy(conf.int = T, conf.level = .95)%>%
  mutate(  across(c(estimate, std.error, conf.low, conf.high),
                  ~ sprintf(paste0("%.", 3, "f"), .x)))
# # A tibble: 2 × 7
# term           estimate std.error statistic p.value conf.low conf.high
# <chr>          <chr>    <chr>         <dbl>   <dbl> <chr>    <chr>    
#   1 (Intercept)    -0.795   0.345         -2.30  0.0260 -1.491   -0.100   
# 2 log(death2021) -0.067   0.035         -1.90  0.0636 -0.139   0.004  

lm.death2022=lm(gap22 ~ log(death2022), data=df.le)
lm.death2022%>% broom::tidy(conf.int = T, conf.level = .95)%>%
  mutate(  across(c(estimate, std.error, conf.low, conf.high),
                  ~ sprintf(paste0("%.", 3, "f"), .x)))
# # A tibble: 2 × 7
# term           estimate std.error statistic p.value conf.low conf.high
# <chr>          <chr>    <chr>         <dbl>   <dbl> <chr>    <chr>    
#   1 (Intercept)    -1.353   0.681         -1.99  0.0530 -2.724   0.019    
# 2 log(death2022) -0.107   0.084         -1.27  0.210  -0.276   0.062  



library(grid);library(gridExtra)
combined_plot <- arrangeGrob(plot_cases, plot_icu, plot_death, ncol=3, 
                             top=textGrob("Life expectancy changes (years) and reported COVID-19 burden by prefecture", 
                                          gp = gpar(fontface = "bold", fontsize = 20)))

ggsave("./figures/fig2_pref_e0_gap_covid_GH.tiff", plot = combined_plot, dpi = 300, width = 15, height = 5, 
       compression = "lzw")

