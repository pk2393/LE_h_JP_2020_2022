rm(list=ls())
library(tidyverse);library(magrittr)
palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

lt_jmd_b<-read.table("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/bltper_1x1.txt", skip=2, header=T)
lt_jmd_f<-read.table("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/fltper_1x1.txt", skip=2, header=T)
lt_jmd_m<-read.table("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/mltper_1x1.txt", skip=2, header=T)


###
### function to calc H and thresholds
### modified from codes  used in Aburto et al. DOI: https://dx.doi.org/10.4054/DemRes.2019.41.4
### https://github.com/jmaburto/The-treshold-age-of-the-lifetable-Entropy/tree/master
###

calc_H <- function(lt){ 
  
  Age <- lt$Age
  Age[Age == "110+"] <- 110
  Age <- as.integer(Age)
  
  lx <- lt$lx / 1e5
  dx <- lt$dx / 1e5
  ex <- lt$ex
  ax <- lt$ax
  
  # e-dagger
  ex_bar <- ex + ax * (ex - c(ex[-1], tail(ex, 1)))
  ex_bar[length(ex_bar)] <- tail(ex, 1)
  e_dag_x <- rev(cumsum(rev(ex_bar * dx))) / lx
  Hx <- e_dag_x / ex        # H(x)
  H0 <- Hx[1]               # Keyfitz H at birth
  e_dag0 <- e_dag_x[1]
  
  # cumulative hazard
  cum_haz <- -log(lx)
  g <- cum_haz + Hx - 1 - H0
  
  # g(a_H)=0
  f_g <- approxfun(Age, g, rule = 2)
  a_H <- uniroot(function(z) f_g(z), c(0, 110))$root

  
  # wx, Wx
  wx <- dx * ex # by definition: discretized
  wxWx <- wx * (g / e_dag0)            
  
  lt_h <- data.frame(
    Age, ex, lx, dx, 
    Hx, cum_haz = cum_haz,
    g, wx, wxWx
  )
  
  threshold <- data.frame(e0 = ex[1], H = H0, H_log = -log(H0),
                          a_H = a_H
                          )

  
  return(list(lt_h = lt_h, threshold = threshold))
}


###
### summarize results
###

res_equality_one_by_one = data.frame(
  year=2000:2022, 
  h_total=NA,   h_male=NA,   h_female=NA, 
  a_h_total=NA, a_h_male=NA, a_h_female=NA, 
  e0_total=NA, e0_male=NA, e0_female=NA
)

res_wxWx <- res_wxWx_female <- res_wxWx_male <-  matrix(NA, nrow=111, ncol=24) %>% data.frame()
colnames(res_wxWx) <- colnames(res_wxWx_female) <- colnames(res_wxWx_male)<-c("age", paste0("wxWx", 2000:2022))
res_wxWx$age<-res_wxWx_male$age<-res_wxWx_female$age<-0:110

for(y in 2000:2022){
  lt_jmd_b%>%
    filter(Year==y)->lt_tmp
  lt_jmd_f%>%
    filter(Year==y)->lt_tmp_f
  lt_jmd_m%>%
    filter(Year==y)->lt_tmp_m

  res_equality_one_by_one$e0_total[y-1999]<-lt_tmp$ex[1]
  res_equality_one_by_one$e0_male[y-1999]<-lt_tmp_m$ex[1]
  res_equality_one_by_one$e0_female[y-1999]<-lt_tmp_f$ex[1]
  
  calc_H(lt_tmp) -> H_res_tmp
  calc_H(lt_tmp_f) -> H_res_tmp_f
  calc_H(lt_tmp_m) -> H_res_tmp_m
  

  res_wxWx[[paste0("wxWx",y)]] <- H_res_tmp$lt_h$wxWx
  res_wxWx_female[[paste0("wxWx",y)]] <- H_res_tmp_f$lt_h$wxWx
  res_wxWx_male[[paste0("wxWx",y)]] <- H_res_tmp_m$lt_h$wxWx
  
  -log(H_res_tmp$lt_h$Hx[1]) -> res_equality_one_by_one$h_total[y-1999]
  -log(H_res_tmp_f$lt_h$Hx[1]) -> res_equality_one_by_one$h_female[y-1999]
  -log(H_res_tmp_m$lt_h$Hx[1]) -> res_equality_one_by_one$h_male[y-1999]

  H_res_tmp$threshold$a_H -> res_equality_one_by_one$a_h_total[y-1999]
  H_res_tmp_f$threshold$a_H -> res_equality_one_by_one$a_h_female[y-1999]
  H_res_tmp_m$threshold$a_H -> res_equality_one_by_one$a_h_male[y-1999]
}


###
### visualise h-e0 and h-year
###

res_equality_one_by_one%<>%
  mutate(lead_e0_total=lead(e0_total), 
         lead_e0_male=lead(e0_male), 
         lead_e0_female=lead(e0_female), 
         lead_h_total=lead(h_total), 
         lead_h_male=lead(h_male), 
         lead_h_female=lead(h_female))%>%
  mutate(arrow_e0_total=.1*e0_total+.9*lead_e0_total, 
         arrow_e0_male=.1*e0_male+.9*lead_e0_male, 
         arrow_e0_female=.1*e0_female+.9*lead_e0_female, 
         arrow_h_total=.1*h_total+.9*lead_h_total, 
         arrow_h_male=.1*h_male+.9*lead_h_male, 
         arrow_h_female=.1*h_female+.9*lead_h_female)%>%
  mutate(year_text_1=case_when(
    year==2011 ~ year, 
    year==2005 ~ year, 
    year==2000 ~ year, 
    T ~ NA
  ))%>%
  mutate(year_text_2=case_when(
    year==2022 ~ year, 
    year==2020 ~ year, 
    year==2019 ~ year, 
    T ~ NA
  ))

highlight_years=c(2000, 2005, 2011, 2019, 2020, 2022)


res_equality_one_by_one %>%
  ggplot(aes(x = e0_total, y = h_total)) +
  geom_point() +  
  geom_point(data = res_equality_one_by_one[res_equality_one_by_one$year %in% highlight_years, ],
             aes(x = e0_total, y = h_total),
             size = 2.5, 
             color=palette[2]) + 
  geom_segment(aes(xend = arrow_e0_total, yend = arrow_h_total), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = palette[12]) +
  theme_minimal() +
  geom_text(aes(label = year_text_1), nudge_x = -0.1, nudge_y = 0.002, size=5) +
  geom_text(aes(label = year_text_2), nudge_x = 0.1, nudge_y = 0.003, size=5) +
  labs(title = "(B) Life Expectancy at Birth and h, 2000-22", 
       x = "Life Expectancy at Birth (years)", 
       y = expression(h))+
  theme(plot.title = element_text(size = 24), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(
    limits=c(81, 85), breaks = scales::pretty_breaks(n = 10))->plot_e0_h_total



plot_e0_h_total

res_equality_one_by_one %>%
  ggplot(aes(x = e0_male, y = h_male)) +
  geom_point() +  
  geom_point(data = res_equality_one_by_one[res_equality_one_by_one$year %in% highlight_years, ],
             aes(x = e0_male, y = h_male),
             size = 2.5, 
             color=palette[2]) + 
  geom_segment(aes(xend = arrow_e0_male, yend = arrow_h_male), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = palette[12]) +
  theme_minimal() +
  geom_text(aes(label = year_text_1), nudge_x = -0.15, nudge_y = 0.002, size=5) +
  geom_text(aes(label = year_text_2), nudge_x = 0.15, nudge_y = 0.003, size=5) +
  labs(title = "(B) Life Expectancy at Birth and h, 2000-22", 
       x = "Life Expectancy at Birth (years)", 
       y = expression(h))+
  theme(plot.title = element_text(size = 24), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(
    limits=c(77.5, 82), breaks = scales::pretty_breaks(n = 10))->plot_e0_h_male

plot_e0_h_male

res_equality_one_by_one %>%
  ggplot(aes(x = e0_female, y = h_female)) +
  geom_point() +  
  geom_point(data = res_equality_one_by_one[res_equality_one_by_one$year %in% highlight_years, ],
             aes(x = e0_female, y = h_female),
             size = 2.5, 
             color=palette[2]) + 
  geom_segment(aes(xend = arrow_e0_female, yend = arrow_h_female), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = palette[12]) +
  theme_minimal() +
  geom_text(aes(label = year_text_1), nudge_x = -0.15, nudge_y = 0.002, size=5) +
  geom_text(aes(label = year_text_2), nudge_x = 0.075, nudge_y = 0.003, size=5) +
  labs(title = "(B) Life Expectancy at Birth and h, 2000-22", 
       x = "Life Expectancy at Birth (years)", 
       y = expression(h))+
  theme(plot.title = element_text(size = 24), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(
    limits=c(84, 88), breaks = scales::pretty_breaks(n = 10))->plot_e0_h_female

plot_e0_h_female


res_equality_one_by_one %>%
  ggplot(aes(x = year, y = h_total)) +
  geom_point() +
  labs(title = "(A) Trend of h, 2000-22", 
       x = "Year", 
       y = expression(h))+
  theme_minimal() +
  theme(plot.title = element_text(size = 24), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits=c(1.99, 2.12))+
  scale_x_continuous(breaks = res_equality_one_by_one$year,
                     labels = as.character(res_equality_one_by_one$year))->plot_year_h_total
plot_year_h_total

res_equality_one_by_one %>%
  ggplot(aes(x = year, y = h_male)) +
  geom_point() +
  labs(title = "(A) Trend of h, 2000-22", 
       x = "Year", 
       y = expression(h))+
  theme_minimal() +
  theme(plot.title = element_text(size = 24), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits=c(1.92, 2.07))+
  scale_x_continuous(breaks = res_equality_one_by_one$year,
                     labels = as.character(res_equality_one_by_one$year))->plot_year_h_male
plot_year_h_male

res_equality_one_by_one %>%
  ggplot(aes(x = year, y = h_female)) +
  geom_point() +
  labs(title = "(A) Trend of h, 2000-22", 
       x = "Year", 
       y = expression(h))+
  theme_minimal() +
  theme(plot.title = element_text(size = 24), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits=c(2.14, 2.265))+
  scale_x_continuous(breaks = res_equality_one_by_one$year,
                     labels = as.character(res_equality_one_by_one$year)) -> plot_year_h_female
plot_year_h_female



library(gridExtra)
combined_plot_total <- arrangeGrob(plot_year_h_total, 
                                   plot_e0_h_total,
                                   ncol=2)
combined_plot_male <- arrangeGrob(plot_year_h_male, 
                                  plot_e0_h_male,
                                  ncol=2)
combined_plot_female <- arrangeGrob(plot_year_h_female, 
                                    plot_e0_h_female,
                                    ncol=2)
plot(combined_plot_total)

ggsave("./figures/fig4_e0_equality_total_GH.tiff", plot = combined_plot_total, dpi = 300, width = 20, height = 10, 
       compression = "lzw")
ggsave("./figures/suppl_fig6_e0_equality_male_GH.tiff", plot = combined_plot_male, dpi = 300, width = 20, height = 10, 
       compression = "lzw")
ggsave("./figures/suppl_fig5_e0_equality_female_GH.tiff", plot = combined_plot_female, dpi = 300, width = 20, height = 10, 
       compression = "lzw")



###
### visualise wxWhx
###
res_wxWx$wxWx2022*-1->wxWhx_2022
res_wxWx$wxWx2021*-1->wxWhx_2021
res_wxWx$wxWx2020*-1->wxWhx_2020

res_wxWhx<-res_wxWx
res_wxWhx[,-1]<-res_wxWx[,-1]*-1
res_wxWhx_male<-res_wxWx_male
res_wxWhx_male[,-1]<-res_wxWx_male[,-1]*-1
res_wxWhx_female<-res_wxWx_female
res_wxWhx_female[,-1]<-res_wxWx_female[,-1]*-1

plot_wxWhx = function(df1, df2, sex){
  
  df1%>%
    ggplot(aes(x = age, y = wxWx2022, color="2022")) +
    geom_line(linewidth=1)+
    geom_line(aes(x = age, y = wxWx2021, color="2021"),linetype="dashed", linewidth=1)+ 
    scale_fill_manual(values = c("2022" = palette[11], "2021" = palette[9])) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    geom_vline(aes(xintercept=df2[[paste0("a_h_", sex)]][df2$year == 2022], color="2022"), linetype="dashed")+
    theme_minimal() +
    labs(title = paste0("(A) w(x)Wh(x), ", stringr::str_to_title(sex), " Population"), 
         x = "Age", 
         y = "w(x)Wh(x)")+
    theme(plot.title = element_text(size = 24), 
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(hjust = 1), 
          legend.title = element_blank(), 
          legend.text=element_text(size=16))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))->plot_wxWhx
  
  return(plot_wxWhx)
}
plot_wxWhx(res_wxWhx, res_equality_one_by_one, "total")
plot_wxWhx(res_wxWhx_male, res_equality_one_by_one, "male")
plot_wxWhx(res_wxWhx_female, res_equality_one_by_one, "female")

#############
log_haz_diff=function(year, lft){

    lft%>%filter(Year==year)->lft_after
  cumhaz_after=-log(lft_after$lx/100000)
  haz_after=c(diff(cumhaz_after), lft_after$mx[111])
  
  lft%>%filter(Year==year-1)->lft_before
  cumhaz_before=-log(lft_before$lx/100000)
  haz_before=c(diff(cumhaz_before), lft_before$mx[111])
  
  log_haz_diff=log(haz_after)-log(haz_before)
  
  return(log_haz_diff)
}


###
### visualise yearly change of loghaz
###


plot_mort_improve=function(lt, df2, sex){
  
  log_haz_diff_df=
    data.frame(age=0:110, 
               lhd2020=log_haz_diff(2020, lt), 
               lhd2021=log_haz_diff(2021, lt), 
               lhd2022=log_haz_diff(2022, lt))
  
  log_haz_diff_df %>%
    ggplot(aes(x = age, y = -lhd2022, color="2021-22")) +
    geom_line(linewidth=1)+
    geom_line(aes(x = age, y = -lhd2021, color="2020-21"), linewidth=1)+ 
    geom_line(aes(x = age, y = -lhd2020, color="2019-20"), linewidth=1)+ 
    scale_fill_manual(values = c("2021-22" = palette[11], "2020-21" = palette[9], "2019-20" = palette[8])) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    geom_vline(aes(xintercept=df2[[paste0("a_h_", sex)]][df2$year == 2022], color="2021-22"), linetype="dashed")+
    theme_minimal() +
    labs(title = paste0("(B) Year-on-year Mortality improvement, ", stringr::str_to_title(sex), " Population"), 
         x = "Age", 
         y = "mortality improvement")+
    theme(plot.title = element_text(size = 24), 
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(hjust = 1), 
          legend.title = element_blank(), 
          legend.text=element_text(size=16))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))->plot_mort_improve
  return(plot_mort_improve)
}

plot_mort_improve(lt_jmd_m, res_equality_one_by_one, "male")


library(gridExtra)

combined_plot_total <- arrangeGrob(plot_wxWhx(res_wxWhx, res_equality_one_by_one, "total"), 
                                   plot_mort_improve(lt_jmd_b, res_equality_one_by_one, "total"),
                                   ncol=2)
combined_plot_male <- arrangeGrob(plot_wxWhx(res_wxWhx_male, res_equality_one_by_one, "male"), 
                                   plot_mort_improve(lt_jmd_m, res_equality_one_by_one, "male"),
                                   ncol=2)
combined_plot_female <- arrangeGrob(plot_wxWhx(res_wxWhx_female, res_equality_one_by_one, "female"), 
                                   plot_mort_improve(lt_jmd_f, res_equality_one_by_one, "female"),
                                   ncol=2)

ggsave("./figures/suppl_fig7_wxWhx_mort_GH.tiff", plot = combined_plot_total, dpi = 300, width = 20, height = 10, 
       compression = "lzw")
ggsave("./figures/suppl_fig8_wxWhx_mort_female_GH.tiff", plot = combined_plot_female, dpi = 300, width = 20, height = 10, 
       compression = "lzw")
ggsave("./figures/suppl_fig9_wxWhx_mort_male_GH.tiff", plot = combined_plot_male, dpi = 300, width = 20, height = 10, 
       compression = "lzw")



###
### output
###

res_equality_one_by_one%>%

  select(year, 
         h_total, h_female, h_male, 
         e0_total, e0_female, e0_male, 
         a_h_total, a_h_female, a_h_male)%>%
  writexl::write_xlsx("./figures/res_h_ah_GH.xlsx")


res_wxWhx%>%
  select(age, wxWx2019, wxWx2020, wxWx2021, wxWx2022)%>%
  rename(wxWhx2019=wxWx2019, 
         wxWhx2020=wxWx2020, 
         wxWhx2021=wxWx2021, 
         wxWhx2022=wxWx2022
         )%>%
  writexl::write_xlsx("./figures/res_wxWhx_GH.xlsx")


log_haz_diff_df=
  data.frame(age=0:110, 
             lhd2020=log_haz_diff(2020, lt_jmd_b), 
             lhd2021=log_haz_diff(2021, lt_jmd_b), 
             lhd2022=log_haz_diff(2022, lt_jmd_b))%>%
  writexl::write_xlsx("./figures/log_haz_diff_GH.xlsx")


