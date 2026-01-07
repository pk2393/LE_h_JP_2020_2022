rm(list=ls())
library(dplyr);library(magrittr)

load(file="life_table_2000_2022.RData")


cause_list=c("c_all", "c02000", "c09000", 
             "c10000", "c20000", "c06000", "c11000", "c14000", "c22000", "c05000")


###
### summarize top 9 + remaining
###
cause_top10_func=function(df, year, sex, cause_list){
  
  sex=paste0("_", sex)
  
  df %<>%
    filter(col2==sex) %>%
    filter(col1 %in% cause_list) %>%
    select(-c(all_age, col2, unknown, ag01, ag02, ag03, ag04)) %>%
    mutate(ag0104=ag0004-ag00, .before=ag0004) %>%
    select(-ag0004)
  
  df %>% colnames() -> ag_name
  ag_name<-ag_name[-1]
  
  df %<>%t()
  nm <- df[1,]; colnames(df) <- nm; 
  df <- df[-1,]
  
  df %<>% data.frame() %<>% mutate(ag=ag_name, .before=c_all)
  rownames(df)<-NULL
  
  df[,c(2:11)] <- lapply(df[,c(2:11)], as.numeric)
  
  df$c_others <- df[,2]-  rowSums(df[,c(3:11)])
  
  
  cause_list <- append(cause_list, "c_others")
  for(i in 2:11){
    df %<>%
      mutate(!!paste0("Rx_", cause_list[i]):= !!sym(paste0(cause_list[i]))/c_all)
  }
  
  colnames(df) <- paste0("y", year, "_", colnames(df))
  
  return(df)
}



###
### decompose by cause
###

cause_decomp_func=function(year, sex, cause_list, lt_list){

  readr::read_csv(paste0("./data_cause/Cx_GH", year-2000, ".", sex, ".csv")) -> Cx_next
  readr::read_csv(paste0("./data_cause/Cx_GH", year-2001, ".", sex, ".csv")) -> Cx_prev
  
  readr::read_csv(file=paste0("./data_cause/cause_", year, "_raw.csv")) -> cause_next
  readr::read_csv(file=paste0("./data_cause/cause_", year-1, "_raw.csv")) -> cause_prev
  
  sex2 = tolower(sex)
  
  cause_top10_func(cause_next, year, sex2, cause_list) -> Cause10_next
  cause_top10_func(cause_prev, year-1, sex2, cause_list) -> Cause10_prev
  
  Cause10_next[, c(1, 13:22)] -> Cause10_next
  Cause10_prev[, c(1, 13:22)] -> Cause10_prev
  
  colnames(Cause10_next) <- c("ag", paste0("Rx", 1:10, "_next"))
  colnames(Cause10_prev) <- c("ag", paste0("Rx", 1:10, "_prev"))
  
  lt_list[[paste0(year-1, ".", sex)]] %>% dplyr::select(mx) %>% unlist() %>% as.vector()->mx_prev
  lt_list[[paste0(year, ".", sex)]] %>% dplyr::select(mx) %>% unlist() %>%as.vector()->mx_next
  
  
  Cx_next%<>%
    cbind(Cause10_next[,-1], Cause10_prev[,-1])%>%
    mutate(mx_prev = mx_prev, 
           mx_next = mx_next)%>%
    mutate(
      nCx = case_when(
        is.na(B)== T ~ A, 
        T ~ A + B
    ))%>%
    mutate(
      cx1 = case_when(
        mx_next == mx_prev ~ 0, 
        T ~ nCx*((Rx1_next*mx_next - Rx1_prev*mx_prev)/(mx_next-mx_prev))
      ), 
      cx2 = case_when(
        mx_next == mx_prev ~ 0, 
        T ~ nCx*((Rx2_next*mx_next - Rx2_prev*mx_prev)/(mx_next-mx_prev))
      ), 
      cx3 = case_when(
        mx_next == mx_prev ~ 0, 
        T ~ nCx*((Rx3_next*mx_next - Rx3_prev*mx_prev)/(mx_next-mx_prev))
      ), 
      cx4=case_when(
        mx_next == mx_prev ~ 0, 
        T ~ nCx*((Rx4_next*mx_next - Rx4_prev*mx_prev)/(mx_next-mx_prev))
      ), 
      cx5=case_when(
        mx_next == mx_prev ~ 0, 
        T ~ nCx*((Rx5_next*mx_next - Rx5_prev*mx_prev)/(mx_next-mx_prev))
      ), 
      cx6=case_when(
        mx_next == mx_prev ~ 0, 
        T ~ nCx*((Rx6_next*mx_next - Rx6_prev*mx_prev)/(mx_next-mx_prev))
      ), 
      cx7=case_when(
        mx_next == mx_prev ~ 0, 
        T ~ nCx*((Rx7_next*mx_next - Rx7_prev*mx_prev)/(mx_next-mx_prev))
      ), 
      cx8=case_when(
        mx_next == mx_prev ~ 0, 
        T ~ nCx*((Rx8_next*mx_next - Rx8_prev*mx_prev)/(mx_next-mx_prev))
      ), 
      cx9=case_when(
        mx_next == mx_prev ~ 0, 
        T ~ nCx*((Rx9_next*mx_next - Rx9_prev*mx_prev)/(mx_next-mx_prev))
      ), 
      cx10=case_when(
        mx_next == mx_prev ~ 0, 
        T ~ nCx*((Rx10_next*mx_next - Rx10_prev*mx_prev)/(mx_next-mx_prev))
      ) 
    )
  
  print(sum(Cx_next$nCx))
  print(Cx_next$Jp_next.Tx[1] - Cx_next$Jp_prev.Tx[1])
  print(sum((Cx_next$cx1+Cx_next$cx2+Cx_next$cx3+Cx_next$cx4+Cx_next$cx5+Cx_next$cx6+Cx_next$cx7+Cx_next$cx8+Cx_next$cx9+Cx_next$cx10), na.rm=T))
  
  
  summ_data = Cx_next %>%
    select(cx1,cx2,cx3,cx4,cx5,cx6,cx7,cx8,cx9,cx10)
  
  Cause_result_total=data.frame(
    Cause=c("Neoplastic", "Cardiovascular","Respiratory", "External", "Nervous", "Digestive", "Genitourinary or Renal", "COVID-19", "Mental", "Others"), 
    Cx=c(sum(summ_data$cx1), sum(summ_data$cx2), sum(summ_data$cx3), sum(summ_data$cx4), 
         sum(summ_data$cx5), sum(summ_data$cx6), sum(summ_data$cx7), sum(summ_data$cx8), 
         sum(summ_data$cx9), sum(summ_data$cx10))
    
  )
  
  
  
  colnames(Cx_next)<-c("Age", 
                    "lx_before", "Lx_before", "Tx_before", 
                    "lx_after", "Lx_after", "Tx_after", 
                    "Direct", "Indirect", "nCx", 
                    "Rx_neoplastic_after", "Rx_cardiovascular_after", "Rx_respiratory_after", "Rx_external_after", "Rx_nervous_after", 
                    "Rx_digestive_after", "Rx_genitourinary_renal_after", "Rx_covid19_after", "Rx_mental_after", "Rx_others_after",
                    "Rx_neoplastic_before", "Rx_cardiovascular_before", "Rx_respiratory_before", "Rx_external_before", "Rx_nervous_before", 
                    "Rx_digestive_before", "Rx_genitourinary_renal_before", "Rx_covid19_before", "Rx_mental_before", "Rx_others_before",
                    "mx_before", "mx_after", 
                    "nCx_neoplastic", "nCx_cardiovascular", "nCx_respiratory", "nCx_external", "nCx_nervous", 
                    "nCx_digestive", "nCx_genitourinary_renal", "nCx_covid19", "nCx_mental", "nCx_others"
  )
  
  res_list=list(Data = summ_data, 
                Cause_result_total = Cause_result_total, 
                Cause_result_detailed = Cx_next)
  
  return(res_list)
}

cause_decomp_func(2022, "Total", cause_list, lt_national_for_arriaga)->res
cause_decomp_func(2021, "Total", cause_list, lt_national_for_arriaga)->res2
cause_decomp_func(2020, "Total", cause_list, lt_national_for_arriaga)->res3


###
### output results
###
library(writexl)
cause_decomp_func(2022, "Total", cause_list, lt_national_for_arriaga)$Cause_result_detailed%>%write_xlsx("./figures/suppl_data_decomp_total_2022_GH.xlsx")
cause_decomp_func(2021, "Total", cause_list, lt_national_for_arriaga)$Cause_result_detailed%>%write_xlsx("./figures/suppl_data_decomp_total_2021_GH.xlsx")
cause_decomp_func(2020, "Total", cause_list, lt_national_for_arriaga)$Cause_result_detailed%>%write_xlsx("./figures/suppl_data_decomp_total_2020_GH.xlsx")

cause_decomp_func(2022, "Male", cause_list, lt_national_for_arriaga)$Cause_result_detailed%>%write_xlsx("./figures/suppl_data_decomp_male_2022_GH.xlsx")
cause_decomp_func(2021, "Male", cause_list, lt_national_for_arriaga)$Cause_result_detailed%>%write_xlsx("./figures/suppl_data_decomp_male_2021_GH.xlsx")
cause_decomp_func(2020, "Male", cause_list, lt_national_for_arriaga)$Cause_result_detailed%>%write_xlsx("./figures/suppl_data_decomp_male_2020_GH.xlsx")

cause_decomp_func(2022, "Female", cause_list, lt_national_for_arriaga)$Cause_result_detailed%>%write_xlsx("./figures/suppl_data_decomp_female_2022_GH.xlsx")
cause_decomp_func(2021, "Female", cause_list, lt_national_for_arriaga)$Cause_result_detailed%>%write_xlsx("./figures/suppl_data_decomp_female_2021_GH.xlsx")
cause_decomp_func(2020, "Female", cause_list, lt_national_for_arriaga)$Cause_result_detailed%>%write_xlsx("./figures/suppl_data_decomp_female_2020_GH.xlsx")


###
### visualise
###

library(reshape2);library(RColorBrewer);library(scales);library(gridExtra);library(grid);library(patchwork) 


plot_func_2=function(year, sex, cause_list, lt_list){

  cause_decomp_func(year, sex, cause_list, lt_list)$Cause_result_detailed->tmp_df
  
  tmp_df%<>%
    select(Age, 
           nCx_neoplastic, nCx_cardiovascular, nCx_respiratory, nCx_external, nCx_nervous,
           nCx_digestive, nCx_genitourinary_renal, nCx_covid19, nCx_mental, nCx_others)
  
  colnames(tmp_df)<-c("Age", 
                      "Neoplastic", "Cardiovascular","Respiratory", "External", "Nervous", 
                      "Digestive", "Genitourinary or Renal", "COVID-19", "Mental", "Others")

  long_data<-melt(tmp_df, id.vars="Age")
  

    ggplot(long_data, aes(fill = variable, y = value, x = factor(Age, levels = unique(Age)))) +
    geom_bar(stat = "identity", position = "stack") + 
    coord_flip() +  
    scale_y_continuous(limits = c(-0.1, 0.1), breaks = pretty_breaks(n = 10)) +
    scale_fill_brewer(palette = "Paired", name = "Cause of Death") +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +  
    labs(x = "Age", 
         y = "Change in Life Expectancy (years)") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.line = element_line(color = "black"), 
          axis.ticks = element_line(color = "black"),
          axis.title.x = element_text(size = 14), 
          axis.title.y = element_text(size = 14),  
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12),  
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 12), 
          plot.title=element_text(size=20, face = "bold")) ->plot
  return(plot)
}


plot_aggr_age_cause = function(year, sex, cause_list, lt_list){
  plot_func_2(year, sex, cause_list, lt_list)->plot_total_2022
  plot_func_2(year-1, sex, cause_list, lt_list)->plot_total_2021
  plot_func_2(year-2, sex, cause_list, lt_list)->plot_total_2020
  p1 <- plot_total_2020+labs(title="(A) 2019-20")+ theme(legend.position = "none")
  p2 <- plot_total_2021+labs(title="(B) 2020-21")+ theme(legend.position = "none")
  p3 <- plot_total_2022+labs(title="(C) 2021-22")+ theme(legend.position = "none")
  combined_plot <- p1 + p2 + p3
  final_plot <- combined_plot + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
  return(final_plot)
}

plot_aggr_age_cause(2022, "Total", cause_list, lt_national_for_arriaga)->plot_age_cause_total
plot_aggr_age_cause(2022, "Male", cause_list, lt_national_for_arriaga)->plot_age_cause_male
plot_aggr_age_cause(2022, "Female", cause_list, lt_national_for_arriaga)->plot_age_cause_female


ggsave("./figures/figxxx_decomp_age_cause_total_GH.tiff", plot = plot_age_cause_total, dpi = 300, width = 30, height = 10, 
       compression = "lzw")
ggsave("./figures/Suppl_figxxx_decomp_age_cause_male_GH.tiff", plot = plot_age_cause_male, dpi = 300, width = 30, height = 10, 
       compression = "lzw")
ggsave("./figures/Suppl_figxxx_decomp_age_cause_female_GH.tiff", plot = plot_age_cause_female, dpi = 300, width = 30, height = 10, 
       compression = "lzw")



