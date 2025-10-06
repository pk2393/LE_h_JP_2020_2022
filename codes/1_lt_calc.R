rm(list=ls())
library(dplyr);library(magrittr)

###
readr::read_csv("pref_info.csv") -> pref.info
pref_list<-pref.info$pref_list

lt_national <- read.table("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/bltper_5x1.txt", skip = 2, header = 2)
death_national <- read.table("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/Deaths_5x1.txt", skip = 2, header = 2)
etr_national <- read.table("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/Exposures_5x1.txt", skip = 2, header = 2)

lt_calc = function(lt_orig, death_vec, etr_vec){
  
  lt_tmp  <- matrix(nrow = dim(lt_orig)[1]-2, 
                    ncol = dim(lt_orig)[2]) %>% data.frame()
  
  colnames(lt_tmp) <- colnames(lt_orig)
  lt_tmp$Year <-lt_orig$Year[1:22]
  lt_tmp$Age  <- c(lt_orig$Age[1:21], "100+")
  if(length(death_vec)>22){
    lt_tmp$mx <- c(death_vec[1:21], sum(death_vec[22:24])) / c(etr_vec[1:21], sum(etr_vec[22:24]))
  }else{
    lt_tmp$mx <- c(death_vec) / c(etr_vec)
  }
  
  lt_tmp$ax[1:21] <- lt_orig$ax[1:21]
  lt_tmp$ax[22] <- 1/lt_tmp$mx[22]
  
  lt_tmp %<>%
    mutate(w = c(1, 4, c(rep(5, 20))), 
           ax = case_when(
             ax > w & Age != "100+" ~ w/2, 
             T ~ ax
           ),
           qx = case_when(
             Age =="100+" ~ 1, 
             T ~ (mx * w) / (1 + (w - ax) * mx) 
           ))
  
  lt_tmp$lx[1] <- 1e5
  for(n in 2:nrow(lt_tmp)){
    lt_tmp$lx[n] <- lt_tmp$lx[n-1] * (1-lt_tmp$qx[n-1])
  }
  
  lt_tmp %<>% 
    mutate(dx = c(-diff(lx), NA))%>%
    mutate(dx = case_when(is.na(dx) == T ~ lx, T ~ dx))
  
  for(n in 1:(nrow(lt_tmp))){
    if(n <nrow(lt_tmp)){
      lt_tmp$Lx[n] = lt_tmp$lx[n+1] * lt_tmp$w[n] + (lt_tmp$ax[n] * lt_tmp$dx[n])
    }else{
      lt_tmp$Lx[n]=(lt_tmp$ax[n] * lt_tmp$dx[n])
    }
  }
  
  Tx_vec = sum(lt_tmp$Lx) - c(0, cumsum(lt_tmp$Lx)[-22])
  
  lt_tmp %>%
    mutate(Tx = Tx_vec, 
           ex = Tx / lx) -> lt_y_shortened
  
  
  return(lt_y_shortened)
}

lt_national_for_arriaga = list()
sex.list=c("b", "m", "f")
sex.list2=c("Total", "Male", "Female")


death_national <- read.table("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/Deaths_5x1.txt", skip = 2, header = 2)
etr_national <- read.table("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/Exposures_5x1.txt", skip = 2, header = 2)

for(s_idx in 1:3){
  sex=sex.list[s_idx]
  lt_national <- read.table(paste0("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/", sex, "ltper_5x1.txt"), skip=2, header=T)
  
  
  for(y in 2000:2022){
    
    lt_y <- lt_national%>%
      filter(Year == y)
    
    death_vec <- (death_national%>%filter(Year == y))[[sex.list2[s_idx]]]
    etr_vec <- (etr_national%>%filter(Year == y))[[sex.list2[s_idx]]]
    
    lt_y_shortened <- lt_calc(lt_y, death_vec, etr_vec)
    
    lt_y_shortened$death <- c(death_vec[1:21], sum(death_vec[22:24]))
    lt_y_shortened$etr <- c(etr_vec[1:21], sum(etr_vec[22:24]))
    
    ### store
    lt_national_for_arriaga[[paste0(y, ".", sex.list2[s_idx])]] <- lt_y_shortened
  }
  
}

save(lt_national_for_arriaga, file="life_table_2000_2022.RData")


#####
##### Prefecture Data
#####

#### pref names and codes ####
pref.info$pref_code<-as.numeric(pref.info$pref_code)

LE.summ=data.frame(
  pref=c(rep(NA, 47)), 
  total2019=c(rep(NA, 47)),
  total2020=c(rep(NA, 47)),
  total2021=c(rep(NA, 47)),
  total2022=c(rep(NA, 47))
)


for(n in 1:47){
  # n=5
  print(n)
  pref=pref.info$pref_list[pref.info$pref_code==n]
  LE.summ[["pref"]][n]<-pref
  
  
  if(n<10){
    number=as.character(paste0("0", n))
  }else{
    number=as.character(n)
  }
  
  etr_jmd=read.table(paste0("https://www.ipss.go.jp/p-toukei/JMD/", number,"/STATS/Exposures_5x1.txt"), skip=2, header=T)
  death_jmd=read.table(paste0("https://www.ipss.go.jp/p-toukei/JMD/", number,"/STATS/Deaths_5x1.txt"), skip=2, header = T)
  
  year.list=c(2019:2022)
  
  lft.japan.arriaga.list=c()
  
  for(s_idx in 1:3){
    
    sex=sex.list[s_idx]
    lt_jmd = read.table(paste0("https://www.ipss.go.jp/p-toukei/JMD/", number, "/STATS/", sex, "ltper_5x1.txt"), skip=2, header=T)
    
    for(year in year.list){
      
      lt_y <- lt_jmd%>%
        filter(Year == year)
      
      for(ax_idx in 1:nrow(lt_y)){
        
        if (!is.finite(lt_y$ax[ax_idx]) || is.na(lt_y$ax[ax_idx]) || lt_y$ax[ax_idx] <= 0) {
          
          ax_pref <- (read.table(paste0("https://www.ipss.go.jp/p-toukei/JMD/", number, "/STATS/", "b", "ltper_5x1.txt"),
                                 skip=2, header=TRUE) %>% filter(Year == year))$ax[ax_idx]
          
          if (!is.finite(ax_pref) || is.na(ax_pref) || ax_pref <= 0) {
            ax_pref <- (read.table(paste0("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/", sex, "ltper_5x1.txt"),
                                   skip=2, header=TRUE) %>% filter(Year == year))$ax[ax_idx]
          }
          lt_y$ax[ax_idx] <- ax_pref
        }
        
      }    
      
      death_vec <- (death_jmd%>%filter(Year == year))[[sex.list2[s_idx]]]
      etr_vec <- (etr_jmd%>%filter(Year == year))[[sex.list2[s_idx]]]
      
      lt_y_shortened <- lt_calc(lt_y, death_vec, etr_vec)
      
      
      lt_y_shortened -> lft.japan.arriaga.list[[paste0(year, ".", sex.list2[s_idx])]]
      
      LE.summ[[paste0(tolower(sex.list2[s_idx]), year)]][n]<-lt_y_shortened$ex[1]
      
    }
    
  }
  
  
  
  save(lft.japan.arriaga.list, file=paste0("./data_prefectures/", "life_table_2000_2022_GH.", pref, ".RData"))
  
}


readr::write_csv(LE.summ, "./data_prefectures/LE_summ_pref_GH.csv")
