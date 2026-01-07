rm(list=ls())
library(dplyr);library(magrittr);library(ggplot2);library(patchwork)

# load(file="life_table_2000_2022.RData")

###
### read life table from JMD
###
lt_jmd_b<-read.table("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/bltper_1x1.txt", skip=2, header=T)
lt_jmd_f<-read.table("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/fltper_1x1.txt", skip=2, header=T)
lt_jmd_m<-read.table("https://www.ipss.go.jp/p-toukei/JMD/00/STATS/mltper_1x1.txt", skip=2, header=T)

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
### function to calc H and thresholds
### modified from codes  used in Aburto et al. DOI: https://dx.doi.org/10.4054/DemRes.2019.41.4
### https://github.com/jmaburto/The-treshold-age-of-the-lifetable-Entropy/tree/master
###

calc_H <- function(lt){ 
  Age <- lt$Age
  Age[Age == "110+"] <- 110
  Age <- as.integer(Age)
  
  lx <- lt$lx / 1e5
  mx <- lt$mx
  dx <- lt$dx / 1e5
  ex <- lt$ex
  ax <- lt$ax
  ax[nrow(lt)] <- 1/mx[nrow(lt)]
  
  # e-dagger
  ex_bar <- ex + ax * (ex - c(ex[-1], tail(ex, 1)))
  ex_bar[length(ex_bar)] <- tail(ex, 1)
  e_dag_x <- rev(cumsum(rev(ex_bar * dx))) / lx
  Hx <- e_dag_x / ex        # H(x)
  H0 <- Hx[1]               # Keyfitz H at birth
  e_dag0 <- e_dag_x[1]
  
  lt_h <- data.frame(
    Age, ex, lx, dx
  )
  
  threshold <- data.frame(e0 = ex[1],
                          e_dag  = e_dag0, 
                          H = H0, 
                          H_log = -log(H0)
  )
  
  
  return(list(lt_h = lt_h, threshold = threshold))
}

###
### calculate life table from mx
###

lt_from_mx <- function(Age, mx, ax) {
  n <- length(Age)
  
  ax[length(ax)] <- 1/mx[length(mx)]
  
  
  qx <- mx / (1 + (1 - ax) * mx)
  qx[n] <- 1
  
  lx <- numeric(n)
  lx[1] <- 1e5
  for (i in 1:(n-1)) {
    lx[i+1] <- lx[i] * (1 - qx[i])
  }
  
  dx <- lx * qx
  
  Lx <- numeric(n)
  for (i in 1:(n-1)) {
    Lx[i] <- lx[i+1] + ax[i] * dx[i]
  }
  Lx[n] <- lx[n] / mx[n] 
  
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  
  data.frame(Age = Age, mx = mx, qx = qx, lx = lx, dx = dx, Lx = Lx, Tx = Tx, ex = ex, ax = ax)
}

### ===========================
### Output as figures
### ===========================

### 111 -> 22 age groups
make_agegrp22 <- function(age_int){
  
  if(age_int == 0) return("0")
  
  if(age_int >= 1 && age_int <= 4) return("1-4")
  
  if(age_int >= 100) return("100+")
  
  lo <- (age_int %/% 5) * 5
  
  hi <- lo + 4
  
  paste0(lo, "-", hi)
}

age_levels22 <- c("0", "1-4",
                  paste0(seq(5, 95, 5), "-", seq(9, 99, 5)),
                  "100+")

abridge_agecause_22x10 <- function(mat_age111_cause10, Age_vec){
  # remove"+" from "110+"
  age_int <- suppressWarnings(as.integer(gsub("\\+", "", as.character(Age_vec))))
  age_grp <- vapply(age_int, make_agegrp22, FUN.VALUE = character(1))
  
  df <- as.data.frame(mat_age111_cause10)
  df$age_grp <- age_grp
  
  out <- aggregate(. ~ age_grp, data = df, FUN = sum)
  
  rownames(out) <- out$age_grp
  out$age_grp <- NULL
  
  out2 <- out[age_levels22, , drop = FALSE]
  out2[is.na(out2)] <- 0
  
  as.matrix(out2)
}



### ================================
### decomposition and visualization
### ================================

h_decomp_res    <- matrix(0, nrow = 3, ncol = 10) %>% as.data.frame()
edag_decomp_res <- matrix(0, nrow = 3, ncol = 10) %>% as.data.frame()

log_e0_decomp_res    <- matrix(0, nrow = 3, ncol = 10) %>% as.data.frame()
log_edag_decomp_res <- matrix(0, nrow = 3, ncol = 10) %>% as.data.frame()

sex_list <- c("Total", "Female", "Male")

h_decomp_list = list()
e_dag_decomp_list = list()
log_e0_decomp_list   <- list()
log_edag_decomp_list <- list()

for(s in 1:3){
  # s = 1
  
  sex <- sex_list[s]
  
  plotdat_h_age   <- list()
  plotdat_h_cause <- list()
  plotdat_h_heat  <- list()
  
  plotdat_edag_age   <- list()
  plotdat_edag_cause <- list()
  plotdat_edag_heat  <- list()
  
  
  for(year in 2020:2022){
    # year = 2020
    cause_next <- readr::read_csv(paste0("./data_cause/cause_", year, "_raw.csv"))
    cause_prev <- readr::read_csv(paste0("./data_cause/cause_", year-1, "_raw.csv"))
    
    sex2 <- tolower(sex)
    
    Cause10_next <- cause_top10_func(cause_next, year,   sex2, cause_list)
    Cause10_prev <- cause_top10_func(cause_prev, year-1, sex2, cause_list)
    
    Cause10_next <- Cause10_next %>% dplyr::select(1, dplyr::contains("Rx_"))
    Cause10_prev <- Cause10_prev %>% dplyr::select(1, dplyr::contains("Rx_"))
    
    colnames(Cause10_next) <- c("ag", paste0("Rx", 1:10, "_next"))
    colnames(Cause10_prev) <- c("ag", paste0("Rx", 1:10, "_prev"))
    
    
    if(s == 1){
      lt_jmd_b%>%
        filter(Year==year)->lt_tmp_next  
      
      lt_jmd_b%>%
        filter(Year==year-1)->lt_tmp_prev
    }
    if(s == 2){
      lt_jmd_f%>%
        filter(Year==year)->lt_tmp_next  
      
      lt_jmd_f%>%
        filter(Year==year-1)->lt_tmp_prev
    }
    if(s == 3){
      lt_jmd_m%>%
        filter(Year==year)->lt_tmp_next  
      
      lt_jmd_m%>%
        filter(Year==year-1)->lt_tmp_prev
    }
    
    
    Age <- lt_tmp_prev$Age
    
    
    mx_age_cause_mat_next <- matrix(nrow = 111, ncol = 10)
    mx_age_cause_mat_prev <- matrix(nrow = 111, ncol = 10)
    
    for(idx_cause in 1:10){
      frac_cause_next_abridged <- Cause10_next[, idx_cause+1]
      frac_cause_next <- c()
      
      frac_cause_prev_abridged <- Cause10_prev[, idx_cause+1]
      frac_cause_prev <- c()
      
      for(ag in 1:22){
        if(ag ==1){
          frac_cause_next <- append(frac_cause_next, frac_cause_next_abridged[ag])
          frac_cause_prev <- append(frac_cause_prev, frac_cause_prev_abridged[ag])
        }
        if(ag==2){
          frac_cause_next <- append(frac_cause_next, rep(frac_cause_next_abridged[ag], 4))
          frac_cause_prev <- append(frac_cause_prev, rep(frac_cause_prev_abridged[ag], 4))
        }
        if(ag >= 3 & ag <=21){
          frac_cause_next <- append(frac_cause_next, rep(frac_cause_next_abridged[ag], 5))
          frac_cause_prev <- append(frac_cause_prev, rep(frac_cause_prev_abridged[ag], 5))
        }
        if(ag ==22){
          frac_cause_next <- append(frac_cause_next, rep(frac_cause_next_abridged[ag], 11))
          frac_cause_prev <- append(frac_cause_prev, rep(frac_cause_prev_abridged[ag], 11))
          
        }
      }
      
      mx_age_cause_mat_next[, idx_cause] <- lt_tmp_next$mx * frac_cause_next
      mx_age_cause_mat_prev[, idx_cause] <- lt_tmp_prev$mx * frac_cause_prev
      
    }
    
    mx_prev <- mx_age_cause_mat_prev
    dmx <- mx_age_cause_mat_next - mx_age_cause_mat_prev
    
    ax_base <- lt_tmp_prev$ax
    n_age <- nrow(lt_tmp_prev)
    
    # move alpha for perturbation
    eval_from_alpha <- function(alpha) {
      
      mx_ac  <- mx_prev + alpha * dmx
      mx_all <- rowSums(mx_ac)
      
      lt <- lt_from_mx(Age, mx_all, ax_base)
      th <- calc_H(lt)$threshold
      
      log_e0   <- log(th$e0)
      log_edag <- log(th$e_dag)
      
      h <- log_e0 - log_edag
      
      c(h = h, e_dag = th$e_dag,
        log_e0   = log_e0,
        log_edag = log_edag)
    }
    
    ###
    ### decomp by age and cause : Horiuchi et al, 2008 ;  https://doi.org/10.1353/dem.0.0033
    ###
    
    decomp_Horiuchi_cause <- function(N = 100) {
      
      n_cause <- ncol(mx_prev)
      step <- 0.5 / N
      
      contrib_edag_mat <- matrix(0, nrow = n_age, ncol = n_cause)
      contrib_loge0_mat    <- matrix(0, nrow = n_age, ncol = n_cause)
      contrib_logedag_mat  <- matrix(0, nrow = n_age, ncol = n_cause)
      contrib_h_mat <- contrib_loge0_mat - contrib_logedag_mat
      
      alpha_tmp <- matrix(0, nrow = n_age, ncol = n_cause)
      
      for (k in 1:N) {
        if(k %% 5 == 0 ) print(k)
        
        mid <- (k - 0.5) / N
        alpha_tmp[] <- mid
        
        for (a in 1:n_age) {
          for (c in 1:n_cause) {
            
            # forward
            alpha_tmp[a, c] <- mid + step
            vp <- eval_from_alpha(alpha_tmp)
            
            # backward
            alpha_tmp[a, c] <- mid - step
            vm <- eval_from_alpha(alpha_tmp)
            
            # reset position
            alpha_tmp[a, c] <- mid
            
            contrib_edag_mat[a, c] <- contrib_edag_mat[a, c] + (vp["e_dag"] - vm["e_dag"])
            contrib_loge0_mat[a, c]   <- contrib_loge0_mat[a, c]   + (vp["log_e0"]   - vm["log_e0"])
            contrib_logedag_mat[a, c] <- contrib_logedag_mat[a, c] + (vp["log_edag"] - vm["log_edag"])
            
            contrib_h_mat[a, c] <- contrib_loge0_mat[a, c]  - contrib_logedag_mat[a, c]
            
          }
        }
      }
      
      rownames(contrib_h_mat) <- Age
      rownames(contrib_edag_mat) <- Age
      rownames(contrib_loge0_mat)   <- Age
      rownames(contrib_logedag_mat) <- Age
      
      colnames(contrib_h_mat) <- paste0("cause", 1:n_cause)
      colnames(contrib_edag_mat) <- paste0("cause", 1:n_cause)
      colnames(contrib_loge0_mat)   <- paste0("cause", 1:n_cause)
      colnames(contrib_logedag_mat) <- paste0("cause", 1:n_cause)
      
      list(
        h_agecause        = contrib_h_mat,
        e_dag_agecause    = contrib_edag_mat,
        log_e0_agecause   = contrib_loge0_mat,
        log_edag_agecause = contrib_logedag_mat,
        
        h_by_age        = rowSums(contrib_h_mat),
        e_dag_by_age    = rowSums(contrib_edag_mat),
        log_e0_by_age   = rowSums(contrib_loge0_mat),
        log_edag_by_age = rowSums(contrib_logedag_mat),
        
        h_by_cause        = colSums(contrib_h_mat),
        e_dag_by_cause    = colSums(contrib_edag_mat),
        log_e0_by_cause   = colSums(contrib_loge0_mat),
        log_edag_by_cause = colSums(contrib_logedag_mat)
      )
    }
    
    contrib_cause_10 <- decomp_Horiuchi_cause(N = 10)
    contrib_cause_20 <- decomp_Horiuchi_cause(N = 20)
    
    rtol <- 1e-5
    
    relF <- function(A_low, A_high) {
      num <- sqrt(sum((A_high - A_low)^2, na.rm=TRUE))
      den <- sqrt(sum((A_high)^2, na.rm=TRUE))
      num / max(den, .Machine$double.eps)
    }
    
    # check convergence
    ok_10_20 <- (
        relF(contrib_cause_10$e_dag_agecause,    contrib_cause_20$e_dag_agecause)    < rtol &&
        relF(contrib_cause_10$log_e0_agecause,   contrib_cause_20$log_e0_agecause)   < rtol &&
        relF(contrib_cause_10$log_edag_agecause, contrib_cause_20$log_edag_agecause) < rtol
    )
    
    if (ok_10_20) {
      contrib_cause <- contrib_cause_20
    } else{
      contrib_cause_40 <- decomp_Horiuchi_cause(N = 40)
      
      ok_20_40 <- (
          relF(contrib_cause_20$e_dag_agecause,    contrib_cause_40$e_dag_agecause)    < rtol &&
          relF(contrib_cause_20$log_e0_agecause,   contrib_cause_40$log_e0_agecause)   < rtol &&
          relF(contrib_cause_20$log_edag_agecause, contrib_cause_40$log_edag_agecause) < rtol
      )
      
      if (ok_20_40) {
        contrib_cause <- contrib_cause_40
      } else {
        contrib_cause <- contrib_cause_40
        message("Warning: not converged by rtol; result from N=40.")
      }
      
    }
    
    ### endpoint balance check
    n_cause <- ncol(mx_prev)
    
    alpha0 <- matrix(0, nrow = n_age, ncol = n_cause)  # prev year
    alpha1 <- matrix(1, nrow = n_age, ncol = n_cause)  # next year
    
    v0 <- eval_from_alpha(alpha0)
    v1 <- eval_from_alpha(alpha1)
    
    delta_edag    <- unname(v1["e_dag"]   - v0["e_dag"])
    delta_loge0   <- unname(v1["log_e0"]  - v0["log_e0"])
    delta_logedag <- unname(v1["log_edag"]- v0["log_edag"])
    delta_h <- delta_loge0 - delta_logedag
    
    sum_edag    <- sum(contrib_cause$e_dag_agecause)
    sum_loge0   <- sum(contrib_cause$log_e0_agecause)
    sum_logedag <- sum(contrib_cause$log_edag_agecause)
    sum_h   <- sum_loge0 - sum_logedag
    
    res_h       <- sum_h       - delta_h
    res_edag    <- sum_edag    - delta_edag
    res_loge0   <- sum_loge0   - delta_loge0
    res_logedag <- sum_logedag - delta_logedag
    
    rel_err <- function(res, delta) res / max(abs(delta), .Machine$double.eps)
    
    re_edag    <- rel_err(res_edag,    delta_edag)
    re_loge0   <- rel_err(res_loge0,   delta_loge0)
    re_logedag <- rel_err(res_logedag, delta_logedag)
    
    tol_balance <- 1e-5
    if (any(abs(c(re_edag, re_loge0, re_logedag)) > tol_balance)) {
      warning(sprintf(
        "Balance check failed (sex=%s, year=%d): relerr(edag)=%.2e, relerr(loge0)=%.2e, relerr(logedag)=%.2e",
        sex, year, re_edag, re_loge0, re_logedag
      ))
    }
    
    h_decomp_res[year-2019, 1:10]        <- contrib_cause$h_by_cause
    edag_decomp_res[year-2019, 1:10]     <- contrib_cause$e_dag_by_cause
    log_e0_decomp_res[year-2019, 1:10]   <- contrib_cause$log_e0_by_cause
    log_edag_decomp_res[year-2019, 1:10] <- contrib_cause$log_edag_by_cause
    
    
    #### label: year-to-year difference
    interval_lab <- paste0(year-1, "-", year)
    
    ### #ag: 111 -> 22
    h_mat_22x10    <- abridge_agecause_22x10(contrib_cause$h_agecause,     Age)
    edag_mat_22x10 <- abridge_agecause_22x10(contrib_cause$e_dag_agecause, Age)
    
    ### give names to cause of death
    cause_names <- c("Neoplastic", "Cardiovascular","Respiratory", "External",
                     "Nervous", "Digestive", "Genitourinary or Renal",
                     "COVID-19", "Mental", "Others")
    
    colnames(h_mat_22x10)    <- cause_names
    colnames(edag_mat_22x10) <- cause_names
    rownames(h_mat_22x10)    <- age_levels22
    rownames(edag_mat_22x10) <- age_levels22
    
    ### aggregate by age and cause
    h_by_age_22    <- rowSums(h_mat_22x10)
    h_by_cause_10  <- colSums(h_mat_22x10)
    
    edag_by_age_22   <- rowSums(edag_mat_22x10)
    edag_by_cause_10 <- colSums(edag_mat_22x10)
    
    ### data for plot
    plotdat_h_age[[interval_lab]] <- data.frame(
      interval = interval_lab,
      age = factor(names(h_by_age_22), levels = age_levels22),
      contrib = as.numeric(h_by_age_22)
    )
    
    plotdat_h_cause[[interval_lab]] <- data.frame(
      interval = interval_lab,
      cause = factor(names(h_by_cause_10), levels = cause_names),
      contrib = as.numeric(h_by_cause_10)
    )
    
    plotdat_edag_age[[interval_lab]] <- data.frame(
      interval = interval_lab,
      age = factor(names(edag_by_age_22), levels = age_levels22),
      contrib = as.numeric(edag_by_age_22)
    )
    
    plotdat_edag_cause[[interval_lab]] <- data.frame(
      interval = interval_lab,
      cause = factor(names(edag_by_cause_10), levels = cause_names),
      contrib = as.numeric(edag_by_cause_10)
    )
    
    # data for heatmaps
    tmp_h_heat <- as.data.frame(h_mat_22x10)
    tmp_h_heat$age <- rownames(tmp_h_heat)
    tmp_h_heat <- tidyr::pivot_longer(tmp_h_heat, -age,
                                      names_to = "cause", values_to = "contrib")
    tmp_h_heat$interval <- interval_lab
    tmp_h_heat$age   <- factor(tmp_h_heat$age, levels = age_levels22)
    tmp_h_heat$cause <- factor(tmp_h_heat$cause, levels = cause_names)
    plotdat_h_heat[[interval_lab]] <- tmp_h_heat
    
    tmp_edag_heat <- as.data.frame(edag_mat_22x10)
    tmp_edag_heat$age <- rownames(tmp_edag_heat)
    tmp_edag_heat <- tidyr::pivot_longer(tmp_edag_heat, -age,
                                         names_to = "cause", values_to = "contrib")
    tmp_edag_heat$interval <- interval_lab
    tmp_edag_heat$age   <- factor(tmp_edag_heat$age, levels = age_levels22)
    tmp_edag_heat$cause <- factor(tmp_edag_heat$cause, levels = cause_names)
    plotdat_edag_heat[[interval_lab]] <- tmp_edag_heat
    
  }
  
  
  interval_levels <- c("2019-2020", "2020-2021", "2021-2022")
  
  h_decomp_list[[sex]]        <- cbind(interval = interval_levels, h_decomp_res)
  e_dag_decomp_list[[sex]]    <- cbind(interval = interval_levels, edag_decomp_res)
  log_e0_decomp_list[[sex]]   <- cbind(interval = interval_levels, log_e0_decomp_res)
  log_edag_decomp_list[[sex]] <- cbind(interval = interval_levels, log_edag_decomp_res)
  
  cause_names <- c("Neoplastic", "Cardiovascular","Respiratory", "External", "Nervous", "Digestive", "Genitourinary or Renal", "COVID-19", "Mental", "Others")
  
  ### aggregate plot data by sex
  df_h_age   <- dplyr::bind_rows(plotdat_h_age)
  df_h_cause <- dplyr::bind_rows(plotdat_h_cause)
  df_h_heat  <- dplyr::bind_rows(plotdat_h_heat)
  
  df_edag_age   <- dplyr::bind_rows(plotdat_edag_age)
  df_edag_cause <- dplyr::bind_rows(plotdat_edag_cause)
  df_edag_heat  <- dplyr::bind_rows(plotdat_edag_heat)
  
  # ---- add sex label ----
  df_h_age       <- df_h_age       %>% dplyr::mutate(sex = sex) %>% dplyr::relocate(sex)
  df_h_cause     <- df_h_cause     %>% dplyr::mutate(sex = sex) %>% dplyr::relocate(sex)
  df_h_heat      <- df_h_heat      %>% dplyr::mutate(sex = sex) %>% dplyr::relocate(sex)
  
  df_edag_age    <- df_edag_age    %>% dplyr::mutate(sex = sex) %>% dplyr::relocate(sex)
  df_edag_cause  <- df_edag_cause  %>% dplyr::mutate(sex = sex) %>% dplyr::relocate(sex)
  df_edag_heat   <- df_edag_heat   %>% dplyr::mutate(sex = sex) %>% dplyr::relocate(sex)
  
  # ---- save decomposition plot data ----
  sex_tag <- gsub("[[:space:]]+", "_", sex)
  
  readr::write_csv(df_h_age,      paste0("./figures/data_h_by_age22_",      sex_tag, ".csv"))
  readr::write_csv(df_h_cause,    paste0("./figures/data_h_by_cause10_",    sex_tag, ".csv"))
  readr::write_csv(df_h_heat,     paste0("./figures/data_h_heat_22x10_",    sex_tag, ".csv"))
  
  readr::write_csv(df_edag_age,   paste0("./figures/data_edag_by_age22_",   sex_tag, ".csv"))
  readr::write_csv(df_edag_cause, paste0("./figures/data_edag_by_cause10_", sex_tag, ".csv"))
  readr::write_csv(df_edag_heat,  paste0("./figures/data_edag_heat_22x10_", sex_tag, ".csv"))
  
  
  ### ordering of plots
  interval_levels <- c("2019-2020", "2020-2021", "2021-2022")
  df_h_age$interval      <- factor(df_h_age$interval, levels = interval_levels)
  df_h_cause$interval    <- factor(df_h_cause$interval, levels = interval_levels)
  df_h_heat$interval     <- factor(df_h_heat$interval, levels = interval_levels)
  df_edag_age$interval   <- factor(df_edag_age$interval, levels = interval_levels)
  df_edag_cause$interval <- factor(df_edag_cause$interval, levels = interval_levels)
  df_edag_heat$interval  <- factor(df_edag_heat$interval, levels = interval_levels)
  
  sex_tag <- gsub("[[:space:]]+", "_", sex)
  
  ###
  ### h: panel A: age / panel B: cause
  ###
  
  p_h_age <- ggplot2::ggplot(df_h_age, ggplot2::aes(x = age, y = contrib)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~interval, ncol = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Panel A: age decomposition (Δh)",
                  x = NULL, y = "Contribution")
  
  p_h_cause <- ggplot2::ggplot(df_h_cause, ggplot2::aes(x = cause, y = contrib)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~interval, ncol = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Panel B: cause decomposition (Δh)",
                  x = NULL, y = "Contribution")
  
  fig_h_AB <- p_h_age + p_h_cause + patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(title = paste0("Horiuchi decomposition of Δh, ", sex))
  
  ggplot2::ggsave(paste0("./figures/fig_h_AB_", sex_tag, ".tiff"),
                  plot = fig_h_AB,
                  dpi = 300, width = 14, height = 10, compression = "lzw")
  
  ###
  ### h: heatmap
  ###
  
  p_h_heat <- ggplot2::ggplot(df_h_heat, ggplot2::aes(x = cause, y = age, fill = contrib)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~interval, nrow = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_fill_gradient2(midpoint = 0) +
    ggplot2::labs(title = paste0("Age × cause contributions to Δh, ", sex),
                  x = NULL, y = NULL, fill = "Contribution")
  
  ggplot2::ggsave(paste0("./figures/fig_h_heat_", sex_tag, ".tiff"),
                  plot = p_h_heat,
                  dpi = 300, width = 15, height = 6, compression = "lzw")
  
  ###
  ### edag: panel A: age / panel B: cause
  ###
  
  p_edag_age <- ggplot2::ggplot(df_edag_age, ggplot2::aes(x = age, y = contrib)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~interval, ncol = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Panel A: age decomposition (Δe-dagger)",
                  x = NULL, y = "Contribution (years)")
  
  p_edag_cause <- ggplot2::ggplot(df_edag_cause, ggplot2::aes(x = cause, y = contrib)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~interval, ncol = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Panel B: cause decomposition (Δe-dagger)",
                  x = NULL, y = "Contribution (years)")
  
  fig_edag_AB <- p_edag_age + p_edag_cause + patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(title = paste0("Horiuchi decomposition of Δe-dagger, ", sex))
  
  ggplot2::ggsave(paste0("./figures/fig_edag_AB_", sex_tag, ".tiff"),
                  plot = fig_edag_AB,
                  dpi = 300, width = 14, height = 10, compression = "lzw")
  
  ###
  ### edag: heatmap
  ###
  p_edag_heat <- ggplot2::ggplot(df_edag_heat, ggplot2::aes(x = cause, y = age, fill = contrib)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~interval, nrow = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_fill_gradient2(midpoint = 0) +
    ggplot2::labs(title = paste0("Age × cause contributions to Δe-dagger, ", sex),
                  x = NULL, y = NULL, fill = "Contribution")
  
  ggplot2::ggsave(paste0("./figures/fig_edag_heat_", sex_tag, ".tiff"),
                  plot = p_edag_heat,
                  dpi = 300, width = 15, height = 6, compression = "lzw")
  

}


df_h_all        <- dplyr::bind_rows(h_decomp_list,        .id = "sex")
df_edag_all     <- dplyr::bind_rows(e_dag_decomp_list,    .id = "sex")
df_log_e0_all   <- dplyr::bind_rows(log_e0_decomp_list,   .id = "sex")
df_log_edag_all <- dplyr::bind_rows(log_edag_decomp_list, .id = "sex")

readr::write_csv(df_h_all,        "./figures/decomp_h_by_cause_sex_2020-2022.csv")
readr::write_csv(df_edag_all,     "./figures/decomp_edag_by_cause_sex_2020-2022.csv")
readr::write_csv(df_log_e0_all,   "./figures/decomp_loge0_by_cause_sex_2020-2022.csv")
readr::write_csv(df_log_edag_all, "./figures/decomp_logedag_by_cause_sex_2020-2022.csv")

