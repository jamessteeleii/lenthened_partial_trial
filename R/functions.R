# Functions for targets ----

# Functions for pre-registration sample size estimation ----
get_prepare_MA_data <- function(file) {
  
  ### Note code copied from original paper 
  
  ##### Read csv as data frame into environment - Note: change source address
  Data <- read.csv(here::here("data","Polito et al. RT Extracted Data.csv"), na.strings=c(""," ","NA"))
  
  # Calculate pre-post SDs from SEs
  Data$RT_pre_sd <- ifelse(is.na(Data$RT_pre_se), Data$RT_pre_sd, Data$RT_pre_se * sqrt(Data$RT_n))
  Data$CON_pre_sd <- ifelse(is.na(Data$CON_pre_se), Data$CON_pre_sd, Data$CON_pre_se * sqrt(Data$CON_n))
  Data$RT_post_sd <- ifelse(is.na(Data$RT_post_se), Data$RT_post_sd, Data$RT_post_se * sqrt(Data$RT_n))
  Data$CON_post_sd <- ifelse(is.na(Data$CON_post_se), Data$CON_post_sd, Data$CON_post_se * sqrt(Data$CON_n))
  
  # Convert p to t (Change scores)
  Data$RT_delta_t_value <- replmiss(Data$RT_delta_t_value, with(Data, qt(RT_delta_p_value/2, df=RT_n-1, lower.tail=FALSE)))
  Data$CON_delta_t_value <- replmiss(Data$CON_delta_t_value, with(Data, qt(CON_delta_p_value/2, df=CON_n-1, lower.tail=FALSE)))
  
  # Convert t to SE (Change scores)
  Data$RT_delta_se <- replmiss(Data$RT_delta_se, with(Data, ifelse(is.na(RT_delta_m), 
                                                                   (RT_post_m - RT_pre_m)/RT_delta_t_value, RT_delta_m/RT_delta_t_value)))
  Data$CON_delta_se <- replmiss(Data$CON_delta_se, with(Data, ifelse(is.na(CON_delta_m), 
                                                                     (CON_post_m - CON_pre_m)/CON_delta_t_value, CON_delta_m/CON_delta_t_value)))
  # Make positive
  Data$RT_delta_se <- ifelse(Data$RT_delta_se < 0, Data$RT_delta_se * -1, Data$RT_delta_se)
  Data$CON_delta_se <- ifelse(Data$CON_delta_se < 0, Data$CON_delta_se * -1, Data$CON_delta_se)
  
  # Convert CI to SE (Change scores)
  Data$RT_delta_se <- replmiss(Data$RT_delta_se, with(Data, (RT_delta_CI_upper - RT_delta_CI_lower)/3.92))
  Data$CON_delta_se <- replmiss(Data$CON_delta_se, with(Data, (CON_delta_CI_upper - CON_delta_CI_lower)/3.92))
  
  # Convert SE to SD (Change scores)
  Data$RT_delta_sd <- replmiss(Data$RT_delta_sd, with(Data, RT_delta_se * sqrt(RT_n)))
  Data$CON_delta_sd <- replmiss(Data$CON_delta_sd, with(Data, CON_delta_se * sqrt(CON_n)))
  
  # Calculate pre-post correlation coefficient for those with pre, post, and delta SDs
  Data$RT_ri <- (Data$RT_pre_sd^2 + Data$RT_post_sd^2 - Data$RT_delta_sd^2)/(2 * Data$RT_pre_sd * Data$RT_post_sd)
  Data$CON_ri <- (Data$CON_pre_sd^2 + Data$CON_post_sd^2 - Data$CON_delta_sd^2)/(2 * Data$CON_pre_sd * Data$CON_post_sd)
  
  # Remove values outside the range of -1 to +1 as they are likely due to misreporting or miscalculations in original studies
  Data$RT_ri <- ifelse(between(Data$RT_ri,-1,1) == FALSE, NA, Data$RT_ri)
  Data$CON_ri <- ifelse(between(Data$CON_ri,-1,1) == FALSE, NA, Data$CON_ri)
  
  # Then we'll convert using Fishers r to z, calculate a meta-analytic point estimate, and impute that across the studies with missing correlations
  Data <- escalc(measure = "ZCOR", ri = RT_ri, ni = RT_n, data = Data)
  
  Meta_RT_ri <- rma.mv(yi, V=vi, data=Data,
                       slab=paste(label),
                       random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                       control=list(optimizer="optim", optmethod="Nelder-Mead"))
  
  RobuEstMeta_RT_ri <- robust(Meta_RT_ri, Data$study)
  
  z2r_RT <- psych::fisherz2r(RobuEstMeta_RT_ri$b[1])
  
  Data$RT_ri <- ifelse(is.na(Data$RT_ri), z2r_RT, Data$RT_ri)
  
  Data <- escalc(measure = "ZCOR", ri = CON_ri, ni = CON_n, data = Data)
  
  ### Note, data is coded with study and arm as having explicit nesting so all random effects are (~ 1 | study, ~ 1 | arm)
  Meta_CON_ri <- rma.mv(yi, V=vi, data=Data,
                        slab=paste(label),
                        random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                        control=list(optimizer="optim", optmethod="Nelder-Mead"))
  
  RobuEstMeta_CON_ri <- robust(Meta_CON_ri, Data$study)
  
  z2r_CON <- psych::fisherz2r(RobuEstMeta_CON_ri$b[1])
  
  Data$CON_ri <- ifelse(is.na(Data$CON_ri), z2r_CON, Data$CON_ri)
  
  # Estimate change score difference SD where only pre-post data available
  Data$RT_delta_sd <- replmiss(Data$RT_delta_sd, with(Data, sqrt(RT_pre_sd^2 + RT_post_sd^2 - (2*RT_ri*RT_pre_sd*RT_post_sd))))
  Data$CON_delta_sd <- replmiss(Data$CON_delta_sd, with(Data, sqrt(CON_pre_sd^2 + CON_post_sd^2 - (2*CON_ri*CON_pre_sd*CON_post_sd))))
  
  # Arm based effect size calculations
  Data_hypertrophy <- Data %>% 
    filter(outcome == "hypertrophy")
  
  Data_SMD_hypertrophy_RT <- Data_hypertrophy |>
    select(study, arm, es, RT_n, RT_pre_m, RT_post_m, RT_pre_sd, RT_post_sd, RT_ri, weeks)
  
  Data_SMD_hypertrophy_RT <-  escalc(measure = "SMCR",
                                     m1i = RT_post_m, m2i = RT_pre_m,
                                     sd1i = RT_post_sd, 
                                     ri = RT_ri, ni = RT_n,
                                     data = Data_SMD_hypertrophy_RT
  )
  
  Data_SMD_hypertrophy_RT <- Data_SMD_hypertrophy_RT %>%
    filter(!is.na(vi) &
             !is.infinite(vi)) %>%
    mutate(se = sqrt(vi),
           wi = 1/sqrt(vi),
           size = 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))) |>
    select(study, arm, es, yi, vi, weeks)   |>
    mutate(arm = as.factor(unclass(factor(unlist(arm)))),
           es = as.factor(unclass(factor(unlist(es)))))
  
  
  Data_SMD_hypertrophy_CON <- Data_hypertrophy |>
    select(study, arm, es, CON_n, CON_pre_m, CON_post_m, CON_pre_sd, CON_post_sd, CON_ri, weeks)
  
  Data_SMD_hypertrophy_CON <-  escalc(measure = "SMCR",
                                      m1i = CON_post_m, m2i = CON_pre_m,
                                      sd1i = CON_post_sd, 
                                      ri = CON_ri, ni = CON_n,
                                      data = Data_SMD_hypertrophy_CON
  )
  
  Data_SMD_hypertrophy_CON <- Data_SMD_hypertrophy_CON %>%
    filter(!is.na(vi) &
             !is.infinite(vi)) %>%
    mutate(se = sqrt(vi),
           wi = 1/sqrt(vi),
           size = 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))) |>
    select(study, arm, es, yi, vi, weeks) %>%
    mutate(weeks = 0)    |>
    distinct(yi, .keep_all = TRUE) |>
    mutate(arm = as.factor(unclass(factor(unlist(arm)))+length(unique(Data_SMD_hypertrophy_RT$arm))),
           es = as.factor(unclass(factor(unlist(es)))+length(unique(Data_SMD_hypertrophy_RT$es))))
  
  
  data <- rbind(Data_SMD_hypertrophy_RT, Data_SMD_hypertrophy_CON)
}

log_growth_MA_model <- function(data) {
  
  
  MultiLevelModel_SMD_hypertrophy_weeks_arm <- rma.mv(yi, V=vi, data=data,
                                                      random = list(~ log1p(weeks) | study, ~ 1 | arm, ~ 1 | es), mods = ~ log1p(weeks),
                                                      method="REML", test="t", struct = "GEN"
  )
  
  RobuEstMultiLevelModel_SMD_hypertrophy_weeks_arm <- robust(MultiLevelModel_SMD_hypertrophy_weeks_arm, data$study)
  
  
}

plot_log_growth_MA_model <- function(model, data) {
  data_log <- cbind(data, predict(model)) %>%
    mutate(wi = 1/sqrt(vi),
           size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi)))
  
  preds_arm <- as.data.frame(predict(model, newmods = c(log1p(0:80)))) |>
    mutate(weeks = 0:80)
  
  preds_arm_12_24_36 <- as.data.frame(predict(model, newmods = log1p(c(12,24,36)))) |>
    mutate(weeks = c(12,24,36))
  
  
  plot_log_growth_MA_model <- ggplot(data_log, aes(x=weeks)) +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_point(aes(y=yi,size = size), alpha = 0.1) +
    geom_ribbon(data=preds_arm,aes(ymax=ci.ub, ymin=ci.lb), alpha = 0.2) +
    geom_line(data=preds_arm, aes(y=pred)) +
    annotate("text", 
             x = 12, y = -1, size = 3,
             label = glue::glue("12wks: {round(preds_arm_12_24_36$pred[1],2)}")) +
    annotate("text", 
             x = 12, y = -1.2, size = 3,
             label = glue::glue("[95%CI: {round(preds_arm_12_24_36$ci.lb[1],2)}, {round(preds_arm_12_24_36$ci.ub[1],2)}]")) +
    annotate("text", 
             x = 24, y = -1, size = 3,
             label = glue::glue("24wks: {round(preds_arm_12_24_36$pred[2],2)}")) +
    annotate("text", 
             x = 24, y = -1.2, size = 3,
             label = glue::glue("[95%CI: {round(preds_arm_12_24_36$ci.lb[2],2)}, {round(preds_arm_12_24_36$ci.ub[2],2)}]")) +
    annotate("text", 
             x = 36, y = -1, size = 3,
             label = glue::glue("36wks: {round(preds_arm_12_24_36$pred[3],2)}")) +
    annotate("text", 
             x = 36, y = -1.2, size = 3,
             label = glue::glue("[95%CI: {round(preds_arm_12_24_36$ci.lb[3],2)}, {round(preds_arm_12_24_36$ci.ub[3],2)}]")) +
    labs(y = "Standardised Mean Change", x = "Weeks") +
    theme_classic() +
    labs(title = "Hypertrophy Outcomes",
         subtitle = "Arm-based model - random intercepts (effect, arm, study) and slopes (study)",
         caption = "Data from Steele et al. (2023), DOI: 10.1080/02640414.2023.2286748") +
    guides(size = "none", fill = "none")
  
  plot_log_growth_MA_model
}

sample_size_simulation <- function(file1, file2) {

  # Using prior data from all Discover Strength studies baseline lean mass (standardised) to get estimates of between study (i.e., site) variance
  
  prior_data <- read_csv(file1) |>
    mutate(
      lean_mass_z = (lean_mass - mean(lean_mass, na.rm=TRUE))/sd(lean_mass, na.rm=TRUE)
    ) 
  
  model_site_variance <- lmer(lean_mass_z ~ 1 + (1 | study),
                              data = prior_data)
  
  site_variance <- VarCorr(model_site_variance)
  
  # Using prior study data to get estimates of measurement error for muscle areas (standardised SEMs)
  
  prior_ma_reliability <- read_csv(file2) |>
    # standardise variables
    mutate(
      arm_ma_z = (arm_ma - mean(arm_ma, na.rm=TRUE))/sd(arm_ma, na.rm=TRUE),
      thigh_ma_z = (thigh_ma - mean(thigh_ma, na.rm=TRUE))/sd(thigh_ma, na.rm=TRUE)
    )
  
  arm_ma_reliability <- SimplyAgree::reli_aov(measure = "arm_ma_z",
                                              item = "measurement_no",
                                              id = "participant",
                                              data = prior_ma_reliability)$SEM[,1]
  
  thigh_ma_reliability <- SimplyAgree::reli_aov(measure = "thigh_ma_z",
                                                item = "measurement_no",
                                                id = "participant",
                                                data = prior_ma_reliability)$SEM[,1]
  # Simulate for power
  
  sim <- function(participant_n = as.double(), site_n = as.double(), time_n = as.double(), measurements_n = as.double(),
                  b0 = as.double(), b_time = as.double(), b_cond = as.double(), b_cond_time = as.double(),         # fixed effects 
                  u_site = as.double(), u_participant_arm = as.double(), u_participant_thigh = as.double(), # random intercepts
                  arm_ma_error = as.double(), thigh_ma_error = as.double(),       # measurement error
                  ... # helps the function work with pmap() below
  ) {
    
    # set up data structure
    dat <- add_random(site = site_n) |>
      add_random(participant = participant_n, .nested_in = "site") |>
      add_between("participant", cond = c("full_ROM", "lengthened_partial")) |>
      add_recode("cond", "cond_dummy", full_ROM = 0, lengthened_partial = 1) |>
      add_within("participant", time = seq(0,1, by = time_n), 
                 measurement = seq(1:measurements_n)) |>
      add_ranef("site", u_site = u_site) |>
      add_ranef("participant", u_participant_arm = u_participant_arm) |>
      add_ranef("participant", u_participant_thigh = u_participant_thigh) |>
      add_ranef(arm_ma_error = arm_ma_error) |>
      add_ranef(thigh_ma_error = thigh_ma_error) |>
      mutate(arm_ma = (b0 + u_site + u_participant_arm) + (b_time * time) + (b_cond * cond_dummy) + (b_cond_time * cond_dummy * time) + arm_ma_error,
             thigh_ma = (b0 + u_site + u_participant_thigh) + (b_time * time) + (b_cond * cond_dummy) + (b_cond_time * cond_dummy * time) + thigh_ma_error) 
    
    # 
    # # run mixed effect model and return relevant values
    # model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (1 | participant),
    #                   data = dat, REML = TRUE)
    # 
    # model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (1 | participant),
    #                     data = dat, REML = TRUE)
    # 
    
    if (time_n == 1 && measurements_n == 1) {
      # run mixed effect model and return relevant values
      model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (1 | participant),
                        data = dat, REML = TRUE)
      
      model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (1 | participant),
                          data = dat, REML = TRUE)
    } else {
      # run mixed effect model and return relevant values
      model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (time | participant),
                        data = dat, REML = TRUE)
      
      model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (time | participant),
                          data = dat, REML = TRUE)
    }
    
    hypothesis_tests_arm <- hypotheses(model_arm, equivalence = c(NA,0.1), vcov = "satterthwaite", df = insight::get_df(model_arm))[4,] |>
      mutate(muscle_site = "arm")
    
    
    hypothesis_tests_thigh <- hypotheses(model_thigh, equivalence = c(NA,0.1), vcov = "satterthwaite", df = insight::get_df(model_thigh))[4,] |>
      mutate(muscle_site = "thigh")
    
    bind_rows(hypothesis_tests_arm, hypothesis_tests_thigh)
  }
  
  
  plan(cluster, workers = 10)
  
  sim <- crossing(
    rep = 1:1000, # number of replicates
    participant_n = seq(5,100, by = 5), # range of participant N
    site_n = 15, # fixed site N
    time_n = c(0.5,1), # to examine inclusion of midpoint test or not
    measurements_n = c(1:3), # range of measurements N
  ) |>
    mutate(
      b0 = 0, b_time = 0.05, b_cond = 0, b_cond_time = 0.025,         # fixed effects 
      u_site = sqrt(site_variance$study)[1],  # random intercept site
      u_participant_arm = (1 - sqrt(site_variance$study)[1] - arm_ma_reliability),   # random intercept participant arm
      u_participant_thigh = (1 - sqrt(site_variance$study)[1] - thigh_ma_reliability),  # random intercept participant thigh
      arm_ma_error = arm_ma_reliability,           # measurement error arm
      thigh_ma_error = thigh_ma_reliability         # measurement error thigh  
    ) %>%
    mutate(analysis = future_pmap(., sim)) %>% # not sure why base pipe doesn't work here
    unnest(analysis)
  
  plan(sequential)
  
  sim
  
}

plot_sample_size_estimates <- function(sim) {
  plot_sim <- sim |>
    filter(term == "time:cond_dummy") |>
    mutate(time_n = case_when(time_n == 0.5 ~ "Include Midpoint",
                              time_n == 1 ~ "No Midpoint",)) |>
    mutate(time_n = factor(time_n, levels = c("No Midpoint", "Include Midpoint"))) |>
    pivot_longer(c("p.value", "p.value.nonsup"), 
                 names_to = "test", 
                 values_to = "p.value") |> 
    mutate(test = case_when(test == "p.value" ~ "Difference",
                            test == "p.value.nonsup" ~ "Non-superiority")) |>
    group_by(measurements_n, participant_n, time_n, test) |>
    summarise(total = sum(p.value < .005),
              power = mean(p.value < .005),
              ci.lower = prop.test(total, 2000)$conf.int[1],
              ci.upper = prop.test(total, 2000)$conf.int[2],
              .groups = "drop") |>
    ggplot(aes(x = participant_n, y = power, color = factor(measurements_n))) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
    facet_grid(time_n ~ test) +
    labs(
      x = "Total Number of Participants per Site (N)",
      y = "Power (1-\u03b2)",
      color = "Number of Measurements per Timepoint (N)",
      title = "Sample Size Estimation (N per Site)",
      subtitle = expression(~italic(y)[ijtk]~" = (\u03b2"[0]~"+ "~u[k]~" + "~u[i]~") + \u03b2"[1]~"condition"[j]~" + \u03b2"[2]~"time"[t]~" + \u03b2"[3]~"condition:time"[jt]~" + \u03f5"[ijtk]),
      caption = "Simulated to detect a standardised effect, or non-superiority, for condition:time with a SESOI of 0.1\nSite number k fixed at 15, \u03b1 = 0.01 (corrected to 0.005 for two primary outcomes i.e., arm and thigh muscle area)\n1000 simulations per combination of participants, number of measurements, and timepoints"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  plot_sim
}

sample_size_simulation_var_effect  <- function(file1, file2) {
  
  # Using prior data from all Discover Strength studies baseline lean mass (standardised) to get estimates of between study (i.e., site) variance
  
  prior_data <- read_csv(file1) |>
    mutate(
      lean_mass_z = (lean_mass - mean(lean_mass, na.rm=TRUE))/sd(lean_mass, na.rm=TRUE)
    ) 
  
  model_site_variance <- lmer(lean_mass_z ~ 1 + (1 | study),
                              data = prior_data)
  
  site_variance <- VarCorr(model_site_variance)
  
  # Using prior study data to get estimates of measurement error for muscle areas (standardised SEMs)
  
  prior_ma_reliability <- read_csv(file2) |>
    # standardise variables
    mutate(
      arm_ma_z = (arm_ma - mean(arm_ma, na.rm=TRUE))/sd(arm_ma, na.rm=TRUE),
      thigh_ma_z = (thigh_ma - mean(thigh_ma, na.rm=TRUE))/sd(thigh_ma, na.rm=TRUE)
    )
  
  arm_ma_reliability <- SimplyAgree::reli_aov(measure = "arm_ma_z",
                                              item = "measurement_no",
                                              id = "participant",
                                              data = prior_ma_reliability)$SEM[,1]
  
  thigh_ma_reliability <- SimplyAgree::reli_aov(measure = "thigh_ma_z",
                                                item = "measurement_no",
                                                id = "participant",
                                                data = prior_ma_reliability)$SEM[,1]
  # Simulate for power
  
  sim <- function(participant_n = as.double(), site_n = as.double(), time_n = as.double(), measurements_n = as.double(),
                  b0 = as.double(), b_time = as.double(), b_cond = as.double(), b_cond_time = as.double(),         # fixed effects 
                  u_site = as.double(), u_participant_arm = as.double(), u_participant_thigh = as.double(), # random intercepts
                  u_time = as.double(),     # random slope for time
                  arm_ma_error = as.double(), thigh_ma_error = as.double(),       # measurement error
                  ... # helps the function work with pmap() below
  ) {
    
    # set up data structure
    dat <- add_random(site = site_n) |>
      add_random(participant = participant_n, .nested_in = "site") |>
      add_between("participant", cond = c("full_ROM", "lengthened_partial")) |>
      add_recode("cond", "cond_dummy", full_ROM = 0, lengthened_partial = 1) |>
      add_within("participant", time = seq(0,1, by = time_n), 
                 measurement = seq(1:measurements_n)) |>
      add_ranef("site", u_site = u_site) |>
      add_ranef("participant", u_participant_arm = u_participant_arm) |>
      add_ranef("participant", u_participant_thigh = u_participant_thigh) |>
      add_ranef("participant", u_time = u_time) |>
      add_ranef(arm_ma_error = arm_ma_error) |>
      add_ranef(thigh_ma_error = thigh_ma_error) |>
      mutate(arm_ma = (b0 + u_site + u_participant_arm) + ((b_time + u_time) * time) + (b_cond * cond_dummy) + (b_cond_time * cond_dummy * time) + arm_ma_error,
             thigh_ma = (b0 + u_site + u_participant_thigh) + ((b_time + u_time) * time) + (b_cond * cond_dummy) + (b_cond_time * cond_dummy * time) + thigh_ma_error) 
    
    if (time_n == 1 && measurements_n == 1) {
      # run mixed effect model and return relevant values
      model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (1 | participant),
                        data = dat, REML = TRUE)
      
      model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (1 | participant),
                          data = dat, REML = TRUE)
    } else {
      # run mixed effect model and return relevant values
      model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (time | participant),
                        data = dat, REML = TRUE)
      
      model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (time | participant),
                          data = dat, REML = TRUE)
    }
    
        hypothesis_tests_arm <- hypotheses(model_arm, equivalence = c(NA,0.1), vcov = "satterthwaite", df = insight::get_df(model_arm))[4,] |>
      mutate(muscle_site = "arm")
    
    
    hypothesis_tests_thigh <- hypotheses(model_thigh, equivalence = c(NA,0.1), vcov = "satterthwaite", df = insight::get_df(model_thigh))[4,] |>
      mutate(muscle_site = "thigh")
    
    bind_rows(hypothesis_tests_arm, hypothesis_tests_thigh)
  }
  
  
  plan(cluster, workers = 10)
  
  sim <- crossing(
    rep = 1:1000, # number of replicates
    participant_n = seq(5,100, by = 5), # range of participant N
    site_n = 15, # fixed site N
    time_n = c(0.5,1), # to examine inclusion of midpoint test or not
    measurements_n = c(1:3), # range of measurements N,
    u_time = c(0.0304, 0.039)       # random slope for time
  ) |>
    mutate(
      b0 = 0, b_time = 0.05, b_cond = 0, b_cond_time = 0.025,         # fixed effects 
      u_site = sqrt(site_variance$study)[1],  # random intercept site
      u_participant_arm = (1 - sqrt(site_variance$study)[1] - arm_ma_reliability),   # random intercept participant arm
      u_participant_thigh = (1 - sqrt(site_variance$study)[1] - thigh_ma_reliability),  # random intercept participant thigh
      arm_ma_error = arm_ma_reliability,           # measurement error arm
      thigh_ma_error = thigh_ma_reliability         # measurement error thigh  
    ) %>%
    mutate(analysis = future_pmap(., sim)) %>% # not sure why base pipe doesn't work here
    unnest(analysis)
  
  plan(sequential)
  
  sim
  
}

plot_sample_size_estimates_var_effect <- function(sim) {
  plot_sim <- sim |>
    filter(term == "time:cond_dummy") |>
    mutate(
      time_n = case_when(
        time_n == 0.5 ~ "Include Midpoint",
        time_n == 1 ~ "No Midpoint",
      )
    ) |>
    mutate(time_n = factor(time_n, levels = c("No Midpoint", "Include Midpoint"))) |>
    mutate(
      u_time = case_when(
        u_time == 0.0304 ~ "5% Negative Slopes",
        u_time == 0.039 ~ "10% Negative Slopes",
      )
    ) |>
    mutate(time_n = factor(time_n, levels = c("No Midpoint", "Include Midpoint"))) |>
    pivot_longer(c("p.value", "p.value.nonsup"), 
                 names_to = "test", 
                 values_to = "p.value") |> 
    mutate(test = case_when(test == "p.value" ~ "Difference",
                            test == "p.value.nonsup" ~ "Non-superiority")) |>
    group_by(measurements_n, participant_n, time_n, u_time, test) |> 
    summarise(total = sum(p.value < .005),
              power = mean(p.value < .005),
              ci.lower = prop.test(total, 2000)$conf.int[1],
              ci.upper = prop.test(total, 2000)$conf.int[2], 
              .groups = "drop") |>
    ggplot(aes(x = participant_n, y = power, color = factor(measurements_n))) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
    facet_nested(time_n + u_time ~ test) +
    labs(
      x = "Total Number of Participants per Site (N)",
      y = "Power (1-\u03b2)",
      color = "Number of Measurements per Timepoint (N)",
      title = "Sample Size Estimation (N per Site)",
      subtitle = expression(~italic(y)[ijtk]~" = (\u03b2"[0]~" + "~u[k]~" + "~u[i]~") + \u03b2"[1]~"condition"[j]~" + (\u03b2"[2]~" + "~u[j]~")time"[t]~" + \u03b2"[3]~"condition:time"[jt]~" + \u03f5"[ijtk]),
      caption = "Simulated to detect a standardised effect, or non-superiority, for condition:time with a SESOI of 0.1\nSite number k fixed at 15, \u03b1 = 0.01 (corrected to 0.005 for two primary outcomes i.e., arm and thigh muscle area)\n1000 simulations per combination of number of participants, number of measurements, random slopes for time, and timepoints"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  plot_sim
}

sample_size_simulation_est_numbers <- function(file1, file2) {
  
  # Using prior data from all Discover Strength studies baseline lean mass (standardised) to get estimates of between study (i.e., site) variance
  
  prior_data <- read_csv(file1) |>
    mutate(
      lean_mass_z = (lean_mass - mean(lean_mass, na.rm=TRUE))/sd(lean_mass, na.rm=TRUE)
    ) 
  
  model_site_variance <- lmer(lean_mass_z ~ 1 + (1 | study),
                              data = prior_data)
  
  site_variance <- VarCorr(model_site_variance)
  
  # Using prior study data to get estimates of measurement error for muscle areas (standardised SEMs)
  
  prior_ma_reliability <- read_csv(file2) |>
    # standardise variables
    mutate(
      arm_ma_z = (arm_ma - mean(arm_ma, na.rm=TRUE))/sd(arm_ma, na.rm=TRUE),
      thigh_ma_z = (thigh_ma - mean(thigh_ma, na.rm=TRUE))/sd(thigh_ma, na.rm=TRUE)
    )
  
  arm_ma_reliability <- SimplyAgree::reli_aov(measure = "arm_ma_z",
                                              item = "measurement_no",
                                              id = "participant",
                                              data = prior_ma_reliability)$SEM[,1]
  
  thigh_ma_reliability <- SimplyAgree::reli_aov(measure = "thigh_ma_z",
                                                item = "measurement_no",
                                                id = "participant",
                                                data = prior_ma_reliability)$SEM[,1]
  # Simulate for power
  
  sim <- function(site_n = as.double(), time_n = as.double(), measurements_n = as.double(),
                  b0 = as.double(), b_time = as.double(), b_cond = as.double(), b_cond_time = as.double(),         # fixed effects 
                  u_site = as.double(), u_participant_arm = as.double(), u_participant_thigh = as.double(), # random intercepts
                  u_time = as.double(),     # random slope for time
                  arm_ma_error = as.double(), thigh_ma_error = as.double(),       # measurement error
                  ... # helps the function work with pmap() below
  ) {
    
    # set up data structure
    dat <- add_random(site = site_n) |>
      
      # Add predicted recruitment numbers directly into sim function
      add_random(participant = c(5,30,15,40,5,25,30,25,25,15,30,10,5,30,5), .nested_in = "site") |>
      
      add_between("participant", cond = c("full_ROM", "lengthened_partial")) |>
      add_recode("cond", "cond_dummy", full_ROM = 0, lengthened_partial = 1) |>
      add_within("participant", time = seq(0,1, by = time_n), 
                 measurement = seq(1:measurements_n)) |>
      add_ranef("site", u_site = u_site) |>
      add_ranef("participant", u_participant_arm = u_participant_arm) |>
      add_ranef("participant", u_participant_thigh = u_participant_thigh) |>
      add_ranef(arm_ma_error = arm_ma_error) |>
      add_ranef(thigh_ma_error = thigh_ma_error) |>
      mutate(arm_ma = (b0 + u_site + u_participant_arm) + (b_time * time) + (b_cond * cond_dummy) + (b_cond_time * cond_dummy * time) + arm_ma_error,
             thigh_ma = (b0 + u_site + u_participant_thigh) + (b_time * time) + (b_cond * cond_dummy) + (b_cond_time * cond_dummy * time) + thigh_ma_error) 
    
    # 
    # # run mixed effect model and return relevant values
    # model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (1 | participant),
    #                   data = dat, REML = TRUE)
    # 
    # model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (1 | participant),
    #                     data = dat, REML = TRUE)
    # 
    
    if (time_n == 1 && measurements_n == 1) {
      # run mixed effect model and return relevant values
      model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (1 | participant),
                        data = dat, REML = TRUE)
      
      model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (1 | participant),
                          data = dat, REML = TRUE)
    } else {
      # run mixed effect model and return relevant values
      model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (time | participant),
                        data = dat, REML = TRUE)
      
      model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (time | participant),
                          data = dat, REML = TRUE)
    }
    
    hypothesis_tests_arm <- hypotheses(model_arm, equivalence = c(NA,0.1), vcov = "satterthwaite", df = insight::get_df(model_arm))[4,] |>
      mutate(muscle_site = "arm")
    
    
    hypothesis_tests_thigh <- hypotheses(model_thigh, equivalence = c(NA,0.1), vcov = "satterthwaite", df = insight::get_df(model_thigh))[4,] |>
      mutate(muscle_site = "thigh")
    
    bind_rows(hypothesis_tests_arm, hypothesis_tests_thigh)
  }
  
  
  plan(cluster, workers = 10)
  
  sim <- crossing(
    rep = 1:1000, # number of replicates
    site_n = 15, # fixed site N
    time_n = c(0.5,1), # to examine inclusion of midpoint test or not
    measurements_n = c(1:3), # range of measurements N
    u_time = c(0.0304, 0.039)       # random slope for time
  ) |>
    mutate(
      b0 = 0, b_time = 0.05, b_cond = 0, b_cond_time = 0.025,         # fixed effects 
      u_site = sqrt(site_variance$study)[1],  # random intercept site
      u_participant_arm = (1 - sqrt(site_variance$study)[1] - arm_ma_reliability),   # random intercept participant arm
      u_participant_thigh = (1 - sqrt(site_variance$study)[1] - thigh_ma_reliability),  # random intercept participant thigh
      arm_ma_error = arm_ma_reliability,           # measurement error arm
      thigh_ma_error = thigh_ma_reliability         # measurement error thigh  
    ) %>%
    mutate(analysis = future_pmap(., sim)) %>% # not sure why base pipe doesn't work here
    unnest(analysis)
  
  plan(sequential)
  
  sim
  
}

plot_sample_size_estimates_est_numbers <- function(sim) {
  plot_sim <- sim |>
    filter(term == "time:cond_dummy") |>
    mutate(
      time_n = case_when(
        time_n == 0.5 ~ "Include Midpoint",
        time_n == 1 ~ "No Midpoint",
      )
    ) |>
    mutate(time_n = factor(time_n, levels = c("No Midpoint", "Include Midpoint"))) |>
    mutate(
      u_time = case_when(
        u_time == 0.0304 ~ "5% Negative Slopes",
        u_time == 0.039 ~ "10% Negative Slopes",
      )
    ) |>
    mutate(time_n = factor(time_n, levels = c("No Midpoint", "Include Midpoint"))) |>
    pivot_longer(c("p.value", "p.value.nonsup"), 
                 names_to = "test", 
                 values_to = "p.value") |> 
    mutate(test = case_when(test == "p.value" ~ "Difference",
                            test == "p.value.nonsup" ~ "Non-superiority")) |>
    group_by(measurements_n, time_n, u_time, test) |> 
    summarise(total = sum(p.value < .005),
              power = mean(p.value < .005),
              ci.lower = prop.test(total, 2000)$conf.int[1],
              ci.upper = prop.test(total, 2000)$conf.int[2], 
              .groups = "drop") |>
    ggplot(aes(x = factor(measurements_n), y = power, color = factor(measurements_n))) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
    facet_nested(time_n + u_time ~ test) +
    labs(
      x = "Number of Measurements per Timepoint (N)",
      y = "Power (1-\u03b2)",
      color = "Number of Measurements per Timepoint (N)",
      title = "Simulated Power with Site Predicted Recruitment Numbers",
      subtitle = expression(~italic(y)[ijtk]~" = (\u03b2"[0]~" + "~u[k]~" + "~u[i]~") + \u03b2"[1]~"condition"[j]~" + (\u03b2"[2]~" + "~u[j]~")time"[t]~" + \u03b2"[3]~"condition:time"[jt]~" + \u03f5"[ijtk]),
      caption = "Simulated to detect a standardised effect, or non-superiority, for condition:time with a SESOI of 0.1\nSite number k fixed at 15, \u03b1 = 0.01 (corrected to 0.005 for two primary outcomes i.e., arm and thigh muscle area)\n1000 simulations per combination of number of participants, number of measurements, random slopes for time, and timepoints"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  plot_sim
}

# Functions for pre-registered primary hypertrophy outcomes ----
prepare_hypertrophy_data <- function(data) {
  data <- data |>
    filter(time != "1") |>
    mutate(
      thigh_ma_1 =(4.68*thigh_circum_1)-(2.09*thigh_skinfold_1)-80.99,
      thigh_ma_2 =(4.68*thigh_circum_2)-(2.09*thigh_skinfold_2)-80.99,
      thigh_ma_3 =(4.68*thigh_circum_3)-(2.09*thigh_skinfold_3)-80.99,
      arm_ma_1 = if_else(sex == "M" ,((arm_circum_1 - (pi * (arm_skinfold_1/10)))^2 / (4*pi))-10, ((arm_circum_1 - (pi * (arm_skinfold_1/10)))^2 / (4*pi))-6.5),
      arm_ma_2 = if_else(sex == "M" ,((arm_circum_2 - (pi * (arm_skinfold_2/10)))^2 / (4*pi))-10, ((arm_circum_2 - (pi * (arm_skinfold_2/10)))^2 / (4*pi))-6.5),
      arm_ma_3 = if_else(sex == "M" ,((arm_circum_3 - (pi * (arm_skinfold_3/10)))^2 / (4*pi))-10, ((arm_circum_3 - (pi * (arm_skinfold_3/10)))^2 / (4*pi))-6.5)
    ) |>
    select(site_id, participant_id, condition, time, contains("ma")) |>
    pivot_longer(5:10,
                 names_to = "muscle",
                 values_to = "estimated_ma") |>
    mutate(
      measurement = case_when(
        str_detect(muscle, "1") == TRUE ~ 1,
        str_detect(muscle, "2") == TRUE ~ 2,
        str_detect(muscle, "3") == TRUE ~ 3
      ),
      muscle = case_when(
        str_detect(muscle, "thigh") == TRUE ~ "thigh",
        str_detect(muscle, "arm") == TRUE ~ "arm"
      )
    ) |>
    mutate(cond_dummy = case_when(
      condition == "fROM" ~ -0.5,
      condition == "lpROM" ~ 0.5
    ),
    timepoint = factor(case_when(
      time == 2 ~ "Pre",
      time == 3 ~ "Post"
    ), levels = c("Pre", "Post"))
    ) |>
    select(site_id, participant_id, condition, cond_dummy, time, timepoint, muscle, estimated_ma)
    

}

fit_hypertrophy_models <- function(data) {
  
  model_results <- tibble(term = as.character(),
                               estimate = as.numeric(),
                               std.error = as.numeric(),
                               statistic = as.numeric(),
                               p.value = as.numeric(),
                               s.value = as.numeric(),
                               conf.low = as.numeric(),
                               conf.high = as.numeric(),
                               df = as.integer(),
                               statistic.noninf = as.numeric(),
                               statistic.nonsup = as.numeric(),
                               p.value.noninf = as.numeric(),
                               p.value.nonsup = as.numeric(),
                               p.value.equiv = as.numeric(),
                               muscle = as.character())
  
  for(i in c("thigh", "arm")) {
    
    data_muscle <- data |>
      filter(muscle == i)
    
    std_lm <- lm(estimated_ma ~ cond_dummy,
                     data = data_muscle)
    
    std <- summary(std_lm)$sigma
    
    data_muscle <- data_muscle |>
      mutate(
        estimated_ma_z = (estimated_ma - mean(estimated_ma, na.rm=TRUE)) / std
      )
    
    model <- lmer(estimated_ma_z ~ time*cond_dummy + (1|site_id) + (time | participant_id),
                      data = data_muscle, REML = TRUE)
    
    tests <- hypotheses(model, equivalence = c(-0.1,0.1), vcov = "satterthwaite", df = insight::get_df(model)) |>
      mutate(muscle = i)
    
    model_results <- bind_rows(model_results,
                              tests)
 
  }
  
  return(model_results)
  
}

plot_hypertrophy_models <- function(data, models) {
  
  raw_plot_arm <- data |>
    filter(muscle == "arm") |>
    ggplot(aes(x=timepoint, y=estimated_ma, color=condition, fill=condition)) +
    geom_line(aes(group=interaction(participant_id)),
              stat = "smooth", se = FALSE, method = "lm",
              alpha = 0.25, size = 0.5) +
    stat_slabinterval(data = filter(data, timepoint == "Pre" & muscle == "arm"),
                      point_interval = "mean_qi",
                      .width = .95,
                      side = c("left"),
                      slab_alpha = 0.25,
                      position = position_dodge(w=-0.1)) +
    stat_slabinterval(data = filter(data, timepoint == "Post" & muscle == "arm"),
                      point_interval = "mean_qi",
                      .width = .95,
                      side = c("right"),
                      slab_alpha = 0.25,
                      position = position_dodge(w=-0.1)) +
    scale_color_manual(values = c("#56B4E9","#E69F00")) +
    scale_fill_manual(values = c("#56B4E9","#E69F00")) +
    labs(
      x = "Timepoint",
      y = expression(paste("Estimated Arm Muscle CSA (",cm^2,")")),
      title = "Arm Muscle",
      color = "Condition",
      fill = "Condition"
    ) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom")
  
  raw_plot_thigh <- data |>
    filter(muscle == "thigh") |>
    ggplot(aes(x=timepoint, y=estimated_ma, color=condition, fill=condition)) +
    geom_line(aes(group=interaction(participant_id)),
              stat = "smooth", se = FALSE, method = "lm",
              alpha = 0.25, size = 0.5) +
    stat_slabinterval(data = filter(data, timepoint == "Pre" & muscle == "thigh"),
                      point_interval = "mean_qi",
                      .width = .95,
                      side = c("left"),
                      slab_alpha = 0.25,
                      position = position_dodge(w=-0.1)) +
    stat_slabinterval(data = filter(data, timepoint == "Post" & muscle == "thigh"),
                      point_interval = "mean_qi",
                      .width = .95,
                      side = c("right"),
                      slab_alpha = 0.25,
                      position = position_dodge(w=-0.1)) +
    scale_color_manual(values = c("#56B4E9","#E69F00")) +
    scale_fill_manual(values = c("#56B4E9","#E69F00")) +
    labs(
      x = "Timepoint",
      y = expression(paste("Estimated Thigh Muscle CSA (",cm^2,")")),
      title = "Thigh Muscle",
      color = "Condition",
      fill = "Condition"
    ) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom")
  
  models <- models |>
    mutate(
      muscle = case_when(
        muscle == "arm" ~ "Arm Muscle",
        muscle == "thigh" ~ "Thigh Muscle"
      )
    )
  
  time_plot <- models |>
    filter(term == "time") |>
    ggplot(aes(x=muscle, y=estimate)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_hline(yintercept = c(-0.1,0.1), linetype = 2) +
    geom_pointinterval(aes(ymin=conf.low, ymax=conf.high)) +
    annotate("text", x = 1.5, y = c(-0.09,0.09),
             label = "Smallest Effect Size of Interest Limits",
             size = 2) +
    labs(
      y = "Standardised Mean Effect",
      title = "Main effect of time"
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  interact_plot <- models |>
    filter(term == "time:cond_dummy") |>
    ggplot(aes(x=muscle, y=estimate)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_hline(yintercept = c(-0.1,0.1), linetype = 2) +
    geom_pointinterval(aes(ymin=conf.low, ymax=conf.high)) +
    geom_label(aes(label = glue("Equivalence Test\np = {scales::pvalue(p.value.equiv, 0.001)}")),
              position = position_nudge(x = c(0.25,-0.25)),
              size = 2.5) + 
    annotate("text", x = 1.5, y = c(-0.09,0.09),
             label = "Smallest Effect Size of Interest Limits",
             size = 2) +
    labs(
      y = "Standardised Mean Effect",
      title = "Condition by time interaction"
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  
  raw_plots <- (raw_plot_arm | raw_plot_thigh) +
    plot_annotation(title = "Primary Pre-registered Hypertrophy Outcomes",
                    caption = "Raw means [95% quantile intervals] and participant level unpooled predictions") +
    plot_layout(guides = "collect", axes = "collect") &
    theme(legend.position = "bottom")
  
  coef_plots <- (time_plot / interact_plot) +
    plot_layout(axes = "collect") +
    plot_annotation(caption = "Note, alpha set at 0.01 and adjusted to 0.005 to account for multiple outcomes\nError bars are 95% confidence intervals")
  
  wrap_elements(raw_plots) / wrap_elements(coef_plots) 
}

# Functions for strength outcomes ----
prepare_strength_data <- function(data) {
  data <- data |>
    filter(time <= 84) |>
    mutate(cond_dummy = case_when(
      condition == "fROM" ~ -0.5,
      condition == "lpROM" ~ 0.5
    ),
    time_scaled = time / 84
    )

}

fit_strength_models <- function(data) {
  
  model_results <- tibble(term = as.character(),
                          estimate = as.numeric(),
                          std.error = as.numeric(),
                          statistic = as.numeric(),
                          p.value = as.numeric(),
                          s.value = as.numeric(),
                          conf.low = as.numeric(),
                          conf.high = as.numeric(),
                          df = as.integer(),
                          statistic.noninf = as.numeric(),
                          statistic.nonsup = as.numeric(),
                          p.value.noninf = as.numeric(),
                          p.value.nonsup = as.numeric(),
                          p.value.equiv = as.numeric(),
                          machine_type = as.character())
  
  for(i in c("cp", "lp", "pd")) {
    
    data_strength <- data |>
      filter(machine_type == i)
    
    std_lm <- lm(estimated_1RM ~ cond_dummy,
                 data = data_strength)
    
    std <- summary(std_lm)$sigma
    
    data_strength <- data_strength |>
      mutate(
        estimated_1RM_z = (estimated_1RM - mean(estimated_1RM, na.rm=TRUE)) / std
      )
    
    model <- lmer(estimated_1RM_z ~ time_scaled*cond_dummy + (time_scaled | participant_id/exercise) + (1 | site_id),
                  data = data_strength, REML = TRUE)
    
    tests <- hypotheses(model, vcov = "satterthwaite", df = insight::get_df(model)) |>
      mutate(machine_type = i)
    
    model_results <- bind_rows(model_results,
                               tests)
    
  }
  
  return(model_results)
  
}

get_prepare_MA_data_strength <- function(file) {
  
  ### Note code copied from original paper 
  
  ##### Read csv as data frame into environment - Note: change source address
  Data <- read.csv(here::here("data","Polito et al. RT Extracted Data.csv"), na.strings=c(""," ","NA"))
  
  # Calculate pre-post SDs from SEs
  Data$RT_pre_sd <- ifelse(is.na(Data$RT_pre_se), Data$RT_pre_sd, Data$RT_pre_se * sqrt(Data$RT_n))
  Data$CON_pre_sd <- ifelse(is.na(Data$CON_pre_se), Data$CON_pre_sd, Data$CON_pre_se * sqrt(Data$CON_n))
  Data$RT_post_sd <- ifelse(is.na(Data$RT_post_se), Data$RT_post_sd, Data$RT_post_se * sqrt(Data$RT_n))
  Data$CON_post_sd <- ifelse(is.na(Data$CON_post_se), Data$CON_post_sd, Data$CON_post_se * sqrt(Data$CON_n))
  
  # Convert p to t (Change scores)
  Data$RT_delta_t_value <- replmiss(Data$RT_delta_t_value, with(Data, qt(RT_delta_p_value/2, df=RT_n-1, lower.tail=FALSE)))
  Data$CON_delta_t_value <- replmiss(Data$CON_delta_t_value, with(Data, qt(CON_delta_p_value/2, df=CON_n-1, lower.tail=FALSE)))
  
  # Convert t to SE (Change scores)
  Data$RT_delta_se <- replmiss(Data$RT_delta_se, with(Data, ifelse(is.na(RT_delta_m), 
                                                                   (RT_post_m - RT_pre_m)/RT_delta_t_value, RT_delta_m/RT_delta_t_value)))
  Data$CON_delta_se <- replmiss(Data$CON_delta_se, with(Data, ifelse(is.na(CON_delta_m), 
                                                                     (CON_post_m - CON_pre_m)/CON_delta_t_value, CON_delta_m/CON_delta_t_value)))
  # Make positive
  Data$RT_delta_se <- ifelse(Data$RT_delta_se < 0, Data$RT_delta_se * -1, Data$RT_delta_se)
  Data$CON_delta_se <- ifelse(Data$CON_delta_se < 0, Data$CON_delta_se * -1, Data$CON_delta_se)
  
  # Convert CI to SE (Change scores)
  Data$RT_delta_se <- replmiss(Data$RT_delta_se, with(Data, (RT_delta_CI_upper - RT_delta_CI_lower)/3.92))
  Data$CON_delta_se <- replmiss(Data$CON_delta_se, with(Data, (CON_delta_CI_upper - CON_delta_CI_lower)/3.92))
  
  # Convert SE to SD (Change scores)
  Data$RT_delta_sd <- replmiss(Data$RT_delta_sd, with(Data, RT_delta_se * sqrt(RT_n)))
  Data$CON_delta_sd <- replmiss(Data$CON_delta_sd, with(Data, CON_delta_se * sqrt(CON_n)))
  
  # Calculate pre-post correlation coefficient for those with pre, post, and delta SDs
  Data$RT_ri <- (Data$RT_pre_sd^2 + Data$RT_post_sd^2 - Data$RT_delta_sd^2)/(2 * Data$RT_pre_sd * Data$RT_post_sd)
  Data$CON_ri <- (Data$CON_pre_sd^2 + Data$CON_post_sd^2 - Data$CON_delta_sd^2)/(2 * Data$CON_pre_sd * Data$CON_post_sd)
  
  # Remove values outside the range of -1 to +1 as they are likely due to misreporting or miscalculations in original studies
  Data$RT_ri <- ifelse(between(Data$RT_ri,-1,1) == FALSE, NA, Data$RT_ri)
  Data$CON_ri <- ifelse(between(Data$CON_ri,-1,1) == FALSE, NA, Data$CON_ri)
  
  # Then we'll convert using Fishers r to z, calculate a meta-analytic point estimate, and impute that across the studies with missing correlations
  Data <- escalc(measure = "ZCOR", ri = RT_ri, ni = RT_n, data = Data)
  
  Meta_RT_ri <- rma.mv(yi, V=vi, data=Data,
                       slab=paste(label),
                       random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                       control=list(optimizer="optim", optmethod="Nelder-Mead"))
  
  RobuEstMeta_RT_ri <- robust(Meta_RT_ri, Data$study)
  
  z2r_RT <- psych::fisherz2r(RobuEstMeta_RT_ri$b[1])
  
  Data$RT_ri <- ifelse(is.na(Data$RT_ri), z2r_RT, Data$RT_ri)
  
  Data <- escalc(measure = "ZCOR", ri = CON_ri, ni = CON_n, data = Data)
  
  ### Note, data is coded with study and arm as having explicit nesting so all random effects are (~ 1 | study, ~ 1 | arm)
  Meta_CON_ri <- rma.mv(yi, V=vi, data=Data,
                        slab=paste(label),
                        random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                        control=list(optimizer="optim", optmethod="Nelder-Mead"))
  
  RobuEstMeta_CON_ri <- robust(Meta_CON_ri, Data$study)
  
  z2r_CON <- psych::fisherz2r(RobuEstMeta_CON_ri$b[1])
  
  Data$CON_ri <- ifelse(is.na(Data$CON_ri), z2r_CON, Data$CON_ri)
  
  # Estimate change score difference SD where only pre-post data available
  Data$RT_delta_sd <- replmiss(Data$RT_delta_sd, with(Data, sqrt(RT_pre_sd^2 + RT_post_sd^2 - (2*RT_ri*RT_pre_sd*RT_post_sd))))
  Data$CON_delta_sd <- replmiss(Data$CON_delta_sd, with(Data, sqrt(CON_pre_sd^2 + CON_post_sd^2 - (2*CON_ri*CON_pre_sd*CON_post_sd))))
  
  # Arm based effect size calculations
  Data_strength <- Data %>% 
    filter(outcome == "strength")
  
  Data_SMD_strength_RT <- Data_strength |>
    select(study, arm, es, RT_n, RT_pre_m, RT_post_m, RT_pre_sd, RT_post_sd, RT_ri, weeks)
  
  Data_SMD_strength_RT <-  escalc(measure = "SMCR",
                                     m1i = RT_post_m, m2i = RT_pre_m,
                                     sd1i = RT_post_sd, 
                                     ri = RT_ri, ni = RT_n,
                                     data = Data_SMD_strength_RT
  )
  
  Data_SMD_strength_RT <- Data_SMD_strength_RT %>%
    filter(!is.na(vi) &
             !is.infinite(vi)) %>%
    mutate(se = sqrt(vi),
           wi = 1/sqrt(vi),
           size = 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))) |>
    select(study, arm, es, yi, vi, weeks)   |>
    mutate(arm = as.factor(unclass(factor(unlist(arm)))),
           es = as.factor(unclass(factor(unlist(es)))))
  
  
  Data_SMD_strength_CON <- Data_strength |>
    select(study, arm, es, CON_n, CON_pre_m, CON_post_m, CON_pre_sd, CON_post_sd, CON_ri, weeks)
  
  Data_SMD_strength_CON <-  escalc(measure = "SMCR",
                                      m1i = CON_post_m, m2i = CON_pre_m,
                                      sd1i = CON_post_sd, 
                                      ri = CON_ri, ni = CON_n,
                                      data = Data_SMD_strength_CON
  )
  
  Data_SMD_strength_CON <- Data_SMD_strength_CON %>%
    filter(!is.na(vi) &
             !is.infinite(vi)) %>%
    mutate(se = sqrt(vi),
           wi = 1/sqrt(vi),
           size = 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))) |>
    select(study, arm, es, yi, vi, weeks) %>%
    mutate(weeks = 0)    |>
    distinct(yi, .keep_all = TRUE) |>
    mutate(arm = as.factor(unclass(factor(unlist(arm)))+length(unique(Data_SMD_strength_RT$arm))),
           es = as.factor(unclass(factor(unlist(es)))+length(unique(Data_SMD_strength_RT$es))))
  
  
  data <- rbind(Data_SMD_strength_RT, Data_SMD_strength_CON)
}

log_growth_MA_model_strength <- function(data) {
  
  
  MultiLevelModel_SMD_strength_weeks_arm <- rma.mv(yi, V=vi, data=data,
                                                      random = list(~ log1p(weeks) | study, ~ 1 | arm, ~ 1 | es), mods = ~ log1p(weeks),
                                                      method="REML", test="t", struct = "GEN"
  )
  
  RobuEstMultiLevelModel_SMD_strength_weeks_arm <- robust(MultiLevelModel_SMD_strength_weeks_arm, data$study)
  
  
}

plot_log_growth_MA_model_strength <- function(model, data) {
  data_log <- cbind(data, predict(model)) %>%
    mutate(wi = 1/sqrt(vi),
           size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi)))
  
  preds_arm <- as.data.frame(predict(model, newmods = c(log1p(0:80)))) |>
    mutate(weeks = 0:80)
  
  preds_arm_12_24_36 <- as.data.frame(predict(model, newmods = log1p(c(12,24,36)))) |>
    mutate(weeks = c(12,24,36))
  
  
  plot_log_growth_MA_model <- ggplot(data_log, aes(x=weeks)) +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_point(aes(y=yi,size = size), alpha = 0.1) +
    geom_ribbon(data=preds_arm,aes(ymax=ci.ub, ymin=ci.lb), alpha = 0.2) +
    geom_line(data=preds_arm, aes(y=pred)) +
    annotate("text", 
             x = 12, y = -1, size = 3,
             label = glue::glue("12wks: {round(preds_arm_12_24_36$pred[1],2)}")) +
    annotate("text", 
             x = 12, y = -1.2, size = 3,
             label = glue::glue("[95%CI: {round(preds_arm_12_24_36$ci.lb[1],2)}, {round(preds_arm_12_24_36$ci.ub[1],2)}]")) +
    annotate("text", 
             x = 24, y = -1, size = 3,
             label = glue::glue("24wks: {round(preds_arm_12_24_36$pred[2],2)}")) +
    annotate("text", 
             x = 24, y = -1.2, size = 3,
             label = glue::glue("[95%CI: {round(preds_arm_12_24_36$ci.lb[2],2)}, {round(preds_arm_12_24_36$ci.ub[2],2)}]")) +
    annotate("text", 
             x = 36, y = -1, size = 3,
             label = glue::glue("36wks: {round(preds_arm_12_24_36$pred[3],2)}")) +
    annotate("text", 
             x = 36, y = -1.2, size = 3,
             label = glue::glue("[95%CI: {round(preds_arm_12_24_36$ci.lb[3],2)}, {round(preds_arm_12_24_36$ci.ub[3],2)}]")) +
    labs(y = "Standardised Mean Change", x = "Weeks") +
    theme_classic() +
    labs(title = "Strength Outcomes",
         subtitle = "Arm-based model - random intercepts (effect, arm, study) and slopes (study)",
         caption = "Data from Steele et al. (2023), DOI: 10.1080/02640414.2023.2286748") +
    guides(size = "none", fill = "none")
  
  plot_log_growth_MA_model
}

  # Note, for some reason this function (or at least it's individual components) will work outside of targets but not inside
  # We get an error in qdrg() - 'data' must be a data.frame, environment, or list
  # So you will not see it in the pipeline, but the interval estimates have been manually added to the exploratory plot
  refit_MA_model <- function(data) {
    
    data <- as.data.frame(data)
    
    MultiLevelModel_SMD_strength_weeks_arm <- rma.mv(yi, V=vi, data=data,
                                                        random = list(~ log1p(weeks) | study, ~ 1 | arm, ~ 1 | es), mods = ~ log1p(weeks),
                                                        method="REML", test="t", struct = "GEN"
    )
    
    RobuEstMultiLevelModel_SMD_strength_weeks_arm <- robust(MultiLevelModel_SMD_strength_weeks_arm, data$study)
    
    ref_grid <- qdrg(object=RobuEstMultiLevelModel_SMD_strength_weeks_arm, data=data, at = list(weeks = c(24, 36)))
    
    emt <- update(contrast(ref_grid, "consec", simple = "weeks"), by = NULL)
    
    slope <- confint(emt)
    
    slope
  }

plot_strength_models <- function(data, models) {
  
  data <- data |>
    filter(!is.na(participant_id))
  
  raw_plot_cp <- ggplot() +
    geom_line(data = filter(data, machine_type == "cp" & condition == "lpROM"),
                            aes(x=time, y=estimated_1RM,
                                group=interaction(participant_id, exercise)),
              stat = "smooth", method = "gam", se = FALSE, color = "#56B4E9",
              alpha = 0.25, size = 0.25) +
    geom_line(data = filter(data, machine_type == "cp" & condition == "fROM"),
              aes(x=time, y=estimated_1RM,
                  group=interaction(participant_id, exercise)),
              stat = "smooth", method = "gam", se = FALSE, color = "#E69F00",
              alpha = 0.5, size = 0.25) +
    labs(
      x = "Time (days)",
      y = "Estimated One Repetition Maximum (kg)",
      title = "Chest Press",
      color = "Condition",
      fill = "Condition"
    ) +
    theme_classic(base_size = 8) +
    theme(legend.position = "bottom")
  
  group_plot_cp <- ggplot() +
    geom_smooth(data = filter(data, machine_type == "cp" & condition == "lpROM"),
                aes(x=time, y=estimated_1RM, group=1,
                    color = "lpROM", fill = "lpROM"),
                method = "gam") +
    geom_smooth(data = filter(data, machine_type == "cp" & condition == "fROM"),
                aes(x=time, y=estimated_1RM, group=1,
                    color = "fROM", fill = "fROM"),
                method = "gam") +
    scale_color_manual(values = c("fROM" = "#56B4E9", "lpROM" = "#E69F00")) +
    scale_fill_manual(values = c("fROM" = "#56B4E9", "lpROM" = "#E69F00")) +
    labs(
      x = "Time (days)",
      y = "Estimated One Repetition Maximum (kg)",
      color = "Condition",
      fill = "Condition"
    ) +
    theme_classic(base_size = 8) +
    theme(legend.position = "bottom")
  
  raw_plot_lp <- ggplot() +
    geom_line(data = filter(data, machine_type == "lp" & condition == "lpROM"),
              aes(x=time, y=estimated_1RM,
                  group=interaction(participant_id, exercise)),
              stat = "smooth", method = "gam", se = FALSE, color = "#56B4E9",
              alpha = 0.25, size = 0.25) +
    geom_line(data = filter(data, machine_type == "lp" & condition == "fROM"),
              aes(x=time, y=estimated_1RM,
                  group=interaction(participant_id, exercise)),
              stat = "smooth", method = "gam", se = FALSE, color = "#E69F00",
              alpha = 0.5, size = 0.25) +
    labs(
      x = "Time (days)",
      y = "Estimated One Repetition Maximum (kg)",
      title = "Leg Press",
      color = "Condition",
      fill = "Condition"
    ) +
    theme_classic(base_size = 8) +
    theme(legend.position = "bottom")
  
  group_plot_lp <- ggplot() +
    geom_smooth(data = filter(data, machine_type == "lp" & condition == "lpROM"),
                aes(x=time, y=estimated_1RM, group=1,
                    color = "lpROM", fill = "lpROM"), 
                method = "gam") +
    geom_smooth(data = filter(data, machine_type == "lp" & condition == "fROM"),
                aes(x=time, y=estimated_1RM, group=1,
                    color = "fROM", fill = "fROM"), 
                method = "gam") +
    scale_color_manual(values = c("fROM" = "#56B4E9", "lpROM" = "#E69F00")) +
    scale_fill_manual(values = c("fROM" = "#56B4E9", "lpROM" = "#E69F00")) +
    labs(
      x = "Time (days)",
      y = "Estimated One Repetition Maximum (kg)",
      color = "Condition",
      fill = "Condition"
    ) +
    theme_classic(base_size = 8) +
    theme(legend.position = "bottom")
  
  raw_plot_pd <- ggplot() +
    geom_line(data = filter(data, machine_type == "pd" & condition == "lpROM"),
              aes(x=time, y=estimated_1RM,
                  group=interaction(participant_id, exercise)),
              stat = "smooth", method = "gam", se = FALSE, color = "#56B4E9",
              alpha = 0.25, size = 0.25) +
    geom_line(data = filter(data, machine_type == "pd" & condition == "fROM"),
              aes(x=time, y=estimated_1RM,
                  group=interaction(participant_id, exercise)),
              stat = "smooth", method = "gam", se = FALSE, color = "#E69F00",
              alpha = 0.5, size = 0.25) +
    labs(
      x = "Time (days)",
      y = "Estimated One Repetition Maximum (kg)",
      title = "Pulldown",
      color = "Condition",
      fill = "Condition"
    ) +
    theme_classic(base_size = 8) +
    theme(legend.position = "bottom")
  
  group_plot_pd <- ggplot() +
    geom_smooth(data = filter(data, machine_type == "pd" & condition == "lpROM"),
                aes(x=time, y=estimated_1RM, group=1,
                    color = "lpROM", fill = "lpROM"), 
                method = "gam") +
    geom_smooth(data = filter(data, machine_type == "pd" & condition == "fROM"),
                aes(x=time, y=estimated_1RM, group=1,
                    color = "fROM", fill = "fROM"), 
                method = "gam") +
    scale_color_manual(values = c("fROM" = "#56B4E9", "lpROM" = "#E69F00")) +
    scale_fill_manual(values = c("fROM" = "#56B4E9", "lpROM" = "#E69F00")) +
    labs(
      x = "Time (days)",
      y = "Estimated One Repetition Maximum (kg)",
      color = "Condition",
      fill = "Condition"
    ) +
    theme_classic(base_size = 8) +
    theme(legend.position = "bottom")
  
  models <- models |>
    mutate(
      machine_type = case_when(
        machine_type == "cp" ~ "Chest Press",
        machine_type == "lp" ~ "Leg Press",
        machine_type == "pd" ~ "Pulldown"
      )
    )
  
  time_plot <- models |>
    filter(term == "time_scaled") |>
    ggplot(aes(x=machine_type, y=estimate)) +
    geom_hline(yintercept = 0, linetype = 3) +
    annotate("rect", 
             xmin = -Inf, xmax = Inf, 
             ymin = 0.105, ymax = 0.129,
             alpha = 0.5) +
    annotate("text", x = 2, y = 0.117,
             label = "Predicted Time Effect Interval Estimate from Meta-Regression",
             size = 2) +
    scale_x_discrete() +
    geom_pointinterval(aes(ymin=conf.low, ymax=conf.high)) +
    labs(
      y = "Standardised Mean Effect",
      title = "Main effect of time"
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  interact_plot <- models |>
    filter(term == "time_scaled:cond_dummy") |>
    ggplot(aes(x=machine_type, y=estimate)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_pointinterval(aes(ymin=conf.low, ymax=conf.high)) +
    labs(
      y = "Standardised Mean Effect",
      title = "Condition by time interaction"
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  raw_plots <- (raw_plot_cp | raw_plot_lp | raw_plot_pd) +
    plot_annotation(title = "Secondary Exploratory Strength Outcomes",
                    caption = "Penalised cubic spline smooths applied to participant and exercise level raw data") +
    plot_layout(axes = "collect")
  
  group_plots <- (group_plot_cp | group_plot_lp | group_plot_pd) +
    plot_annotation(caption = "Penalised cubic spline smooths applied to group level raw data") +
    plot_layout(guides = "collect", axes = "collect") &
    theme(legend.position = "bottom")
  
  coef_plots <- (time_plot / interact_plot) +
    plot_annotation(caption = "Error bars are 95% confidence intervals") +
    plot_layout(axes = "collect")
  
  wrap_elements(raw_plots) / wrap_elements(group_plots) /  wrap_elements(coef_plots)
  
}


# Additional exploratory models for hypertrophy outcomes ----
# Here we drop the random slopes, utilise a multilevel ANCOVA, and include all baseline data
prepare_explore_hypertrophy_data <- function(data) {
  data <- data |>
    mutate(
      thigh_ma_1 =(4.68*thigh_circum_1)-(2.09*thigh_skinfold_1)-80.99,
      thigh_ma_2 =(4.68*thigh_circum_2)-(2.09*thigh_skinfold_2)-80.99,
      thigh_ma_3 =(4.68*thigh_circum_3)-(2.09*thigh_skinfold_3)-80.99,
      arm_ma_1 = if_else(sex == "M" ,((arm_circum_1 - (pi * (arm_skinfold_1/10)))^2 / (4*pi))-10, ((arm_circum_1 - (pi * (arm_skinfold_1/10)))^2 / (4*pi))-6.5),
      arm_ma_2 = if_else(sex == "M" ,((arm_circum_2 - (pi * (arm_skinfold_2/10)))^2 / (4*pi))-10, ((arm_circum_2 - (pi * (arm_skinfold_2/10)))^2 / (4*pi))-6.5),
      arm_ma_3 = if_else(sex == "M" ,((arm_circum_3 - (pi * (arm_skinfold_3/10)))^2 / (4*pi))-10, ((arm_circum_3 - (pi * (arm_skinfold_3/10)))^2 / (4*pi))-6.5)
    ) |>
    select(site_id, participant_id, condition, time, contains("ma")) |>
    pivot_longer(5:10,
                 names_to = "muscle",
                 values_to = "estimated_ma") |>
    mutate(
      measurement = case_when(
        str_detect(muscle, "1") == TRUE ~ 1,
        str_detect(muscle, "2") == TRUE ~ 2,
        str_detect(muscle, "3") == TRUE ~ 3
      ),
      muscle = case_when(
        str_detect(muscle, "thigh") == TRUE ~ "thigh",
        str_detect(muscle, "arm") == TRUE ~ "arm"
      )
    ) |>
    mutate(cond_dummy = case_when(
      condition == "fROM" ~ -0.5,
      condition == "lpROM" ~ 0.5
    ),
    timepoint = factor(case_when(
      time == 1 ~ "Pre",
      time == 2 ~ "Pre",
      time == 3 ~ "Post"
    ), levels = c("Pre", "Post"))
    ) |>
    select(site_id, participant_id, condition, cond_dummy, time, timepoint, muscle, estimated_ma)
  
  
}

fit_explore_hypertrophy_models <- function(data) {
  
  model_results <- tibble(term = as.character(),
                          estimate = as.numeric(),
                          std.error = as.numeric(),
                          statistic = as.numeric(),
                          p.value = as.numeric(),
                          s.value = as.numeric(),
                          conf.low = as.numeric(),
                          conf.high = as.numeric(),
                          df = as.integer(),
                          statistic.noninf = as.numeric(),
                          statistic.nonsup = as.numeric(),
                          p.value.noninf = as.numeric(),
                          p.value.nonsup = as.numeric(),
                          p.value.equiv = as.numeric(),
                          muscle = as.character())
  
  for(i in c("thigh", "arm")) {
    
    data_muscle <- data |>
      filter(muscle == i)
    
    std_lm <- lm(estimated_ma ~ cond_dummy,
                 data = data_muscle)
    
    std <- summary(std_lm)$sigma
    
    data_muscle <- data_muscle |>
      mutate(
        estimated_ma_z = (estimated_ma - mean(estimated_ma, na.rm=TRUE)) / std
      )
    
    model <- lmer(estimated_ma_z ~ time + time:cond_dummy + (1|site_id) + (1 | participant_id),
                  data = data_muscle, REML = TRUE)
    
    tests <- hypotheses(model, equivalence = c(-0.1,0.1), vcov = "satterthwaite", df = insight::get_df(model)) |>
      mutate(muscle = i)
    
    model_results <- bind_rows(model_results,
                               tests)
    
  }
  
  return(model_results)
  
}

  # Note, for some reason this function (or at least it's individual components) will work outside of targets but not inside
  # We get an error in qdrg() - 'data' must be a data.frame, environment, or list
  # So you will not see it in the pipeline, but the interval estimates have been manually added to the exploratory plot
  refit_MA_model <- function(data) {
    
    data <- as.data.frame(data)
    
    MultiLevelModel_SMD_hypertrophy_weeks_arm <- rma.mv(yi, V=vi, data=data,
                                                        random = list(~ log1p(weeks) | study, ~ 1 | arm, ~ 1 | es), mods = ~ log1p(weeks),
                                                        method="REML", test="t", struct = "GEN"
    )
    
    RobuEstMultiLevelModel_SMD_hypertrophy_weeks_arm <- robust(MultiLevelModel_SMD_hypertrophy_weeks_arm, data$study)
    
    ref_grid <- qdrg(object=RobuEstMultiLevelModel_SMD_hypertrophy_weeks_arm, data=data, at = list(weeks = c(24, 36)))
  
    emt <- update(contrast(ref_grid, "consec", simple = "weeks"), by = NULL)
  
    slope <- confint(emt)
  
    slope
  }

plot_explore_hypertrophy_models <- function(models) {
  
  models <- models |>
    mutate(
      muscle = case_when(
        muscle == "arm" ~ "Arm Muscle",
        muscle == "thigh" ~ "Thigh Muscle"
      )
    )
  
  time_plot <- models |>
    filter(term == "time") |>
    ggplot(aes(x=muscle, y=estimate)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_hline(yintercept = c(-0.1,0.1), linetype = 2) +
    annotate("rect", 
             xmin = -Inf, xmax = Inf, 
             ymin = 0.04176119, ymax = 0.05676972,
             alpha = 0.5) +
    annotate("text", x = 1.5, y = 0.0493,
             label = "Predicted Time Effect Interval Estimate from Meta-Regression",
             size = 2) +
    scale_x_discrete() +
    geom_pointinterval(aes(ymin=conf.low, ymax=conf.high)) +
    annotate("text", x = 1.5, y = c(-0.09,0.09),
             label = "Smallest Effect Size of Interest Limits",
             size = 2) +
    labs(
      y = "Standardised Mean Effect",
      title = "Main effect of time"
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  interact_plot <- models |>
    filter(term == "time:cond_dummy") |>
    ggplot(aes(x=muscle, y=estimate)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_hline(yintercept = c(-0.1,0.1), linetype = 2) +
    geom_pointinterval(aes(ymin=conf.low, ymax=conf.high)) +
    annotate("text", x = 1.5, y = c(-0.09,0.09),
             label = "Smallest Effect Size of Interest Limits",
             size = 2) +
    labs(
      y = "Standardised Mean Effect",
      title = "Condition by time interaction"
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  coef_plots <- (time_plot / interact_plot) +
    plot_annotation(title = "Exploratory Multilevel ANCOVA Models for Hypertrophy Outcomes",
                    caption = "Error bars are 95% confidence intervals") +
    plot_layout(axes = "collect") 
  
  wrap_elements(coef_plots) 
}
# Quick look at Wolf et al. data https://osf.io/q9djw 

# Misc -----
calculate_wolf_effect_sizes <- function() {
  
  wolf_data <- read_xlsx("data/wolf_data.xlsx") |>
    janitor::clean_names() |>
    slice(1:50)
  
  wolf_summary <- wolf_data |>
    group_by(group) |>
    summarise(
      bb55_premean = mean(bb55_pre),
      bb45_premean = mean(bb45_pre),
      tb55_premean = mean(tb55_pre),
      tb45_premean = mean(tb45_pre),
      bb55_postmean = mean(bb55_post),
      bb45_postmean = mean(bb45_post),
      tb55_postmean = mean(tb55_post),
      tb45_postmean = mean(tb45_post),
      bb55_presd = sd(bb55_pre),
      bb45_presd = sd(bb45_pre),
      tb55_presd = sd(tb55_pre),
      tb45_presd = sd(tb45_pre),
      bb55_postsd = sd(bb55_post),
      bb45_postsd = sd(bb45_post),
      tb55_postsd = sd(tb55_post),
      tb45_postsd = sd(tb45_post),
      bb55_cor = cor(bb55_pre, bb55_post),
      bb45_cor = cor(bb45_pre, bb45_post),
      tb55_cor = cor(tb55_pre, tb55_post),
      tb45_cor = cor(tb45_pre, tb45_post),
      n = n()
    ) |>
    pivot_longer(c(-group, -n), 
                 names_to = c("Var", ".value"), 
                 names_sep="_" )
  
  wolf_summary <-  escalc(measure = "SMCR",
                          m1i = postmean, m2i = premean,
                          sd1i = presd, 
                          ri = cor, ni = n,
                          data = wolf_summary
  )
  
  ppc_effects <- wolf_summary |>
    select(group, Var, yi, vi) |>
    pivot_wider(id_cols = "Var",
                names_from = "group",
                values_from = c("yi", "vi")) |>
    mutate(smd_ppc = yi_2 - yi_1, 
           smd_ppc_lower = smd_ppc - sqrt((vi_1 + vi_2)*1.96),
           smd_ppc_higher = smd_ppc + sqrt((vi_1 + vi_2)*1.96)
    )
}