##### Notes from Luke, James and I meeting about power, alpha, SESOIs, site recruitment etc.

# test for equivalence
# SESOI at 0.1?
# number of locations - 15 
# ballpark - 1/2 25-35, and 1/2 5-15 participants
# duration of recruitment - 4-6 weeks of recruitment
# dates ?? 


# Functions for targets

sample_size_simulation<- function(file1, file2) {

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
      mutate(muscle_site = "arm")
    
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
    summarise(power = mean(p.value < .025),
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
      caption = "Simulated to detect a standardised effect, or non-superiority, for condition:time with a SESOI of 0.1\nSite number k fixed at 15, \u03b1 = 0.05 (corrected to 0.025 for two primary outcomes i.e., arm and thigh muscle area)\n1000 simulations per combination of participants, number of measurements, and timepoints"
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
      mutate(muscle_site = "arm")
    
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
    summarise(power = mean(p.value < .005), 
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
      caption = "Simulated to detect a standardised effect, or non-superiority, for condition:time with a SESOI of 0.1\nSite number k fixed at 15, \u03b1 = 0.05 (corrected to 0.025 for two primary outcomes i.e., arm and thigh muscle area)\n1000 simulations per combination of participants, number of measurements, and timepoints"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  plot_sim
}
