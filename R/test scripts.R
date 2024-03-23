library(tidyverse)
library(lme4)
library(lmerTest)
library(faux)
library(broom.mixed)
library(furrr)

# Using prior data from all Discover Strength studies baseline lean mass (standardised) to get estimates of between study (i.e., site) variance

prior_data <- read_csv("data/prior_studies_data.csv") |>
  mutate(
    lean_mass_z = (lean_mass - mean(lean_mass, na.rm=TRUE))/sd(lean_mass, na.rm=TRUE)
    ) 

model_site_variance <- lmer(lean_mass_z ~ 1 + (1 | study),
                         data = prior_data)

site_variance <- VarCorr(model_site_variance)

# # Using prior data from body composition challenge lean mass (standardised) to get estimates of between participant variance
# 
# body_comp_challenge_data <- read.csv("data/body_comp_challenge_data.csv")
# 
# body_comp_challenge_data <- body_comp_challenge_data %>%
#   pivot_longer(9:44,
#                names_to = "names",
#                values_to = "values") %>%
#   separate(names, c("condition","time", "outcome", "variable"), "_") %>%
#   unite("outcome", outcome:variable) |>
#   filter(outcome == "bodpod_leanmass") |>
#   mutate(
#     lean_mass_z = (values - mean(values, na.rm=TRUE))/sd(values, na.rm=TRUE)
#   ) 
# 
# model_participant_variance <- lmer(lean_mass_z ~ 1 + time*condition + (1 | participant),
#                             data = body_comp_challenge_data)
# 
# participant_variance <- VarCorr(model_participant_variance)


# Using prior study data to get estimates of measurement error for muscle areas (standardised SEMs)

prior_ma_reliability <- read_csv("data/prior_ma_reliability.csv") |>
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

# Simulation parameters
sites <- 15
participants <- 50
measurements <- 3
b0 <- 0
b_time <- 0.15 
b_cond <- 0
b_cond_time <- 0.1 
u_site <- sqrt(site_variance$study)[1] 
u_participant <- (1 - sqrt(site_variance$study)[1] - arm_ma_reliability)
sigma <- arm_ma_reliability
m_error <- arm_ma_reliability


dat <- add_random(site = sites) |>
  add_random(participant = 10, .nested_in = "site") |>
  add_between("participant", cond = c("full_ROM", "lengthened_partial")) |>
  add_recode("cond", "cond_dummy", full_ROM = 0, lengthened_partial = 1) |>
  add_within("participant", time = c("pre", "mid", "post")) |>
  add_recode("time", "time_cont", pre = 0, mid = 0.5, post = 1) |>
  add_ranef("site", u_site = u_site) |>
  add_ranef("participant", u_participant = u_participant) |>
  add_ranef(sigma = sigma) |>
  mutate(m_error = m_error) |>
  mutate(arm_ma = (b0 + u_site + u_participant) + (b_time * time_cont) + (b_cond * cond_dummy) + (b_cond_time * cond_dummy * time_cont) + sigma) 

dat |>
  filter(time != "mid") |>
  pivot_wider(id_cols = c("site", "participant", "cond"), names_from = "time",
              values_from = "arm_ma") |>
  # ggplot(aes(x=pre,y=post)) +
  # geom_point() 
  summarise(r = cor(pre,post))

dat |>
  ggplot(aes(x = time_cont, y = arm_ma, color = cond)) +
  # geom_point(position = position_jitter(w=0.01), alpha = 0.25) +
  # geom_smooth(aes(group = participant), se = FALSE, method = "lm", alpha = 0.5, size = 0.5) +
  geom_smooth(method = "lm")

model <- lmer(arm_ma ~ time_cont * cond_dummy + (1 | site) + (1 | participant),
     data = dat)

summary(model)

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
    # add_recode("time", "time_cont", pre = 0, mid = 0.5, post = 1) |>
    add_ranef("site", u_site = u_site) |>
    add_ranef("participant", u_participant_arm = u_participant_arm) |>
    add_ranef("participant", u_participant_thigh = u_participant_thigh) |>
    add_ranef(arm_ma_error = arm_ma_error) |>
    add_ranef(thigh_ma_error = thigh_ma_error) |>
    mutate(arm_ma = (b0 + u_site + u_participant_arm) + (b_time * time) + (b_cond * cond_dummy) + (b_cond_time * cond_dummy * time) + arm_ma_error,
           thigh_ma = (b0 + u_site + u_participant_thigh) + (b_time * time) + (b_cond * cond_dummy) + (b_cond_time * cond_dummy * time) + thigh_ma_error) 
  
  
  # run mixed effect model and return relevant values
  model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (1 | participant),
            data = dat, REML = TRUE)
  
  model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (1 | participant),
                    data = dat, REML = TRUE)
  
  tidy_model_arm <- broom.mixed::tidy(model_arm) |>
    mutate(muscle_site = "arm")
  
  
  tidy_model_thigh <- broom.mixed::tidy(model_thigh) |>
    mutate(muscle_site = "thigh")
  
  bind_rows(tidy_model_arm, tidy_model_thigh)
}


sim_0_2 <- crossing(
  rep = 1:1000, # number of replicates
  participant_n = seq(2,20, by = 2), # range of participant N
  site_n = 15, # fixed site N
  time_n = c(0.5,1), # to examine inclusion of midpoint test or not
  measurements_n = c(1:3), # range of measurements N
) |>
  mutate(
  b0 = 0, b_time = 0.3, b_cond = 0, b_cond_time = 0.2,         # fixed effects 
  u_site = sqrt(site_variance$study)[1],  # random intercept site
  u_participant_arm = (1 - sqrt(site_variance$study)[1] - arm_ma_reliability),   # random intercept participant arm
  u_participant_thigh = (1 - sqrt(site_variance$study)[1] - thigh_ma_reliability),  # random intercetp participant thigh
  arm_ma_error = arm_ma_reliability,           # measurement error arm
  thigh_ma_error = thigh_ma_reliability         # measurement error thigh  
  ) %>%
  mutate(analysis = pmap(., sim)) %>% # not sure why base pipe doesn't work here
  unnest(analysis)


plot_0_2 <- sim_0_2 |>
  filter(effect == "fixed", term == "time:cond_dummy") |>
  mutate(
    time_n = case_when(
      time_n == 0.5 ~ "Include Midpoint",
      time_n == 1 ~ "No Midpoint",
    )
  ) |>
  mutate(time_n = factor(time_n, levels = c("No Midpoint", "Include Midpoint"))) |>
  group_by(measurements_n, participant_n, time_n) |> 
  summarise(power = mean(p.value < .025), 
            .groups = "drop") |>
  ggplot(aes(x = participant_n, y = power, color = factor(measurements_n))) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point() +
  geom_smooth(se = FALSE) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
  facet_wrap("time_n") +
  labs(
    x = "Total Number of Participants per Site (N)",
    y = "Power (1-\u03b2)",
    color = "Number of Measurements per Timepoint (N)",
    title = "Sample Size Estimation (N per Site)",
    subtitle = "Simulated to detect a standardised effect for condition:time of 0.1\nSite number fixed at 15, \u03b1 = 0.05 (corrected to 0.025 for two primary outcomes)",
    caption = "1000 simulations per combination of participants, number of measurements, and timepoints"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

plot_0_2

plan(cluster, workers = 10)

tic()
sim_0_1 <- crossing(
  rep = 1:5, # number of replicates
  participant_n = seq(2,20, by = 2), # range of participant N
  site_n = 15, # fixed site N
  time_n = c(0.5,1), # to examine inclusion of midpoint test or not
  measurements_n = c(1:3), # range of measurements N
) |>
  mutate(
    b0 = 0, b_time = 0.15, b_cond = 0, b_cond_time = 0.1,         # fixed effects 
    u_site = sqrt(site_variance$study)[1],  # random intercept site
    u_participant_arm = (1 - sqrt(site_variance$study)[1] - arm_ma_reliability),   # random intercept participant arm
    u_participant_thigh = (1 - sqrt(site_variance$study)[1] - thigh_ma_reliability),  # random intercetp participant thigh
    arm_ma_error = arm_ma_reliability,           # measurement error arm
    thigh_ma_error = thigh_ma_reliability         # measurement error thigh  
  ) %>%
  mutate(analysis = future_pmap(., sim)) %>% # not sure why base pipe doesn't work here
  unnest(analysis)
toc()

sim_0_1 <- crossing(
  rep = 1:5, # number of replicates
  participant_n = seq(2,20, by = 2), # range of participant N
  site_n = 15, # fixed site N
  time_n = c(0.5,1), # to examine inclusion of midpoint test or not
  measurements_n = c(1:3), # range of measurements N
) |>
  mutate(
    b0 = 0, b_time = 0.15, b_cond = 0, b_cond_time = 0.1,         # fixed effects 
    u_site = sqrt(site_variance$study)[1],  # random intercept site
    u_participant_arm = (1 - sqrt(site_variance$study)[1] - arm_ma_reliability),   # random intercept participant arm
    u_participant_thigh = (1 - sqrt(site_variance$study)[1] - thigh_ma_reliability),  # random intercetp participant thigh
    arm_ma_error = arm_ma_reliability,           # measurement error arm
    thigh_ma_error = thigh_ma_reliability         # measurement error thigh  
  ) %>%
  mutate(analysis = future_pmap(., sim)) %>% # not sure why base pipe doesn't work here
  unnest(analysis)
toc()

plot_0_1 <- sim_0_1 |>
  filter(effect == "fixed", term == "time:cond_dummy") |>
  mutate(
    time_n = case_when(
      time_n == 0.5 ~ "Include Midpoint",
      time_n == 1 ~ "No Midpoint",
    )
  ) |>
  mutate(time_n = factor(time_n, levels = c("No Midpoint", "Include Midpoint"))) |>
  group_by(measurements_n, participant_n, time_n) |> 
  summarise(power = mean(p.value < .025), 
            .groups = "drop") |>
  ggplot(aes(x = participant_n, y = power, color = factor(measurements_n))) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point() +
  geom_smooth(se = FALSE) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
  facet_wrap("time_n") +
  labs(
    x = "Total Number of Participants per Site (N)",
    y = "Power (1-\u03b2)",
    color = "Number of Measurements per Timepoint (N)",
    title = "Sample Size Estimation (N per Site)",
    subtitle = "Simulated to detect a standardised effect for condition:time of 0.1\nSite number fixed at 15, \u03b1 = 0.05 (corrected to 0.025 for two primary outcomes)",
    caption = "1000 simulations per combination of participants, number of measurements, and timepoints"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

plot_0_1
