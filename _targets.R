# _targets.R file
library(targets)
library(tarchetypes)
library(crew)
source("R/functions.R")
tar_option_set(
  controller = crew_controller_local(workers = 2),
  packages = c(
    "here",
    # "readxl",
    # "janitor",
    "tidyverse",
    "base",
    "metafor",
    "lme4",
    "lmerTest",
    "faux",
    "furrr",
    # "scales",
    # "ggtext",
    "ggh4x",
    # "zoo",
    # "performance",
    # "see",
    # "rstan",
    # "brms",
    # "bayesplot",
    "marginaleffects",
    "broom.mixed"
    # "patchwork",
    # "kableExtra",
    # "knitr",
    # "quarto",
    # "officer",
    # "officedown",
  )
)


list(
  # Load and prepare data from previous meta-analysis
  tar_target(MA_data_file, here("data","Polito et al. RT Extracted Data.csv"), format = "file"),
  tar_target(MA_data, get_prepare_MA_data(MA_data_file)),
  
  # Fit and plot log growth MA model for hypertrophy
  tar_target(MA_model, log_growth_MA_model(MA_data)),
  tar_target(MA_plot, plot_log_growth_MA_model(MA_model,MA_data)),
  tar_target(
    MA_plot_png,
    ggsave(
      plot = MA_plot,
      filename = "sample_estimates/MA_plot.png",
      device = "png",
      dpi = 300,
      w = 10,
      h = 5
    )
  ),
  
  
  # Load in data for sample simulations
  tar_target(prior_studies_data, here("data","prior_studies_data.csv"), format = "file"),
  tar_target(prior_ma_reliability_data, here("data","prior_ma_reliability.csv"), format = "file"),
  
  # Run sample size simulations for power
  tar_target(sim_constant_effect, sample_size_simulation(prior_studies_data,prior_ma_reliability_data)),
  tar_target(sim_var_effect, sample_size_simulation_var_effect(prior_studies_data,prior_ma_reliability_data)),
  tar_target(sim_est_numbers, sample_size_simulation_est_numbers(prior_studies_data,prior_ma_reliability_data)),
  
  # Plot sample size simulations
  tar_target(plot_sim_constant_effect, plot_sample_size_estimates(sim_constant_effect)),
  tar_target(
    plot_sim_constant_effect_png,
    ggsave(
      plot = plot_sim_constant_effect,
      filename = "sample_estimates/plot_sim_constant_effect.png",
      device = "png",
      dpi = 300,
      w = 10,
      h = 5
    )
  ),
  tar_target(plot_sim_var_effect, plot_sample_size_estimates_var_effect(sim_var_effect)),
  tar_target(
    plot_sim_var_effect_png,
    ggsave(
      plot = plot_sim_var_effect,
      filename = "sample_estimates/plot_sim_var_effect.png",
      device = "png",
      dpi = 300,
      w = 10,
      h = 10
    )
  ),
  tar_target(plot_sim_est_numbers, plot_sample_size_estimates_est_numbers(sim_est_numbers)),
  tar_target(
    plot_sim_est_numbers_png,
    ggsave(
      plot = plot_sim_est_numbers,
      filename = "sample_estimates/plot_sim_est_numbers.png",
      device = "png",
      dpi = 300,
      w = 10,
      h = 10
    )
  )
  
)
