# _targets.R file
library(targets)
library(tarchetypes)
library(crew)
source("R/functions.R")
tar_option_set(
  controller = crew_controller_local(workers = 2),
  packages = c(
    "here",
    "readxl",
    "janitor",
    "tidyverse",
    "base",
    "metafor",
    "lme4",
    "lmerTest",
    "faux",
    "furrr",
    "ggh4x",
    "marginaleffects",
    "emmeans",
    "broom.mixed",
    "ggdist",
    "glue",
    "patchwork",
    "grateful",
    "osfr"
    # "kableExtra",
    # "knitr",
    # "quarto",
    # "officer",
    # "officedown",
  )
)

list(
  #####
  # Pre-regsitration targets
  # Load and prepare data from previous meta-analysis
  tar_target(
    MA_data_file,
    here("data", "Polito et al. RT Extracted Data.csv"),
    format = "file"
  ),
  tar_target(MA_data, get_prepare_MA_data(MA_data_file)),
  
  # Fit and plot log growth MA model for hypertrophy
  tar_target(MA_model, log_growth_MA_model(MA_data)),
  tar_target(MA_plot, plot_log_growth_MA_model(MA_model, MA_data)),
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
  tar_target(
    prior_studies_data,
    here("data", "prior_studies_data.csv"),
    format = "file"
  ),
  tar_target(
    prior_ma_reliability_data,
    here("data", "prior_ma_reliability.csv"),
    format = "file"
  ),
  
  # Run sample size simulations for power
  tar_target(
    sim_constant_effect,
    sample_size_simulation(prior_studies_data, prior_ma_reliability_data)
  ),
  tar_target(
    sim_var_effect,
    sample_size_simulation_var_effect(prior_studies_data, prior_ma_reliability_data)
  ),
  tar_target(
    sim_est_numbers,
    sample_size_simulation_est_numbers(prior_studies_data, prior_ma_reliability_data)
  ),
  
  # Plot sample size simulations
  tar_target(
    plot_sim_constant_effect,
    plot_sample_size_estimates(sim_constant_effect)
  ),
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
  tar_target(
    plot_sim_var_effect,
    plot_sample_size_estimates_var_effect(sim_var_effect)
  ),
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
  tar_target(
    plot_sim_est_numbers,
    plot_sample_size_estimates_est_numbers(sim_est_numbers)
  ),
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
  ),
  
  #####
  # Study results targets
  
  # Load in data
  tar_target(
    file_data_hypertrophy,
    here("data", "PD_data_hypertrophy.csv"),
    format = "file"
  ),
  tar_target(data_hypertrophy, read_csv(file_data_hypertrophy)),
  
  tar_target(
    file_data_strength,
    here("data", "PD_data_strength.csv"),
    format = "file"
  ),
  tar_target(data_strength, read_csv(file_data_strength)),
  
  # Primary pre-registered hypertrophy outcomes
  tar_target(
    data_hypertrophy_ready,
    prepare_hypertrophy_data(data_hypertrophy)
  ),
  tar_target(
    results_hypertrophy,
    fit_hypertrophy_models(data_hypertrophy_ready)
  ),
  tar_target(
    plots_hypertrophy,
    plot_hypertrophy_models(data_hypertrophy_ready, results_hypertrophy)
  ),
  tar_target(
    plots_hypertrophy_tiff,
    ggsave(
      plot = plots_hypertrophy,
      filename = "plots/hypertrophy_plot.tiff",
      device = "tiff",
      dpi = 300,
      width = 6,
      height = 9
    )
  ),
  
  # Strength outcomes
  tar_target(data_strength_ready, prepare_strength_data(data_strength)),
  tar_target(results_strength, fit_strength_models(data_strength_ready)),
  
    # Fit and plot log growth MA model for hypertrophy
    tar_target(MA_data_strength, get_prepare_MA_data_strength(MA_data_file)),
    tar_target(MA_model_strength, log_growth_MA_model_strength(MA_data_strength)),
    tar_target(MA_plot_strength, plot_log_growth_MA_model_strength(MA_model_strength, MA_data_strength)),
    tar_target(
      MA_plot_strength_png,
      ggsave(
        plot = MA_plot_strength,
        filename = "plots/MA_plot_strength.png",
        device = "png",
        dpi = 300,
        w = 10,
        h = 5
      )
    ),
  
  tar_target(
    plots_strength,
    plot_strength_models(data_strength_ready, results_strength)
  ),
  tar_target(
    plots_strength_tiff,
    ggsave(
      plot = plots_strength,
      filename = "plots/strength_plot.tiff",
      device = "tiff",
      dpi = 300,
      width = 6,
      height = 9
    )
  ),
  
  # Exploratory hypertrophy outcomes
  tar_target(
    data_explore_hypertrophy_ready,
    prepare_explore_hypertrophy_data(data_hypertrophy)
  ),
  tar_target(
    results_explore_hypertrophy,
    fit_explore_hypertrophy_models(data_explore_hypertrophy_ready)
  ),
  
  # Note, for some reason this function will work outside of targets but not inside
  # So you will not see it in the pipeline, but the interval estimates have been manually added to the exploratory plot
  # tar_target(
  #   MA_slope,
  #   refit_MA_model(MA_data)
  # ),
  
  tar_target(
    plots_explore_hypertrophy,
    plot_explore_hypertrophy_models(results_explore_hypertrophy)
  ),
  tar_target(
    plots_explore_hypertrophy_tiff,
    ggsave(
      plot = plots_explore_hypertrophy,
      filename = "plots/explore_hypertrophy_plot.tiff",
      device = "tiff",
      dpi = 300,
      width = 6,
      height = 4.5
    )
  ),
  
  # Manuscript
  tar_target(grateful_report, cite_packages(out.dir = "pre_print", cite.tidyverse = TRUE, out.format = "pdf")),
  tar_target(wolf_effect_sizes, calculate_wolf_effect_sizes())
  
)
