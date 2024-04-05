install.packages("ggdag")

y = (b0 + u_site + u_participant) + (b_time * time_cont) + (b_cond * cond_dummy) + (b_cond_time * cond_dummy * time_cont) + sigma

library(tidyverse)
library(ggdag)

dagified <- dagify(
  y_1 ~ u_site + u_participant + sigma, 
  y_2 ~ y_1 + time_cont +  cond_dummy + cond_dummy_time_cont + sigma,
  exposure = "cond_dummy",
  outcome = "y_2"
)

bigger_dag <- dag_paths(dagified)

ggdag_parents(bigger_dag, "cond_dummy") +
  theme_dag_blank()
