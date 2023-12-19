# helper function to read in all the results for the timed outbreak
readtimedoutbreak <- function(Reff, kpar, phi, symp, IQlevel, prtrace, prselftrace, appuse, genT,
                              prapptrace, delaysi, delaystt, delaytt, delayta, dtest, dir_output_path = DIR_PATH_OUTBREAK) {
  # reconstruct filenames of results
  result_fname1 <- paste("to", Reff, kpar, genT, symp, sep = "_")
  genfile <- glue("{result_fname1}_rep")
  genfile <- file.path(dir_output_path, genfile)
  log_info("Scenario result file: {genfile}")

  # read in the results
  allresults <- list()
  for (i in 0:9) {
    fname <- paste0(genfile, i)
    if (file.exists(fname)) {
      allresults <- c(
        allresults,
        readRDS(fname)
      )
    }
  }
  return(allresults)
}

# path where result files are stored
DIR_PATH_OUTBREAK <- "results"

# source files to get parameter values
source("R/functions/simulate_utils.R")
source("R/functions/analyse_controls_utils.R")
source("R/functions/models.R")

# default scenario parameter values
scenario_default <- list(
  Reff = "R1-3",
  kpar = "k0-1",
  phi = "phi0-2",
  IQlevel = "IQ0-1",
  prtrace = PR_TRACE,
  prselftrace = PR_SELF_TRACE,
  appuse = "app0-16",
  genT = "gTm4s4U",
  symp = "symp0-7",
  prapptrace = PR_APP_TRACE,
  dtest = TEST_CONTACT_DELAY,
  delaytt = TRACE_DELAY,
  delayta = APP_TRACE_DELAY,
  delaystt = SELF_TRACE_DELAY,
  delaysi = ISO_DELAY
)
# name of scenario
scenario_list_name <- "scenario_default"

# read result files for timed outbreak of default scenario
res <- do.call(readtimedoutbreak, get(scenario_list_name))

# gather infection times from the result files
df <- purrr::map(1:length(res), \(i) {
  df_tmp <- tibble(
    generation = res[[i]]$generations,
    infectiontime = res[[i]]$infectiontimes,
  ) |>
    mutate(sim = i) # simulation
}) |>
  bind_rows() |>
  mutate( # day of infection
    inf_day = infectiontime |> floor()
  )

# prepare data for plotting
df_plt <- df |>
  filter(generation > 0) |> # remove index cases from generation 0
  group_by(generation, inf_day) |>
  summarise(cases = n())

# plot course of infection for largest outbreak
p <- df_plt |>
  ggplot(aes(x = inf_day, y = cases, col = generation |> as.factor())) +
  geom_line() +
  labs(x = "time", y = "Number of cases", col = "generation") +
  theme_minimal()

# save resulting plot
ggsave("figures/supplement_cases_over_time_per_generation.png",
  plot = p,
  bg = "white", width = 20, height = 20, units = "cm", dpi = 1200
)
