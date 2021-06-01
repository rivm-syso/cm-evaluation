# Copyright 2021 Rijksinstituut voor Volksgezondheid en Milieu (RIVM).
#
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU Affero General Public License as published by the Free 
# Software Foundation, either version 3 of the License, or any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program.  If not, see <https://www.gnu.org/licenses/>.”
#
###############################################################################
### Make all simulations for the initalization of the model ###
# - transmission chains
# - clustering
# - times
# - behaviours
#
# parallelization using forking:
################################################################################

# load model functions
source("R/functions/simulate_utils.R")
source("R/functions/models.R")
source("R/functions/models_utils.R")

# load libraries
library(parallel)
library(logger)
log_layout(layout_glue_colors)
library(tidyverse)

# Epidemic outbreak output directory
RESULTS_DIR_PATH <- "results"

# check that directory results exists on disk, otherwise create it
if (!dir.exists(RESULTS_DIR_PATH)) {
  dir.create(RESULTS_DIR_PATH)
}

# number of cores
num_cores <- 16


############ First: transmission chains #######################################
# parameter values to iterate over note genmaxs is a function of R_eff for
# computational reasons (high R_eff blows up epi over number of generations)
genmaxs <- c(10, 20)
Reffs <- c(1.3, 1.05)
kpars <- c(0.5, 0.05, 0.1)

stopifnot(length(genmaxs) == length(Reffs))

num_sim <- length(genmaxs)

for (i in 1:num_sim) {
  for (kpar in kpars) {
    log_info("Apply function make1000tcs with R_eff = {Reffs[i]}, genmax = {genmaxs[i]}, kpar = {kpar} and save to {RESULTS_DIR_PATH}")
    mclapply(0:9, function(x) {make1000tcs(Reffs[i], genmaxs[i], kpar, x)}, mc.cores = num_cores)
  }
}


############################ Second: clustering ################################
# list all transmission chain files to simulate the network for
tc_files <- list.files(path = RESULTS_DIR_PATH, pattern = "^tc")
# clustering coefficients to simulate
clustering <- 0.2
mclapply(tc_files, function(x) make1000cls(clustering, x), mc.cores = num_cores)

# for other clustering coefficients, only simulate on specific files, otherwise takes really long
tc_files <- list.files(path = RESULTS_DIR_PATH, pattern = "^tc_R1-3_k0-1")
# clustering coefficients to simulate
clustering <- 0.3
mclapply(tc_files, function(x) make1000cls(clustering, x), mc.cores = num_cores)
clustering <- 0
mclapply(tc_files, function(x) make1000cls(clustering, x), mc.cores = num_cores)


##################### Third: times #############################################
# set parameters for timed outbreak: each index represents one set of parameters for a timed outbreak
correlateds <- c(F, F, F, F)
genTshapes <- c(4, 4, 4, 4)
genTmeans <- c(4, 5, 4, 5)
prsymp <- c(PR_SYMPTOMATIC, PR_SYMPTOMATIC, 0.3, 0.3)

stopifnot( # all vectors should be the same length
  length(correlateds) == length(genTshapes),
  length(genTshapes) == length(genTmeans),
  length(genTmeans) == length(prsymp)
  )
num_sims <- length(correlateds)

tc_files <- list.files(path = RESULTS_DIR_PATH, pattern = "^tc")
for (i in 1:num_sims) {
  log_info("Simulate timed outbreak: set of parameters {i} of {num_sims}")
  mclapply(tc_files, function(x) make1000tos(correlateds[i], genTmeans[i], genTshapes[i], prsymp[i], x), mc.cores = num_cores)
}


##################### Fourth: behaviours #######################################
nothings <- c(FRAC_NOTHING, FRAC_NOTHING, FRAC_NOTHING, FRAC_NOTHING, FRAC_NOTHING, FRAC_NOTHING, 0.05, 0.15, 0.05, 0.15)
use_tracing_app <- c(PR_USE_TRACINGAPP, 0.2, 0.3, 0.4, 0.6, 0.8, 0.24, 0.08, PR_USE_TRACINGAPP, PR_USE_TRACINGAPP)

stopifnot(length(nothings) == length(use_tracing_app))

num_models <- length(nothings)

tc_files <- list.files(path = RESULTS_DIR_PATH, pattern = "^tc")

for (i in 1:num_models) { 
  log_info("Simulate behaviours: {i} out of {num_models}")
  mclapply(tc_files, function(x) make1000ics(nothings[i], use_tracing_app[i], x), mc.cores = num_cores)
}

