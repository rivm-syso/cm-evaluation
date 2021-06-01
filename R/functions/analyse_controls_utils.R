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

################################################################################
#
# Util functions for analysing scenario results
#
# Read scenario result
# Summarise scenario result
#
################################################################################


#' Helper function to read all result files for given parameter values
readscenarioresults <- function(Reff, kpar, phi, symp, IQlevel, prtrace, prselftrace, appuse, genT,
                                prapptrace, delaysi, delaystt, delaytt, delayta, dtest, dir_output_path = DIR_PATH_OUTBREAK) {
  # reconstruct filenames of results
  result_fname1 <- paste("res", Reff, kpar, phi, genT, symp, IQlevel, appuse, sep = "_")
  genfile <- glue("{result_fname1}_pt{num_to_fname(prtrace)}_pat{num_to_fname(prapptrace)}_pst{num_to_fname(prselftrace)}_dsi{num_to_fname(delaysi)}_dstt{num_to_fname(delaystt)}_dtt{num_to_fname(delaytt)}_dta{num_to_fname(delayta)}_dtest{num_to_fname(dtest)}_rep")
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


#' Put summary of results into a tibble
#' 
#' @param summarytibble tibble - name of tibble to which summary of results should be appended
#' @param scenario_list_name str - name of list containing the scenario
#' @param scenario_description str - human readable name of scenario
summarise_scenarioresult <- function(scenario_description, scenario_list_name, summarytibble = tibble()) {
  # read in results
  res <- do.call(readscenarioresults, get(scenario_list_name))
  # calculate reproduction numbers and numbers tested through different routes and directly prevented
  # set number of generations
  Reff <- get(scenario_list_name)$Reff
  Rlength <- case_when(Reff == "R1-3" ~ 9, Reff == "R1-05" ~ 19)
  num_gen <- Rlength - 2
  
  summarytibble <<-
    bind_rows(
      summarytibble,
      tibble(
        scenario = scenario_description,
        R0 = do.call(`/`, as.list(rowSums(sapply(res, function(x) x$nocontrol))[c(Rlength, 2)]))^(1 / num_gen), 
        Riso = do.call(`/`, as.list(rowSums(sapply(res, function(x) x$testiso))[c(Rlength, 2)]))^(1 / num_gen),
        RSelfT = do.call(`/`, as.list(rowSums(sapply(res, function(x) x$selftracetest))[c(Rlength, 2)]))^(1 / num_gen),
        RBCOT = do.call(`/`, as.list(rowSums(sapply(res, function(x) x$tracetest))[c(Rlength, 2)]))^(1 / num_gen),
        RAppT = do.call(`/`, as.list(rowSums(sapply(res, function(x) x$apptracetest))[c(Rlength, 2)]))^(1 / num_gen),
        n_test_approute = sapply(res, function(x) x$n_test_approute) %>% unlist %>% sum(na.rm = TRUE),
        n_test_bcoroute = sapply(res, function(x) x$n_test_bcoroute) %>% unlist %>% sum(na.rm = TRUE),
        n_test_selfroute = sapply(res, function(x) x$n_test_selfroute) %>% unlist %>% sum(na.rm = TRUE),
        test_approute_prevented = sapply(res, function(x) x$direct_prev_approute) %>% unlist %>% sum(na.rm = TRUE)
      )
    )
}
