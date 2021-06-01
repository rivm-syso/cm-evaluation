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
#################################################
#
# Helper functions for simulation of one partition of 
# the model building blocks
# 1. transmission chain
# 2. timed outbreak
# 3. transmission network
# 4. initialize controls (initialize adherence)
#
#################################################

source("R/functions/models.R")
library(glue)

#' Helper to check if simulation results are already present on disk
#' 
#' @param file_path path to file
#' Stops function with error message if file already exists 
check_results_on_disk <- function(file_path) {
  if (file.exists(file_path)) {
    stop("Simulation results already present in file {file_path}, 
          delete or rename if you want to rerun simulations with the same settings")
  }
}


#' Helper function to replace dots in doubles with dashes to use in filenames
#' @param num numeric object
#' @return numeric converted to string with dot replaced by dash
num_to_fname <- function(num) {
  str_replace(num, "\\.", "-")
}


#' Helper to simulate transmission chains for one partition
#' 
#' @param Reff double - effective reproduction number
#' @param genmax int - maximum number of generations in the outbreak
#' @param kpar double - dispersion parameter
#' @param repl int - partition id number
#' @param result_dir_path str - path to directory where results are written to
#' 1000 simulations of transmission chains are written away into one RDS file
make1000tcs <- function(Reff, genmax, kpar, repl, results_dir_path = RESULTS_DIR_PATH) {
  file_path <- glue("tc_R{num_to_fname(Reff)}_k{num_to_fname(kpar)}_rep{repl}")
  check_results_on_disk(file.path(results_dir_path, file_path))
  sims <- sapply(1:1000, function(x) {
    list(
      transchain(
        genmax = genmax,
        Reff = Reff,
        kpar = kpar,
        seedstart = x + repl * 1000
      )
    )
  })
  saveRDS(sims, file = file.path(results_dir_path, file_path))
}


#' Helper to simulate transmission networks for one partition
#' 
#' @param clust double - clustering coefficient
#' @param tc_file str - filename containing transmission chain results to simulate network for
#' @param results_dir_path str - directory to write results to
#' 1000 simulations of transmission networks for given transmission chains are 
#' written away into one RDS file
make1000cls <- function(clust, tc_file, results_dir_path = RESULTS_DIR_PATH) {
  str_tn_params <- glue("phi{num_to_fname(clust)}")
  result_fname <- str_replace(tc_file, "(tc)(.+)(_rep\\d)", glue("tn\\2_{str_tn_params}\\3"))
  check_results_on_disk(file.path(results_dir_path, result_fname))
  sims <- readRDS(file.path(results_dir_path, tc_file))
  for (i in 1:1000) {
    sims[[i]] <- transnetwork(tcdata = sims[[i]], prtriangle = clust)
  }
  saveRDS(sims, file = file.path(results_dir_path, result_fname))
}


#' Helper to simulate timing of outbreak events for one partition
#' 
#' @param corr_infectivity bool - whether infectivities are correlated
#' @param genTmean double
#' @param genTshape
#' @param prsymp 
#' @param tc_file str - filename containing transmission chain results to simulate event times for
#' @param results_dir_path - str - directory to write results to
#' 1000 simulations of timed outbreaks for given transmission chains are 
#' written away into one RDS file
make1000tos <- function(corr_infectivity, genTmean, genTshape, prsymp, tc_file, results_dir_path = RESULTS_DIR_PATH) {
  str_to_params <- glue("gTm{num_to_fname(genTmean)}s{num_to_fname(genTshape)}{c(\"U\", \"C\")[1 + corr_infectivity]}_symp{num_to_fname(prsymp)}")
  result_fname <- str_replace(tc_file, "(tc)(.+)(_rep\\d)", glue("to\\2_{str_to_params}\\3"))
  check_results_on_disk(file.path(results_dir_path, result_fname))
  sims <- readRDS(file.path(results_dir_path, tc_file))
  for (j in 1:1000) {
    sims[[j]] <- timedoutbreak(tcdata = sims[[j]],
      correlatedinfectivity = corr_infectivity,
      genTmean = genTmean,
      genTGashape = genTshape,
      probsymptomatic = prsymp
    )
  }
  saveRDS(sims, file = file.path(results_dir_path, result_fname))
}


#' Helper to simulate adherences for one partition
#' 
#' @param frac_nothing
#' @param frac_appuse
#' @param tc_file
#' @param results_dir_path
#' 
make1000ics <- function(frac_nothing, frac_appuse, tc_file, results_dir_path = RESULTS_DIR_PATH) {
  str_ic_params <- glue("IQ{num_to_fname(frac_nothing)}_app{num_to_fname(frac_appuse)}")
  result_fname <- str_replace(tc_file, "(tc)(.+)(_rep\\d)", glue("ic\\2_{str_ic_params}\\3"))
  check_results_on_disk(result_fname)
  sims <- readRDS(file.path(results_dir_path, tc_file))
  for (j in 1:1000) {
      sims[[j]] <- initializecontrol(tcdata = sims[[j]],
                                     pr_nothing = frac_nothing,
                                     prusetracingapp = frac_appuse
      )
  }
  saveRDS(sims, file = file.path(results_dir_path, result_fname))
}


#' Simulate controls on epidemic 1000 times - one partition and does bookkeeping analysis of number of
#' infectives under the control measures for each simulation round. Saves file with results.
#'
#' @param tnfile str - name of result file of the transmission network simulations without last part "_rep<n>"
#' @param tofile str - name of result file of the timed outbreak simulations without last part "_rep<n>"
#' @param icfile str - name of result file of the initialize controls simulations without last part "_rep<n>"
#' @param repl int - identifier of partition
#' @param prselftrace double - probability of self tracing given a link
#' @param prtrace double - probability of bco tracing given a link
#' @param prapptrace double - probability of app tracing given a link and both use the app
#' @param dselfiso double - delay in self isolation after showing symptoms
#' @param dselftracetest double - delay in contact (found in informal tracing) asking for test after index tests positive
#' @param dtracetest double - delay in contact (found through bco tracing) asking for test after index tests positive
#' @param dtestapp double - delay in contact (found through app tracing) asking for test after index tests positive
#' @param dtest double - delay in asking for test till positive result comes back for a symptomatic case
#' @param dir_path_outbreak str - path to directory that contains the results of model building blocks tn, to, ic
#' @return None - file with results is saved to disk
simallcontrol <- function(tnfile, tofile, icfile, repl,
  prselftrace, prtrace, prapptrace,
  dselfiso, dselftracetest, dtracetest, dtestapp, dtest, dir_path_outbreak = DIR_PATH_OUTBREAK) {
  # check if directory to write results to already exists and create new dir if not
  if (!dir.exists(dir_path_outbreak)) {
    stop("dir_path for base simulation files does not exist. Run the simulations or point to the correct directory")
  }
  # construct filenames and check if results exist
  file_path_tn <- file.path(dir_path_outbreak, glue("{tnfile}_rep{repl}"))
  file_path_to <- file.path(dir_path_outbreak, glue("{tofile}_rep{repl}"))
  file_path_ic <- file.path(dir_path_outbreak, glue("{icfile}_rep{repl}"))
  group1 <- "(tn|to|ic)" # first group "tn", "to" or "ic"
  group2 <- "_(R(.*?)_k(.*?))_" # second group "R<x1>_k<x2>"
  str_pattern <- paste0(group1, group2) # str_remove to obtain third group "remainder"
  result_fname1 <- glue("res{str_extract(tnfile, group2)}{str_remove(tnfile, str_pattern)}_{str_remove(tofile, str_pattern)}_{str_remove(icfile, str_pattern)}")
  result_fname2 <- glue("{result_fname1}_pt{num_to_fname(prtrace)}_pat{num_to_fname(prapptrace)}_pst{num_to_fname(prselftrace)}_dsi{num_to_fname(dselfiso)}_dstt{num_to_fname(dselftracetest)}_dtt{num_to_fname(dtracetest)}_dta{num_to_fname(dtestapp)}_dtest{num_to_fname(dtest)}_rep{repl}")
  result_fpath <- file.path(dir_path_outbreak, result_fname2)
  if (file.exists(result_fpath)) {
    stop("Simulation results for the controls already exist at {result_fpath}")
  }
  if (c(file_path_tn, file_path_to, file_path_ic) %>% sapply(X = ., FUN = file.exists) %>% all()) {
    log_info("Read in model building blocks to simulate controls on")
    log_info("tn-file = {file_path_tn}, to-file = {file_path_to}, ic-file = {file_path_ic}")
    simstn <- readRDS(file_path_tn)
    simsto <- readRDS(file_path_to)
    simsic <- readRDS(file_path_ic)
  } else {
    log_info(file_path_ic, file_path_to, file_path_tn)
    stop("Some result files are missing. Run the simulations for the epidemic outbreak")
  }
  # Start simulation of controls
  reslist <- list()
  # Simulate test + trace controls on the input epidemic outbreaks (1000 per partition)
  for (i in 1:1000) {
    inputlist <- c(simstn[[i]], simsto[[i]], simsic[[i]])
    # self isolation
    tn_iso <- applyselfisolation(
      toicdata = inputlist,
      testcontactdelay = dtest,
      isodelay = dselfiso
    )
    # self isolation & informal trace and test
    tn_selftracetest <- applyselftracingtesting(
      titndata = tn_iso,
      prtrace = prselftrace,
      tracedelay = dselftracetest,
      testcontactdelay = dtest
    )
    # self isolation & informal trace+test & bco trace+test
    tn_tracetest <- applytracingtesting(
      titndata = tn_iso,
      tracedelay = dtracetest,
      testcontactdelay = dtest,
      prtrace = prtrace,
      prselftrace = prselftrace,
      selftracedelay = dselftracetest
    )
    # self isolation & informal trace+test & bco trace+test & app trace+test
    tn_apptest <- applyapptesting(
      titndata = tn_iso,
      prtrace = prtrace,
      prapptrace = prapptrace,
      prselftrace = prselftrace,
      testcontactdelay = dtest,
      apptracedelay = dtestapp,
      tracedelay = dtracetest,
      selftracedelay = dselftracetest
    )
    # Number of cases directly prevented by app: in three steps
    # 1. collect the cases that have been traced by app and ask for test - 
    # note: last generation can't prevent any cases, discard them
    cases_test_approute <- which(tn_apptest$apptestrouteapp & inputlist$generations != GEN_MAX)
    # 2. collect the infectees of those cases - list of vectors
    infectees_test_approute <- sapply(cases_test_approute, function(infector) which(inputlist$infectors == infector))
    # 3. for each infectee check if he/she is a case without app control AND prevented with app control
    direct_prev_approute <- sapply(infectees_test_approute, function(x) !tn_tracetest$isoquartestprevented[x] & tn_apptest$appisoquartestprevented[x])
    n_direct_prev_approute <- direct_prev_approute %>% unlist %>% sum(na.rm = TRUE)
    # Number of cases traced through app, bco, informal
    n_test_approute <- sum(tn_apptest$apptestrouteapp & inputlist$generations != GEN_MAX)
    n_test_bcoroute <- sum(tn_apptest$apptestroutebco & inputlist$generations != GEN_MAX)
    n_test_selfroute <- sum(tn_apptest$apptestrouteself & inputlist$generations != GEN_MAX)
    
    # put results per outbreak from controls that get written away into ctrlres list
    ctrlres <- list(
      networkseed = inputlist$networkseed, # seed to associate controls to outbreak
      direct_prev_approute = n_direct_prev_approute, # number of cases directly prevented through app
      n_test_approute = n_test_approute, # number of cases traced+test through app under app control
      n_test_bcoroute = n_test_bcoroute, # number of cases traced+test through bco under app control
      n_test_selfroute = n_test_selfroute, # number of cases traced+test through informal under app control
      # number of cases per generation without any controls
      nocontrol = tabulate(inputlist$generations, nbins = GEN_MAX), 
      # number of cases per generation with self isolation
      testiso = tabulate(inputlist$generations[!tn_iso$isoprevented], nbins = GEN_MAX), 
      # number of cases per generation with informal bco
      selftracetest = tabulate(inputlist$generations[!tn_selftracetest$selfisoquartestprevented], nbins = GEN_MAX), 
      # number of cases per generation with bco
      tracetest = tabulate(inputlist$generations[!tn_tracetest$isoquartestprevented], nbins = GEN_MAX), 
      # number of cases per generation with app
      apptracetest = tabulate(inputlist$generations[!tn_apptest$appisoquartestprevented], nbins = GEN_MAX) 
    )
    # append to results
    reslist <- c(reslist, list(ctrlres))
  }
  # save all results of controls to file
  saveRDS(reslist, file = result_fpath)
}


#' Wrapper around simallcontrols to obtain all 10 partitioned simulation results 
#' of the controls on the epidemic for a given scenario 
#'
#' @param Reff string - coding of Reff in results filename, e.g. "R1-3" for R_eff=1.3
#' @param kpar string - coding of k_parr in results filename, e.g. "k0-1" for k_par=0.1
#' @param phi string = coding for phi in results filename, e.g. "phi0-2"
#' @param genT - string - coding for generation interval, e.g. "gTm4s5U"
#' @param symp - string - coding for fraction symptomatic, e.g. "symp0-6"
#' @param IQlevel string - coding for frac of population that doesn't get tested
#' when symptomatic and traced through bco, e.g. "IQ0-7"
#' @param appuse string - coding for frac of population that uses the app, e.g. "app0-16"
#' @param prtrace double - probability to trace through bco tracing given a link
#' @param prselftrace double - probability to trace through informal tracing given a link
#' @param prapptrace double - probability to trace through app tracing given a link and both have the app
#' @param dtest double - delay in asking for test and outcome for symptomatic
#' @param delaytt double - delay in bco tracing - test
#' @param delayta double - delay in app tracing - test
#' @param delaystt double - delay in informal tracing - test
#' @param dselfiso double - delay in getting testing due to symptoms
#' @param dir_path_outbreak str - directory path pointing to where simulations of outbreak are stored
#' @return None - results are saved in RDS format in given directory
analysescenario <- function(Reff, kpar, phi, IQlevel, prtrace, prselftrace, appuse,
  genT, symp, prapptrace, dtest, delaysi, delaystt, delaytt,
  delayta, dir_path_outbreak = DIR_PATH_OUTBREAK) {
  genfile <- paste(Reff, kpar, sep = "_") # encoding the transmission chain results for specified parameters
  tnfile <- paste("tn", genfile, phi, sep = "_") # transmission chain results for specified parameters (w/o partitioning)
  tofile <- paste("to", genfile, genT, symp, sep = "_") # timed outbreak results for specified parameters (w/o partitioning)
  icfile <- paste("ic", genfile, IQlevel, appuse, sep = "_") # initialize control results for specified parameters (w/o partitioning)
  
  log_info("Reading from {dir_path_outbreak} the files: tnfile = {tnfile}, tofile = {tofile}, icfile = {icfile}")
  log_info("Results will be saved to {dir_path_outbreak}")
  mclapply(0:9,
    function(x) {
      simallcontrol(tnfile, tofile, icfile, x,
        prtrace = prtrace,
        prselftrace = prselftrace,
        prapptrace = prapptrace,
        dselfiso = delaysi,
        dselftracetest = delaystt,
        dtracetest = delaytt,
        dtestapp = delayta,
        dtest = dtest,
        dir_path_outbreak = dir_path_outbreak
      )
    },
    mc.cores = num_cores
  )
}
