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
######################################################
#
# Model building blocks
#
# 1. Contact network and transmission timings
#
#     - transchain: transmission tree (inital outbreak: branching process)
#     - transnetwork: add links between individuals in a (transmission) tree
#     - timedoutbreak: add time of events (infection times and symptom onset
#         times) to a transmission tree
#     - initializecontrol: individual compliance to control measures
#
# 2. Control measures
#
#     - applyselfisolation: self isolation due to symptoms (pre and post positive test result)
#     - applyselftracingtesting: testing due to tracing
#     - applytracingtesting: testing due to tracing
#     - applyapptesting: testing due to app notification
#
######################################################

# load scripts and libraries
source("R/functions/models_utils.R")
source("R/functions/models_controls_utils.R")

library(logger)
log_layout(layout_glue_colors)
library(tidyverse)


# GLOBAL PARAMETERS - defaults in function calls
SEED_START <- 1 # make sure simulations are reproducible by setting seed in simulations
# assumptions
GEN_MAX <- 10 # number of generations to simulate
R_EFF <- 1.3 # reproduction number
K_PARR <- 0.1 # dispersion coefficient
PR_TRIANGLE <- 0.2 # clustering coefficient
CONNECTION_MAX <- 3 # maximum number of connections in creating triangles
MEAN_INC_TIME <- 5 # mean incubation period
PR_SYMPTOMATIC <- 0.7 # probability that case is symptomatic
CORR_INFECTIVITY <- FALSE # correlation between infectivitiy
QUARANTINE_DURATION <- 10000 # number of days in quarantine: set high to make this parameter irrelevant

# tracing probabilities
PR_USE_TRACINGAPP <- 0.16 # fraction of the population with the app
PR_TRACE <- 0.40 # probability to trace the link with bco
PR_SELF_TRACE <- 0.32 # ASSUMPTION probability to self trace the link, lower than bco
PR_APP_TRACE <- 0.75 # probability to trace the link with the app (if both would have app)

# delays (in days)
ISO_DELAY <- 1.6 # delay from sympton onset to asking for test
TEST_CONTACT_DELAY <- 1.3 # delay from asking for test to positive test outcome with symptoms
TRACE_DELAY <- 1.5 # delay in tracing contacts and them asking for test via bco after positive test
ASYM_TEST_DELAY <- 2.3 # delay from asking for test to positive test outcome for asymptomatic case (has to wait number of days before it can test)

APP_TRACE_DELAY <- 1 # delay in tracing contacts and them asking for test via app after positive test
SELF_TRACE_DELAY <- 0.17 # ASSUMPTION: delay in tracing contacts and them asking for test via self after positive test: 4hours

# behaviours
FRAC_NOTHING <- 0.1 # ASSUMPTION: fraction of the population that won't participate in any controls
PR_SYM_BCO_ASKTEST <- 1 # default should be 1: already normalized by frac_nothing - minimum behaviour - probability to ask for test when symptomatic and traced through bco
PR_ASYM_BCO_ASKTEST <- 1 # probability to ask for test when asymptomatic at moment of being traced through bco, relative to pr_sym_bco_asktest
PR_SYM_APP_ASKTEST <- 0.9 # probability to ask for test when symptomatic at moment of being traced through app, relative to pr_sym_bco_asktest
PR_ASYM_APP_ASKTEST <- 0.9 # probability to ask for test when asymptomatic at moment of being traced through app, relative to pr_sym_bco_asktest
PR_SYM_SELF_ASKTEST <- 0.56 # probability to ask for test due to symptoms: relative to pr_sym_bco_asktest - %pop that ask for test multiplied by 1-frac_nothing
PR_SYM_SELFBCO_ASKTEST <- 0.75 # probablility to ask for test when symptomatic at moment of being traced through informal bco - relative to pr_sym_bco_asktest

PR_POSTTEST <- 0.9 # probability to isolate after positive testoutcome: is applied only after positive testoutcome
PR_SYM_POSTTEST <- PR_POSTTEST
PR_ASYM_POSTTEST <- PR_POSTTEST

PR_SYM_BCO_PRETEST <- 0.75 # probability to quarantine when symptomatic at the moment of being traced with bco - given that he/she asks for test
PR_SYM_APP_PRETEST <- 0.75 # probability to quarantine when symptomatic at the moment of being traced with app - given that he/she asks for test
PR_SYM_SELFBCO_PRETEST <- 0.75 # probability to quarantine when symptomatic at the moment of being traced through informal bco - given that he/she asks for test
PR_SYM_SELF_PRETEST <- 0.5 # probability to quarantine when test reason is symptoms - given that he/she asks for test
PR_ASYM_BCO_PRETEST <- 0.5 # probability to quarantine when asymptomatic at the moment of being traced with bco - given that he/she asks for test
PR_ASYM_APP_PRETEST <- 0.5 # probability to quarantine when asymptomatic at the moment of being traced with app - given that he/she asks for test


################ Constructing "toicdata" #######################################
################ contact network, transmission dynamics and behaviours ##########


#' Simulates uncontrolled transmission blueprint of 1 to 'genmax' generations
#' after index case (no times), including contact network and quantiles for
#' incubation periods, generation intervals, and behavioural probabilities
#'
#' @param genmax maximum number of generations to simulate
#' @param Reff effective reproduction number
#' @param kpar dispersion coefficient
#' @param seedstart set seed for reproducibility
#' @return list of transmission chain and reproducibility of model components on top of transmissin chain:
#' - infectors = vector - each element j with index i represents the infector j of case i,
#' - generations = vector - each element j with index i represents the generation j of case i
#' - seedstart = int - reproducibility of transmission chain
#' - networkseed = int - reproducibility of network
#' - incubationseed = int - reproducibility of incubation times
#' - generationseed = int - reproducibility of generation intervals
#' - symptomaticseed= int - reproducibility of symptom onset times
#' - tracingseed = int - reproducibility of nco/informal bco tracing links
#' - appcontactseed = int - reproducibility of app tracing links
#' - asktest_seed = int - reproducibility of behaviour with regards to asking for test
#' - isoquar_seed = int - reproducibility of behaviour with regards to isolation and quarantine
#' - tracingappuse_seed = int - reproducibility of behaviour with regards to using the app
#' @examples
#' # Using default parameters
#' tcdata <- transchain()
#' # Set parameters manually
#' tcdata <- transchain(genmax = 12, Reff = 1.2, kpar = 0.5, seedstart = 888888)
transchain <- function(genmax = GEN_MAX, Reff = R_EFF, kpar = K_PARR, seedstart = SEED_START, ...) {
  set.seed(seedstart)
  # offspring from index case: at least 1
  while ((indexoffspring <- rnbinom(1, kpar, mu = Reff)) == 0) {}

  # start bookkeeping, index case is case 1 and has infector 0
  allinfectors <- c(0, rep(1, indexoffspring))
  allgenerations <- allinfectors

  # case by case expansion, generation wise, iterate through all cases until
  # genmax is reached
  nextcase <- 2
  while (
    length(allinfectors) >= nextcase
  ) {
    if (allgenerations[nextcase] < genmax) {
      offspring <- rnbinom(1, kpar, mu = Reff)
      allinfectors <- c(allinfectors, rep(nextcase, offspring))
      allgenerations <- c(allgenerations, rep(allgenerations[nextcase] + 1, offspring))
    }
    nextcase <- nextcase + 1
  }

  # assign fixed seeds to model components for reproducibility
  # transmission network
  networkseed <- seedstart + 10000

  # distribution quantiles
  incubationseed <- seedstart + 20000
  generationseed <- seedstart + 30000
  symptomaticseed <- seedstart + 40000
  tracingseed <- seedstart + 50000
  appcontactseed <- seedstart + 60000

  # behaviour quantiles: willingness to follow control measures
  asktest_seed <- seedstart + 7000 # ALT
  isoquar_seed <- seedstart + 8000 # ALT
  tracingappuse_seed <- seedstart + 23000

  # complete all output
  return(list(
    infectors = allinfectors,
    generations = allgenerations,
    seedstart = seedstart,
    networkseed = networkseed,
    incubationseed = incubationseed,
    generationseed = generationseed,
    symptomaticseed = symptomaticseed,
    tracingseed = tracingseed,
    appcontactseed = appcontactseed,
    asktest_seed = asktest_seed,
    isoquar_seed = isoquar_seed,
    tracingappuse_seed = tracingappuse_seed
  ))
}


#' Add clustering to the transmission chain: create some triangles in
#' transmission chain
#'
#' @param tcdata transmission chain data - list with elements as returned by the function transchain
#' @param prtriangle double - probability that a randomly chosen triple forms a triangle
#' @param connection_max int - max depth of clustering
#' @return list of
#' - tcdata - list containing transmission chain data
#' - networkmatrix matrix - adjacency matrix with extra triangles included
#' @examples
#' # Using defaults
#' tndata <- transnetwork()
#' # Manual settings
#' tcdata <- transchain(genmax = 12, Reff = 1.2, kpar = 0.5, seedstart = 888888)
#' tndata <- transnetwork(tcdata = tcdata, prtriangle = 0.5)
transnetwork <- function(tcdata = NULL, prtriangle = PR_TRIANGLE, connection_max = CONNECTION_MAX, ...) {
  if (is.null(tcdata)) {
    tcdata <- transchain(...)
    log_warn("No transmission chain data provided in transnetwork, using default")
  }

  # number of cases
  pop_size <- length(tcdata$infectors)
  # generate triangle probabilities
  triangleprobs <- generate_adjacency_probs(pop_size, tcdata$networkseed)

  # bluematrix -> links between infector and cases
  bluematrix <- initialize_bluematrix(pop_size, tcdata$infectors)

  # red matrix -> confirmed non-links
  # initialize, individual can't be linked to itself
  redmatrix <- diag(TRUE, nrow = pop_size)

  # propose links based on existing links and non-links
  # identify all possible links that create triangles
  candidates <- (crossprod(bluematrix) & !bluematrix) & !redmatrix
  connectionorder <- 1
  # conditions 1: candidate links, 2: total number of possibilities = pop_size^2, 3: order of connection depth
  while (any(candidates) & (sum(bluematrix) + sum(redmatrix) < pop_size^2) & connectionorder <= connection_max) {
    # accept allowed links
    bluematrix <- bluematrix | (candidates & (triangleprobs < prtriangle^connectionorder))
    # reject unallowed links
    redmatrix <- redmatrix | (candidates & (triangleprobs > prtriangle^connectionorder))
    # propose new links
    candidates <- (crossprod(bluematrix) & !bluematrix) & !redmatrix
    connectionorder <- connectionorder + 1
  }
  diag(bluematrix) <- FALSE

  return(c(
    tcdata,
    list(networkmatrix = bluematrix)
  ))
}


#' Assigns times of events (infection and symptom onset) for a given
#' transmission chain
#'
#' @param tcdata transmission chain data - list with elements as returned by the function transchain
#' @param incTmean double -  mean incubation time
#' @param correlatedinfectivity bool - wheteher infectivity of infector and
#'   offspring are correlated
#' @param probsymptomatic double - probability that an individual is symptomatic
#' @return list containing
#' - tcdata - transmission chain data used to generate timed outbreak
#' - infectiontimes = vector - element t with index i is the infection time of case i
#' - symptomonsettimes = vector - element t with index i is the symptom onset time of case i.
#' Note: symptom onset time is Inf if no symptoms develop
#' @examples
#' # default settings
#' to <- timedoutbreak()
#' # Manual settings
#' tcdata <- transchain(genmax = 12, Reff = 1.2, kpar = 0.5, seedstart = 888888)
#' todata <- transnetwork(tcdata = tcdata, probsymptomatic = 0.3)
timedoutbreak <- function(tcdata = NULL,
                          incTmean = MEAN_INC_TIME,
                          correlatedinfectivity = FALSE,
                          probsymptomatic = PR_SYMPTOMATIC, ...) {
  if (is.null(tcdata)) {
    tcdata <- transchain(...)
    log_warn("No transmission chain provided in timed outbreak, default is used")
  }
  # total number of cases - infected population size
  pop_size <- length(tcdata$infectors)

  # sample distribution quantiles
  set.seed(tcdata$incubationseed)
  incubationtimesQ <- runif(pop_size)
  set.seed(tcdata$generationseed)
  generationtimesQ <- runif(pop_size)
  set.seed(tcdata$symptomaticseed)
  symptomprobabilitiesQ <- runif(pop_size)

  incubationtimes <- time_incub(distQ = incubationtimesQ, incTmean = incTmean, ...)

  if (correlatedinfectivity) {
    generationtimes <- time_gen(
      distQ = generationtimesQ,
      relincubation = c(1, incubationtimes[tcdata$infectors]) / incTmean, ...
    )
  } else {
    generationtimes <- time_gen(distQ = generationtimesQ, ...)
  }

  # assign infection times
  infectiontimes <- c(0) # index case infected at time 0
  for (nextcase in 2:pop_size) { # all cases minus index
    # definition of generation time
    infector <- tcdata$infectors[nextcase]
    infectiontime_nextcase <- infectiontimes[infector] + generationtimes[nextcase]
    infectiontimes <- c(infectiontimes, infectiontime_nextcase) # append
  }
  # assign symptom onset times, definition of incubation time
  symptomonsettimes <- infectiontimes + incubationtimes
  # those not displaying symptoms
  symptomonsettimes[symptomprobabilitiesQ > probsymptomatic] <- Inf

  return(c(
    tcdata,
    list(
      infectiontimes = infectiontimes,
      symptomonsettimes = symptomonsettimes
    )
  ))
}


#' Initialize behaviour with respect to control measures for each individual in
#' the transmission chain. Note: behaviours are nested
#'
#' @param tcdata transmission chain data - list with elements as returned by the function transchain
#' @param pr_nothing double - fraction of the population that won't participate in any controls
#' @param pr_sym_bco_asktest double - default should be 1: already normalized by
#' frac_nothing - minimum behaviour - probability to ask for test when symptomatic
#' and traced through bco
#' @param pr_sym_posttest double - probability to isolate after positive testoutcome when symptomatic:
#' is applied only after positive testoutcome
#' @param pr_sym_bco_pretest double - probability to quarantine when symptomatic
#' at the moment of being traced with bco - given that he/she asks for test
#' @param pr_sym_app_asktest double - probability to ask for test when symptomatic
#' at moment of being traced through app, relative to pr_sym_bco_asktest
#' @param pr_sym_app_pretest double - probability to quarantine when symptomatic
#' at the moment of being traced with app - given that he/she asks for test
#' @param pr_sym_self_asktest double - probability to ask for test due to symptoms:
#' relative to pr_sym_bco_asktest - %pop that ask for test multiplied by 1-frac_nothing
#' @param pr_sym_self_pretest double - probability to quarantine when test reason
#' is symptoms - given that he/she asks for test
#' @param pr_asym_bco_asktest double - probability to ask for test when asymptomatic
#' at moment of being traced through bco, relative to pr_sym_bco_asktest
#' @param pr_asym_posttest double - probability to isolate after positive testoutcome
#' when asymptomatic: is applied only after positive testoutcome
#' @param pr_asym_bco_pretest double - probability to quarantine when asymptomatic
#' at the moment of being traced with bco - given that he/she asks for test
#' @param pr_asym_app_asktest double - probability to ask for test when asymptomatic
#' at moment of being traced through app, relative to pr_sym_bco_asktest
#' @param pr_asym_app_pretest double - probability to quarantine when asymptomatic
#' at the moment of being traced with app - given that he/she asks for test
#' @param prusetracingapp double - fraction of the population that uses the app
#' @param pr_sym_selfbco_asktest double - probablility to ask for test when symptomatic
#' at moment of being traced through informal bco - relative to pr_sym_bco_asktest
#' @param pr_sym_selfbco_pretest double - probability to quarantine when symptomatic
#' at the moment of being traced through informal bco - given that he/she asks for test
#' @return list containing the following elements
#' - tcdata - list containing transmission chain data
#' - sym_app_asktest - vector of booleans
#' - sym_app_pretest - vector of booleans
#' - sym_bco_asktest  - vector of booleans
#' - sym_posttest - vector of booleans
#' - sym_bco_pretest  - vector of booleans
#' - sym_self_asktestst - vector of booleans
#' - sym_self_pretestst - vector of booleans
#' - asym_app_asktestst - vector of booleans
#' - asym_app_pretestst - vector of booleans
#' - asym_bco_asktestst - vector of booleans
#' - asym_posttest - vector of booleans
#' - asym_bco_pretestst - vector of booleans
#' - asym_selfbco_asktest - vector of booleans
#' - asym_selfbco_pretest - vector of booleans
#' - tracingappuse - vector of booleans
#' - sym_selfbco_asktestktest - vector of booleans
#' - sym_selfbco_pretestetest - vector of booleans
#' @examples
#' # Default transmission chain and behaviour parameters
#' icdata <- initializecontrol()
initializecontrol <- function(tcdata = NULL,
                              pr_nothing = FRAC_NOTHING,
                              pr_sym_bco_asktest = PR_SYM_BCO_ASKTEST,
                              pr_sym_posttest = PR_SYM_POSTTEST,
                              pr_sym_bco_pretest = PR_SYM_BCO_PRETEST,
                              pr_sym_app_asktest = PR_SYM_APP_ASKTEST,
                              pr_sym_app_pretest = PR_SYM_APP_PRETEST,
                              pr_sym_self_asktest = PR_SYM_SELF_ASKTEST,
                              pr_sym_self_pretest = PR_SYM_SELF_PRETEST,
                              pr_asym_bco_asktest = PR_ASYM_BCO_ASKTEST,
                              pr_asym_posttest = PR_ASYM_POSTTEST,
                              pr_asym_bco_pretest = PR_ASYM_BCO_PRETEST,
                              pr_asym_app_asktest = PR_ASYM_APP_ASKTEST,
                              pr_asym_app_pretest = PR_ASYM_APP_PRETEST,
                              prusetracingapp = PR_USE_TRACINGAPP,
                              pr_sym_selfbco_asktest = PR_SYM_SELFBCO_ASKTEST,
                              pr_sym_selfbco_pretest = PR_SYM_SELFBCO_PRETEST,
                              ...) {
  # check that probabilities are consistent: throw warning/error if not
  # These consistencies need to hold for the current model assumptions,
  # In general different choices can be made, change them in case
  stopifnot(exprs = {
    pr_sym_bco_pretest / pr_sym_posttest <= 1
    pr_sym_app_pretest / pr_sym_posttest <= 1
    pr_sym_self_pretest / pr_sym_posttest <= 1
    pr_sym_selfbco_pretest / pr_sym_posttest <= 1
    pr_asym_bco_pretest / pr_asym_posttest <= 1
    pr_asym_app_pretest / pr_asym_posttest <= 1
    pr_sym_app_asktest / pr_sym_bco_asktest <= 1
    pr_sym_app_pretest / pr_sym_bco_pretest <= 1
    pr_sym_selfbco_asktest / pr_sym_bco_asktest <= 1
    pr_sym_selfbco_pretest / pr_sym_bco_pretest <= 1
    pr_asym_bco_asktest / pr_sym_bco_asktest <= 1
    pr_asym_bco_pretest / pr_sym_bco_pretest <= 1
    pr_sym_self_asktest / pr_sym_bco_asktest <= 1
    pr_sym_self_pretest / pr_sym_bco_pretest <= 1
    pr_sym_selfbco_asktest / pr_sym_bco_asktest <= 1
    pr_sym_selfbco_pretest / pr_sym_bco_pretest <= 1
    prusetracingapp / (1 - pr_nothing) <= 1
  })

  if (is.null(tcdata)) {
    tcdata <- transchain(...)
    log_warn("No transmission chain provided in initialize control, default is used")
  }

  # size of infected population
  pop_size <- length(tcdata$infectors)

  # sample distribution quantiles (low values = more compliant)
  set.seed(tcdata$asktest_seed)
  asktest_Q <- runif(pop_size)
  set.seed(tcdata$isoquar_seed)
  isoquar_Q <- runif(pop_size)
  set.seed(tcdata$tracingappuse_seed)
  tracingappuse_Q <- runif(pop_size)

  # behaviour: willingness to follow control measures, TRUE or FALSE
  # Note: we assume the behaviours to be nested
  nothing <- asktest_Q > 1 - pr_nothing

  sym_bco_asktest <- asktest_Q < pr_sym_bco_asktest * (1 - pr_nothing)
  sym_posttest <- sym_bco_asktest & isoquar_Q < pr_sym_posttest
  sym_bco_pretest <- sym_bco_asktest & isoquar_Q < pr_sym_bco_pretest

  asym_posttest <- sym_posttest

  sym_app_asktest <- asktest_Q < pr_sym_app_asktest * (1 - pr_nothing)
  sym_app_pretest <- sym_app_asktest & isoquar_Q < pr_sym_app_pretest

  asym_bco_asktest <- asktest_Q < pr_asym_bco_asktest * (1 - pr_nothing)
  asym_bco_pretest <- asym_bco_asktest & isoquar_Q < pr_asym_bco_pretest

  asym_app_asktest <- asktest_Q < pr_asym_app_asktest * (1 - pr_nothing)
  asym_app_pretest <- asym_app_asktest & isoquar_Q < pr_asym_app_pretest

  sym_self_asktest <- asktest_Q < pr_sym_self_asktest * (1 - pr_nothing)
  sym_self_pretest <- sym_self_asktest & isoquar_Q < pr_sym_self_pretest

  sym_selfbco_asktest <- asktest_Q < pr_sym_selfbco_asktest * (1 - pr_nothing)
  sym_selfbco_pretest <- sym_selfbco_asktest & isoquar_Q < pr_sym_selfbco_pretest

  tracingappuse <- !nothing & (tracingappuse_Q < prusetracingapp / (1 - pr_nothing))

  return(c(
    tcdata,
    list(
      sym_app_asktest = sym_app_asktest,
      sym_app_pretest = sym_app_pretest,
      sym_bco_asktest = sym_bco_asktest,
      sym_posttest = sym_posttest,
      sym_bco_pretest = sym_bco_pretest,
      sym_self_asktest = sym_self_asktest,
      sym_self_pretest = sym_self_pretest,
      asym_app_asktest = asym_app_asktest,
      asym_app_pretest = asym_app_pretest,
      asym_bco_asktest = asym_bco_asktest,
      asym_posttest = asym_posttest,
      asym_bco_pretest = asym_bco_pretest,
      asym_selfbco_asktest = rep(FALSE, pop_size), # asymptomatics without official tracing can't test
      asym_selfbco_pretest = rep(FALSE, pop_size), # asymptomatics that are not officially traced won't quarantine
      tracingappuse = tracingappuse,
      sym_selfbco_asktest = sym_selfbco_asktest,
      sym_selfbco_pretest = sym_selfbco_pretest
    )
  ))
}


######################## control measures: test and trace ######################


#' Testing, quarantine and isolation due to showing symptoms
#'
#' @param toicdata list containing all elements corresponding to
#' - transmission chain (function transchain),
#' - network (function transnetwork),
#' - timed outbreak (function timedoutbreak),
#' - behaviours (function initializecontrols)
#' @param isodelay double - delay from sympton onset to asking for test
#' @param testcontactdelay double - delay from asking for test to positive test outcome with symptoms
#' @param quarantineduration double - number of days in quarantine
#' @return list:
#' - toicdata
#' - isoprevented - vector of booleans: TRUE if case i is prevented by selfisolation
#' - isotimes - vector of doubles: t at index i is time of iso/quar for case i, Inf for no isolation
#' - testisotimes - vecotr of doubles: t at index i is time of testing positive for case i
#' Inf when this never happens
#' @examples
#' tcdata <- transchain(seedstart = 1)
#' todata <- timedoutbreak(tcdata = tcdata)
#' tndata <- transnetwork(tcdata = todata)
#' toicdata <- initializecontrol(tcdata = tndata)
#' titndata <- applyselfisolation(toicdata = tcdata)
#' @export
applyselfisolation <- function(toicdata = NULL, isodelay = ISO_DELAY,
                               testcontactdelay = TEST_CONTACT_DELAY, quarantineduration = QUARANTINE_DURATION, ...) {
  if (is.null(toicdata)) {
    log_warn("No timed-outbreak-initalized-control data provided in applyselfisolation, using default")
    toicdata <- transchain(...) %>%
      timedoutbreak(tcdata = ., ...) %>%
      initializecontrol(tcdata = ., ...)
  }
  # times of selfisolation for those who are willing to isolate and develop symptoms
  # note that those without symptoms have Inf symptom onset time
  isotimes <- rep(Inf, length(toicdata$infectors))
  # quarantine immediately after calling for test
  isotimes[toicdata$sym_self_asktest & toicdata$sym_self_pretest] <- toicdata$symptomonsettimes[toicdata$sym_self_asktest & toicdata$sym_self_pretest] + isodelay
  # isolate after positive test result only: delay to ask test and delay to positive test outcome after asking for test
  isotimes[toicdata$sym_self_asktest & !toicdata$sym_self_pretest & toicdata$sym_posttest] <- toicdata$symptomonsettimes[toicdata$sym_self_asktest & !toicdata$sym_self_pretest & toicdata$sym_posttest] + isodelay + testcontactdelay

  # time of positive test outcome
  testisotimes <- rep(Inf, length(toicdata$infectors))
  testisotimes[toicdata$sym_self_asktest] <- toicdata$symptomonsettimes[toicdata$sym_self_asktest] + isodelay + testcontactdelay

  # Iteratively identify cases that are prevented due to timely selfisolation of
  # infector: initalize: infection time of case is larger than isolation time of
  # infector, index case can't be prevented
  prevented <- c(
    FALSE,
    toicdata$infectiontimes[2:length(isotimes)] > isotimes[toicdata$infectors[2:length(isotimes)]]
  )

  # Identify new prevented case iteratively:
  # prevented | (toicdata$infectors %in% which(prevented)): either individual
  # itself is prevented OR when indivual is prevented, it prevents all cases
  # that would otherwise be its offspring, reassign result to prevented vector
  # (prevented <- ...) and check that update of prevented cases is made, stop
  # iteration if not  oldprevented <- prevented
  oldprevented <- prevented
  while (any((prevented <- prevented | (toicdata$infectors %in% which(prevented))) != oldprevented)) {
    oldprevented <- prevented
  }
  # prevented cases do not isolate nor do they test positive
  isotimes[prevented] <- Inf
  testisotimes[prevented] <- Inf

  # complete all output
  return(c(toicdata, list(
    isoprevented = prevented,
    isotimes = isotimes,
    testisotimes = testisotimes
  )))
}


#' Selfisolation + informal tracing + testing of contacts combined with 
#' quarantine (before test) and isolation (after positive test outcome)
#'
#' @param titndata list - epidemic outbreak data including results from the function applyselfisolation
#' @param prtrace double - probability that an existing link is traced (through informal tracing)
#' @param tracedelay double - delay in tracing contacts and them asking for test via bco after positive test
#' @param testcontactdelay double - delay from asking for test to positive test outcome with symptoms
#' @param asym_testdelay double - delay from asking for test to positive test 
#' outcome for asymptomatic case (has to wait number of days before it can test)
#' @param quarantineduration double - duration of the quarantine period
#' @return list with the following elements
#' - selfisoquarprevented: vector of booleans - element i is whether case i is prevented under applyselftracingtesting control
#' - selfisoquartesttimes: vector of doubles - element i is the time at which case i goes into iso/quar
#' - selftestpositivetimes: vector of doubles - element i is the time at which case i tests positive
applyselftracingtesting <- function(titndata = NULL,
                                    prtrace = PR_SELF_TRACE,
                                    tracedelay = SELF_TRACE_DELAY,
                                    testcontactdelay = TEST_CONTACT_DELAY,
                                    asym_testdelay = ASYM_TEST_DELAY,
                                    quarantineduration = QUARANTINE_DURATION, ...) {
  if (is.null(titndata)) {
    titndata <- applyselfisolation(...) %>% transnetwork(tcdata = ., ...)
    log_warn("No test-isolation-transmission-network provided in applytracingtesting, using default")
  }
  pop_size <- length(titndata$infectors)
  # setup tracing
  tracelinks <- trace_fun(titndata$tracingseed, pop_size, prtrace, titndata$networkmatrix)

  # start with test-positive times and prevented cases through selfisolation
  prevented <- titndata$isoprevented
  testpositivetimes <- titndata$testisotimes
  isoquartimes <- titndata$isotimes

  anothertracinground <- TRUE
  while (anothertracinground) {
    # determine tracing times of contacts
    # determine self tracing times of contacts: bco contacts
    tracingtimes <- tracingtimes_fun(
      testpositivetimes,
      tracelinks,
      tracedelay,
      titndata$infectiontimes
    )
    # determine new testpositivetimes
    proposedtestpositivetimes <- testpositivetimes_tracetesting(
      tracingtimes,
      testpositivetimes,
      titndata$symptomonsettimes,
      testcontactdelay,
      asym_testdelay,
      titndata$sym_selfbco_asktest,
      titndata$asym_selfbco_asktest
    )

    # determine isoquartimes
    selfquartimes <- quartimes_fun(
      tracingtimes,
      titndata$symptomonsettimes,
      titndata$sym_selfbco_pretest,
      titndata$asym_selfbco_pretest
    )
    isoquartimes <- pmin(isoquartimes, selfquartimes)
    # those cases prevented by isolation or quarantine (infector iso/quar or itself quarantine)
    prevented <- isoquar_prevented(
      titndata$infectors,
      titndata$infectiontimes,
      prevented,
      isoquartimes,
      quarantineduration
    )
    newtestpositivetimes <- pmin(
      testpositivetimes,
      proposedtestpositivetimes
    )
    newtestpositivetimes[prevented] <- Inf

    # those cases that go in isolation after positive test
    isoquartimes <- isotimes_fun(
      isoquartimes,
      newtestpositivetimes,
      titndata$symptomonsettimes,
      titndata$sym_posttest,
      titndata$asym_posttest
    )
    # those additionally prevented by those that isolate after positive test
    prevented <- isoquar_prevented(
      titndata$infectors,
      titndata$infectiontimes,
      prevented,
      isoquartimes,
      quarantineduration
    )
    newtestpositivetimes[prevented] <- Inf

    anothertracinground <- any(newtestpositivetimes != testpositivetimes)
    testpositivetimes <- newtestpositivetimes
  }

  return(list(
    selfisoquartestprevented = prevented,
    selfisoquartesttimes = isoquartimes,
    selftestpositivetimes = testpositivetimes
  ))
}


#' Selfisolation + informal tracing+testing + tracing+testing via formal BCO 
#' combined with quarantine (before test) and isolation (after positive test outcome)
#' 
#' @param titndata test isolation, transmission network data
#' @param prtrace double - probability that an existing link is traced (through formal bco tracing)
#' @param pr_selftrace double - probability that an existing link is traced (through informal tracing)
#' @param tracedelay double - delay in tracing contacts and them asking for test via bco after positive test
#' @param selftracedelay double - delay in tracing contacts and them asking for test via informal trace after positive test
#' @param testcontactdelay double - delay from asking for test to positive test outcome with symptoms
#' @param quarantineduration double - number of days in quarantine
#' @return list with the following elements
#' - selftracingtimes - vector of doubles: element i is the time at which case i is traced through the informal route
#' - bcotracingtimes - vector of doubles: element i is the time at which case i is traced throught bco
#' - isoquartestroutebco - vector of booleans: element i is whether case i is traced through bco
#' - isoquartestrouteself - vector of booleans: element i is whether case i is traced through informal route
#' - isoquartestprevented - vector of booleans: element i is whether case i is prevented under applytracingtesting control
#' - isoquartesttimes - vecotr of doubles: element i is time at which case i goes into iso/quar
#' - isoquartestpositivetimes - vector of doubles: element i is time at which case i tests positive
applytracingtesting <- function(titndata = NULL,
                                prtrace = PR_TRACE,
                                prselftrace = PR_SELF_TRACE,
                                tracedelay = TRACE_DELAY,
                                selftracedelay = SELF_TRACE_DELAY,
                                testcontactdelay = TEST_CONTACT_DELAY,
                                asym_testdelay = ASYM_TEST_DELAY,
                                quarantineduration = QUARANTINE_DURATION, ...) {
  if (is.null(titndata)) {
    titndata <- applyselfisolation(...) %>% transnetwork(tcdata = ., ...)
    log_warn("No test-isolation-transmission-network provided in applytracingtesting, using default")
  }
  pop_size <- length(titndata$infectors)
  # setup tracing
  tracelinks <- trace_fun(titndata$tracingseed, pop_size, prtrace, titndata$networkmatrix)
  selftracelinks <- trace_fun(titndata$tracingseed, pop_size, prselftrace, titndata$networkmatrix)

  # start with test-positive times and prevented cases through selfisolation
  prevented <- titndata$isoprevented
  testpositivetimes <- titndata$testisotimes
  isoquartimes <- titndata$isotimes

  anothertracinground <- TRUE
  while (anothertracinground) {
    # determine tracing times of contacts
    bcotracingtimes <- tracingtimes_fun(
      testpositivetimes,
      tracelinks,
      tracedelay,
      titndata$infectiontimes
    )
    # determine new testpositivetimes
    proposedbcotestpositivetimes <- testpositivetimes_tracetesting(
      bcotracingtimes,
      testpositivetimes,
      titndata$symptomonsettimes,
      testcontactdelay,
      asym_testdelay,
      titndata$sym_bco_asktest,
      titndata$asym_bco_asktest
    )
    bcoisoquartimes <- quartimes_fun(
      bcotracingtimes,
      titndata$symptomonsettimes,
      titndata$sym_bco_pretest,
      titndata$asym_bco_pretest
    )
    # determine self tracing times of contacts: bco contacts
    selftracingtimes <- tracingtimes_fun(
      testpositivetimes,
      selftracelinks,
      selftracedelay,
      titndata$infectiontimes
    )
    proposedselftestpositivetimes <- testpositivetimes_tracetesting(
      selftracingtimes,
      testpositivetimes,
      titndata$symptomonsettimes,
      testcontactdelay,
      asym_testdelay,
      titndata$sym_selfbco_asktest,
      titndata$asym_selfbco_asktest
    )
    # determine isoquartimes and prevented
    selfisoquartimes <- quartimes_fun(
      selftracingtimes,
      titndata$symptomonsettimes,
      titndata$sym_selfbco_pretest,
      titndata$asym_selfbco_pretest
    )
    isoquartimes <- pmin(isoquartimes, bcoisoquartimes, selfisoquartimes)
    # those cases prevented by isolation or quarantine (infector iso/quar or itself quarantine)
    prevented <- isoquar_prevented(
      titndata$infectors,
      titndata$infectiontimes,
      prevented,
      isoquartimes,
      quarantineduration
    )
    newtestpositivetimes <- pmin(
      testpositivetimes,
      proposedbcotestpositivetimes,
      proposedselftestpositivetimes
    )
    newtestpositivetimes[prevented] <- Inf
    # those that isolate after positive test
    isoquartimes <- isotimes_fun(
      isoquartimes,
      newtestpositivetimes,
      titndata$symptomonsettimes,
      titndata$sym_posttest,
      titndata$asym_posttest
    )
    # additionally prevented
    prevented <- isoquar_prevented(
      titndata$infectors,
      titndata$infectiontimes,
      prevented,
      isoquartimes,
      quarantineduration
    )
    newtestpositivetimes[prevented] <- Inf
    # check if there is anything left to trace
    anothertracinground <- any(newtestpositivetimes != testpositivetimes)
    testpositivetimes <- newtestpositivetimes
  }
  # Which tracing route lead to testing positive:
  # who gets positive testresult exactly after being traced + delay between being traced and positive testresult
  # calculate both for symptomatic and asymptomatic delay
  testroute_bco_symp <- (bcotracingtimes + testcontactdelay <= testpositivetimes) & titndata$sym_bco_asktest
  testroute_bco_asymp <- (bcotracingtimes + asym_testdelay <= testpositivetimes) & titndata$asym_bco_asktest
  testroute_bco <- testroute_bco_symp | testroute_bco_asymp
  testroute_bco[testpositivetimes == Inf] <- FALSE

  testroute_self <- (selftracingtimes + testcontactdelay <= testpositivetimes) & titndata$sym_self_asktest
  testroute_self[testpositivetimes == Inf] <- FALSE

  return(list(
    selftracingtimes = selftracingtimes,
    bcotracingtimes = bcotracingtimes,
    isoquartestroutebco = testroute_bco,
    isoquartestrouteself = testroute_self,
    isoquartestprevented = prevented,
    isoquartesttimes = isoquartimes,
    isoquartestpositivetimes = testpositivetimes
  ))
}


#' Selfisolation + informal tracing+testing + tracing+testing via formal BCO + 
#' tracing+testing via app combined with quarantine (before test) and isolation 
#' (after positive test outcome)
#'
#' @param titndata
#' @param prapptrace double - probability that a link is traced through the app 
#' given that both contacts in the link have the app
#' @param prtrace double - probability that a link is traced through formal bco tracing
#' @param prselftrace double - probability that a link is traced through informal tracing
#' @param apptracedelay double - delay in tracing contacts and them asking for test 
#' via app tracing route after positive test
#' @param tracedelay double - delay in tracing contacts and them asking for test 
#' via bco tracing route after positive test
#' @param selftracedelay double - delay in tracing contacts and them asking for test 
#' via informal tracing route after positive test
#' @param testcontactdelay double - delay from asking for test to positive test outcome with symptoms
#' @param asym_testdelay double - delay from asking for test to positive test 
#' outcome for asymptomatic case (has to wait number of days before it can test)
#' @param quarantineduration double - duration of the quarantine period
#' @return list with elements
#' - bcotracingtimes - vector of doubles: element i is time at which case i is traced through bco
#' - selftracingtimes - vector of doubles: element i is time at which case i is traced through informal bco
#' - apptracingtimes - vector of doubles: element i is time at whcih case i is traced through app
#' - apptestrouteapp - vector of booleans: element i is whether case i asks for test through app route 
#' - apptestroutebco - vector of booleans: element i is whether case i asks for test through bco route 
#' - apptestrouteself - vector of booleans: element i is whether case i asks for test through informal route 
#' - appisoquartestprevented - vector of booleans: element i is whether case i is prevented applyapptesting route
#' - appisoquartesttimes - vector of doubles: element i is time at which case i goes into iso/quar
#' - apptestpositivetimes - vector of doubles: element i is time at which case i tests positive
applyapptesting <- function(titndata = NULL,
                            prapptrace = PR_APP_TRACE,
                            prtrace = PR_TRACE,
                            prselftrace = PR_SELF_TRACE,
                            apptracedelay = APP_TRACE_DELAY,
                            tracedelay = TRACE_DELAY,
                            selftracedelay = SELF_TRACE_DELAY,
                            testcontactdelay = TEST_CONTACT_DELAY,
                            asym_testdelay = ASYM_TEST_DELAY,
                            quarantineduration = QUARANTINE_DURATION, ...) {
  if (is.null(titndata)) {
    titndata <- applyselfisolation(...) %>% transnetwork(tcdata = ., ...)
  }
  # set up tracelinks for different tracing methods
  pop_size <- length(titndata$infectors)
  appusematrix <- generate_appusematrix(pop_size, titndata$tracingappuse)

  # identify contacts that will be mentioned (one way: column mentioned by row)
  tracelinks <- trace_fun(titndata$tracingseed, pop_size, prtrace, titndata$networkmatrix)
  selftracelinks <- trace_fun(titndata$tracingseed, pop_size, prselftrace, titndata$networkmatrix)
  apptracelinks <- trace_fun(titndata$appcontactseed, pop_size, prapptrace, titndata$networkmatrix) & appusematrix

  # start with test-positive times and prevented cases through selfisolation
  prevented <- titndata$isoprevented
  testpositivetimes <- titndata$testisotimes
  isoquartimes <- titndata$isotimes

  anothertracinground <- TRUE
  while (anothertracinground) {
    # determine bco tracing times of contacts
    bcotracingtimes <- tracingtimes_fun(
      testpositivetimes,
      tracelinks,
      tracedelay,
      titndata$infectiontimes
    )
    # determine new testpositivetimes
    proposedbcotestpositivetimes <- testpositivetimes_tracetesting(
      bcotracingtimes,
      testpositivetimes,
      titndata$symptomonsettimes,
      testcontactdelay,
      asym_testdelay,
      titndata$sym_bco_asktest,
      titndata$asym_bco_asktest
    )
    bcoisoquartimes <- quartimes_fun(
      bcotracingtimes,
      titndata$symptomonsettimes,
      titndata$sym_bco_pretest,
      titndata$asym_bco_pretest
    )
    # determine app tracing times of contacts
    apptracingtimes <- tracingtimes_fun(
      testpositivetimes,
      apptracelinks,
      apptracedelay,
      titndata$infectiontimes
    )
    proposedapptestpositivetimes <- testpositivetimes_tracetesting(
      apptracingtimes,
      testpositivetimes,
      titndata$symptomonsettimes,
      testcontactdelay,
      asym_testdelay,
      titndata$sym_app_asktest,
      titndata$asym_app_asktest
    )
    appisoquartimes <- quartimes_fun(
      apptracingtimes,
      titndata$symptomonsettimes,
      titndata$sym_app_pretest,
      titndata$asym_app_pretest
    )
    # determine self tracing times of contacts: bco contacts
    selftracingtimes <- tracingtimes_fun(
      testpositivetimes,
      selftracelinks,
      selftracedelay,
      titndata$infectiontimes
    )
    proposedselftestpositivetimes <- testpositivetimes_tracetesting(
      selftracingtimes,
      testpositivetimes,
      titndata$symptomonsettimes,
      testcontactdelay,
      asym_testdelay,
      titndata$sym_selfbco_asktest,
      titndata$asym_selfbco_asktest
    )
    # determine isoquartimes and prevented
    selfisoquartimes <- quartimes_fun(
      selftracingtimes,
      titndata$symptomonsettimes,
      titndata$sym_selfbco_pretest,
      titndata$asym_selfbco_pretest
    )
    # # determine isoquartimes and prevented
    # isoquar through tracing is only valid if before pure selfisolation time
    isoquartimes <- pmin(isoquartimes, bcoisoquartimes, appisoquartimes, selfisoquartimes)
    prevented <- isoquar_prevented(
      titndata$infectors,
      titndata$infectiontimes,
      prevented,
      isoquartimes,
      quarantineduration
    )
    new_testpositivetimes <- pmin(
      testpositivetimes,
      proposedbcotestpositivetimes,
      proposedapptestpositivetimes,
      proposedselftestpositivetimes
    )
    new_testpositivetimes[prevented] <- Inf
    # those that isolate after positive test
    isoquartimes <- isotimes_fun(
      isoquartimes,
      new_testpositivetimes,
      titndata$symptomonsettimes,
      titndata$sym_posttest,
      titndata$asym_posttest
    )
    # additionally prevented
    prevented <- isoquar_prevented(
      titndata$infectors,
      titndata$infectiontimes,
      prevented,
      isoquartimes,
      quarantineduration
    )
    new_testpositivetimes[prevented] <- Inf
    anothertracinground <- any(new_testpositivetimes != testpositivetimes)
    testpositivetimes <- new_testpositivetimes
  }
  # Which tracing route lead to testing positive:
  # who gets positive testresult exactly after being traced + delay between being traced and positive testresult
  # calculate both for symptomatic and asymptomatic delay
  testroute_app_symp <- (apptracingtimes + testcontactdelay <= testpositivetimes) & titndata$sym_app_asktest
  testroute_app_asymp <- (apptracingtimes + asym_testdelay <= testpositivetimes) & titndata$asym_app_asktest
  testroute_app <- testroute_app_symp | testroute_app_asymp
  testroute_app[testpositivetimes == Inf] <- FALSE

  testroute_bco_symp <- (bcotracingtimes + testcontactdelay <= testpositivetimes) & titndata$sym_bco_asktest
  testroute_bco_asymp <- (bcotracingtimes + asym_testdelay <= testpositivetimes) & titndata$asym_bco_asktest
  testroute_bco <- testroute_bco_symp | testroute_bco_asymp
  testroute_bco[testpositivetimes == Inf] <- FALSE

  testroute_self <- (selftracingtimes + testcontactdelay <= testpositivetimes) & titndata$sym_selfbco_asktest
  testroute_self[testpositivetimes == Inf] <- FALSE

  return(list(
    bcotracingtimes = bcotracingtimes,
    selftracingtimes = selftracingtimes,
    apptracingtimes = apptracingtimes,
    apptestrouteapp = testroute_app,
    apptestroutebco = testroute_bco,
    apptestrouteself = testroute_self,
    appisoquartestprevented = prevented,
    appisoquartesttimes = isoquartimes,
    apptestpositivetimes = testpositivetimes
  ))
}
