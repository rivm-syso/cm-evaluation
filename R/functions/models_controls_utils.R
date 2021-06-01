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
#####################################################
#
# Util functions for applying controls to epidemic
#
#####################################################


#' Cases that are prevented because of iso/quar
#'
#' @description
#' An individual is prevented from getting infected if
#' 1. its infector is prevented OR
#' 2. its infector is in quarantine when the individual would
#' otherwise have gotten infected OR
#' 3. individual itself is in quarantine in otherwise infection time.
#' Note 1 and 2 cover both cases that infector isolates or quarantines separately
#'
#' @param infectors vector of ints - element i is the infector of case i
#' @param infectiontimes vector of doubles- element i is the infection time of case i
#' @param prevented vector of bool - element i is whether case i is prevented
#' @param isoquartimes vector of doubles - element i is whether case is goes into iso/quarantine
#' @param quarantineduration duration of quarantine in days
#' @return prevented - vector of bool: element i is whether case i is prevented
isoquar_prevented <- function(infectors,
                              infectiontimes,
                              prevented,
                              isoquartimes,
                              quarantineduration) {
  for (nextcase in order(infectiontimes)[-1]) {
    inf_time <- infectiontimes[nextcase] # infection time of nextcase
    infector <- infectors[nextcase] # infector of nextcase
    # whether infector is in quarantine
    infector_quarantine <- (inf_time %>%
      between(isoquartimes[infector], isoquartimes[infector] + quarantineduration))
    # whether nextcase is in quarantine
    individual_quarantine <- inf_time %>%
      between(isoquartimes[nextcase], isoquartimes[nextcase] + quarantineduration)
    if (prevented[infector] | infector_quarantine | individual_quarantine) {
      prevented[nextcase] <- TRUE
    }
  }
  return(prevented)
}


#' Calculate tracing times of individuals
#'
#' @param testpositivetimes vector of doubles - element i is time at which individual i tests positive
#' @param tracelinks matrix of booleans- element m_(i,j) is whether there is a link between cases i and j through tracingroute of interest
#' @param tracedelay double - delay in case getting traced and his/her contacts asking for test as a result - tracingroute dependent
#' @param infectiontimes vector of doubles - times at which case i gets infected
#' @return tracingtimes - vector of doubles - element i is time at which case i gets traced
tracingtimes_fun <- function(testpositivetimes,
                             tracelinks,
                             tracedelay,
                             infectiontimes) {
  tracingtimes <- sapply(1:length(testpositivetimes), function(x) {
    candidatetimes <- c(Inf, testpositivetimes[tracelinks[, x]] + tracedelay)
    return(min(candidatetimes[candidatetimes > infectiontimes[x]]))
  })
  return(tracingtimes)
}


#TODO delete if deleting all quarantine functions
#' Calculate proposed times at which individuals can get tested positive due to
#' trace and quarantine of individuals. Model policy is that you can only test after
#' developing symptoms
#'
#' @param tracingtimes time at which contact gets traced
#' @param symptomonsettimes time at which case develops symptoms (and is eligble for testing)
#' @param testcontactdelay delay in testing of contact (in days)
#' @param sym_asktest_complaince
#' @return proposedtestpositivetimes
#'
testpositivetimes_tracequarantine <- function(tracingtimes,
                                              symptomonsettimes,
                                              testcontactdelay,
                                              sym_asktest_compliance) {
  # if symptomonsettimes is before tracingtime: traced individual will get tested at tracingtime,
  # otherwise possibly at symptomonsettime
  # include testcontactdelay
  proposedtestpositivetimes <- pmax(tracingtimes, symptomonsettimes) + testcontactdelay
  # remove those that are unwilling to ask for test
  proposedtestpositivetimes[!sym_asktest_compliance] <- Inf
  return(proposedtestpositivetimes)
}


#' Calculation of testpositivetimes due to one round of tracing and testing
#'
#' @param tracingtimes vector of doubles - element i is time at which case i is traced
#' @param testpositivetimes vector of doubles - element i is time at which 
#' case i tests positive in previous round
#' @param symptomonsettimes vector of doubles - element i is time at which case i 
#' shows symptoms - Inf if he/she stays asymptomatic
#' @param testcontactdelay double - delay from asking for test to positive test outcome with symptoms
#' @param asym_testdelay double - delay from asking for test to positive test 
#' outcome for asymptomatic case (has to wait number of days before it can test)
#' @param sym_asktest_compliance vector of bools - element i is whether case i asks for test if symptomatic
#' @param asym_asktest_compliance vector of bools - element i is whether case i asks for test if asymptomatic
#' @return vector of doubles - element i is testpositivetime of case i due to 
#' one round of tracing and testing 
testpositivetimes_tracetesting <- function(tracingtimes,
                                           testpositivetimes,
                                           symptomonsettimes,
                                           testcontactdelay,
                                           asym_test_delay,
                                           sym_asktest_compliance,
                                           asym_asktest_compliance) {
  # for those showing symptoms
  proposedtestpositivetimes <- pmax(tracingtimes, symptomonsettimes) + testcontactdelay
  # those that do not ask for test will not test positive
  proposedtestpositivetimes[!sym_asktest_compliance] <- Inf
  # asymptomatic testing: only if symptomonsettimes doesn't lead to an earlier testpositivetime
  asym_testing_cond <- (symptomonsettimes > tracingtimes + asym_test_delay)
  # need to condition on both asym testing condition and adherence to asym testing
  proposedtestpositivetimes[asym_testing_cond & asym_asktest_compliance] <- tracingtimes[asym_testing_cond & asym_asktest_compliance] + asym_test_delay
  # new testpositivetime only holds if it is before testpositivetime of previous trace step
  proposedtestpositivetimes <- pmin(proposedtestpositivetimes, testpositivetimes)
  return(proposedtestpositivetimes)
}


#' Times at which individuals go into quarantine when being traced
#'
#' @param tracingtimes vector of doubles - element i is time at which case i is being traced
#' @param symptomonsettimes vector of doubles - element i is time at which case i shows symptoms
#' @param sym_pretest_compliance vector of bools - element i is whether case i 
#' goes into quarantine when being traced and shows symptoms at tracingtime
#' @param asym_pretest_compliance vector of bools - element i is whether case i 
#' goes into quarantine when being traced and shows no symptoms at tracingtime
#' @return vector of doubles - element i is time at which case i goes into quarantine
quartimes_fun <- function(tracingtimes,
                          symptomonsettimes,
                          sym_pretest_compliance,
                          asym_pretest_compliance) {
  # quarantine at moment of being traced
  quartimes <- tracingtimes
  # remove individuals unwilling to quarantine pretest: 
  # consider separately symptoms/no symptoms at moment of tracing
  has_symptoms_pretest <- symptomonsettimes < tracingtimes & tracingtimes < Inf
  quartimes[has_symptoms_pretest & !sym_pretest_compliance] <- Inf
  quartimes[!has_symptoms_pretest & !asym_pretest_compliance] <- Inf
  return(quartimes)
}


#' Times at which individuals go into isolation 
#' after new testpositivetime if not already in iso/quarantine or 
#' if new isotime is before current isoquartime
#' 
#' @param isoquartimes vector of doubles - element i is time at which case i currently is in iso/quar
#' @param testpositivetimes vector of doubles - element i is time at which case i tests positive
#' @param symptomonsettimes vector of doubles - element i is time at which case i shows symptoms
#' @param sym_posttest_compliance vector of bools - element i is whether case i adheres to isolation after positive test when showing symptoms
#' @param asym_posttest_compliance vector of bools - element i is whether case i adheres to isolation after positive test when asymptomatic
#' @return vector of doubles - element i is time at which case i goes into iso/quarantine
isotimes_fun <- function(isoquartimes,
                         testpositivetimes,
                         symptomonsettimes,
                         sym_posttest_compliance,
                         asym_posttest_compliance) {
  # isolate after testing positive if not already in quarantine
  isoquartimes <- pmin(isoquartimes, testpositivetimes)
  # remove individuals unwilling to isolate posttest
  # consider separately those that have symptoms/no symptoms at moment of testing positive
  has_symptoms_posttest <- symptomonsettimes < testpositivetimes & testpositivetimes < Inf
  isoquartimes[has_symptoms_posttest & !sym_posttest_compliance] <- Inf
  isoquartimes[!has_symptoms_posttest & !asym_posttest_compliance] <- Inf
  return(isoquartimes)
}


#' Calculate the links through which contacts can be traced
#'
#' @param seed int - seed number for reproducibility
#' @param pop_size int - size of population
#' @param pr_for_trace double - probability that a contact is traced, dependent on trace route
#' @param networkmatrix matrix of bools - adjacency matrix m[i,j] = m[j,i] whether 
#' there is a contact link between individuals i and j
#' @return matrix of bools - tracelinks[i,j] = tracelinks[j,i] is
#'   whether contact j mentions contact i in tracing process and vice versa
trace_fun <- function(seed, pop_size, pr_for_trace, networkmatrix) {
  # probability matrix p_ij is the probability that contact j is mentioned by case i
  tracingmatrixQ <- generate_adjacency_probs(pop_size, seed)
  # identify contacts that will be mentioned
  # include contact network (note: tracelinks is symmetric matrix)
  tracelinks <- (tracingmatrixQ > 1 - pr_for_trace) & networkmatrix
  return(tracelinks)
}
