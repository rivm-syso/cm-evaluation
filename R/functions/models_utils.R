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
#############################################
#
# Helper functions for simulation models
#
#############################################


#' Generate incubation times from gamma distribution
#'
#' @param distQ vector of probabilities randomly generated 
#' @param incTGashape double - shape of incubation time distribution
#' @param incTmean double - mean incubation time
#' @return vector of doubles - element i is the incubation time of case i
time_incub <- function(distQ, incTGashape = 4, incTmean = 5, ...) {
  return(qgamma(distQ, shape = incTGashape, scale = incTmean / incTGashape))
}


#' Offspring generation interval, distribution may depend on incubation time
#'
#' @param distQ vector of probabilities randomly generated
#' @param genTGashape double - shape of generation interval distribution
#' @param genTmean double - mean generation time
#' @param relincubation double 
#' @return vector of doubles - gamma distributed offspring interval time
time_gen <- function(distQ, genTGashape = 4, genTmean = 4, relincubation = 1, ...) {
  return(qgamma(
    distQ,
    shape = genTGashape, scale = genTmean * relincubation / genTGashape
  ))
}


#' Probabilities matrix of adjacencies
#' 
#' Element p_ij = p_ji = probability that case i and
#' j form a triangle if in a triple, p_ii = 0
#'
#' @param pop_size integer, population size of cases
#' @param matrix_seed integer, seed to ensure reproducibility
#' @return symmetric matrix of doubles with dimension (pop_size, pop_size), zeros on diagonal
generate_adjacency_probs <- function(pop_size, matrix_seed) {
  set.seed(matrix_seed)
  adjacency_probs <- matrix(runif(pop_size^2),
    ncol = pop_size,
    nrow = pop_size
  )
  # make the matrix of probabilities symmetric
  adjacency_probs[lower.tri(adjacency_probs, diag = T)] <- 0
  adjacency_probs <- adjacency_probs + t(adjacency_probs)
  return(adjacency_probs)
}


#' Initialize adjacency matrix between infectors and offspring from transmission chain data
#'
#' @param pop_size integer, population size of cases
#' @param infectors vector, element i at index j is the infector i of individual j
#' @return symmetric bool matrix of size (pop_size, pop_size), FALSE on diagona. 
#' m[i, j] = m[j, i] = TRUE if there is a link between cases i and j and FALSE if not
initialize_bluematrix <- function(pop_size, infectors) {
  bluematrix <- diag(FALSE, nrow = pop_size) # no connection with oneself
  bluematrix[cbind(1:pop_size, infectors)] <- TRUE # connections from case to its infector
  bluematrix <- bluematrix | t(bluematrix) # make symmetric: connection from infector to case
  return(bluematrix)
}


#' Helper function to generate app use matrix 
#' 
#' @param pop_size int - size of the population
#' @param traceappuse vector of bools - element i is whether case i uses the app
#' @return matrix of bools - m[i,j] = m[j,i] = TRUE if i and j both use the app 
#' and FALSE if either i or j don't use the app
generate_appusematrix <- function(pop_size, traceappuse){
  appusematrix <- matrix(
    rep(traceappuse, pop_size) * rep(traceappuse, each = pop_size), 
    ncol = pop_size
  )
  return(appusematrix)
}
