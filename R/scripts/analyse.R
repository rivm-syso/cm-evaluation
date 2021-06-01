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
################################################################################
#
# Define, simulate and analyse the scenarios
# 
# Note this script assumes that model building blocks are present in specified directory
# Run the simulations in simulate.R
################################################################################

source("R/functions/simulate_utils.R")
source("R/functions/analyse_controls_utils.R")
source("R/functions/models.R")

library(parallel)
library(glue)


# Use current date as identifyer for analysis
RESULT_ID <- as.character(strftime(Sys.Date(), "%Y%m%d"))

# dir path where the simulated outbreak lives from simulations.R
DIR_PATH_OUTBREAK <- "results_final" 


################################### define scenarios #######################################
# default scenario
scenario_default <- list(Reff = "R1-3", 
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
                         delaysi = ISO_DELAY)

####### appgebruik tov de rest default setting
scenario_default_20 <- scenario_default
scenario_default_20["appuse"] <- "app0-2"

scenario_default_30 <- scenario_default
scenario_default_30["appuse"] <- "app0-3"

scenario_default_40 <- scenario_default
scenario_default_40["appuse"] <- "app0-4"

scenario_default_60 <- scenario_default
scenario_default_60["appuse"] <- "app0-6"

scenario_default_80 <- scenario_default
scenario_default_80["appuse"] <- "app0-8"

###### shorter app delay
## for varying app uses
scenario_appdelay16 <- scenario_default
scenario_appdelay16["delayta"] <- 0.33

scenario_appdelay30 <- scenario_appdelay16
scenario_appdelay30["appuse"] <- "app0-3"

scenario_appdelay20 <- scenario_appdelay16
scenario_appdelay20["appuse"] <- "app0-2"

scenario_appdelay40 <- scenario_appdelay16
scenario_appdelay40["appuse"] <- "app0-4"

scenario_appdelay60 <- scenario_appdelay16
scenario_appdelay60["appuse"] <- "app0-6"

scenario_appdelay80 <- scenario_appdelay16
scenario_appdelay80["appuse"] <- "app0-8"


######## zonder bco
## for varying app uses
scenario_future1 <- scenario_default
scenario_future1["prtrace"] <- 0
scenario_future1["prselftrace"] <- 0.16

scenario_future1_20 <- scenario_future1
scenario_future1_20["appuse"] <- "app0-2"

scenario_future1_30 <- scenario_future1
scenario_future1_30["appuse"] <- "app0-3"

scenario_future1_40 <- scenario_future1
scenario_future1_40["appuse"] <- "app0-4"

scenario_future1_60 <- scenario_future1
scenario_future1_60["appuse"] <- "app0-6"

scenario_future1_80 <- scenario_future1
scenario_future1_80["appuse"] <- "app0-8"

########## afgeschaald bco
## for varying app use
scenario_future2 <- scenario_default
scenario_future2["prtrace"] <- 0.20
scenario_future2["prselftrace"] <- 0.16

scenario_future2_20 <- scenario_future2
scenario_future2_20["appuse"] <- "app0-2"

scenario_future2_30 <- scenario_future2
scenario_future2_30["appuse"] <- "app0-3"

scenario_future2_40 <- scenario_future2
scenario_future2_40["appuse"] <- "app0-4"

scenario_future2_60 <- scenario_future2
scenario_future2_60["appuse"] <- "app0-6"

scenario_future2_80 <- scenario_future2
scenario_future2_80["appuse"] <- "app0-8"

########### Optimistisch toekomst scenarios
#### Korte app delay
## Geen BCO
###### Optimaal scenario1
scenario_optimaal1 <- scenario_default
scenario_optimaal1["appuse"] <- "app0-3"
scenario_optimaal1["delayta"] <- 0.33
scenario_optimaal1["prtrace"] <- 0
scenario_optimaal1["prselftrace"] <- 0.16 

# # varying app use
scenario_optimaal1_16 <- scenario_optimaal1
scenario_optimaal1_16["appuse"] <- "app0-16"

scenario_optimaal1_20 <- scenario_optimaal1
scenario_optimaal1_20["appuse"] <- "app0-2"

scenario_optimaal1_40 <- scenario_optimaal1
scenario_optimaal1_40["appuse"] <- "app0-4"

scenario_optimaal1_60 <- scenario_optimaal1
scenario_optimaal1_60["appuse"] <- "app0-6"

scenario_optimaal1_80 <- scenario_optimaal1
scenario_optimaal1_80["appuse"] <- "app0-8"

### Afgeschaald BCO
# Optimaal scenario2
scenario_optimaal2 <- scenario_optimaal1
scenario_optimaal2["prtrace"] <- 0.20 

scenario_optimaal2_16 <- scenario_optimaal2
scenario_optimaal2_16["appuse"] <- "app0-16"

scenario_optimaal2_20 <- scenario_optimaal2
scenario_optimaal2_20["appuse"] <- "app0-2"

scenario_optimaal2_40 <- scenario_optimaal2
scenario_optimaal2_40["appuse"] <- "app0-4"

scenario_optimaal2_60 <- scenario_optimaal2
scenario_optimaal2_60["appuse"] <- "app0-6"

scenario_optimaal2_80 <- scenario_optimaal2
scenario_optimaal2_80["appuse"] <- "app0-8"


###### niets doen
## clustering van app use en populatie dat niks doet
scenario_appclusterplus <- scenario_default
scenario_appclusterplus["IQlevel"] <- "IQ0-05"
scenario_appclusterplus["appuse"] <- "app0-24"

scenario_appclustermin <- scenario_default
scenario_appclustermin["IQlevel"] <- "IQ0-15"
scenario_appclustermin["appuse"] <- "app0-08"

# clustering van populatie dat niets doen
scenario_nothing10 <- scenario_default
scenario_nothing10["IQlevel"] <- "IQ0-05"

scenario_nothing30 <- scenario_default
scenario_nothing30["IQlevel"] <- "IQ0-15"


########### sensitivity analysis
# Sensitivity: pr_trace from coronit
scenario_trace <- scenario_default
scenario_trace["prtrace"] <- 0.31
scenario_trace["prselftrace"] <- 0.25

########## infection parameters
### R_eff
# optimaal scenario met R_eff = 1.05
scenario_R_default <- scenario_default
scenario_R_default["Reff"] <- "R1-05"

### Generatieinterval
# default
scenario_gen5 <- scenario_default
scenario_gen5["genT"] <- "gTm5s4U"

### dispersionparam
scenario_kpar_large <- scenario_default
scenario_kpar_large["kpar"] <- "k0-5"

scenario_kpar_small <- scenario_default
scenario_kpar_small["kpar"] <- "k0-05"

### clustercoefficient
scenario_cluster <- scenario_default
scenario_cluster["phi"] <- "phi0-3"

# Sensitivity: zero clustering
scenario_nocluster <- scenario_default
scenario_nocluster["phi"] <- "phi0"

# No cluster, higher app use
scenario_nocluster40 <- scenario_nocluster
scenario_nocluster40["appuse"] <- "app0-4"

### pr_symptomatic
scenario_symp <- scenario_default
scenario_symp["symp"] <- "symp0-3"

############################## save all scenarios with parameters 
list_analysis <- ls()[grepl("^scenario_", ls())]

scenario_row <- rbind(get(list_analysis[1])) %>% as.data.frame()
df_scenarios <- scenario_row %>% mutate(scenario = list_analysis[1])

for (name in list_analysis) {
  scenario_row <- rbind(get(name)) %>% as.data.frame()
  scenario_row <- scenario_row %>% mutate(scenario = name)
  df_scenarios <- df_scenarios %>% add_row(scenario_row)
}
df_scenarios <- df_scenarios %>% mutate(
  pr_sym_bco_asktest = PR_SYM_BCO_ASKTEST,
  pr_asym_bco_asktest = PR_ASYM_BCO_ASKTEST,
  pr_sym_app_asktest = PR_SYM_APP_ASKTEST,
  pr_asym_app_asktest = PR_ASYM_APP_ASKTEST,
  pr_sym_self_asktest = PR_SYM_SELF_ASKTEST,
  pr_sym_selfbco_asktest = PR_SYM_SELFBCO_ASKTEST,
  pr_sym_posttest = PR_SYM_POSTTEST,
  pr_asym_posttest = PR_ASYM_POSTTEST,
  pr_sym_bco_pretest = PR_SYM_BCO_PRETEST, 
  pr_sym_app_pretest = PR_SYM_APP_PRETEST, 
  pr_sym_selfbco_pretest = PR_SYM_SELFBCO_PRETEST,
  pr_sym_self_pretest = PR_SYM_SELF_PRETEST,
  pr_asym_bco_pretest = PR_ASYM_BCO_PRETEST,
  pr_asym_app_pretest = PR_ASYM_APP_PRETEST)

df_scenarios %>% saveRDS(file.path(DIR_PATH_OUTBREAK, glue("scenario_parameters_{RESULT_ID}.rds")))


######################## simulate scenarios ##########################################
do.call(analysescenario, scenario_default)
do.call(analysescenario, scenario_default_20)
do.call(analysescenario, scenario_default_30)
do.call(analysescenario, scenario_default_40)
do.call(analysescenario, scenario_default_60)
do.call(analysescenario, scenario_default_80)
do.call(analysescenario, scenario_appdelay16)
do.call(analysescenario, scenario_appdelay20)
do.call(analysescenario, scenario_appdelay30)
do.call(analysescenario, scenario_appdelay40)
do.call(analysescenario, scenario_appdelay60)
do.call(analysescenario, scenario_appdelay80)
do.call(analysescenario, scenario_future1)
do.call(analysescenario, scenario_future1_20)
do.call(analysescenario, scenario_future1_30)
do.call(analysescenario, scenario_future1_40)
do.call(analysescenario, scenario_future1_60)
do.call(analysescenario, scenario_future1_80)
do.call(analysescenario, scenario_future2)
do.call(analysescenario, scenario_future2_20)
do.call(analysescenario, scenario_future2_30)
do.call(analysescenario, scenario_future2_40)
do.call(analysescenario, scenario_future2_60)
do.call(analysescenario, scenario_future2_80)
do.call(analysescenario, scenario_optimaal1)
do.call(analysescenario, scenario_optimaal1_16)
do.call(analysescenario, scenario_optimaal1_20)
do.call(analysescenario, scenario_optimaal1_40)
do.call(analysescenario, scenario_optimaal1_60)
do.call(analysescenario, scenario_optimaal1_80)
do.call(analysescenario, scenario_optimaal2)
do.call(analysescenario, scenario_optimaal2_16)
do.call(analysescenario, scenario_optimaal2_20)
do.call(analysescenario, scenario_optimaal2_40)
do.call(analysescenario, scenario_optimaal2_60)
do.call(analysescenario, scenario_optimaal2_80)

########### sensitivity compared to default
do.call(analysescenario, scenario_trace)
do.call(analysescenario, scenario_appclustermin)
do.call(analysescenario, scenario_appclusterplus)
do.call(analysescenario, scenario_nocluster)
do.call(analysescenario, scenario_nocluster40)
do.call(analysescenario, scenario_kpar_large)
do.call(analysescenario, scenario_kpar_small)
do.call(analysescenario, scenario_cluster)
do.call(analysescenario, scenario_symp)
do.call(analysescenario, scenario_nothing10)
do.call(analysescenario, scenario_nothing30)
do.call(analysescenario, scenario_R_default)
do.call(analysescenario, scenario_gen5)


####################### summarise results #######################################
# append summary of results into one tibble
summarytibble = tibble()

summarise_scenarioresult(scenario_description = "default", scenario_list_name = "scenario_default", summarytibble = summarytibble)

summarise_scenarioresult(scenario_name = "default_20", scenario_list_name = "scenario_default_20", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "default_30", scenario_list_name = "scenario_default_30", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "default_40", scenario_list_name = "scenario_default_40", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "default_60", scenario_list_name = "scenario_default_60", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "default_80", scenario_list_name = "scenario_default_80", summarytibble = summarytibble)

summarise_scenarioresult(scenario_name = "default_appdelay_16", scenario_list_name = "scenario_appdelay16", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "default_appdelay_20", scenario_list_name = "scenario_appdelay20", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "default_appdelay_30", scenario_list_name = "scenario_appdelay30", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "default_appdelay_40", scenario_list_name = "scenario_appdelay40", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "default_appdelay_60", scenario_list_name = "scenario_appdelay60", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "default_appdelay_80", scenario_list_name = "scenario_appdelay80", summarytibble = summarytibble)
 
summarise_scenarioresult(scenario_name = "future B app 16", scenario_list_name = "scenario_future1", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "future B app 20", scenario_list_name = "scenario_future1_20", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "future B app 30", scenario_list_name = "scenario_future1_30", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "future B app 40", scenario_list_name = "scenario_future1_40", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "future B app 60", scenario_list_name = "scenario_future1_60", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "future B app 80", scenario_list_name = "scenario_future1_80", summarytibble = summarytibble)
 
summarise_scenarioresult(scenario_name = "future A app 16", scenario_list_name = "scenario_future2", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "future A app 20", scenario_list_name = "scenario_future2_20", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "future A app 30", scenario_list_name = "scenario_future2_30", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "future A app 40", scenario_list_name = "scenario_future2_40", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "future A app 60", scenario_list_name = "scenario_future2_60", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "future A app 80", scenario_list_name = "scenario_future2_80", summarytibble = summarytibble)

summarise_scenarioresult(scenario_name = "app 16 w/o bco", scenario_list_name = "scenario_optimaal1_16", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "app 20 w/o bco", scenario_list_name = "scenario_optimaal1_20", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "app 30 w/o bco", scenario_list_name = "scenario_optimaal1", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "app 40 w/o bco", scenario_list_name = "scenario_optimaal1_40", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "app 60 w/o bco", scenario_list_name = "scenario_optimaal1_60", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "app 80 w/o bco", scenario_list_name = "scenario_optimaal1_80", summarytibble = summarytibble)
 
summarise_scenarioresult(scenario_name = "app 16 w/ bco", scenario_list_name = "scenario_optimaal2_16", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "app 20 w/ bco", scenario_list_name = "scenario_optimaal2_20", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "app 30 w/ bco", scenario_list_name = "scenario_optimaal2", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "app 40 w/ bco", scenario_list_name = "scenario_optimaal2_40", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "app 60 w/ bco", scenario_list_name = "scenario_optimaal2_60", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "app 80 w/ bco", scenario_list_name = "scenario_optimaal2_80", summarytibble = summarytibble)

summarise_scenarioresult(scenario_name = "R_eff 1.05 default", scenario_list_name = "scenario_R_default", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "gen 5 days", scenario_list_name = "scenario_gen5", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "kpar 0.5", scenario_list_name = "scenario_kpar_large", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "nothing 10", summarytibble = summarytibble, scenario_list_name = "scenario_nothing10")
summarise_scenarioresult(scenario_name = "nothing 30", summarytibble = summarytibble, scenario_list_name = "scenario_nothing30")
summarise_scenarioresult(scenario_name = "kpar 0.05", summarytibble = summarytibble, scenario_list_name = "scenario_kpar_small")
summarise_scenarioresult(scenario_name = "cluster 0.3", summarytibble = summarytibble, scenario_list_name = "scenario_cluster")
summarise_scenarioresult(scenario_name = "no_cluster", scenario_list_name = "scenario_nocluster", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "pr_symptomatic 60", scenario_list_name = "scenario_symp", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "pr_trace", scenario_list_name = "scenario_trace", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "app_nothing_optimistic_cluster", scenario_list_name = "scenario_appclusterplus", summarytibble = summarytibble)
summarise_scenarioresult(scenario_name = "app_nothing_negative_cluster", scenario_list_name = "scenario_appclustermin", summarytibble = summarytibble)


###### create summary tables of scenario results for final analysis

# Decrease in R relative to R0
summarytibble <- summarytibble %>% mutate(
  dRSelfTR0 = (RSelfT - R0) / R0 * 100,
  dRBCOTR0 = (RBCOT - R0) / R0 * 100,
  dRAppTR0 = (RAppT - R0) / R0 * 100
) 

sm_tbl <- summarytibble %>%
  select(scenario, R0, Riso, RSelfT, RBCOT, RAppT, dRSelfTR0, dRBCOTR0, dRAppTR0)

# Fraction of app notificiations that directly prevent cases
prevented_by_app <- summarytibble %>%
  select(scenario, n_test_approute, test_approute_prevented, n_test_bcoroute, n_test_selfroute) %>%
  mutate(frac = test_approute_prevented / n_test_approute)

## Save results into rds files
sm_tbl %>% saveRDS(file.path(DIR_OUTPUT_PATH, glue("summarytable_reduction_R_{RESULT_ID}.rds")))
prevented_by_app %>% saveRDS(file.path(DIR_OUTPUT_PATH, glue("summarytable_fraction_prevented_by_app_{RESULT_ID}.rds")))

