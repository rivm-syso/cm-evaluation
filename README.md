# Coronamelder evaluation 

This repository contains the code developed and used for the evaluation of the Coronamelder (CM), the Dutch covid digital contact tracing app. The public report (in Dutch) that the codebase is part of is RIVM (2021) CoronaMelder: modelstudie naar effectiviteit - Digitaal contactonderzoek in de bestrijding van COVID-19 RIVM-briefrapport 2021-0092 and can be found [here](https://doi.org/10.21945/RIVM-2021-0092). 

The repository contains 

a. the simulation results that were used to produce all the numbers and figures as presented in the report (`data`)

b. the scripts to reproduce the simulation results, including all the chosen parameter values

c. the implementation of the contact tracing model 

All code has been written in the programming language [R](https://www.r-project.org/about.html); see [Requirements](#requirements) for detailed specification.

## <a name = "overview"></a> Overview

The codebase consists of two parts. The main part is the functions that encodes the different components of the model; see the report for a detailed description of the model

a. Epidemic outbreak

- transmission chain (`tc`)
- timing of events (`to`)
- network of contacts between the cases in the outbreak (`tn`)

b. Adherence of individuals in the outbreak to testing and tracing (`ic`)

c. Testing and tracing controls: these controls are applied on an epidemic outbreak with predetermined adherences of individuals to all testing and tracing interventions. 

These functions are the basis for analysing several scenarios. A scenario is described by a set of model parameters. For each scenario 10.000 simulations are run: for each run a different outbreak is simulated. The results over all these simulations are then averaged for the scenario. 

The codebase consists scripts to reproduce the results in the report, that is, the scenarios as they are specified in the report are simulated and the results are analysed, and visualised. One can adjust these scripts or write customized scripts to consider different scenarios. In the remainder of the README the organization of the codebase is explained.

## <a name = "requirements"></a> Requirements

The code has been developed and runs under the RIVM R-Studio servers, and makes use of parallelization through forking.

```
R version 3.6.0 (2019-04-26)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: Red Hat Cloud Infrastructure
```

Next to the R base packages, the following packages need to be installed

```
glue
logger
tidyverse
ggforce
lubridate
```

The code has not been optimized for speed. To give an indication, one can expect for one scenario that it takes about one hour on an RIVM R-studio server. To reproduce all the results, one should estimate an overnight job.

## Usage

### Reproducing the results

To reproduce the figures from the report, the results from the simulations are included under `publicdata`, together with the standalone script `finalanalysis.R` that does the final analysis and produces the figures from the report in `R/scripts`. 

To reproduce the results in the report including the simulation results, two additional scripts are included under `R/scripts`, and should be run in order. Please be aware that the time it takes to execute these scripts is substantial, see also [Requirements](#requirements).

```
simulate.R
analyse.R
```

1. `simulate.R` takes care of simulating a) the epidemic outbreak and b) the adherences of individuals to testing and tracing (see [Overview](#overview))

2. `analyse.R` i) defines the scenarios in the report ii) runs the simulations for the scenarios iii) processes the results from the simulations for the scenarios and iv) saves the results for final analysis and visualization

Then 

3. `finalanalysis.R` does the final analysis and produces the figures as presented in the report. One can either run this script independently from `simulate.R` and `analysis.R` using the public data that is included (`publicdata`), or use their own results e.g. produced by executing the scripts `simulate.R` and `analysis.R`. In the latter case, one needs to modify the paths in `finalanalysis.R`. Note that, to obtain the output plots in the console, one should source the script with echo, i.e.
    

      > source('R/scripts/finalanalysis.R', echo=TRUE)

Note that all parameter values that were used in the report are listed in the report together with the motivation for the values, and are provided in the code repository, so all results can be reproduced without any additional inputs. 

Furthermore, note that `simulate.R` and `analyse.R` simulate *all* the scenarios used in the report, and as said, this takes up quite some time. If the user is interested in simulating a specific scenario, then this is also possible, see [Detailed usage](#details) for details. Essentially, scenarios, such as the default scenario in the report, are simulated in `analyse.R` with the call

    > do.call(analysescenario, scenario_default)

which means that the list of parameter values `scenario_default` needs to be defined, and the files from the simulations of the epidemic outbreaks need to be present on disk (`simulate.R` takes care of all the epidemic outbreak files necessary to simulate *all* scenarios - one could either choose to rewrite `simulate.R` to suit the needs for a specific scenario or run the entire script and then simulate the desired scenario(s)). We refer the user to [Detailed usage](#details) as well as the documentation in the scripts themselves.


### <a name = "details"></a> Detailed usage

The working directory is the project root `cm-evaluation`, and is set by default if one uses Rstudio. 

The three scripts `simulate.R`, `analyse.R`, and `finalanalysis.R` are custom written around the model for simulating the scenarios as presented in the report. To help guide the user to do his/her own custom analysis using the model, in this section we explain the essential parts of the codebase, its main functions with example usage and expected output. These functions are located in `R/functions`, source the directory to have access to all functions of the code base. Documentation is found with each function, access them in the source code. 

To source all functions at once, one can e.g. run
    
    files.sources = list.files("R/functions")
    sapply(files.sources, function(x) source(file.path("R", "functions", x)))
    
Or source each function file separately. In the usage below it is assumed the user has sourced the functions and therefore can run all examples.

Example usage to simulate **once** the full set of controls of self isolation, informal trace+test, bco trace+test and app trace+test. In the console run:
    
    # simulate model building blocks
    > tc <- transchain(seedstart = 2)
    > to <- timedoutbreak(tcdata = tc)
    > ic <- initializecontrol(tcdata = tc)
    > tn <- transnetwork(tcdata = tc)
    > inputlist <- c(tn, to, ic)
    # outbreak with the control of cases testing and isolating due to symptoms
    > testres <- applyselfisolation(inputlist)
    # testres are the cases that can trigger any tracing and testing
    # apply the full set of controls including the app
    > apptracingtesting <- applyapptesting(titndata = testres)

This set of function calls yields the resulting list `apptracingtesting` that includes among others the vector of booleans named `appisoquartestprevented` which indicates for each case (case in the epidemic outbreak without any controls) whether he/she gets prevented to become a case under the full set of controls.

Note that all model parameters that are unspecified in the function calls take on the default parameter values as specified in `R/functions/models.R`.

The script `simulate.R` is built up as follows. Simulation results are written into RDS files, where the runs are partitioned so that each RDS file contains the results of 1000 runs, and each partition has an identifier so that each model part can be linked across partitions (partitions are indicated in each filename by `rep<n>` where `<n>` is an integer to identify the partition). For example `tc_<XXX>_rep0.rds` contain the simulation results of the transmission chains of partition 0 and the corresponding timing of events in those outbreaks are written into `to_<XXX>_rep0.rds`

The output filenames are built up into three parts:

1. First part (`tc`, `to`, `tn`, `ic`, `res`) indicates which model part has been simulated

2. Third part (`rep<n>`) indicates the partition id of the total number of simulations. Note we count partitions from 0

3. Second part (everything between first `_` and last `_`) indicate the parameters used for those simulations

The parameters that make up a scenario:

*The parameters are described in more detail in the accompanying report*

1. `R<x>` for the effective reproduction number
2. `k<x>` for the dispersion parameter
3. `phi<x>` for the clustering coefficient
4. `gTm<x1>s<x2><Y>` for the generationinterval where x1 is the mean, x2 the shape and Y = U, C whether correlated or uncorrelated infectivity 
5. `IQ<x>` for the fraction of the population that will not ask for a test when traced through bco
6. `app<x>` for the fraction of the population that have the app in use
7. `pt<x>` for the bco tracing probability
8. `pat<x>` for the app tracing probability
9. `pst<x>` for the self tracing probability
10. `dsi<x>` for the delay in self isolation
11. `dstt<x>`for the delay in self tracing
12. `dtt<x>` for the delay in bco tracing
13. `dta<x>` for the delay in app tracing
14. `dtest<x>` for the delay in asking for test and positive test outcome for symptomatic individuals

Note that there are remaining model parameters that are in principle not varied throughout scenarios and set to default in `R/functions/models.R`.

Unless otherwise indicated `<x>` denotes a numeric object converted to a string where the dot is replaced by a dash, e.g. an app use of 0.16 is found back in the filenames as `app0-16`.

The last part of the simulations consists of simulating the controls on the epidemic using the model building blocks that are simulated in the previous steps. The key function to simulate controls and save the results is in `R/functions/simulate_utils.R`. Again, the simulations are partioned in the same way as before into 10 RDS files of 1000 simulations each. The output filename is `res_XXX_rep<x>` built up as before into the three parts, only now the first part is `res` to indicate the results of simulating the controls (as opposed to `tn`, `to`, `ic`). 

*Description of how controls are applied to the epidemic outbreak are described in the accompanying report*

Results that are saved are

- the networkseed to associate the controls to the outbreak
- the number of cases per generation under the selfisolation control
- the number of cases per generation under the selfisolation & informal bco control
- the number of cases per generation under the selfisolation, informal bco and bco control
- the number of cases per generation under the selfisolation, informal bco, bco, and app control
- the number of cases that get tested and are traced through informal bco
- the number of cases that get tested and are traced through bco
- the number of cases that get tested and are traced through the app
- the number of cases that are directly prevented through the app

A scenario is described by the set of model parameters (the set that one varies and the set that is set to default). One can simulate 10.000 runs of a scenario by calling the function `analysescenario` (from `R/functions/simulate_controls.R`). Note that `analysescenario` requires the files from the simulations of the epidemic outbreak to be in place (`tn`, `to`, `ic`, see the details above), otherwise an error is thrown 

Example usage

    > analysescenario(
        Reff = "R1-3", 
        kpar = "k0-1", 
        phi = "phi0-2", 
        genT = "gTm4s4U",
        symp = "symp0-7",
        IQlevel = "IQ0-1", 
        appuse = "app0-16",
        prtrace = 0.4, 
        prselftrace = 0.32, 
        prapptrace = 0.75, 
        dtest = 1.3, 
        delaysi = 1.6, 
        delaystt = 0.17, 
        delaytt = 1.5,
        delayta = 1,
        dir_path_outbreak = DIR_PATH_OUTBREAK
    )

This function call applies the controls to 10.000 simulations of an outbreak and saves these simulations into 10 files `res_R1-3_k0-1_phi0-2_gTm4s5U_symp0-7_IQ0-1_app0-16_pt0-4_pat0-75_pst0-32_dsi1-6_dstt0-17_dtt1-5_dta1_dtest1-3_rep<x>`, where `x = 0,...,9`. Here `dir_path_outbreak` is the path to the directory that contains the simulation files for the model building blocks (`tn`, `to`, and `ic`) as described above and in which the results for the scenario will be written into. Note that the model building blocks need to be present for the given scenario, otherwise an error will be thrown. Likewise, when the results for the scenario already exist, an error is thrown.

Once a scenario is simulated and files with the results are written one can read in the results and perform analysis. Helper functions to read in and process the results are in `R/functions/analyse_controls_utils.R`. 

Example usage:

    > res <- readscenarioresults(
        Reff = "R1-3", 
        kpar = "k0-1", 
        phi = "0-2", 
        genT = "gT4m4sU",
        symp = "symp0-7",
        IQlevel = "IQ0-1", 
        appuse = "app0-16",
        prtrace = 0.4, 
        prselftrace = 0.32, 
        prapptrace = 0.75, 
        dtest = 1.3, 
        delaysi = 1.6, 
        delaystt = 0.17, 
        delaytt = 1.5,
        delayta = 1, 
        dir_output_path = DIR_PATH_OUTBREAK
    )

where `res` is a list. Each list contains the results of one simulation run, where the results of one simulation run is again a list containing the results that are saved of the controls as described earlier (e.g. the number of cases per generation without any controls). 

A wrapper function `summarise_scenarioresult` in `R/functions/analyse_controls_utils.R` is available that 

1. reads the results for a given scenario
2. processes the results and append them to a given tibble:

- Reproduction numbers without controls and under the controls `R0`, `Riso`, `RSelfT`, `RBCOT`, `RAppT` for no controls, self test+isolation, informal bco, bco, and app, respectively
- Numbers of cases that got tested through tracing due to informal route, bco, app route (`n_test_selfroute`, `n_test_bcoroute`, and `n_test_approute`, respectively) over all simulation runs of given scenario
- Number of cases directly prevented through the app `test_approute_prevented` over all simulation runs of given scenario

Example usage
    
    # define the scenario as a list
    > scenario_default <- list(Reff = "R1-3", 
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

    # read and process results for scenario into empty tibble named summarytibble
    > summarise_scenarioresult(scenario_description = "default", 
                             scenario_list_name = "scenario_default", 
                             summarytibble = tibble())
    # display the resulting tibble                         
    > summarytibble
        # A tibble: 1 x 10
        scenario    R0  Riso RSelfT RBCOT RAppT n_test_approute n_test_bcoroute n_test_selfroute test_approute_pr…
        <chr>    <dbl> <dbl>  <dbl> <dbl> <dbl>           <int>           <int>            <int>             <int>
    1 default   1.30  1.24   1.23  1.14  1.14           22247          365525           264034              2728

Here `scenario_description` is the human readable scenario name that one can specify. One can proceed with additional processing using the summarytibble.


## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.


## Feedback

If you encounter a clear bug, please file an issue with a minimal reproducible example on GitHub.