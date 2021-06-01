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
#################################################################################
#
# stand-alone code to make the figures for the CoronaMelder report, and to 
# calculate the number of prevented cases
# imported data: results generated in "analyse", and manually copied into "data"
#   - cleaned_results = R values, and reductions of R values, for a set of named scenarios
#   - frac_notification_prevented = numbers of people tested through different routes, 
#     and cases prevented after app notification, for a set of named scenarios
#
#################################################################################

library(tidyverse)
library(lubridate)
library(ggforce)

# results from "R/scripts/analyse.R"
notificationprevented <- readRDS("publicdata/summarytable_fraction_prevented_by_app_20210423.rds")
reductionR <- readRDS("publicdata/summarytable_reduction_R_20210423.rds")

# params
date_start <- "2020-12-01" %>% ymd
date_end <- "2021-03-31" %>% ymd

################
### Figure 1 ###
################
### summarising the raw data from test streets (not publicly available) ###
# # read raw data, filter and select
# teststraten <- readRDS("data/Teststraten_data_20210413_0914.rds") %>% 
#   # filter: correct period, and date_appointment <= date_sampling <= date_result
#   filter(
#     Afspraak_aangemaakt %>% as.Date() %>% between(date_start, date_end),
#     Datum_uitslag %>% as.Date() %>% between(date_start, date_end),
#     (Afspraak_aangemaakt %>% floor_date("hour")) <= (Datum_monsterafname %>% floor_date("hour")),
#     (Datum_monsterafname %>% floor_date("hour")) <= (Datum_uitslag %>% floor_date("hour"))
#   ) %>%
#   # select relevant variables
#   select(Afspraak_aangemaakt, Datum_monsterafname, Datum_uitslag, BCONummer, Coronamelder, Uitslag)
# # make time series: number of CM cases per day
# CMseries <- teststraten %>%
#   filter(Uitslag == "Positief") %>%
#   mutate(date = as.Date(Datum_monsterafname)) %>%
#   group_by(date, Coronamelder) %>%
#   summarise(CMincidence = n()) %>%
#   ungroup()
### read the time series ###
CMseries <- readRDS("publicdata/teststreetsummary_20210413.rds")

### make figure 1 ###
CMseries %>%
  group_by(date) %>%
  mutate(CMincidence = if_else(Coronamelder == "Ja", CMincidence, sum(CMincidence)),
         Coronamelder = if_else(Coronamelder == "Ja", "Testen na notificatie CoronaMelder", "Alle testen")) %>%
  ungroup() %>%
  ggplot(aes(x = date, y = CMincidence, fill = Coronamelder)) +
  geom_col() + 
  labs(title = "Aantal positieve testen in GGD-teststraten",
       x = "Dag van testafname", y = "Aantal positieve testen",
       fill = NULL) +
  guides(fill = FALSE) +
  theme_light() +
  facet_wrap(~Coronamelder, scales = "free_y")


####################################
### Figure 2 and prevented cases ###
####################################

# calculate counterfactual cases (as in Fig S4 of Wymant et al)
### simple discretized version. The counterfactual (prevented) cases on "day_now" are the non-realized progeny
### of historical cases prior to "day" that were directly prevented by CoronaMelder. 
### Each CM-notification that happened on day "day_hist" prior to "day_now" contributes as follows:
### (1) multiplication with the number of directly prevented cases per notification. This is
###     a result from the model simulations, now in the object 'notificationprevented'
### (2) multiplication with the ratio incidence(day_now)/incidence(day_hist), under the assumption that
###     the counterfactual epidemic growth would have been the same as the realized epidemic growth
### (3) multiplication with 0.25, to account for a generation interval of 4 days, so that the epidemic
###     growth ratio of step (2) is applied only once per generation interval

### for step (2), we need the daily incidence. This is calculated for the corona dashboard of the
### ministry of Health by RIVM on a regular basis. On the dashboard it is presented as prevalence
### of infectious cases, here we use incidence of symptomatic cases:
incidencecurve <- readRDS("publicdata/incidence_2021-04-12.rds")

counterfactualcases <- function(day_now, iniday, directlyprevented) {
  # for step (2): incidence(day_now)
  inc_ondaynow <- incidencecurve %>%
    filter(date == day_now) %>% pull(incidence_mean)
  
  # for step (2): incidence(day_hist), the inverse of which is multiplied by daily notifications
  inverse_inc_history <- incidencecurve %>%
    # leave out last two days, because of generation interval of 4 days
    filter(date >= iniday & date <= day_now - 3) %>%
    left_join(CMseries %>% filter(Coronamelder == "Ja"), by = "date") %>% 
    mutate(CMincidence = if_else(is.na(CMincidence), 0L, CMincidence)) %>% 
    mutate(dayweight = CMincidence / incidence_mean) %>%
    pull(dayweight) %>% sum()
  
  return(0.25 * directlyprevented * inc_ondaynow * inverse_inc_history)
}

### calculate counterfactual curves ###
CFcurve <- sapply(seq(date_start, date_end, 1), counterfactualcases, iniday = date_start, 
                  directlyprevented = notificationprevented$frac[notificationprevented$scenario == "default"])
CFcurvemin <- sapply(seq(date_start, date_end, 1), counterfactualcases, iniday = date_start, 
                     directlyprevented = notificationprevented$frac[notificationprevented$scenario == "pr_symptomatic 60"])
CFcurvemax <- sapply(seq(date_start, date_end, 1), counterfactualcases, iniday = date_start, 
                     directlyprevented = notificationprevented$frac[notificationprevented$scenario == "gen 5 days"])

### make figure 2 ###
incidencecurve %>% 
  filter(date >= date_start, date <= date_end) %>%
  select(date, incidence_mean) %>%
  mutate(cfinc = CFcurve) %>%
  pivot_longer(2:3, names_to = "whichcases", values_to = "incidence") %>%
  mutate(whichcases = if_else(whichcases == "cfinc",
                              "Door Coronamelder voorkomen\n infecties, volgens model",
                              " Verloop van epidemie geschat\n o.b.v. ziekenhuisopnames")) %>%
  ggplot(aes(x = date, y = incidence, fill = whichcases)) +
  geom_col() +
  guides(fill = F) +
  theme_light() +
  facet_wrap(~whichcases, scales = "free_y") +
  labs(title = "Impact CoronaMelder op verloop epidemie", x = "Dag van infectie", y = "Aantal nieuwe infecties", fill = NULL)

### total sum of prevented cases
sum(CFcurve)
sum(CFcurvemin)
sum(CFcurvemax)


########################
### Figures 3  and 4 ###
########################
plotinput <- reductionR %>%
  # positions of pie segments on unit circle (x1pos and x2pos) and input for labels (dR)
  mutate(x1pos_1 = 0,
         x1pos_2 = - 0.02 * pi * dRSelfTR0,
         x1pos_3 = - 0.02 * pi * dRBCOTR0,
         x1pos_4 = - 0.02 * pi * dRAppTR0,
         x2pos_1 = x1pos_2,
         x2pos_2 = x1pos_3,
         x2pos_3 = x1pos_4,
         x2pos_4 = 2 * pi,
         dR_1 = -dRSelfTR0,
         dR_2 = -dRBCOTR0 - dR_1,
         dR_3 = -dRAppTR0 - dR_1 - dR_2,
         dR_4 = 100 - dR_1 - dR_2 - dR_3) %>%
  select(scenario, x1pos_1:dR_4) %>%
  # create variable 'wh' indicating which fraction of R (1=self, 2=bco, 3=app, 4=rest)
  pivot_longer(x1pos_1:dR_4, names_to = c("var", "wh"), names_sep = "_", values_to = "val") %>%
  # create variables for positions and labels
  pivot_wider(names_from = "var", values_from = "val") %>%
  # positions (xlab,ylab) and text (dR) of labels, keep label value in dRval
  mutate(xlab = if_else(wh == "4", 0, 1.2 * sin((x1pos + x2pos)/2)),
         ylab = if_else(wh == "4", -0.1, 1.1 * cos((x1pos + x2pos)/2)),
         dRval = dR,
         dR = paste0(format(round(dR, 1), nsmall = 1), "%")) %>%
  # move xlabs to prevent too much overlapping labels
  group_by(scenario) %>%
  mutate(xlab = if_else(wh == "2" & ylab > ylab[wh == "1"] - 0.05, pmax(xlab, xlab[wh == "1"] + 0.3), xlab),
         xlab = if_else(wh == "3" & ylab > ylab[wh == "2"] - 0.05, pmax(xlab, xlab[wh == "2"] + 0.3), xlab)) %>%
  ungroup() %>%
  # change scenario names so that they are coded by app notifier (who enters the code) and app usage level
  mutate(scenario = recode(scenario, "default" = "current_appGGD_16", "default_20" = "current_appGGD_20", 
                           "default_30" = "current_appGGD_30", "default_40" = "current_appGGD_40",
                           "default_60" = "current_appGGD_60", "default_80" = "current_appGGD_80",
                           "default_appdelay_16" = "current_appself_16", "default_appdelay_20"  = "current_appself_20",
                           "default_appdelay_30" = "current_appself_30", "default_appdelay_40" = "current_appself_40",
                           "default_appdelay_60" = "current_appself_60", "default_appdelay_80" = "current_appself_80",
                           "future B" = "futurenoBCO_appGGD_16", "future B app 20"= "futurenoBCO_appGGD_20",
                           "future B app 30" = "futurenoBCO_appGGD_30", "future B app 40" = "futurenoBCO_appGGD_40",
                           "future B app 60" = "futurenoBCO_appGGD_60", "future B app 80" = "futurenoBCO_appGGD_80",
                           "future A" = "futureBCO_appGGD_16", "future A app 20"= "futureBCO_appGGD_20",
                           "future A app 30" = "futureBCO_appGGD_30", "future A app 40" = "futureBCO_appGGD_40",
                           "future A app 60" = "futureBCO_appGGD_60", "future A app 80" = "futureBCO_appGGD_80",
                           "app 16 w/o bco" = "futurenoBCO_appself_16", "app 20 w/o bco" = "futurenoBCO_appself_20",
                           "app 30 w/o bco" = "futurenoBCO_appself_30", "app 40 w/o bco" = "futurenoBCO_appself_40",
                           "app 60 w/o bco" = "futurenoBCO_appself_60", "app 80 w/o bco" = "futurenoBCO_appself_80",
                           "app 16 w/ bco" = "futureBCO_appself_16", "app 20 w/ bco" = "futureBCO_appself_20",
                           "app 30 w/ bco" = "futureBCO_appself_30", "app 40 w/ bco" = "futureBCO_appself_40",
                           "app 60 w/ bco" = "futureBCO_appself_60", "app 80 w/ bco" = "futureBCO_appself_80",
                           "R_eff 1.05 default" = "sensRlow_appGGD_16", "gen 5 days" = "sensgenT_appGGD_16",
                           "kpar 0.5" = "senslowsuperspr_appGGD_16", "nothing 10" = "senshighadh_appGGD_16",
                           "nothing 30" = "senslowadh_appGGD_16", "kpar 0.05" = "senshighsuperspr_appGGD_16",
                           "cluster 0.3" = "senshighclust_appGGD_16", "no_cluster" = "sensnoclust_appGGD_16", 
                           "pr_symptomatic 60" = "sensprsymp_appGGD_16", "pr_trace" = "sensprtrace_appGGD_16", 
                           "old settings" = "toremove",
                           "app_nothing_optimistic_cluster" = "compliant_appGGD_16", 
                           "app_nothing_negative_cluster" = "noncompliant_appGGD_16")) %>%
  filter(scenario != "toremove") %>%
  # create variable names for app notifier and app usage level
  separate(scenario, into = c("scenario", "appnotifier", "appusage")) 

### make figure 3A ###
plotinput %>%
  filter(scenario == "current" & appnotifier == "appGGD" & appusage == "16") %>%
  mutate(scenario = recode(scenario, "current" = "basismodel (tabellen 1 t/m 4)")) %>%
  ggplot() +
  geom_arc_bar(aes(x0=0, y0=0, r0 = 0, r = c(1,1,1,.3)[as.numeric(wh)], 
                   start=x1pos, end = x2pos, fill = wh), color = NA) +
  geom_blank(aes(x = xlab * 1.1)) +
  scale_fill_manual(values = c("blue", "orange", "red", "grey"),
                    name = "Reproductiegetal",
                    labels = c("afname door testen", "afname door BCO", "afname door CM", "blijft over")) +
  scale_color_manual(values = c("blue", "orange", "red", "black"),
                     name = "Reproductiegetal",
                     labels = c("afname door testen", "afname door BCO", "afname door CM", "blijft over")) +
  geom_text(aes(x = xlab, 
                y=  ylab, 
                label = dR, color = wh), show.legend = F) +
  coord_fixed() +
  theme_void() +
  facet_wrap(~scenario)

### make figure 3B ###
plotinput %>%
  filter(scenario == "compliant" | scenario == "noncompliant") %>%
  mutate(scenario = recode(scenario, "compliant" = "sociale groepen met\nhoge adherentie",
                           "noncompliant" = "sociale groepen met\nlage adherentie")) %>%
  ggplot() +
  geom_arc_bar(aes(x0=0, y0=0, r0 = 0, r = c(1,1,1,.3)[as.numeric(wh)], 
                   start=x1pos, end = x2pos, fill = wh), color = NA) +
  geom_blank(aes(x = xlab * 1.1)) +
  scale_fill_manual(values = c("blue", "orange", "red", "grey")) +
  scale_color_manual(values = c("blue", "orange", "red", "black")) +
  geom_text(aes(x = xlab, 
                y=  ylab, 
                label = dR, color = wh)) +
  coord_fixed() +
  guides(fill = FALSE, color = FALSE) +
  theme_void() +
  facet_wrap(~scenario, nrow = 2)

### make figure 3C
plotinput %>%
  filter(scenario %in% c("sensRlow", "senslowsuperspr", "sensnoclust", "sensprsymp", "senslowadh",
                         "sensgenT", "senshighsuperspr", "senshighclust", "sensprtrace", "senshighadh")) %>%
  mutate(scenario = factor(scenario,
                           levels = c("sensRlow", "senslowsuperspr", "sensnoclust", "sensprsymp", "senslowadh",
                                      "sensgenT", "senshighsuperspr", "senshighclust", "sensprtrace", "senshighadh"))) %>%
  mutate(scenario = recode(scenario, "sensRlow" = "lager reproductiegetal", 
                           "senshighsuperspr" = "meer superspreading", 
                           "senshighclust" = "meer contacten in\ngeclusterd sociaal netwerk", 
                           "senshighadh" = "minder mensen die zich\nnooit laten testen",
                           "sensgenT" = "langer generatie-interval", 
                           "senslowsuperspr" = "minder superspreading", 
                           "sensprsymp" = "minder symptomatische\ninfecties", 
                           "senslowadh" = "meer mensen die zich\nnooit laten testen",
                           "sensprtrace" = "lagere traceringskans\n(BCO en zelf informeren)",
                           "sensnoclust" = "geen extra contacten in\ngeclusterd sociaal netwerk")) %>%
  ggplot() +
  geom_arc_bar(aes(x0=0, y0=0, r0 = 0, r = c(1,1,1,.3)[as.numeric(wh)], 
                   start=x1pos, end = x2pos, fill = wh), color = NA) +
  geom_blank(aes(x = xlab * 1.1)) +
  scale_fill_manual(values = c("blue", "orange", "red", "grey")) +
  scale_color_manual(values = c("blue", "orange", "red", "black"),
                     name = "Reproductiegetal",
                     labels = c("afname door testen", "afname door BCO", "afname door CM", "blijft over")) +
  geom_text(aes(x = xlab, 
                y=  ylab, 
                label = dR, color = wh), show.legend = F) +
  guides(fill = FALSE)+
  coord_fixed() +
  theme_void() +
  facet_wrap(~scenario, nrow = 2) +
  theme(strip.text = element_text(size = 11))

### make figure 4A ###
plotinput %>%
  filter(scenario %in% c("current", "futureBCO", "futurenoBCO"), appusage == "16") %>%
  mutate(scenario = recode(scenario, 
                           "current" = "basisscenario",
                           "futureBCO" = "toekomstscenario met BCO",
                           "futurenoBCO" = "toekomstscenario zonder BCO"),
         appnotifier = factor(appnotifier, levels = c("appGGD", "appself")),
         appnotifier = recode(appnotifier,
                              "appGGD" = "CM-notificatie via GGD\n(huidige situatie)",
                              "appself" = "CM-notificatie direct\ndoor gebruiker")) %>%
  ggplot() +
  geom_arc_bar(aes(x0=0, y0=0, r0 = 0, r = c(1,1,1,.3)[as.numeric(wh)], 
                   start=x1pos, end = x2pos, fill = wh), color = NA) +
  geom_blank(aes(x = xlab * 1.1)) +
  scale_fill_manual(values = c("blue", "orange", "red", "grey")) +
  scale_color_manual(values = c("blue", "orange", "red", "black")) +
  geom_text(aes(x = xlab, 
                y=  ylab, 
                label = dR, color = wh)) +
  coord_fixed() +
  guides(fill = FALSE, color = FALSE) +
  theme_void() +
  facet_grid(appnotifier~scenario, switch = "y") +
  theme(strip.text.y.left = element_text(angle = 0))

# make figure 4B
plotinput %>%
  filter(scenario %in% c("current", "futureBCO", "futurenoBCO"), wh != "4", appusage != "50") %>%
  mutate(scenario = recode(scenario, 
                           "current" = "basisscenario",
                           "futureBCO" = "toekomstscenario met BCO",
                           "futurenoBCO" = "toekomstscenario zonder BCO"),
         appnotifier = factor(appnotifier, levels = c("appGGD", "appself")),
         appnotifier = recode(appnotifier,
                              "appGGD" = "CM-notificatie via GGD\n(huidige situatie)",
                              "appself" = "CM-notificatie direct\ndoor gebruiker")) %>%
  mutate(wh = factor(wh, levels = c("3", "2", "1")),
         labeltext = if_else(wh == "3", dR, "")) %>%
  group_by(scenario, appnotifier, appusage) %>%
  mutate(labely = sum(dRval) + 1) %>%
  ungroup() %>%
  ggplot(aes(x = appusage, y = dRval, fill = wh)) +
  geom_col() +
  geom_text(aes(y = labely, label = labeltext, color = wh), size = 3) +
  theme_light() +
  scale_fill_manual(values = c("red", "orange", "blue", "grey"),
                    name = "Reproductiegetal",
                    labels = c("afname door CM", "afname door BCO", "afname door testen", "blijft over")) +
  scale_color_manual(values = c("red", "orange", "blue", "black"),
                     name = "Reproductiegetal",
                     labels = c("afname door CM", "afname door BCO", "afname door testen", "blijft over")) +
  labs(x = "CM-gebruik als percentage van de bevolking",
       y = "Afname in reproductiegetal (in procenten)") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 11)) +
  facet_grid(appnotifier~scenario) 
