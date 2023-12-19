##### from finalanalysis.R, with adaptations of figures to English for publication in scientific journal

### stand-alone code to make the figures for the CoronaMelder report, and to calculate the number of prevented cases
### imported data: results generated in "analyse", and manually copied into "data"
###   - cleaned_results = R values, and reductions of R values, for a set of named scenarios
###   - frac_notification_prevented = numbers of people tested through different routes, and cases prevented after app notification, for a set of named scenarios

library(tidyverse)
library(lubridate)
library(ggforce)
library(cowplot)

# results from "R/scripts/analyse.R"
notificationprevented <- readRDS("publicdata/summarytable_fraction_prevented_by_app_20210423.rds")
reductionR <- readRDS("publicdata/summarytable_reduction_R_20210423.rds")

# params
date_start <- "2020-12-01" %>% ymd
date_end <- "2021-03-31" %>% ymd

########################
### Figures 1 and 2 ###
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
  # mutate(xlab = if_else(wh == "2" & ylab > ylab[wh == "1"] - 0.05, pmax(xlab, xlab[wh == "1"] + 0.3), xlab),
  #        xlab = if_else(wh == "3" & ylab > ylab[wh == "2"] - 0.05, pmax(xlab, xlab[wh == "2"] + 0.3), xlab)) %>%
  mutate(xlab = if_else(wh == "2" & ylab > ylab[wh == "1"] - 0.1, pmax(xlab, xlab[wh == "1"] + 0.3), xlab),
         ylab = if_else(wh == "2" & ylab > ylab[wh == "1"] - 0.1, pmin(ylab, ylab[wh == "1"] - 0.15), ylab),
         xlab = if_else(wh == "3" & ylab > ylab[wh == "2"] - 0.1, pmax(xlab, xlab[wh == "2"] + 0.3), xlab),
         ylab = if_else(wh == "3" & ylab > ylab[wh == "2"] - 0.1, pmin(ylab, ylab[wh == "2"] - 0.15), ylab)) %>%
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

### make figure 3A ### => adjusted for Epidemics poster
p1A <- plotinput %>%
  filter(scenario == "current" & appnotifier == "appGGD" & appusage == "16") %>%
  mutate(scenario = recode(scenario, "current" = "Baseline parameter set\n (tables 1-4)")) %>%
  ggplot() +
  geom_arc_bar(aes(x0=0, y0=0, r0 = 0, r = c(1,1,1,.3)[as.numeric(wh)], 
                   start=x1pos, end = x2pos, fill = wh), color = NA) +
  geom_blank(aes(x = xlab * 1.1)) +
  scale_fill_manual(values = c("blue", "orange", "red", "grey"),
                    name = "Reduction in R by...",
                    labels = c("...testing and\n   informal tracing", "...manual tracing", "...tracing app", "Remaining")) +
  scale_color_manual(values = c("blue", "orange", "red", "black"),
                     name = "Reduction in R by...",
                     labels = c("...testing and\n   informal tracing", "...manual tracing", "...tracing app", "Remaining")) +
  geom_text(aes(x = xlab, 
                y=  ylab, 
                label = dR, color = wh), show.legend = F) +
  coord_fixed() +
  theme_void() +
  facet_wrap(~scenario)
# ggsave("results/figures_revised/Article_F1A.png", width = 7, height = 5, units = "cm", dpi = 1200)

### make figure 3B ###
p1B <- plotinput %>%
  filter(scenario == "compliant" | scenario == "noncompliant") %>%
  mutate(scenario = recode(scenario, "compliant" = "Baseline parameter set\n(high adherence group)",
                           "noncompliant" = "Baseline parameter set\n(low adherence group)")) %>%
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
  guides(fill = "none", color = "none") +
  theme_void() +
  facet_wrap(~scenario, nrow = 1)
# ggsave("results/figures_revised/Article_F1B.png", width = 9, height = 5, units = "cm", dpi = 1200)

### make figure 3C => adjusted for Epidemics poster
p1C <- plotinput %>%
  # mutate(xlab = case_when(wh == "1" ~ xlab - 0.1, 
  #                         wh == "3" ~ xlab + 0.1,
  #                         TRUE ~ xlab),
  #        ylab = case_when(wh == "1" ~ ylab - 0.1, 
  #                         wh == "2" ~ ylab - 0.1,
  #                         wh == "3" ~ ylab - 0.2,
  #                         TRUE ~ ylab)) %>%
  filter(scenario %in% c("sensnoclust", "sensRlow", "senslowsuperspr", "sensprsymp", "senslowadh",
                         "senshighclust", "sensgenT", "senshighsuperspr", "sensprtrace", "senshighadh")) %>%
  mutate(scenario = factor(scenario,
                           levels = c("sensnoclust", "sensRlow", "senslowsuperspr", "sensprsymp", "senslowadh",
                                      "senshighclust", "sensgenT", "senshighsuperspr", "sensprtrace", "senshighadh"))) %>%
  mutate(scenario = recode(scenario, "sensRlow" = "lower\nreproduction number", 
                           "senshighsuperspr" = "more superspreading", 
                           "senshighclust" = "more clustered\ncontacts", 
                           "senshighadh" = "fewer non-compliers",
                           "sensgenT" = "longer\ngeneration interval", 
                           "senslowsuperspr" = "less superspreading", 
                           "sensprsymp" = "fewer symptomatic\ninfections", 
                           "senslowadh" = "more non-compliers",
                           "sensprtrace" = "lower tracing probability\n(informal and manual)",
                           "sensnoclust" = "no clustered\ncontacts")) %>%
  ggplot() +
  geom_arc_bar(aes(x0=0, y0=0, r0 = 0, r = c(.9,.9,.9,.3)[as.numeric(wh)], 
                   start=x1pos, end = x2pos, fill = wh), color = NA) +
  geom_blank(aes(x = xlab * 1.1)) +
  scale_fill_manual(values = c("blue", "orange", "red", "grey")) +
  scale_color_manual(values = c("blue", "orange", "red", "black"),
                     name = "Reproductiegetal",
                     labels = c("afname door testen", "afname door BCO", "afname door CM", "blijft over")) +
  geom_text(aes(x = xlab, 
                y=  ylab, 
                label = dR, color = wh), show.legend = F, size = 4) +
  guides(fill = "none")+
  coord_fixed() +
  theme_void() +
  facet_wrap(~scenario, nrow = 2) +
  theme(strip.text = element_text(size = 9))
# ggsave("results/figures_revised/Article_F1C.png", width = 16, height = 7, units = "cm", dpi = 1200)


plot_grid(
  plot_grid(p1A, p1B, 
          nrow = 1,
          rel_widths = c(7, 9),
          labels = c("a", "b")), p1C,
  nrow = 2, 
  rel_heights = c(5, 7),
  labels = c("", "c")
)
ggsave("results/figures_revised//Fig1.png", bg = "white", width = 20, height = 15, units = "cm", dpi = 1200)
ggsave("results/figures_revised/Fig1.jpg", bg = "white", width = 20, height = 15, units = "cm", dpi = 1200)



### make figure 4A ###
p2A <- plotinput %>%
  filter(scenario %in% c("current", "futureBCO", "futurenoBCO"), appusage == "16") %>%
  mutate(scenario = recode(scenario, 
                           "current" = "Baseline parameter set",
                           "futureBCO" = "Scenario: less effective\nmanual tracing",
                           "futurenoBCO" = "Scenario: no\nmanual tracing"),
         appnotifier = factor(appnotifier, levels = c("appGGD", "appself")),
         appnotifier = recode(appnotifier,
                              "appGGD" = "app-notification\nby authorities\n",
                              "appself" = "app-notification\ndirectly by user\n")) %>%
  ggplot() +
  geom_arc_bar(aes(x0=0, y0=0, r0 = 0, r = c(.9,.9,.9,.3)[as.numeric(wh)], 
                   start=x1pos, end = x2pos, fill = wh), color = NA) +
  geom_blank(aes(x = xlab * 1.1)) +
  scale_fill_manual(values = c("blue", "orange", "red", "grey"),
                    name = "Reduction in R by...",
                    labels = c("...testing and\n   informal tracing", "...manual tracing", "...tracing app", "Remaining")) +
  scale_color_manual(values = c("blue", "orange", "red", "black"),
                     name = "Reduction in R by...",
                     labels = c("...testing and\n   informal tracing", "...manual tracing", "...tracing app", "Remaining")) +
  geom_text(aes(x = xlab, 
                y=  ylab, 
                label = dR, color = wh)) +
  coord_fixed() +
  guides(color = "none") +
  theme_void() +
  facet_grid(appnotifier~scenario, switch = "y") +
  theme(strip.text.y.left = element_text(angle = 90))
# ggsave("results/figures_revised/Article_F2A.png", width = 16, height = 7, units = "cm", dpi = 1200)

# make figure 4B => adjusted for Epidemics poster
p2B <- plotinput %>%
  filter(scenario %in% c("current", "futureBCO", "futurenoBCO"), wh != "4", appusage != "50") %>%
  mutate(scenario = recode(scenario, 
                           "current" = "Baseline parameter set",
                           "futureBCO" = "Scenario: less effective\nmanual tracing",
                           "futurenoBCO" = "Scenario: no\nmanual tracing"),
         appnotifier = factor(appnotifier, levels = c("appGGD", "appself")),
         appnotifier = recode(appnotifier,
                              "appGGD" = "app-notification\nby authorities",
                              "appself" = "app-notification\ndirectly by user")) %>%
  mutate(wh = factor(wh, levels = c("3", "2", "1")),
         labeltext = if_else(wh == "3", dR, "")) %>%
  group_by(scenario, appnotifier, appusage) %>%
  mutate(labely = sum(dRval) + 2) %>%
  ungroup() %>%
  # filter(scenario == "Base parameter set", appnotifier == "app-notification\nby authorities") %>%
  ggplot(aes(x = appusage, y = dRval, fill = wh)) +
  geom_col() +
  geom_text(aes(y = labely, label = labeltext, color = wh), size = 3) +
  theme_light() +
  scale_fill_manual(values = c("red", "orange", "blue", "grey"),
                    name = "Reduction by...",
                    labels = c("app-based,", "manual, and", "informal tracing", "blijft over")) +
  scale_color_manual(values = c("red", "orange", "blue", "black"),
                     name = "Reduction by...",
                     labels = c("app-based,", "manual, and", "informal tracing", "blijft over")) +
  labs(x = "Percentage of population using the Coronamelder app (baseline: 16%)",
       y = "Percent reduction in R") +
  theme(legend.position = "none",
        strip.text = element_text(size = 9)) +
  facet_grid(appnotifier ~ scenario) 
# ggsave("results/figures_revised/Article_F2B.png", width = 16, height = 9, units = "cm", dpi = 1200)


plot_grid(
  p2A, p2B,
  nrow = 2, 
  rel_heights = c(7, 9),
  labels = c("a", "b")
)
ggsave("results/figures_revised/Fig2.png", bg = "white", width = 20, height = 20, units = "cm", dpi = 1200)
ggsave("results/figures_revised/Fig2.jpg", bg = "white", width = 20, height = 20, units = "cm", dpi = 1200)
