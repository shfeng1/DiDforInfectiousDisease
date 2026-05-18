rm(list=ls())
# Lin to the paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8223099/
here::i_am("2b_Kansas_Masking/0_Kansas_Clean_Data.R")
source("./global_options.R")
source("./1a_Scripts/0_Estimate_Rt.R")

covidestim <- vroom::vroom("./0_Data/covidestim-daily-fips.csv.xz") %>% 
  dplyr::select(fips, date, cases.fitted, infections, immune, deaths, cum.incidence, Rt)
df <- haven::read_dta("./0_Data/mask_ordinances_JAMA_letter_v5.dta")
county.trt <- sort(unique(df$ncounty[df$did21==1])) # list of treated county IDs

# combine Kansas case data with COVIDestim
df.model <- df %>% filter(date >= "2020-03-13") %>% # to make sure the treatment (July 24) falls on the start of a week not the middle of a week
  merge(covidestim, by.x = c("ncofips", "date"), by.y = c("fips", "date"), all.x = T) %>%
  group_by(ncounty) %>% arrange(date) %>%
  mutate(time = 1:n(), growth = infections / lag(infections, 7),
         E_lag = lag(infections, 1),
         sus_frac = (coestpop2019 - cum.incidence) / coestpop2019)
################################################################################################################################
# Compute prevalence from the incidence time-series by convolving the residence times in both 
# Exposed and Infectious compartments according to a Geometric distribution
inf_days <- 5; inf_days_std <- 2; delta <- 3 # Cori et al. instantaneous Rt --> we call it prevalence Rt
df.model$prevalence <- compute_prevalence(inf_mean=inf_days, ID=df.model$ncounty, inc=df.model$infections, 
                                          time=df.model$time, Ttot=max(df.model$time))
df.model$infected_est <- compute_infected(delta=delta, ID=df.model$ncounty, inc=df.model$infections, 
                                          time=df.model$time, Ttot=max(df.model$time)-1)
df.model$I_est <- compute_prevalence(inf_mean=inf_days, ID=df.model$ncounty, inc=df.model$E_lag / delta, 
                                     time=df.model$time, Ttot=max(df.model$time))
saveRDS(df.model, "./0_Data/Kansas.rds")
################################################################################################################################
df.first <- df.model %>% filter(dayssincefirstcase == 1)
# clean up the data set; calculate growth, Rt, beta; and aggregate to weekly level
df.clean <- df.model %>%
  group_by(ncounty) %>%  arrange(date) %>%
  mutate(week = ceiling(time / 7),
         stnnewcases7davg = ifelse(stnnewcases7davg<0, 0, stnnewcases7davg),
         prevalence_lag = lag(prevalence, 1),
         I_est_lag = lag(I_est, 1),
         ncounty = as.numeric(haven::as_factor(ncounty))) %>%
  group_by(ncounty, week) %>%
  summarise(start_date = min(date), dayssincefirstcase = min(dayssincefirstcase),
            sus_frac = mean(sus_frac), coestpop2019 = mean(coestpop2019),
            stnnewcases7davg = mean(stnnewcases7davg), growth = mean(growth),
            Rt = mean(Rt), infections = mean(infections), E_lag = mean(E_lag),
            infected_est = mean(infected_est, na.rm = T),
            prevalence_lag = mean(prevalence_lag, na.rm = T),
            Rt_est = sum(infected_est) / sum(prevalence_lag),
            Rt_exposure = sum(infections) / sum(I_est_lag)) %>%
  dplyr::select(ncounty, week, start_date, dayssincefirstcase, coestpop2019, sus_frac, stnnewcases7davg, 
                infections, growth, infected_est, prevalence_lag, E_lag, Rt, Rt_est, Rt_exposure) %>%
  ungroup() %>%
  mutate(Rt_est = ifelse(prevalence_lag==0, NA, Rt_est),
         beta_exposure = Rt_exposure / sus_frac, # instantaneous Rt on the exposure
         beta_est = Rt_est / sus_frac, # Rt prevalence-based
         beta = Rt / sus_frac, # Rt cohort
         trt.time = (start_date >= "2020-07-24"),
         trt.unit = (ncounty %in% county.trt),
         trt_post = (trt.time & trt.unit),
         ncounty = relevel(factor(ncounty), ref = "104")) %>% 
  filter(! ncounty %in% df.first$ncounty[df.first$date >= "2020-06-24"])
saveRDS(df.clean %>% filter(start_date >= "2020-06-05", start_date <= "2020-12-04"), "./0_Data/Kansas_Cleaned.rds")
