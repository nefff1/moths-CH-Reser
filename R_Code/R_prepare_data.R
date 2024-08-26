# SETUP ------------------------------- ########################################
################################################################################.

# ... packages #################################################################
################################################################################.

library(tidyverse); theme_set(theme_classic())

# ... read data ################################################################
################################################################################.

# moth records data, available from GBIF (downloaded as Darwin Core Archive)
# https://doi.org/10.15468/dl.gcagva
d_moths <- data.table::fread("Data/GBIF/occurrence.txt", header = T)
d_moths <- d_moths |> 
  # correct few details that differ from the dataset used for the analyses
  mutate(locality = ifelse(gbifID %in% c(3002611147, 3003135068, 3003121096),
                           "Meride, Fontana", locality),
         year = ifelse(gbifID == 3003041056, 2010, year),
         year = ifelse(gbifID == 4923547328, 2016, year),
         day = ifelse(gbifID == 3002595452, 9, day),
         taxonID = ifelse(gbifID %in% c(4906151261, 4906136165), 
                          "infospecies.ch:infofauna:32538", taxonID),
         taxonID = ifelse(gbifID == 3122744392, 
                          "infospecies.ch:infofauna:32183", taxonID),
         samplingEffort = ifelse(locality == "Ins, Landwirtschaftliche Schule Seeland" &
                                   year == 1979, "LF /80 MLL /1979:6.5.-24.11 /1 /80MLL /1979_6 /NA",
                                 samplingEffort),
         individualCount = ifelse(gbifID %in% c(3003091312, 3002794229),
                                  1, individualCount)) |> 
  select(A = year, M = month, J = day, active_hours = sampleSizeValue,
         DETAILS = samplingEffort, LOC = locality, taxonID, individualCount)

# taxonomy used in the analyses, available from the GitHub repository
d_std <- read.table("Data/d_taxonomy.txt", header = T)

# sampling details available from the GitHub repository (site x year level)
d_samplings <- read.table("Data/d_samplings.txt",
                          header = T)

# nights in which nothing was caught (missing from GBIF data)
d_nullnights <- read.table("Data/d_nullnights.txt", header = T) |> 
  mutate(Samplingdate = as.Date(Samplingdate))

# temperature and precipitation data per night and LOC
# based on data from https://www.meteoswiss.admin.ch
# variables 'TabsD' and 'RhiresD' on WGS84 grid
d_weather <- read.table("Data/d_weather.txt", header = T) |> 
  mutate(Samplingdate = as.Date(Samplingdate,
                                format = "%Y-%m-%d"))
# contains columns LOC (locality ID), Samplingdate, 
# T_2day (average temperature over two days), P_2day (precipitation sum over two days)

# trait group data available from Table S1 in manuscript
# "Changes in moth communities across 50 years depend on elevation"
d_traits <- read.table("Data/d_traits.txt", header = T, sep = "\t")
names(d_traits) <- c("Name_std", "mass_cat", "Tavg_mean_cat", 
                     "Spec", "overwintering_stage")

# estimated body mass (species-level), available from the GitHub repository
d_mass <- read.table("Data/d_mass.txt", header = T)

# ... define functions #########################################################
################################################################################.

f_dates <- function(x, year){
  x_split <- strsplit(x, ";")[[1]]
  x_split <- trimws(x_split)
  x_split <- x_split[grepl(year, x_split)]
  x_split <- gsub(year, "", x_split)
  x_split <- gsub(":", "", x_split)
  
  year_split <- strsplit(x_split, "\\+")[[1]]
  year_split <- trimws(year_split)
  
  out <- data.frame()
  for (i in seq_len(length(year_split))){
    start <- strsplit(year_split[i], "-")[[1]][1]
    start <- paste0(start, ".", year)
    start <- gsub("\\.\\.", "\\.", start)
    
    end <- strsplit(year_split[i], "-")[[1]][2]
    end <- paste0(end, ".", year)
    end <- gsub("\\.\\.", "\\.", end)
    
    
    out_sub <- data.frame(as.Date(start, format = c("%d.%m.%Y")),
                          as.Date(end, year, format = c("%d.%m.%Y")))
    names(out_sub) <- paste0(c("start", "end"), i)
    
    if (i == 1) {
      out <- out_sub
    } else {
      out <- out %>% 
        bind_cols(out_sub)
    }
    
  }
  out
}

f_scale <- function(x, sub) {
  x %>% 
    left_join(d_scalings %>% filter(data == sub), by = "var") %>% 
    mutate(out = (orig - mean) / sd) %>% 
    select(out) %>% 
    deframe()
}

# EDIT DATA --------------------------- ########################################
################################################################################.

d_moths <- d_moths |> 
  separate(DETAILS, sep = " /",
           into = c("traptype", "bulb_detail", "Dates_active",
                    "n_trap", "bulbtype", "spattemp_cluster", "samplingpair")) |>
  mutate(across(c(traptype, bulb_detail, Dates_active, n_trap, bulbtype),
                ~ trimws(.))) |>
  select(-c(spattemp_cluster, samplingpair)) |> # incomplete, complete data in d_samplings
  left_join(d_std, by = "taxonID") |>
  select(-taxonID) |>
  mutate(Samplingdate = as.Date(paste(J, M, A, sep = "-"),
                                format = "%d-%m-%Y"),
         yday = yday(Samplingdate))
  
# ... data per site ############################################################
################################################################################.

d_sites <- d_samplings %>%
  select(LOC, height) %>%
  distinct() %>%
  mutate(height_cat = ifelse(height >= 1000, "high", "low"),
         height_cat = factor(height_cat, levels = c("low", "high")))

d_samplings <- d_samplings |> 
  select(-height)

# ... data per year and site ###################################################
################################################################################.

d_samplings <- d_samplings |> 
  left_join(d_moths |> 
              select(LOC, A, traptype, Dates_active, 
                     n_trap, bulbtype) |> 
              distinct(),
            by = c("LOC", "A"),
            relationship = "one-to-one") |> 
  mutate(n_trap = as.ordered(n_trap),
         traptype2 = ifelse(traptype == "LF-Changins", "LF", traptype),
         traptype2 = as.factor(traptype2),
         traptype = as.factor(traptype),
         bulbtype = as.factor(bulbtype),
         Dates_active = ifelse(Dates_active == "NA", NA, Dates_active)) |> 
  # correct an error:
  mutate(Dates_active = ifelse(LOC == "Ins, Landwirtschaftliche Schule Seeland" &
                                 A == 1979,
                               gsub("1978", "1979", Dates_active), Dates_active)) |> 
  # special cases (manually operated fixed trap):
  mutate(Dates_active = ifelse((LOC == "Gotthard-Hospiz" &
                                 A == 1983) | LOC == "Braunwald, Gumenbahn Talstation",
                               NA, Dates_active))

# add boolean whether sampling duration information is available or not --------.

d_samplings <- d_samplings |> 
  left_join(d_moths |> 
              group_by(LOC, A) |> 
              summarise(hours_data = any(!is.na(active_hours))))

# ... data per night ###########################################################
################################################################################.

# Extract dates ----------------------------------------------------------------.

d_samplingperiods <- data.frame()
for (i in seq_len(nrow(d_samplings))) {
  if (!is.na(d_samplings$Dates_active[i])){
    d_samplingperiods <- f_dates(d_samplings$Dates_active[i], 
                                 d_samplings$A[i]) %>% 
      mutate(LOC = d_samplings$LOC[i],
             A = d_samplings$A[i],
             Dates_active_or = d_samplings$Dates_active[i]) %>% 
      bind_rows(d_samplingperiods, .)
  }
}

d_samplingdates <- data.frame()

# from fixed traps:
for (i in seq_len(nrow(d_samplingperiods))){
  
  d_samplingdates <- data.frame(LOC = d_samplingperiods$LOC[i],
                                A = d_samplingperiods$A[i],
                                Samplingdate = seq.Date(d_samplingperiods$start1[i], 
                                                        d_samplingperiods$end1[i], 1)) %>% 
    bind_rows(d_samplingdates, .)
  
  if (!is.na(d_samplingperiods$start2[i])){
    d_samplingdates <- data.frame(LOC = d_samplingperiods$LOC[i],
                                  A = d_samplingperiods$A[i],
                                  Samplingdate = seq.Date(d_samplingperiods$start2[i], 
                                                          d_samplingperiods$end2[i], 1)) %>% 
      bind_rows(d_samplingdates, .)
  }
}

# add  manually operated fixed traps (special cases)
d_samplingdates <- d_moths %>% 
  filter((LOC == "Gotthard-Hospiz" & 
            A == 1983) | 
           LOC == "Braunwald, Gumenbahn Talstation") %>% 
  select(LOC, A, Samplingdate) %>% 
  distinct() %>% 
  bind_rows(d_samplingdates, .)

# from manual traps:

for (i in seq_len(nrow(d_samplings))) {
  if (d_samplings$traptype[i] == "p") {
    d_samplingdates <- d_moths %>% 
      filter(LOC == d_samplings$LOC[i],
             A == d_samplings$A[i]) %>% 
      select(LOC, J, M, A) %>% 
      distinct() %>% 
      mutate(Samplingdate = as.Date(paste(J, M, A, sep = "-"),
                                    format = "%d-%m-%Y")) %>% 
      select(-c(J, M)) %>% 
      arrange(Samplingdate) %>% 
      bind_rows(d_samplingdates, .)
  }
}


# add null nights (nights in which nothing was caught) -------------------------.
d_samplingdates <- d_samplingdates |> 
  bind_rows(d_nullnights |> select(-active_hours))

d_samplingdates <- d_samplingdates %>% 
  mutate(yday = yday(Samplingdate))

# calculate stretches of zeronights --------------------------------------------.

d_samplingdates <- d_moths %>% 
  group_by(LOC, A, yday) %>% 
  summarise(n_rec = n(),
            .groups = "drop") %>% 
  right_join(d_samplingdates, by = c("LOC", "A", "yday")) %>% 
  mutate_at(vars(n_rec), ~ifelse(is.na(.), 0, .)) 

d_tmp <- data.frame()
for (loc_i in unique(d_samplingdates$LOC)){
  d_target <- d_samplingdates %>% 
    filter(LOC == loc_i) %>% 
    arrange(Samplingdate)
  d_target$period_zero <- NA
  counter <- 0
  for (i in seq_len(nrow(d_target))){
    if (d_target$n_rec[i] == 0){
      if (i == 1){
        counter <- counter + 1
      } else if (d_target$n_rec[i-1] != 0 | d_target$Samplingdate[i-1] != d_target$Samplingdate[i] - 1){
        counter <- counter + 1
      }
      d_target$period_zero[i] <- counter
    }
  }
  d_tmp <- d_tmp %>% 
    bind_rows(d_target)
}

d_samplingdates <- d_tmp %>% 
  group_by(LOC, period_zero) %>% 
  mutate(n_zeronights = n()) %>% 
  ungroup() %>% 
  mutate(n_zeronights = ifelse(is.na(period_zero), 0, n_zeronights))

# add sampling duration information --------------------------------------------.

d_samplingdates <- d_samplingdates |> 
  left_join(d_moths |> 
              select(LOC, Samplingdate, active_hours) |> 
              distinct() |> 
              group_by(LOC, Samplingdate) |> 
              filter(!(n() > 1 & any(!is.na(active_hours)) & is.na(active_hours))) |> 
              ungroup() |> 
              bind_rows(d_nullnights |> 
                          select(LOC, Samplingdate, active_hours)),
            by = c("LOC", "Samplingdate"),
            relationship = c("one-to-one")) 

# ... combined dataset for modelling ###########################################
################################################################################.

d_mod <-
  d_moths %>%
  mutate(visit_ID = paste(LOC, A, yday)) %>%
  left_join(d_mass, by = "Name_std",
            relationship = "many-to-one") |> 
  mutate(mass_sum = individualCount * mass / 1000) |> 
  group_by(LOC, visit_ID, A, yday) %>%
  summarise(abu_tot = sum(individualCount, na.rm = T),
            sric = length(unique(Name_std)),
            mass_tot = sum(mass_sum),
            .groups = "drop") %>%
  right_join(d_samplingdates %>%
               filter(n_zeronights < 10),  # exclude stretches of 10 an more zeronights
             by = c("LOC", "A", "yday")) %>%
  mutate(visit_ID = paste(LOC, A, yday)) %>%
  mutate_at(vars(abu_tot, sric, mass_tot), ~ ifelse(is.na(.), 0, .)) %>%
  mutate(Date_prev = Samplingdate - 1,
         visit_ID_prev = paste(LOC, year(Date_prev), yday(Date_prev)),
         sample_previous = ifelse(visit_ID_prev %in% visit_ID, "yes", "no")) %>%
  select(-visit_ID_prev) %>%
  left_join(d_weather, by = c("LOC", "Samplingdate")) %>%
  mutate(trap_ID_A = paste(LOC, A),
         trap_ID_A  = as.factor(trap_ID_A)) %>%
  left_join(d_samplings, by = c("LOC", "A"), relationship = "many-to-one") %>%
  mutate(night_ID = paste(Samplingdate, ifelse(is.na(samplingpair), 
                                               LOC, samplingpair), sep = "|"),
         sample_previous = as.factor(sample_previous),
         A_id = A) %>% # not to be transformed
  left_join(d_sites %>% select(LOC, height, height_cat), 
            by = "LOC", relationship = "many-to-one")

d_mod_ne <- d_mod %>%
  filter(!(!hours_data & traptype == "p"))

d_mod_lf <- d_mod %>%
  filter(traptype2 == "LF")

d_mod_p <- d_mod %>%
  filter(traptype == "p")

# ... scaling dataset ##########################################################
################################################################################.

# site x year level ------------------------------------------------------------.

d_scalings <- d_samplings %>%
  select(-c(A, LOC, traptype, traptype2, Dates_active,
            n_trap, bulbtype, hours_data, spattemp_cluster, samplingpair)) %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  group_by(var) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            .groups = "drop") %>%
  mutate(data = "full")

d_scalings <- d_samplings %>%
  filter(!(traptype == "p" & !hours_data)) %>%
  select(-c(A, LOC, traptype, traptype2, Dates_active,
            n_trap, bulbtype, hours_data, spattemp_cluster, samplingpair)) %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  group_by(var) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            .groups = "drop") %>%
  mutate(data = "no_estimates") %>%
  bind_rows(d_scalings, .)


d_scalings <- d_samplings %>%
  filter(traptype2 == "LF") %>%
  select(-c(A, LOC, traptype, traptype2, Dates_active,
            n_trap, bulbtype, hours_data, spattemp_cluster, samplingpair)) %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  group_by(var) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            .groups = "drop") %>%
  mutate(data = "LF") %>%
  bind_rows(d_scalings, .)

d_scalings <- d_samplings %>%
  filter(traptype == "p") %>%
  select(-c(A, LOC, traptype, traptype2, Dates_active,
            n_trap, bulbtype, hours_data, spattemp_cluster, samplingpair)) %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  group_by(var) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            .groups = "drop") %>%
  mutate(data = "p") %>%
  bind_rows(d_scalings, .)

# site level -------------------------------------------------------------------.

d_scalings <- d_sites %>%
  summarise(mean = mean(height),
            sd = sd(height),
            .groups = "drop") %>%
  mutate(var = "height",
         data = "full") %>%
  bind_rows(d_scalings, .)

d_scalings <- d_sites %>%
  filter(LOC %in% d_samplings$LOC[!(d_samplings$traptype == "p" & 
                                      !d_samplings$hours_data)]) %>%
  summarise(mean = mean(height),
            sd = sd(height),
            .groups = "drop") %>%
  mutate(var = "height",
         data = "no_estimates") %>%
  bind_rows(d_scalings, .)

d_scalings <- d_sites %>%
  filter(LOC %in% d_samplings$LOC[d_samplings$traptype2 == "LF"]) %>%
  summarise(mean = mean(height),
            sd = sd(height),
            .groups = "drop") %>%
  mutate(var = "height",
         data = "LF") %>%
  bind_rows(d_scalings, .)

d_scalings <- d_sites %>%
  filter(LOC %in% d_samplings$LOC[d_samplings$traptype == "p"]) %>%
  summarise(mean = mean(height),
            sd = sd(height),
            .groups = "drop") %>%
  mutate(var = "height",
         data = "p") %>%
  bind_rows(d_scalings, .)

# night level ------------------------------------------------------------------.

d_scalings <- d_mod %>%
  select(A, yday, T_2day, P_2day) %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  group_by(var) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            .groups = "drop") %>%
  mutate(data = "full") %>%
  bind_rows(d_scalings, .)

d_scalings <- d_moths %>%
  filter(!is.na(active_hours)) |> 
  select(LOC, J, M, A, active_hours) |> 
  distinct() |> 
 summarise(mean = mean(active_hours),
           sd = sd(active_hours)) %>%
  mutate(var = "active_hours",
         data = "full") %>%
  bind_rows(d_scalings, .)

d_scalings <- d_mod_ne %>%
  select(A, yday, T_2day, P_2day) %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  group_by(var) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            .groups = "drop") %>%
  mutate(data = "no_estimates") %>%
  bind_rows(d_scalings, .)

d_scalings <- d_moths %>%
  filter(!is.na(active_hours)) |> 
  select(LOC, J, M, A, active_hours) |> 
  distinct() %>%
  summarise(mean = mean(active_hours),
            sd = sd(active_hours)) %>%
  mutate(var = "active_hours",
         data = "no_estimates") %>%
  bind_rows(d_scalings, .)

d_scalings <- d_mod_lf %>%
  select(A, yday, T_2day, P_2day) %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  group_by(var) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            .groups = "drop") %>%
  mutate(data = "LF") %>%
  bind_rows(d_scalings, .)

d_scalings <- d_mod_p %>%
  select(A, yday, T_2day, P_2day) %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  group_by(var) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            .groups = "drop") %>%
  mutate(data = "p") %>%
  bind_rows(d_scalings, .)

d_scalings <- d_moths %>%
  filter(!is.na(active_hours)) |> 
  select(LOC, J, M, A, active_hours) |> 
  distinct() |> 
  summarise(mean = mean(active_hours),
            sd = sd(active_hours)) %>%
  mutate(var = "active_hours",
         data = "p") %>%
  bind_rows(d_scalings, .)

# ... scaled data ##############################################################
################################################################################.

sel_pars <- unique(d_scalings$var)

d_mod_z <- d_mod %>%
  mutate(across(c(!! enquo(sel_pars)),
                   ~ f_scale(data.frame(orig = ., var = cur_column()), "full")))

d_mod_ne_z <- d_mod_ne %>%
  mutate(across(c(!! enquo(sel_pars)),
                ~ f_scale(data.frame(orig = ., var = cur_column()), "no_estimates")))

d_mod_lf_z <- d_mod_lf %>%
  mutate(across(c(!! enquo(sel_pars)),
                ~ f_scale(data.frame(orig = ., var = cur_column()), "LF")))

d_mod_p_z <- d_mod_p %>%
  mutate(across(c(!! enquo(sel_pars)),
                ~ f_scale(data.frame(orig = ., var = cur_column()), "p")))

# ... datasets for traits ######################################################
################################################################################.

# Mass groups ------------------------------------------------------------------.

d_mod_mass <-
d_moths %>%
  mutate(visit_ID = paste(LOC, A, yday)) %>%
  left_join(d_mass, by = "Name_std",
            relationship = "many-to-one") |> 
  mutate(mass_sum = individualCount * mass) |> 
  left_join(d_traits, by = "Name_std") |> 
  group_by(LOC, visit_ID, A, yday, mass_cat) %>%
  summarise(abu_tot = sum(individualCount, na.rm = T),
            sric = length(unique(Name_std)),
            mass_tot = sum(mass_sum),
            .groups = "drop") %>%
  pivot_wider(everything(), names_from = "mass_cat",
              values_from = c("abu_tot", "sric", "mass_tot"),
              values_fill = 0, names_sep = ".") %>%
  pivot_longer(-c(LOC, visit_ID, A, yday),
               names_to = c(".value", "mass_cat"), names_sep = "\\.") %>%
  right_join(d_samplingdates %>%
               filter(n_zeronights < 10),  # exclude stretches of 10 an more zeronights
             by = c("LOC", "A", "yday")) %>%
  mutate(visit_ID = paste(LOC, A, yday)) %>%
  mutate_at(vars(abu_tot, sric, mass_tot), ~ ifelse(is.na(.), 0, .)) %>%
  mutate(Date_prev = Samplingdate - 1,
         visit_ID_prev = paste(LOC, year(Date_prev), yday(Date_prev)),
         sample_previous = ifelse(visit_ID_prev %in% visit_ID, "yes", "no")) %>%
  select(-visit_ID_prev) %>%
  left_join(d_weather, by = c("LOC", "Samplingdate")) %>%
  mutate(trap_ID_A = paste(LOC, A),
         trap_ID_A  = as.factor(trap_ID_A)) %>%
  left_join(d_samplings, by = c("LOC", "A"), relationship = "many-to-one") %>%
  mutate(night_ID = paste(Samplingdate, ifelse(is.na(samplingpair), 
                                               LOC, samplingpair), sep = "|"),
         sample_previous = as.factor(sample_previous),
         A_id = A) %>% # not to be transformed
  left_join(d_sites %>% select(LOC, height, height_cat), by = "LOC",
            relationship = "many-to-one")

d_mod_mass_z <- d_mod_mass %>%
  mutate(across(c(!! enquo(sel_pars)),
                ~ f_scale(data.frame(orig = ., var = cur_column()), "full")))

# T niche groups ---------------------------------------------------------------.

d_mod_Tavg <-
  d_moths %>%
  mutate(visit_ID = paste(LOC, A, yday)) %>%
  left_join(d_mass, by = "Name_std",
            relationship = "many-to-one") |> 
  mutate(mass_sum = individualCount * mass) |> 
  left_join(d_traits, by = "Name_std") |> 
  group_by(LOC, visit_ID, A, yday, Tavg_mean_cat) %>%
  summarise(abu_tot = sum(individualCount, na.rm = T),
            sric = length(unique(Name_std)),
            mass_tot = sum(mass_sum),
            .groups = "drop") %>%
  pivot_wider(everything(), names_from = "Tavg_mean_cat",
              values_from = c("abu_tot", "sric", "mass_tot"),
              values_fill = 0, names_sep = ".") %>%
  pivot_longer(-c(LOC, visit_ID, A, yday),
               names_to = c(".value", "Tavg_mean_cat"), names_sep = "\\.") %>%
  right_join(d_samplingdates %>%
               filter(n_zeronights < 10),  # exclude stretches of 10 an more zeronights
             by = c("LOC", "A", "yday")) %>%
  mutate(visit_ID = paste(LOC, A, yday)) %>%
  mutate_at(vars(abu_tot, sric, mass_tot), ~ ifelse(is.na(.), 0, .)) %>%
  mutate(Date_prev = Samplingdate - 1,
         visit_ID_prev = paste(LOC, year(Date_prev), yday(Date_prev)),
         sample_previous = ifelse(visit_ID_prev %in% visit_ID, "yes", "no")) %>%
  select(-visit_ID_prev) %>%
  left_join(d_weather, by = c("LOC", "Samplingdate")) %>%
  mutate(trap_ID_A = paste(LOC, A),
         trap_ID_A  = as.factor(trap_ID_A)) %>%
  left_join(d_samplings, by = c("LOC", "A"), relationship = "many-to-one") %>%
  mutate(night_ID = paste(Samplingdate, ifelse(is.na(samplingpair), 
                                               LOC, samplingpair), sep = "|"),
         sample_previous = as.factor(sample_previous),
         A_id = A) %>% # not to be transformed
  left_join(d_sites %>% select(LOC, height, height_cat), by = "LOC",
            relationship = "many-to-one")

d_mod_Tavg_z <- d_mod_Tavg %>%
  mutate(across(c(!! enquo(sel_pars)),
                ~ f_scale(data.frame(orig = ., var = cur_column()), "full")))


#  overwintering stage groups --------------------------------------------------.

d_mod_hib <-
  d_moths %>%
  mutate(visit_ID = paste(LOC, A, yday)) %>%
  left_join(d_mass, by = "Name_std",
            relationship = "many-to-one") |> 
  mutate(mass_sum = individualCount * mass) |> 
  left_join(d_traits, by = "Name_std") |> 
  group_by(LOC, visit_ID, A, yday, overwintering_stage) %>%
  summarise(abu_tot = sum(individualCount, na.rm = T),
            sric = length(unique(Name_std)),
            mass_tot = sum(mass_sum),
            .groups = "drop") %>%
  pivot_wider(everything(), names_from = "overwintering_stage",
              values_from = c("abu_tot", "sric", "mass_tot"),
              values_fill = 0, names_sep = ".") %>%
  pivot_longer(-c(LOC, visit_ID, A, yday),
               names_to = c(".value", "overwintering_stage"), names_sep = "\\.") %>%
  right_join(d_samplingdates %>%
               filter(n_zeronights < 10),  # exclude stretches of 10 an more zeronights
             by = c("LOC", "A", "yday")) %>%
  mutate(visit_ID = paste(LOC, A, yday)) %>%
  mutate_at(vars(abu_tot, sric, mass_tot), ~ ifelse(is.na(.), 0, .)) %>%
  mutate(Date_prev = Samplingdate - 1,
         visit_ID_prev = paste(LOC, year(Date_prev), yday(Date_prev)),
         sample_previous = ifelse(visit_ID_prev %in% visit_ID, "yes", "no")) %>%
  select(-visit_ID_prev) %>%
  left_join(d_weather, by = c("LOC", "Samplingdate")) %>%
  mutate(trap_ID_A = paste(LOC, A),
         trap_ID_A  = as.factor(trap_ID_A)) %>%
  left_join(d_samplings, by = c("LOC", "A"), relationship = "many-to-one") %>%
  mutate(night_ID = paste(Samplingdate, ifelse(is.na(samplingpair), 
                                               LOC, samplingpair), sep = "|"),
         sample_previous = as.factor(sample_previous),
         A_id = A) %>% # not to be transformed
  left_join(d_sites %>% select(LOC, height, height_cat), by = "LOC",
            relationship = "many-to-one")

d_mod_hib_z <- d_mod_hib %>%
  mutate(across(c(!! enquo(sel_pars)),
                ~ f_scale(data.frame(orig = ., var = cur_column()), "full")))

# food specialisation ----------------------------------------------------------.

d_mod_spec <-
  d_moths %>%
  mutate(visit_ID = paste(LOC, A, yday)) %>%
  left_join(d_mass, by = "Name_std",
            relationship = "many-to-one") |> 
  mutate(mass_sum = individualCount * mass) |> 
  left_join(d_traits, by = "Name_std") |> 
  group_by(LOC, visit_ID, A, yday, Spec) %>%
  summarise(abu_tot = sum(individualCount, na.rm = T),
            sric = length(unique(Name_std)),
            mass_tot = sum(mass_sum),
            .groups = "drop") %>%
  pivot_wider(everything(), names_from = "Spec",
              values_from = c("abu_tot", "sric", "mass_tot"),
              values_fill = 0, names_sep = ".") %>%
  pivot_longer(-c(LOC, visit_ID, A, yday),
               names_to = c(".value", "Spec"), names_sep = "\\.") %>%
  right_join(d_samplingdates %>%
               filter(n_zeronights < 10),  # exclude stretches of 10 an more zeronights
             by = c("LOC", "A", "yday")) %>%
  mutate(visit_ID = paste(LOC, A, yday)) %>%
  mutate_at(vars(abu_tot, sric, mass_tot), ~ ifelse(is.na(.), 0, .)) %>%
  mutate(Date_prev = Samplingdate - 1,
         visit_ID_prev = paste(LOC, year(Date_prev), yday(Date_prev)),
         sample_previous = ifelse(visit_ID_prev %in% visit_ID, "yes", "no")) %>%
  select(-visit_ID_prev) %>%
  left_join(d_weather, by = c("LOC", "Samplingdate")) %>%
  mutate(trap_ID_A = paste(LOC, A),
         trap_ID_A  = as.factor(trap_ID_A)) %>%
  left_join(d_samplings, by = c("LOC", "A"), relationship = "many-to-one") %>%
  mutate(night_ID = paste(Samplingdate, ifelse(is.na(samplingpair), 
                                               LOC, samplingpair), sep = "|"),
         sample_previous = as.factor(sample_previous),
         A_id = A) %>% # not to be transformed
  left_join(d_sites %>% select(LOC, height, height_cat), 
            by = "LOC", relationship = "many-to-one")

d_mod_spec_z <- d_mod_spec %>%
  mutate(across(c(!! enquo(sel_pars)),
                ~ f_scale(data.frame(orig = ., var = cur_column()), "full")))


