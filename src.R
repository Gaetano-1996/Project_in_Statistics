#==============================================================================#
############################### PROJECT FUNCTIONS ##############################
#==============================================================================#

#==== LIBRARY IMPORTS ======================================================####
library(tidyverse)
library(pcalg)
library(mice)
library(ggmice)
library(micd)
library(tpc)
library(visdat)
library(future)
library(parallel)
library(data.table)
library(haven)
library(skimr)

#==== UTILITY FUNCTIONS ====================================================####

create_directories = function() {
  dirs = c("plots2", "plots2/eda", "plots2/chains", "plots2/graphs", "summaries")
  sapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)
  cat("Directory setup complete.\n")
}

remove_incomplete_observations = function(data, 
                                           threshold = 0.25, 
                                           verbose = TRUE) {
  
  # Store original dimensions
  n_obs_original = nrow(data)
  n_vars = ncol(data)
  
  # Calculate proportion missing per observation
  miss_per_obs = rowSums(is.na(data)) / n_vars
  
  # Identify observations to remove
  to_remove = miss_per_obs > threshold
  
  # Create filtered dataset
  data_filtered = data[!to_remove, ]
  
  # Calculate statistics
  n_removed = sum(to_remove)
  pct_removed = (n_removed / n_obs_original) * 100
  
  # With multiple variables, threshold means >threshold_vars variables missing
  threshold_vars = ceiling(threshold * n_vars)
  
  if(verbose) {
    cat("=== OBSERVATION REMOVAL SUMMARY ===\n")
    cat(sprintf("Total variables: %d\n", n_vars))
    cat(sprintf("Threshold: %.1f%% (>%d variables missing)\n", 
                threshold * 100, threshold_vars))
    cat(sprintf("Original observations: %d\n", n_obs_original))
    cat(sprintf("Observations removed: %d (%.1f%%)\n", n_removed, pct_removed))
    cat(sprintf("Observations retained: %d (%.1f%%)\n", 
                nrow(data_filtered), 100 - pct_removed))
    cat("\n")
    
    # Distribution of missingness in removed observations
    if(n_removed > 0) {
      removed_miss = miss_per_obs[to_remove]
      cat("Missingness in removed observations:\n")
      cat(sprintf("  Min: %.1f%% (%d vars)\n", 
                  min(removed_miss) * 100, 
                  round(min(removed_miss) * n_vars)))
      cat(sprintf("  Mean: %.1f%% (%d vars)\n", 
                  mean(removed_miss) * 100, 
                  round(mean(removed_miss) * n_vars)))
      cat(sprintf("  Max: %.1f%% (%d vars)\n", 
                  max(removed_miss) * 100, 
                  round(max(removed_miss) * n_vars)))
    }
  }
  
  # Return filtered data
  return(data_filtered)
}

#=== DATA PREPARATION FUNCTIONS ============================================####

prepare_ggs_data = function(input_path,
                            output_path) {
  cat("=== PREPARING GGS DATA ===\n")
  
  # Set working directory for data files
  old_wd = getwd()
  setwd(input_path)
  
  # DATA IMPORT
  cat("Importing data files...\n")
  dw1 = fread("GGS_Wave1.csv")
  dw2 = fread("GGS_Wave2.csv")
  df.stata = read_dta("readyforMACHINE.dta")
  
  # DATA PREP
  cat("Processing Wave 1 data...\n")
  
  # FILTERING BY AGE
  dw1 = dw1[aage >= 18 & aage <= 45]
  
  # Removing problematic countries
  dw1 = dw1[!(acountry %in% c(17))]
  dw2 = dw2[!(bcountry %in% c(17))]
  
  # Removing useless birth years
  dw2 = dw2 %>% 
    select(-matches("^b253y_1[5-9]$")) 
  
  # Adjusting numbers for Czech Republic
  dw2 = dw2 %>% 
    mutate(bnkids = if_else(bcountry == "28", bnkids - 10, bnkids))
  
  # Creating unique ID
  dw1$ID = as.numeric(dw1$acountry) * 1000000000 + as.numeric(dw1$arid)
  dw2$ID = as.numeric(dw2$bcountry) * 1000000000 + as.numeric(dw2$brid)
  
  cat("Creating outcome variable...\n")
  
  # OUTCOME VARIABLE CREATION
  bcols = c("ID", "bnkids", grep("^b253y_", names(dw2), value = TRUE))
  dwu = merge(dw1[, c("ID","ayear",'amonth', "acountry",'ankids')], 
              dw2[, ..bcols], 
              by = "ID")
  
  by_cols = grep("^b253y_", names(dw2), value = TRUE)
  patterns = "^$|^\\.a$|^\\.b$|^\\.c$"
  dwu = dwu %>%
    mutate(
      across(all_of(by_cols),
             ~ replace(.x, str_detect(.x, patterns), NA)
      )
    )
  
  # OUTCOME DEFINITION
  dwu = dwu %>% 
    mutate(
      morekids = bnkids - ankids,
      birth_year = rowSums(!is.na(across(all_of(by_cols)))) > 0,
      birth_in_range = rowSums(
        across(all_of(by_cols), ~ . >= ayear & . <= ayear + 3), na.rm = TRUE),
      newchild = case_when(
        morekids > 0 & birth_in_range ~ 1,
        morekids <= 0 ~ 0,
        morekids > 0 & !birth_year ~ NA,
        TRUE ~ 0
      )
    )
  
  dwu = dwu %>% 
    mutate(
      new_child = factor(
        newchild,
        levels = c(0,1),
        labels = c("no","yes")
      )
    )
  
  # Drop NA in newchild
  dwu = dwu %>% 
    filter(!is.na(new_child))
  
  cat("Processing harmonization and variable creation...\n")
  
  # Clean missing patterns
  patterns = "^$|^\\.a$|^\\.b$|^\\.c$|^\\.d$"
  dw1 = dw1 %>% 
    mutate(across(everything(), ~ if_else(str_detect(., patterns), NA, .)))
  
  # HARMONIZATION AND VARIABLE CREATION (Full processing from data_prep.R)
  
  # country: country of interview 
  dw1[,
      country := factor(
        acountry,
        levels = c(11, 12, 13, 14, 15,16,18,19, 21,22,23, 25,26, 28,29),
        labels = c("Bulgaria", "Russia", "Georgia", "Germany",
                   "France","Hungary","Netherlands","Romania",
                   "Austria", "Estonia",'Belgium',"Lithuania",'Poland', "CzechRepubl",
                   'Sweden')
      )]
  
  # hh_disability: handicap or limited member of the hh
  discols = grep("^ahg9_", names(dw1), value = TRUE)
  dw1 = dw1 %>% 
    mutate(
      any1    = rowSums(across(all_of(discols), ~ . == 1), na.rm = TRUE) > 0,
      all_na  = rowSums(is.na(across(all_of(discols)))) == length(discols),
      dishh   = ifelse(all_na, NA_integer_, as.integer(any1))
    ) %>% 
    mutate(
      hh_disability = factor(dishh,
                             levels = c(0,1),
                             labels = c("no","yes"))
    )
  
  # migrant: born in country of interview
  dw1 = dw1 %>% 
    mutate(
      migrant = factor(
        a105,
        levels = c(1,2),
        labels = c("yes","no")
      )
    )
  
  # p_migrant: migrant status  partner
  dw1 = dw1 %>% 
    mutate(
      p_migrant = factor(
        a374a,
        levels = c(1,2),
        labels = c("yes","no")
      )
    )
  
  # nr_rooms: number of rooms of R's dwelling
  dw1 = dw1 %>% 
    mutate(
      nr_rooms_cat = cut(
        as.numeric(dw1$a119),
        breaks = c(1,2,3,4,5,Inf),
        right = F,
        include.lowest = T,
        labels = c("1","2","3","4",">=5"),
        ordered = T
      ),
      nr_rooms_num = as.numeric(dw1$a119))
  
  # dwell_owner: dwelling ownership
  dw1 = dw1 %>% 
    mutate(
      dwell_ownership_harm = factor(case_when(
        a122 %in% c("1202", "1203", "2", "1801") ~ "not owner",
        a122 %in% c("1204", "1") ~ "owner",
        TRUE ~ "other"
      ), levels = c("owner", "not owner", "other")),
      dwell_ownership = factor(
        as.numeric(a122),
        levels = c(1,2,3,4,1202,1203,1204,1801),
        labels = c("owner", "paying rent",
                   "accommodation rent-free", 'other',
                   "rent from firm or private person",
                   "rent and for communal services",
                   "pay only for communal services",
                   "not owner")
      )
    )
  
  # dwell_sat: dwelling satisfaction
  dw1 = dw1 %>% 
    mutate(a145 = if_else(as.numeric(a145)>10,NA,a145),
           dwell_sat = factor(
             as.numeric(a145),
             levels = seq(1,10),
             labels = c('1','2','3','4','5','6','7','8','9','10'),
             ordered = T
           ))
  
  # education: consolidated highest reached education level
  dw1 = dw1 %>% 
    mutate(
      aeduc2 = case_when(
        a148 == "1501" ~ 0,
        a148 == "1502" ~ 2,
        a148 %in% c("1503","1504","1505") ~ 3,
        a148 %in% c("1506","1507")  ~ 5,
        aeduc == "6" ~ 5,
        T ~ as.numeric(aeduc)
      ),
      education = cut(
        aeduc2,
        breaks = c(-Inf, 0, 1, 2, 3, 4, 5),
        right = T,
        include.lowest = T,
        labels = c("pre-primary education",
                   "primary level",
                   "lower secondary level",
                   "upper secondary level",
                   "post secondary non-tertiary",
                   "tertiary or higher"),
        ordered = T
      )
    )
  
  # help_childcare: regular help with childcare 
  dw1 = dw1 %>% 
    mutate(
      help_childcare = factor(
        a203a,
        levels = c(1,2),
        labels = c("yes","no")
      )
    )
  
  # gave_birth: have given birth  
  dw1 = dw1 %>% 
    mutate(
      gave_birth = factor(
        a209,
        levels = c(1,2),
        labels = c("yes","no")
      )
    )
  
  # dead_child: R mention at least 1 dead child
  deadchild_cols = grep("^a254_", names(dw1), value = TRUE)[1:13]
  
  dw1 = dw1 %>%
    mutate(
      n_na  = rowSums(is.na(across(all_of(deadchild_cols)))),
      has1  = rowSums(across(all_of(deadchild_cols), ~ .==1), na.rm=TRUE) > 0,
      all2  = rowSums(across(
        all_of(deadchild_cols), ~ .==2), na.rm=TRUE) == (length(deadchild_cols) - n_na),
      dead_child = as.factor(case_when(
        n_na  == length(deadchild_cols) ~ NA_character_,
        has1                   ~ "yes",
        all2                   ~ "no",
        TRUE                   ~ NA_character_
      ))
    ) %>%
    select(-n_na, -has1, -all2)
  
  # yr_coliving: year starting living together
  dw1 = dw1 %>% 
    mutate(
      yrs_coliv = as.numeric(ayear) - as.numeric(a301y),
      yrs_coliv = if_else(yrs_coliv <= 0, NA, yrs_coliv),
      yrs_coliv_cat = cut(
        yrs_coliv,
        breaks = quantile(yrs_coliv, na.rm = T),
        right = T,
        include.lowest = T,
        ordered = T
      ),
      yrs_coliv_num = as.numeric(yrs_coliv)
    )
  
  # elapsed_marriage: Time elapsed between start of relationship and marriage
  dw1 = dw1 %>% 
    mutate(
      elapsed_marriage =  as.numeric(a372bTdiff),
      elapsed_marriage_cat = cut(
        elapsed_marriage,
        breaks = c(-Inf, 0, 1, 2, 3,4, Inf),
        right = T,
        include.lowest = T,
        labels = c("0","1","2","3","4",">=5"),
        ordered = T
      ),
      elapsed_marriage_num = as.numeric(elapsed_marriage)
    )
  
  # hh_tasks_sat: satisfaction with tasks division in hh
  dw1 = dw1 %>% 
    mutate(
      hh_tasks_sat = factor(
        a402,
        levels = c(0,1,2,3,4,5,6,7,8,9,10),
        labels = c("0","1","2","3","4","5","6","7","8","9","10"),
        ordered = T
      )
    )
  
  # partner_sat: satisfaction with partner 
  dw1 = dw1 %>% 
    mutate(
      partner_sat_harm = factor(case_when(
        a407 %in% c('0','1','2','1805') ~ "not at all satisfied",
        a407 %in% c('3','4','1804') ~ "less satisfied",
        a407 %in% c('5','6','1803') ~ 'satisfied',
        a407 %in% c('7','8','1802') ~ 'mostly satisfied',
        a407 %in% c('9','10','1801') ~ 'completely satisfied',
      ), levels = c("not at all satisfied","less satisfied",'satisfied',
                    'mostly satisfied','completely satisfied'),
      ordered = T),
      partner_sat = factor(
        as.numeric(a407),
        levels = c(0:10,
                   1801:1805),
        labels = c('0','1','2','3','4','5','6','7','8','9','10',
                   'completely satisfied', 'mostly satisfied', 'satisfied',
                   'less satisfied', 'not at all satisfied'),
      )
    )
  
  # disability_parents: disability status of the parents
  dw1 = dw1 %>% 
    mutate(
      disability_mom = if_else(a518 == "1", 1, 0),
      disability_dad = if_else(a538 == "1", 1, 0),
      disability_sum = rowSums(across(c(disability_mom, disability_dad)), na.rm = TRUE),
      both_parents_na = is.na(disability_mom) & is.na(disability_dad),
      disability_parents = case_when(
        both_parents_na ~ NA_character_,
        disability_sum == 2 ~ "both",
        disability_sum == 1 ~ "one", 
        disability_sum == 0 ~ "none"
      ),
      disability_parents = factor(disability_parents, 
                                  levels = c("none", "one", "both"),
                                  labels = c("none", "one parent disabled", 
                                             "both parents disabled")
    )) %>%
    select(-disability_sum, -both_parents_na)
  
  # divorced_parents: whether R's parents are divorced
  dw1 = dw1 %>% 
    mutate(
      a550 = case_when(
        a550 == "1" ~ 1,
        a550 == "2" ~ 2,
        a550 %in% c("3","4","5") ~ 3,
        T ~ as.numeric(a550)
      ),
      divorced_parents = factor(
        a550,
        levels = c(1,2,3),
        labels = c("yes","no", "unknown")
      )
    )
  
  # nr_brothers: number of brothers
  dw1 = dw1 %>% 
    mutate(
      nr_bros_cat = cut(
        as.numeric(a5106b_b),
        breaks = c(-Inf, 0, 1, 2, 3, 4, Inf),
        right = T,
        include.lowest = T,
        labels = c("0","1","2","3","4",">=5"),
        ordered = T
      ),
      nr_bros_num = as.numeric(a5106b_b)
    )
  
  # nr_sisters: number of sisters
  dw1 = dw1 %>% 
    mutate(
      nr_sis_cat = cut(
        as.numeric(a5106b_s),
        breaks = c(-Inf, 0, 1, 2, 3, 4, Inf),
        right = T,
        include.lowest = T,
        labels = c("0","1","2","3","4",">=5"),
        ordered = T
      ),
      nr_sis_num = as.numeric(a5106b_s)
    )
  
  # nr_siblings: number of siblings
  dw1 = dw1 %>%
    mutate(
      a5106b_b = as.numeric(a5106b_b),
      a5106b_s = as.numeric(a5106b_s),
      sibling_sum = rowSums(across(c(a5106b_b, a5106b_s)), na.rm = TRUE),
      both_na     = is.na(a5106b_b) & is.na(a5106b_s),
      nr_siblings_num = if_else(both_na, NA_integer_, sibling_sum)
    ) %>%
    mutate(
      nr_siblings_cat = cut(
        nr_siblings_num,
        breaks       = c(-Inf, 0, 1, 2, 3, 4, Inf),
        right        = TRUE,
        include.lowest = TRUE,
        labels       = c("0","1","2","3","4",">=5"),
        ordered      = TRUE
      )
    ) %>% select(-sibling_sum, -both_na)
  
  # granpas_alive: number of grandparents alive
  dw1 = dw1 %>% 
    mutate(
      granpas_alive = as.factor(if_else(
        as.numeric(a5107) > 0, "yes","no"
      ))
    )
  
  # current_preg: R's current pregnancy status
  dw1 = dw1 %>% 
    mutate(
      current_preg = factor(
        a602,
        levels = c(1,2,3),
        labels = c("yes","no","maybe")
      )
    )
  
  # int_now: R's intention to have a child now
  dw1 = dw1 %>% 
    mutate(
      int_now = factor(case_when(
        a611 == "1" ~ "yes",
        is.na(a611) ~ NA,
        T ~ "no"
      ), 
      levels = c("yes","no"),
      labels = c("yes","no")
      )
    )
  
  # p_int_now: partner's intention to have a child now
  dw1 = dw1 %>% 
    mutate(
      p_int_now = factor(
        a615,
        levels = c(1,2,3),
        labels = c("yes","no", "not sure")
      ) 
    )
  
  # possibility_child: physically possible to have another baby?
  dw1 = dw1 %>% 
    mutate(
      possibility_child_harm = as.factor(
        case_when(
          a612 %in% c(4, 1602)                ~ "yes",
          a612 %in% c(3, 2)                   ~ "not sure",
          a612 %in% c(1, 1501, 1601)          ~ "no",
          TRUE                                ~ NA_character_
        )
      ),
      possibility_child = factor(
        as.numeric(a612),
        levels = c(1,2,3,4,1501,1601,1602),
        labels = c("definitely not",
                   "probably not",
                   "probably yes",
                   "definitely yes",
                   "yes - health risk",
                   "no",
                   "yes")
      )
    )
  
  # p_possibility_child: physically possible to have another baby (partner)?
  labs = names(attr(df.stata$var_616, "labels")[1:4])
  dw1 = dw1 %>% 
    mutate(
      p_possibility_child = factor(
        a616,
        levels = c(1,2,3,4),
        labels = labs
      )
    )
  
  # int_couple: couple intentions
  labs = names(attr(df.stata$var_ertintent, "labels"))
  dw1 = dw1 %>% 
    mutate(
      int_couple = factor(
        fertintent,
        levels = c(1,2,3,4),
        labels = labs
      )
    )
  
  # int_3yr: intention to have another baby in the future
  dw1 = dw1 %>% 
    mutate(
      int_3yrs_harm = case_when(
        a622 %in% c("1","1602","1801") ~ "yes",
        a622 %in% c("2","3","4","5","1601","1802","1803") ~ "no"
      ),
      int_3yrs_harm = factor(
        int_3yrs_harm,
        levels = c("yes","no"),
        labels = c("yes","no")
      ),
      int_3yrs = factor(
        a622,
        levels = c(1,2,3,4,1601,1602,1801,1802,1803),
        labels = c("definitely not",
                   "probably not",
                   "probably yes",
                   "definitely yes",
                   "no",
                   "yes",
                   "want to have a child within 3 years",
                   "wants to have a child in > 3 years",
                   "does not want to have (more) children"
        )
      )
    )
  
  # int_atall: intention to have a child at all
  labs = names(attr(df.stata$var_624, "labels")[1:4])
  dw1 = dw1 %>% 
    mutate(
      int_atall_harm = if_else(a624 == "1601", 1, 
                               if_else(a624 == "1602", 4, as.numeric(a624))),
      int_atall_harm = factor(
        int_atall_harm,
        levels = c(1,2,3,4),
        labels = labs,
        ordered = T
      ),
      int_atall = factor(
        as.numeric(a624),
        levels = c(1,2,3,4,1601,1602),
        labels = c("definitely not",
                   "probably not",
                   "probably yes",
                   "definitely yes",
                   "no",
                   "yes"
        )
      )
    )
  
  # child_sex: preferred child sex
  labs = names(attr(df.stata$var_625, "labels")[1:3])
  dw1 = dw1 %>% 
    mutate(
      child_sex = factor(
        as.numeric(a625),
        levels = c(1,2,3),
        labels = labs
      )
    )
  
  # how_many_more: how many more children
  dw1 = dw1 %>% 
    mutate(
      a626 = if_else(as.numeric(a626) > 23, NA, a626),
      how_many_more = cut(
        as.numeric(a626),
        breaks = c(-Inf, 0, 1, 2, Inf),
        right = T,
        include.lowest = T,
        labels = c("0","1","2",">=3"),
        ordered = T
      ))
  
  # general_health: perceived general health status
  labs = names(attr(df.stata$var_701, "labels")[1:5])
  
  dw1 = dw1 %>% 
    mutate(
      general_health = factor(
        a701,
        levels = c(1:5),
        labels = labs,
        ordered = T
      )
    )
  
  # chronic: chronic or longstanding illness
  dw1 = dw1 %>% 
    mutate(
      chronic = factor(
        a702a,
        levels = 1:2,
        labels = c("yes","no")
      )
    )
  
  # health_limitation: Health related limitation or disability
  dw1 = dw1 %>% 
    mutate(
      health_limitation = factor(
        a703a,
        levels = 1:2,
        labels = c("yes","no")
      )
    )
  
  # work_type: part/full time 
  dw1 = dw1 %>% 
    mutate(
      work_type = factor(
        a834,
        levels = 1:2,
        labels = c("full-time","part-time")
      )
    )
  
  # p_work_type: partner's work type
  dw1 = dw1 %>% 
    mutate(
      p_work_type = factor(
        a922,
        levels = 1:2,
        labels = c("full-time","part-time")
      )
    )
  
  # work_sat: work satisfaction
  dw1 = dw1 %>% 
    mutate(
      work_sat = factor(
        as.numeric(a839),
        levels = seq(1,10),
        labels = c('1','2','3','4','5','6','7','8','9','10'),
        ordered = T
      )
    ) 
  
  # contract_type: type of work contract
  dw1 = dw1 %>% 
    mutate(
      contract_type = factor(
        as.numeric(a845),
        levels = seq(1:4),
        labels = c('permanent','fixed-term','temporary','no written contract')
      )
    ) 
  
  # add_income: additional income
  dw1 = dw1 %>% 
    mutate(
      add_income = factor(
        as.numeric(a860),
        levels = seq(1:2),
        labels = c('yes','no')
      )
    ) 
  
  # financial: financial situation
  dw1 = dw1 %>% 
    mutate(
      financial = factor(case_when(
        a1002 %in% c("1","1901") ~ "very difficult",
        a1002 %in% c("2","1902") ~ "difficult",
        a1002 %in% c("3","1903") ~ "somewhat difficult",
        a1002 %in% c("4","1904","1905") ~ "fairly easy",
        a1002 %in% c("5","1906") ~ "easy",
        a1002 %in% c("6","1907") ~ "very easy",
      ), 
      levels = c("very difficult","difficult",
                 "somewhat difficult","fairly easy",
                 "easy","very easy"),
      labels = c("very difficult","difficult",
                 "somewhat difficult","fairly easy",
                 "easy","very easy"),
      ordered = T)
    )
  
  # Partner variables: sex, education, age
  part_cols = c("femage",'maleage','femeduc','maleeduc')
  
  c1 = (!is.na(dw1$femage) | !is.na(dw1$femeduc)) & 
    (is.na(dw1$maleage) & is.na(dw1$maleeduc))
  
  c2 = is.na(dw1$femage) & is.na(dw1$femeduc) & 
    is.na(dw1$maleage) & is.na(dw1$maleeduc)
  
  dw1 = dw1 %>% 
    mutate(
      p_sex = if_else(c1, "female", if_else(c2, NA_character_, "male")),
      p_education = if_else(p_sex=="female", femeduc, maleeduc, missing = NA_character_),
      p_age_num = as.numeric(if_else(p_sex=="female", femage, maleage, missing = NA_character_))
    ) %>% 
    mutate(
      p_sex = as.factor(p_sex),
      p_education = factor(
        case_when(
          p_education == "0" ~ "pre-primary education",
          p_education == "1" ~ "primary level",
          p_education == "2" ~ "lower secondary level",
          p_education == "3" ~ "upper secondary level",
          p_education == "4" ~ "post secondary non-tertiary",
          p_education %in% c("5", "6") ~ "tertiary or higher",
          TRUE ~ NA_character_
        ),
        ordered = TRUE
      ),
      p_age_cat = cut(
        as.numeric(p_age_num),
        breaks = c(0, 26, 33, 40, Inf),
        right = TRUE,
        include.lowest = TRUE,
        labels = c("<=26",
                   "(26-33]",
                   "(33-40]",
                   ">40"),  
        ordered = TRUE
      )
    ) %>% 
    select(-all_of(part_cols))
  
  # sex: R's sex
  dw1 = dw1 %>% 
    mutate(sex =factor(
      asex,
      levels = c(1,2),
      labels = c("male","female")
    ))
  
  # age: R's age
  dw1 = dw1 %>% 
    mutate(
      age_num = as.numeric(aage),
      age_cat = cut(
        as.numeric(aage),
        breaks = c(-Inf, 18, 26,33,40,Inf),
        right = TRUE,
        include.lowest = TRUE,
        labels = c("<=18",
                   "(18-26]",
                   "(26-33]",
                   "(33-40]",
                   ">40"),
        ordered = TRUE
      )
    )
  
  # hh_type: type of household
  labs = names(attr(df.stata$var_hhtype, "labels"))[1:9]
  dw1 = dw1 %>% 
    mutate(
      hh_type = factor(
        ahhtype,
        levels = c(1,2,3,4,5,6,7,8,9),
        labels = labs
      )
    )
  
  # hh_size: number of people in the household
  dw1 = dw1 %>% 
    mutate(
      hh_size_num = as.numeric(ahhsize),
      hh_size_cat = cut(
        as.numeric(ahhsize),
        breaks = c(-Inf, 1, 2, 3, 4, Inf),
        right = TRUE,
        include.lowest = TRUE,
        labels = c("1","2","3","4",">=5"),
        ordered = TRUE
      )
    )
  
  # activity_status: R's activity status
  labs = names(attr(df.stata$var_actstat, "labels"))[1:10]
  dw1 = dw1 %>% 
    mutate(
      activity_status = factor(
        aactstat,
        levels = c(1,2,3,4,5,6,7,8,9,10),
        labels = labs
      )
    )
  
  # marital_status: marital status of R
  labs = names(attr(df.stata$var_marstat, "labels"))[1:4]
  dw1 = dw1 %>% 
    mutate(
      marital_status = factor(
        amarstat,
        levels = c(1,2,3,4),
        labels = labs
      )
    )
  
  # nr_kids: number of R's kids
  dw1 = dw1 %>% 
    mutate(
      nr_kids_num = as.numeric(ankids),
      nr_kids_cat = cut(
        as.numeric(ankids),
        breaks = c(-Inf, 0, 1, 2, Inf),
        right = TRUE,
        include.lowest = TRUE,
        labels = c("0","1","2",">=3"),
        ordered = TRUE
      )
    )
  
  # nr_partners: number of R's partners
  dw1 = dw1 %>% 
    mutate(
      nr_partners_num = as.numeric(anpartner),
      nr_partners_cat = cut(
        as.numeric(anpartner),
        breaks =  c(-Inf, 0, 1, Inf),
        right = TRUE,
        include.lowest = TRUE,
        labels = c("0","1",">=2"),
        ordered = TRUE
      )
    )
  
  # partner_status: status of partner
  labs = names(attr(df.stata$var_parstat, "labels"))[1:3]
  dw1 = dw1 %>% 
    mutate(
      partner_status = factor(
        aparstat,
        levels = c(1,2,3),
        labels = labs
      )
    )
  
  # sett_type: settlement type
  dw1 = dw1 %>% 
    mutate(
      sett_type_harm = case_when(
        atype %in% c("2","1201","1202","1203","1301",
                    "1401","1402","1403","1404",
                    "1405","1406","1502","1503",
                    "1504","1801","1802","1803",
                    "2504","2505") ~ "urban",
        atype %in% c("1","1204","1407","1408","1409",
                    "1410","1501","1505",
                    "1804","1805","2501","2502",
                    "2503") ~ "rural",
        T~NA
      ),
      sett_type_harm = factor(
        sett_type_harm,
        levels = c("urban","rural"),
        labels = c("urban","rural")
      )
      )
  
  # nr_dissol: number of R's dissolved partnerships
  dw1 = dw1 %>% 
    mutate(
      nr_dissol_num = as.numeric(numdissol),
      nr_dissol_cat = cut(
        as.numeric(numdissol),
        breaks = c(-Inf, 0,1, Inf),
        right = T,
        include.lowest = T,
        labels = c("0","1",">=2"),
        ordered = T
      )
    )
  
  # nr_marriage: number of R's marriages
  dw1 = dw1 %>% 
    mutate(
      nr_marriage_num = as.numeric(nummarriage),
      nr_marriage_cat = cut(
        as.numeric(nummarriage),
        breaks = c(-Inf, 0, 1, Inf),
        right = T,
        include.lowest = T,
        labels = c("0","1",">=2"),
        ordered = T
      )
    )
  
  # child_prev: any children from previous partnerships?
  dw1 = dw1 %>% 
    mutate(
      child_prev = factor(
        childprevp,
        levels = c(1,2),
        labels = c("yes","no")
      )
    )
  
  # age_youngest: age of youngest child
  dw1 = dw1 %>% 
    mutate(
      age_youngest_num = as.numeric(ageyoungest),
      age_youngest_cat = cut(
        ageyoungest,
        breaks =quantile(as.numeric(ageyoungest),na.rm =T),
        include.lowest = T,
        ordered = T)
    )
  
  #age_oldest: age of oldest child
  dw1 = dw1 %>% 
    mutate(
      age_oldest_num = as.numeric(ageoldest),
      age_oldest_cat = cut(
        ageoldest,
        breaks =quantile(as.numeric(ageoldest),na.rm =T),
        include.lowest = T,
        ordered = T)
    )
  
  # cores_child: coresident children
  dw1 = dw1 %>% 
    mutate(
      cores_child = factor(
        coreschild,
        levels = c(1,2),
        labels = c("yes","no")
      )
    )
  
  # cores_parent: coresident parents
  dw1 = dw1 %>% 
    mutate(
      cores_parent = factor(
        coresparen,
        levels = c(1,2),
        labels = c("yes","no")
      )
    )
  
  # cores_grandp: coresident grandparents
  dw1 = dw1 %>% 
    mutate(
      cores_grandp = factor(
        coresgrandp,
        levels = c(1,2),
        labels = c("yes","no")
      )
    )
  
  # cores_sibl: coresident siblings
  dw1 = dw1 %>% 
    mutate(
      cores_sibl = factor(
        coressibl,
        levels = c(1,2),
        labels = c("yes","no")
      )
    )
  
  # nonres_child: non-resident children
  dw1 = dw1 %>% 
    mutate(
      nonres_child = factor(
        if_else(numnonres > 0, 1, 0),
        levels  = c(1, 0),
        labels  = c("yes", "no")
      )
    )
  
  # Partner variables handling
  partner_vars = c("p_migrant", "partner_sat_harm", "partner_sat", "p_int_now", 
                   "int_couple", "p_possibility_child", "p_work_type", "p_sex",
                   "p_education", "p_age_cat", "elapsed_marriage_cat", "yrs_coliv_cat")
  
  dw1 = dw1 %>%
    mutate(
      across(
        all_of(partner_vars),
        ~ {
          orig_levels = levels(.)
          is_ordered = is.ordered(.)
          new_levels = if("no partner" %in% orig_levels) orig_levels else c(orig_levels, "no partner")
          result = case_when(
            partner_status == "no partner" ~ "no partner",
            TRUE ~ as.character(.)
          )
          if(is_ordered) {
            ordered(result, levels = new_levels)
          } else {
            factor(result, levels = new_levels)
          }
        }
      )
    )
  
  # Merging intentions
  dw1 = dw1 %>%
    mutate(
      intentions = case_when(
        is.na(int_now) & is.na(int_atall_harm) & is.na(int_3yrs_harm) ~ NA_character_,
        int_now == "yes" | 
          int_atall_harm %in% c("probably yes", "definitely yes") | 
          int_3yrs_harm == "yes" ~ "yes",
        TRUE ~ "no"
      ),
      intentions = factor(intentions, levels = c("yes", "no"))
    )
  
  # Final data merge and selection
  df = merge(dw1,
             dwu[,c("ID", "new_child")],
             by = "ID")
  
  df = df %>%
    filter(!is.na(new_child))
  
  # Clean unused levels
  df = df %>% 
    mutate(across(where(is.factor), fct_drop))
  
  # DATASET DEFINITION
  df_mix = df %>% 
    filter(!country %in% c("Hungary","Bulgaria","CzechRepubl")) %>% 
    select(
      new_child, 
      country, 
      sex, 
      partner_sat_harm,
      migrant, 
      education, 
      age_num, 
      marital_status,
      partner_status,
      nr_kids_num, 
      gave_birth, 
      dead_child,
      age_youngest_num,
      activity_status,
      dwell_ownership_harm, 
      general_health, 
      hh_type, 
      sett_type_harm,
      intentions
    ) %>% 
    rename(
      partner_sat = partner_sat_harm,
      age = age_num,
      nr_kids = nr_kids_num,
      age_youngest = age_youngest_num,
      dwell_ownership = dwell_ownership_harm,
      sett_type = sett_type_harm
    ) 
  
  df_mix = df_mix %>% 
    mutate(across(where(is.factor), fct_drop))
  
  # Apply missing data filtering (as done in eda.R)
  cat("Filtering observations with >25% missing values...\n")
  df_final = remove_incomplete_observations(df_mix, threshold = 0.25, verbose = TRUE)
  
  # Save datasets
  setwd(old_wd)
  cat("Saving datasets...\n")
  saveRDS(df_mix, file = paste0(output_path,"/df_mix.rds"))
  saveRDS(df_final, file =  paste0(output_path,"/df_final.rds"))
  
  cat("Data preparation completed successfully.\n")
  cat("Datasets saved: df_mix.rds, df_final.rds\n")
  cat("Final dataset dimensions:", nrow(df_final), "x", ncol(df_final), "\n")
  
  return(df_final)
}

#=== EDA FUNCTIONS =========================================================####

perform_eda = function(data, eda_dir = "plots/eda", summaries_dir = "summaries") {
  cat("=== PERFORMING EXPLORATORY DATA ANALYSIS ===\n")
  
  # Overall summary
  cat("Creating overall data summary...\n")
  skim_na = skim_with(
    factor = sfl(
      prop_na = ~ round(mean(is.na(.)), 3),
      prop_0 = ~ round(mean(. == "no", na.rm = TRUE), 2),
      prop_1 = ~ round(mean(. == "yes", na.rm = TRUE), 2),
      miss_40 = ~ ifelse(mean(is.na(.)) > 0.40, T, F),
      miss_50 = ~ ifelse(mean(is.na(.)) > 0.50, T, F),
      miss_70 = ~ ifelse(mean(is.na(.)) > 0.70, T, F),
      miss_90 = ~ ifelse(mean(is.na(.)) > 0.90, T, F)
    ),
    numeric = sfl(
      prop_na = ~ round(mean(is.na(.)), 3),
      miss_50 = ~ ifelse(mean(is.na(.)) > 0.50, T, F),
      miss_70 = ~ ifelse(mean(is.na(.)) > 0.70, T, F),
      miss_90 = ~ ifelse(mean(is.na(.)) > 0.90, T, F)
    ),
    append = T
  )
  
  overall_summary = skim_na(data)
  
  # Missing data patterns (whole dataset)
  cat("Analyzing missing data patterns...\n")
  missing_plot = vis_miss(data,
                           warn_large_data = F,
                           sort_miss = T) + 
    ggtitle("Missing Data Patterns") + 
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 20, r = 50, b = 40, l = 20, unit = "pt"),
      axis.text.x = element_text(size = 8)
    )
  
  # Save overall plots
  ggsave(file.path(eda_dir, "missing_patterns.png"), missing_plot, 
         width = 12, height = 9, dpi = 300)
  
  # Missing data patterns (by country)
  missing_plot_country = vis_miss(data, 
                                  warn_large_data = F,
                                  sort_miss = T,
                                  facet = country) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(size = 14)
    )
  
  ggsave(file.path(eda_dir, "missing_patterns_bycountry.png"), missing_plot_country, 
         width = 12, height = 9, dpi = 300)
  
  # Outcome distribution
  outcome_plot = data %>%
    ggplot(aes(x = new_child, fill = new_child)) +
    geom_bar() +
    geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5) +
    labs(title = "Distribution of New Child Outcome",
         x = "New Child", y = "Count") +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  ggsave(file.path(eda_dir, "outcome_distribution.png"), outcome_plot,
         width = 8, height = 6, dpi = 300)
  
  # Save overall summary
  write_csv(as.data.frame(overall_summary), file.path(summaries_dir, "overall_summary.csv"))
  
  cat("EDA completed successfully.\n")
  
  return(list(
    overall_summary = overall_summary,
    missing_plot = missing_plot,
    outcome_plot = outcome_plot
  ))
}

generate_country_summaries = function(data, eda_dir = "plots/eda", summaries_dir = "summaries") {
  cat("=== GENERATING COUNTRY-SPECIFIC SUMMARIES ===\n")
  
  country_stats = list()
  
  for(country in unique(data$country)) {
    cat("Processing country:", country, "\n")
    
    country_data = data %>% filter(country == !!country)
    
    # Country-specific summary
    country_summary = skim(country_data)
    
    # Outcome distribution by country
    outcome_country = country_data %>%
      ggplot(aes(x = new_child, fill = new_child)) +
      geom_bar() +
      geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5) +
      labs(title = paste("New Child Distribution -", country),
           x = "New Child", y = "Count") +
      theme_minimal() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    # Save country-specific files
    ggsave(file.path(eda_dir, paste0("outcome_", country, ".png")), outcome_country,
           width = 8, height = 6, dpi = 300)
    
    write_csv(as.data.frame(country_summary), 
             file.path(summaries_dir, paste0("summary_", country, ".csv")))
    
    country_stats[[country]] = list(
      summary = country_summary,
      plot = outcome_country,
      n_obs = nrow(country_data)
    )
  }
  
  # Overall country comparison
  country_comparison = data %>%
    group_by(country, new_child) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(country) %>%
    mutate(prop = count / sum(count)) %>%
    filter(new_child == "yes") %>%
    ggplot(aes(x = reorder(country, prop), y = prop)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1)), 
              hjust = -0.1, size = 3) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Proportion of New Children by Country",
         x = "Country", y = "Proportion") +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 20, r = 50, b = 20, l = 20, unit = "pt")
    )
  
  ggsave(file.path(eda_dir, "country_comparison.png"), country_comparison,
         width = 10, height = 6, dpi = 300)
  
  cat("Country summaries completed successfully.\n")
  
  return(country_stats)
}

# === TIER CREATION FUNCTIONS ==============================================####

create_tiers_vector = function(df_names) {
  
  # Define tier variables (from main.R)
  t0 = c("sex", "age")
  t1 = c("education", "migrant")
  t2 = c("gave_birth", 'dead_child')
  t3 = c('nr_kids', 'age_youngest',"general_health", 'marital_status', 'partner_status', 'activity_status',
         'dwell_ownership', 'sett_type', 'partner_sat', 'hh_type')
  t4 = c("intentions")
  t5 = c("new_child")
  
  # Initialize vector with NA values
  tiers_vector = rep(NA, length(df_names))
  names(tiers_vector) = df_names
  
  # Assign tier numbers to matching variables
  tiers_vector[names(tiers_vector) %in% t0] = 0
  tiers_vector[names(tiers_vector) %in% t1] = 1
  tiers_vector[names(tiers_vector) %in% t2] = 2  
  tiers_vector[names(tiers_vector) %in% t3] = 3
  tiers_vector[names(tiers_vector) %in% t4] = 4
  tiers_vector[names(tiers_vector) %in% t5] = 5
  
  # Check for unmatched variables
  unmatched = names(tiers_vector)[is.na(tiers_vector)]
  if(length(unmatched) > 0) {
    warning("Unmatched variables found: ", paste(unmatched, collapse = ", "))
  }
  
  # Remove unmatched variables (NA values)
  tiers_vector = tiers_vector[!is.na(tiers_vector)]
  
  return(tiers_vector)
}

# === DISCOVERY FUNCTIONS ==================================================####

ggs_discovery = function(
    data,
    tiers = NULL,
    alpha = 0.1,
    context.var = NULL,
    imput.meth = "rf",
    maxit = 10,
    m = 25,
    mincor = 0.2,
    n.core = detectCores(),
    plot.graph = F,
    save.graph = T,
    plot.chains = T,
    save.chains = T,
    neighborhood = c(1,2),
    verbose = 1,
    chains_path = "plots/chains/",
    graphs_path = "plots/graphs/",
    seed = 123){
  
  cat("==== MI - tPC ====\n")
  cat("\n__ Imputation Params ___\n",
      "N. imputed dataset:",m,"\n",
      "Max N. iteration MICE:",maxit,"\n",
      "Imputation method:",imput.meth,"\n",
      "mincor for predictors selection:",mincor,"\n")
  cat("\n__ tPC Params ___\n",
      "Sig. lvl:",alpha,"\n",
      "N. tiers:",length(unique(tiers)),"\n\n")
  
  results = list()
  
  for (i in unique(data$country)){
    cat("=== Country: ",i,"===\n")
    # subset data
    df = data %>% 
      filter(country == i) %>% 
      select(-country)
    
    cat("N. observations:",nrow(df),"\n")
    
    # predictors for imputation model
    pred = quickpred(
      df,
      mincor = mincor
    )
    
    # multiple imputation 
    n.core = min(detectCores()-1,n.core)
    cat("Running mice...\n")
    mi_res = futuremice(
      df,
      m=m,
      maxit = maxit,
      predictorMatrix = pred,
      method=imput.meth,
      parallelseed = seed,
      n.core = n.core,
      ntree=100
    )
    
    p_chains = plot_trace(mi_res) +
      ggtitle(
        paste0(
          "MICE - m:",m,
          ", maxit:",maxit,
          ", mincor:",mincor,
          ", method:",imput.meth)) +
      theme(legend.position = "none")
    
    if (plot.chains){
      plot(p_chains)
    }
    
    if (save.chains){
      ggsave(paste0(chains_path,"/mice-chains-",i,".png"),
             p_chains, width = 8, height = 10, dpi = 300)
    }
    
    mi_data = complete(mi_res, action = "all")
    cat("Running tPC...\n")
    disco = tpc::tpc(suffStat=mi_data, 
                     indepTest=micd::mixMItest, 
                     skel.method = "stable.parallel",
                     labels=colnames(df),
                     alpha=alpha, 
                     tiers=tiers, 
                     context.all = context.var,
                     verbose=if_else(verbose>1,T,F),
                     numCores = n.core)
    
    if (plot.graph){
      plot(disco, main = "")
    }
    
    # Save full causal graph (new feature)
    if (save.graph) {
      # Save standard plot
      full_graph_filename = paste0(graphs_path,"/MI-tPC-",i,"-full.png")
      png(full_graph_filename, width = 2000, height = 1500, res = 500)
      plot(disco, main = paste0("MI-tPC-",i))
      dev.off()
    }
    
    # Save neighborhood graphs
    if (save.graph & !is.null(neighborhood)){
      for (j in neighborhood){
        
        disco_graph = as(disco@graph, "graphNEL")
        disco_graph@nodes = seq(1:ncol(df))
        names(disco_graph@edgeL) = seq(1:ncol(df))
        var_names = names(df)
        node_labels = setNames(var_names, as.character(1:ncol(df)))
        nAttrs = list(label = node_labels)
        
        filename= paste0(graphs_path,paste0("MI-tPC-",i,"-",j,".png"))
        if (j == 1){
          res = 300
        } else{
          res = 500
        }
        png(filename, width = 2000, height = 1500, res = res)
        plotSG(graphObj = disco_graph,
               y = 1,
               dist = j,
               directed = TRUE,
               nodeAttrs = nAttrs,
               main = "")
        dev.off()
      }
    }
    
    results[[i]] = list(
      discovery = disco,
      imputation = mi_res,
      chains_plot = p_chains
    )
  }
  
  cat("Discovery analysis completed successfully.\n")
  return(results)
}