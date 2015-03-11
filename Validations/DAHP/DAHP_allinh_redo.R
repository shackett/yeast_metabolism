setwd("~/Desktop/Rabinowitz/FBA_SRH/Validations/DAHP")

options(stringsAsFactors = F)

library(reshape2)
library(data.table)
library(scales)
library(dplyr)
library(ggplot2)
library(grid)

#library(devtools)
#install_github("dgrtwo/broom")
library(broom)
library(tidyr)

ARO3 <- read.delim('DAHP_allInh_ARO3.txt', fill = NA, header = F) 
ARO4 <- read.delim('DAHP_allInh_ARO4.txt', fill = NA, header = F) 

DAHP_plate_names <- c("ARO3", "ARO4")

DAHP_list <- list()
for(DAHP_plate in DAHP_plate_names){
  
  kineticData <- get(DAHP_plate)
  
  plate_starts <- grep('^Read', kineticData[,1])
  plates <- list()
  
  for(a_plate in 1:length(plate_starts)){
    
    if(a_plate == length(plate_starts)){
      plate_data <- kineticData[plate_starts[a_plate]:nrow(kineticData),]
    }else{
      plate_data <- kineticData[plate_starts[a_plate]:(plate_starts[a_plate + 1]-1),]
    }
    
    header_rows <- !(apply(!is.na(plate_data), 1, sum) > 2)
    
    plate_matrix <- plate_data[!header_rows,1:13]
    plate_matrix[1,1] <- "Row"
    colnames(plate_matrix) <- plate_matrix[1,]; plate_matrix <- plate_matrix[-1,]
    
    if(nrow(plate_matrix) == 0){next}
    
    plate_summary <- melt(plate_matrix, variable.name = "Col", value.name = "OD", id.vars = "Row")
    plate_summary <- data.frame(plate_ID = last(plate_data[header_rows,1]), plate_summary)
    
    plates[[a_plate]] <- plate_summary
  }
  
  plates_aggregate <- do.call("rbind", plates)
  plates_aggregate <- plates_aggregate[!is.na(plates_aggregate$OD),]
  
  DAHP_list[[DAHP_plate]] <- plates_aggregate
  
}

DAHP_plates <- tbl_df(melt(DAHP_list, value.name = "Absorbance"))
DAHP_plates <- DAHP_plates %>% dplyr::select(plate_ID, Row, Col, Absorbance, plate = L1)

plate_IDs <- DAHP_plates %>% group_by(plate_ID, plate) %>% summarize(type = "Kinetic", time = NA, cumTime = NA)
plate_IDs$type[grep('Pathlength', plate_IDs$plate_ID)] <- "PL"
plate_IDs$time[plate_IDs$type == "Kinetic"] <- regmatches(plate_IDs$plate_ID, regexpr('[0-9]{1}:[0-9]{2}:[0-9]{2}', plate_IDs$plate_ID))

formatted_time <- as.POSIXlt(plate_IDs$time[plate_IDs$type == "Kinetic"], format = "%H:%M:%S")

plate_IDs$cumTime[plate_IDs$type == "Kinetic"] <- formatted_time - formatted_time[1]

DAHP_plates <- tbl_df(merge(DAHP_plates, plate_IDs))

### Define experimental design of each plate ###

# INHinteract - lower concentrations of carbamoyl-phosphate

INH_ID_mat <- matrix(rep(rep(c("Tyrosine", "Tyrosol", "Phenylalanine", "Phenylpyruvate", "1-PhenylEthanol", "2-PhenylEthanol"), each = 2), 8), ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = T)
INH_ID_mat[8,] <- "None"
INH_ID <- melt(INH_ID_mat, value.name = "INH_ID")

INH_CONC <- melt(matrix(rep(c(10, 3, 1, 0.3, 0.1, 0.03, 0.01, 0), 12), ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = F), value.name = "INH_CONC")

ENZYME_CONC <- melt(matrix(1, ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = F), value.name = "ENZ_CONC")
ENZYME_CONC$ENZ_CONC[ENZYME_CONC$Row == "H" & ENZYME_CONC$Col %in% c(1:6)] <- 0

ARO3_design <- data.frame(merge(merge(INH_ID, INH_CONC), ENZYME_CONC), enzyme = "ARO3", plate = "ARO3")
ARO4_design <- data.frame(merge(merge(INH_ID, INH_CONC), ENZYME_CONC), enzyme = "ARO4", plate = "ARO4")

DAHP_plate_design <- rbind(ARO3_design, ARO4_design)

####### Plot OD(t) labelling #########

DAHP_plates_timecourse <- tbl_df(merge(DAHP_plates, DAHP_plate_design))
DAHP_plates_timecourse <- DAHP_plates_timecourse %>% filter(type == "Kinetic") %>% mutate(wellID = paste(plate, Row, Col, sep = "_"))
DAHP_plates_timecourse <- DAHP_plates_timecourse %>% group_by(wellID) %>% mutate(A_change = Absorbance - Absorbance[cumTime == 0])

scatter_facet_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                       panel.background = element_rect(fill = "gray80"), legend.position = "top", 
                       axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "black"),
                       axis.text = element_text(color = "black"), axis.text.x = element_text(size = 10, angle = 90, hjust = 0, vjust = 0.5),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), strip.background = element_rect(color = "black", fill = "burlywood1", size = 1.2)
                       )


# INHinteract
ggplot(DAHP_plates_timecourse, aes(x = cumTime, y = A_change, color = factor(INH_CONC), group = wellID)) + 
  geom_line(size = 1) + scatter_facet_theme + facet_grid(INH_ID ~ enzyme) +
  scale_y_continuous(expression(Delta ~ "Absorbance")) +
  scale_x_continuous("Time elapsed (s)")



#### Estimate intial rate for each well ####
  
# filter high OD concentrations of inhibitors because their kinetics cannot be assessed
DAHP_plates_summarize <- DAHP_plates_timecourse %>% filter(
  !(INH_ID %in% c("Tyrosine", "Tyrosol", "Phenylpyruvate") & INH_CONC %in% c(3,10))
  )


# take the slope by condition

#DAHP_plates_rates <- DAHP_plates_summarize %>% group_by(wellID) %>% do(tidy(lm(A_change ~ cumTime, data=.))) %>%
#  filter(term == "cumTime") %>% select(wellID, rate = estimate)

DAHP_plates_rates <- DAHP_plates_summarize %>% group_by(wellID) %>% do(tidy(lm(A_change ~ cumTime + I(cumTime^2), data=.))) %>%
  filter(term == "cumTime") %>% dplyr::select(wellID, rate = estimate)

DAHP_plates_PL <- DAHP_plates %>% filter(type == "PL") %>% select(plate, Row, Col, PL = Absorbance) %>% mutate(wellID = paste(plate, Row, Col, sep = "_"))

DAHP_plate_summaries <- inner_join(DAHP_plates_rates, DAHP_plates_PL)

DAHP_merged <- tbl_df(merge(DAHP_plate_summaries, DAHP_plate_design)) %>%
  mutate(corrected_rate = -1*rate / PL) %>% group_by(enzyme) %>% mutate(relative_rate = corrected_rate / mean(corrected_rate[INH_ID == "None" & ENZ_CONC == 1]))
DAHP_merged <- DAHP_merged %>% group_by(INH_ID, INH_CONC, ENZ_CONC, enzyme) %>% mutate(Rep = 1:length(enzyme)) 

barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 20, face = "bold"), 
                       panel.background = element_rect(fill = "gray80"), legend.position = "right", 
                       axis.ticks = element_line(color = "black", size = 1.5),
                       axis.text = element_text(color = "black"), axis.text.x = element_text(size = 20, angle = 90, hjust = 1, vjust = 0.5),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1.5),
                       axis.line = element_line(color = "black", size = 1), strip.background = element_rect(fill = "cadetblue2"),
                       panel.margin = unit(0.6, "lines")
                       )

ggplot(DAHP_merged %>% filter(ENZ_CONC == 1 & INH_CONC != 0), aes(x = factor(INH_CONC), y = relative_rate, group = Rep)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "cornflowerblue", color = "black") + 
  barplot_theme + facet_grid(INH_ID ~ enzyme) +
  scale_x_discrete("Inhibitior Concentration (mM)") + scale_y_continuous("Fractional Activity", label = percent_format())
ggsave("ARO3_4_allInh.pdf", height = 14, width = 7)

### All proposed inhibitors that were actually inert retained ~ 80% of total activity - use these inert samples to get the fractional activity
# for inhibitors with a real effect

DAHP_merged <- tbl_df(merge(DAHP_plate_summaries, DAHP_plate_design)) %>%
  mutate(corrected_rate = -1*rate / PL) 

DAHP_merged <- DAHP_merged %>% group_by(enzyme) %>% mutate(relative_rate = corrected_rate / mean(corrected_rate[INH_ID %in% c("1-PhenylEthanol", "2-Phenylethanol", "Tyrosol", "Phenylpyruvate") & INH_CONC == 0.01 & ENZ_CONC == 1]))

DAHP_merged <- DAHP_merged %>% group_by(INH_ID, INH_CONC, ENZ_CONC, enzyme) %>% mutate(Rep = 1:length(enzyme)) %>%
  ungroup() %>% filter(INH_ID %in% c("Phenylalanine", "Tyrosine", "Phenylpyruvate")) %>%
  mutate(INH_ID = factor(INH_ID, levels = c("Phenylalanine", "Tyrosine", "Phenylpyruvate")))


ggplot(DAHP_merged %>% filter(ENZ_CONC == 1 & INH_CONC != 0), aes(x = factor(INH_CONC), y = relative_rate, group = Rep)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "cornflowerblue", color = "black") + 
  barplot_theme + facet_grid(INH_ID ~ enzyme) +
  scale_x_discrete("Inhibitior Concentration (mM)") + scale_y_continuous("Fractional Activity", label = percent_format(), breaks = c(0, 0.25, 0.5, 0.75, 1))

# Estimate inhibitory constants

inh_eval_const <- function(inh_k, INH, RATE, inh_hill = 1, frac_const = 0){
  # with supplied vector of INH concnetrations and same length vector of RATE [0,1] evaluate the fit of proposed kinetic constant
  fitted_RATE <- inh_k ^ inh_hill / (inh_k ^ inh_hill + INH ^ inh_hill)
  fitted_RATE <- frac_const + (1-frac_const)*fitted_RATE
  sum((RATE - fitted_RATE)^2)
  }

# consider first just straight MM-kinetics

inhibited_pairs <- DAHP_merged %>% filter((plate == "ARO3" & INH_ID == "Phenylalanine") | (plate == "ARO4" & INH_ID %in% c("Phenylalanine", "Tyrosine")))
inhibited_pairs <- inhibited_pairs %>% mutate(pair = paste(enzyme, INH_ID, sep = "-")) %>% dplyr::select(pair, INH_CONC, relative_rate)

test_params <- 10^(seq(-3, 3, by = 0.01))

mm_param_eval <- expand.grid(Ki = test_params, pair = unique(inhibited_pairs$pair))
mm_param_eval <- mm_param_eval %>% left_join(inhibited_pairs) %>% tbl_df() %>% group_by(Ki, pair) %>%
  summarize(SS = inh_eval_const(Ki, INH = INH_CONC, RATE = relative_rate))

MM_inhibitory_constants <- mm_param_eval %>% group_by(pair) %>% filter(SS == min(SS))

inhibited_pairs <- inhibited_pairs %>% left_join(MM_inhibitory_constants %>% dplyr::select(-SS)) %>% mutate(fitted_rate = Ki / (Ki + INH_CONC)) %>% separate(pair, into = c("enzyme", "INH_ID"), remove = F) %>%
  mutate(INH_ID = factor(INH_ID, levels = levels(DAHP_merged$INH_ID)))

ggplot() + geom_bar(data = DAHP_merged %>% filter(ENZ_CONC == 1 & INH_CONC != 0), aes(x = factor(INH_CONC), y = relative_rate, group = Rep),
                    stat = "identity", position = "dodge", fill = "cornflowerblue", color = "black") + 
  geom_line(data = inhibited_pairs, aes(x = factor(INH_CONC), y = fitted_rate, group = pair), size = 1.3) + 
  geom_point(data = inhibited_pairs, aes(x = factor(INH_CONC), y = fitted_rate), fill = "coral", shape = 21, size = 3) +
  barplot_theme + facet_grid(INH_ID ~ enzyme) +
  scale_x_discrete("Inhibitior Concentration (mM)") + scale_y_continuous("Fractional Activity", label = percent_format(), breaks = c(0, 0.25, 0.5, 0.75, 1))
ggsave("ARO3_4_subsetInh.pdf", height = 7, width = 9)
ggsave("ARO3_4_subsetInh.eps", height = 7, width = 9)

