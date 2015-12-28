setwd("~/Desktop/Rabinowitz/FBA_SRH/Validations/PYK")

options(stringsAsFactors = F)

library(reshape2)
# library(data.table)
library(scales)
library(dplyr)
library(ggplot2)
library(grid)

#library(devtools)
#install_github("dgrtwo/broom")
library(broom)
library(tidyr)

PYK_dualReg <- read.delim('2015.11.24-pyk_dualReg3.txt', fill = NA, header = F)
PYK_FBP_citrate <- read.delim('2015.11.24-pyk_fbp_citrate.txt', fill = NA, header = F)
PYK_FBP_AMP <- read.delim('2015.11.24-pyk_fbp_amp.txt', fill = NA, header = F)

all_plate_names <- c("PYK_dualReg", "PYK_FBP_citrate", "PYK_FBP_AMP")

all_plate_list <- list()
for(one_plate in all_plate_names){
  
  kineticData <- get(one_plate)
  
  plate_starts <- grep('^Read', kineticData[,1])
  plates <- list()
  
  for(a_plate in 1:length(plate_starts)){
    
    if(a_plate == length(plate_starts)){
      plate_data <- kineticData[plate_starts[a_plate]:nrow(kineticData),]
    }else{
      plate_data <- kineticData[plate_starts[a_plate]:(plate_starts[a_plate + 1]-1),]
    }
    
    header_rows <- !(apply(!is.na(plate_data), 1, sum) > 4)
    header_rows[header_rows] <- !(plate_data[header_rows,1] %in% LETTERS[1:8])
    
    plate_matrix <- plate_data[!header_rows,1:(ncol(plate_data)-1)]
    plate_matrix[1,1] <- "Row"
    colnames(plate_matrix) <- plate_matrix[1,]; plate_matrix <- plate_matrix[-1,]
    
    if(nrow(plate_matrix) == 0){next}
    
    plate_summary <- melt(plate_matrix, variable.name = "Col", value.name = "OD", id.vars = "Row")
    plate_summary <- data.frame(plate_ID = last(plate_data[header_rows,1]), plate_summary)
    
    plates[[a_plate]] <- plate_summary
  }
  
  plates_aggregate <- do.call("rbind", plates)
  plates_aggregate <- plates_aggregate[!is.na(plates_aggregate$OD),]
  
  all_plate_list[[one_plate]] <- plates_aggregate %>% mutate(plate = one_plate)
  
}

all_plate_data <- do.call("rbind", all_plate_list) %>% tbl_df()

# Reformat plate time and cumulative time
plate_IDs <- all_plate_data %>% group_by(plate_ID, plate) %>% summarize(type = "Kinetic", time = NA, cumTime = NA)
plate_IDs$type[grep('Pathlength', plate_IDs$plate_ID)] <- "PL"
plate_IDs$time[plate_IDs$type == "Kinetic"] <- regmatches(plate_IDs$plate_ID, regexpr('[0-9]{1}:[0-9]{2}:[0-9]{2}', plate_IDs$plate_ID))

formatted_time <- as.POSIXlt(plate_IDs$time[plate_IDs$type == "Kinetic"], format = "%H:%M:%S")

plate_IDs$cumTime[plate_IDs$type == "Kinetic"] <- formatted_time - formatted_time[1]

all_plate_data <- all_plate_data %>% left_join(plate_IDs, by = c("plate_ID", "plate")) %>%
  mutate(Col = as.numeric(Col))


### Define experimental design of each plate ###

plate_design_track <- NULL

# PYK_dualReg

INH_ID_mat <- matrix(rep(c(rep("Citrate", 5), rep("AMP", 5), "None", "None"), times = 8),
                     ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = T)
INH_ID <- melt(INH_ID_mat, value.name = "INH_ID")

INH_CONC <- melt(matrix(rep(c(rep(c(2, 6, 10, 15, 20), times = 2), 0, 0)*10/12, times = 8), 
                        ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = T), value.name = "INH_CONC")

FBP_CONC <- melt(matrix(c(rep(0.5, 96/2), rep(5, 96/2))*10/12,
                        ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = T), value.name = "FBP_CONC")

PEP_CONC <- melt(matrix(1 * 10/12, ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = T), value.name = "PEP_CONC")

ENZYME_CONC <- melt(matrix(c(rep(1, 88), rep(0, 8)), ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = F), value.name = "ENZ_CONC")

plate_design_track <- rbind(plate_design_track, INH_ID %>%
                              left_join(INH_CONC, by = c("Row", "Col")) %>%
                              left_join(FBP_CONC, by = c("Row", "Col")) %>%
                              left_join(PEP_CONC, by = c("Row", "Col")) %>%
                              left_join(ENZYME_CONC, by = c("Row", "Col")) %>%
                              mutate(plate = "PYK_dualReg"))

# PYK_FBP_citrate

INH_ID_mat <- matrix(c(rep("Citrate", 12*6), rep("None", 12*2)),
                     ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = T)
INH_ID <- melt(INH_ID_mat, value.name = "INH_ID")

INH_CONC <- melt(matrix(rep(c(1, 2, 4, 6, 12, 20, 0, 0), 12)*10/13,
                        ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = F), value.name = "INH_CONC")

FBP_CONC <- melt(matrix(rep(c(0.5, 1, 2, 3, 6, 10), times = 16)*10/13,
                        ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = T), value.name = "FBP_CONC")

PEP_CONC <- melt(matrix(rep(c(0.2, 1), each = 48) * 10/13,
                        ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = F), value.name = "PEP_CONC")

ENZYME_CONC <- melt(matrix(c(rep(1, 84), rep(0, 12)), ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = T), value.name = "ENZ_CONC")

plate_design_track <- rbind(plate_design_track, INH_ID %>%
                              left_join(INH_CONC, by = c("Row", "Col")) %>%
                              left_join(FBP_CONC, by = c("Row", "Col")) %>%
                              left_join(PEP_CONC, by = c("Row", "Col")) %>%
                              left_join(ENZYME_CONC, by = c("Row", "Col")) %>%
                              mutate(plate = "PYK_FBP_citrate"))

# PYK_FBP_AMP

INH_ID <- INH_ID %>% mutate(INH_ID = gsub('Citrate', 'AMP', INH_ID))

plate_design_track <- rbind(plate_design_track, INH_ID %>%
                              left_join(INH_CONC, by = c("Row", "Col")) %>%
                              left_join(FBP_CONC, by = c("Row", "Col")) %>%
                              left_join(PEP_CONC, by = c("Row", "Col")) %>%
                              left_join(ENZYME_CONC, by = c("Row", "Col")) %>%
                              mutate(plate = "PYK_FBP_AMP"))

# Clean-up values slightly (rounding)

round2 <- function(x){round(x,2)}
plate_design_track <- plate_design_track %>% mutate_each(funs(round2), INH_CONC:PEP_CONC)


# Add experimental design to kinetic data

plate_design_track <- plate_design_track %>% mutate(Row = as.character(Row))

all_plate_data <- all_plate_data %>% left_join(plate_design_track, c("Row", "Col", "plate"))

scatter_data <- all_plate_data %>% filter(type == "Kinetic") %>% mutate(wellID = paste(plate, Row, Col, sep = "_"))

#### Filter time as needed #### 
# inspect raw OD plots and adjust filtering parameters on an assay-by-assay basis as needed

filter_spec <- function(OD, cumTime, is_decreasing = F, frac_filter_beginning = 0.1, frac_filter_end = 0.5){
  # for a trajectory o\t : OD ~ cumTime
  # filter beginning and end of the trajectory (max defined by frac_filter_beginning and frac_filter_end)
  # to minimize the RMSE
  # is_decreasing defines the expected direction of OD change
  
  # return:
  # RMSE - for subsequent outlier filtering
  
  #Filter then end to remove large pathologies like running out of substrate, then remove irregularities at the beginning
  
  # Filter the end:
  
  if(frac_filter_end != 0){
    
    possible_end_points <- seq(ceiling(length(cumTime) * (1-frac_filter_end)), length(cumTime))
    
    end_fit <- lapply(possible_end_points, function(x){
      model <- lm(OD[1:x] ~ cumTime[1:x])
      data.frame(
        end_point = x,
        slope = coef(model)[2],
        RMSE = sqrt(deviance(model)/(df.residual(model))),
        df = df.residual(model)
      )})
    end_fit <- do.call("rbind", end_fit)
    
    # filter to slopes in the appropriate direction (if any sequence exists)
    if(is_decreasing == T){
      dir_filter <- end_fit %>% filter(slope <= 0)
    }else{
      dir_filter <- end_fit %>% filter(slope >= 0)
    }
    if(length(dir_filter) == 0){
      dir_filter <- end_fit 
    }
    
  end_point = dir_filter$end_point[which.min(dir_filter$RMSE * (dir_filter$df/(dir_filter$df-10)))]
    
  }else{
    end_point = length(cumTime)
  }
  
  # Filter the beginning:
  
  if(frac_filter_beginning != 0){
    
    possible_start_points <- seq(1, ceiling(length(cumTime) * frac_filter_beginning))
    
    start_fit <- lapply(possible_start_points, function(x){
      model <- lm(OD[x:end_point] ~ cumTime[x:end_point])
      data.frame(
        start_point = x,
        slope = coef(model)[2],
        RMSE = sqrt(deviance(model)/(df.residual(model))),
        df = df.residual(model)
      )})
    start_fit <- do.call("rbind", start_fit)
    
    # filter to slopes in the appropriate direction (if any sequence exists)
    if(is_decreasing == T){
      dir_filter <- start_fit %>% filter(slope <= 0)
    }else{
      dir_filter <- start_fit %>% filter(slope >= 0)
    }
    if(length(dir_filter) == 0){
      dir_filter <- start_fit
    }
    
    start_point = dir_filter$start_point[which.min(dir_filter$RMSE *(dir_filter$df/(dir_filter$df-10)))]
    
  }else{
    start_point = 1
  }
  
  return(data.frame(start = start_point, end = end_point,
                    slope = dir_filter$slope[dir_filter$start_point == start_point],
                    RMSE = dir_filter$RMSE[dir_filter$start_point == start_point]
  ))
  
}

filter_setup <- data.frame(wellID = unique(scatter_data$wellID)) %>%
  mutate(is_decreasing = T,
         frac_filter_beginning = 0.05,
         frac_filter_end = 0.6)

filter_bounds <- lapply(1:nrow(filter_setup), function(i){
  plate_data <- scatter_data %>% filter(wellID == filter_setup$wellID[i])
  filter_spec(OD = plate_data$OD, cumTime = plate_data$cumTime,
              is_decreasing = filter_setup$is_decreasing[i],
              frac_filter_beginning = filter_setup$frac_filter_beginning[i],
              frac_filter_end = filter_setup$frac_filter_end[i])
})
filter_bounds <- do.call("rbind", filter_bounds)

filter_bounds <- data.frame(wellID = filter_setup$wellID, filter_bounds)

filtered_scatter_data <- scatter_data %>% filter(type == "Kinetic") %>% 
  group_by(wellID) %>% mutate(timePoint = 1:n()) %>%
  left_join(filter_bounds, by = "wellID")

# show included and excluded data from time zero

# either filter or not (comment out one line)
filtered_scatter_data <- filtered_scatter_data %>% mutate(filtered = ifelse(timePoint < start | timePoint > end, T, F))
#filtered_scatter_data <- filtered_scatter_data %>% mutate(filtered = F)

# fit a quadratic 
quad_fit <- filtered_scatter_data %>% filter(!filtered) %>% group_by(wellID) %>% mutate(cumTime = cumTime - min(cumTime)) %>%
  do(tidy(lm(OD ~ cumTime + I(cumTime^2), data = .))) %>%
  filter(term == "cumTime") %>% dplyr::select(wellID, quad_fit = estimate)



scatter_facet_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                       panel.background = element_rect(fill = "gray80"), legend.position = "top", 
                       axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "black"),
                       axis.text = element_text(color = "black"), axis.text.x = element_text(size = 10, angle = 90, hjust = 0, vjust = 0.5),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), strip.background = element_rect(color = "black", fill = "burlywood1", size = 1.2)
                       )


for(one_plate in all_plate_names){
  
  plate_data <- filtered_scatter_data %>% filter(plate == one_plate) %>% 
    group_by(wellID) %>% mutate(A_change = OD - OD[cumTime == 0])
  
  # all points showing filtering
  print(
    ggplot(plate_data, aes(x = cumTime, y = A_change, color = factor(INH_CONC), group = wellID, alpha = -1*filtered)) + 
      geom_line(size = 1) + scatter_facet_theme + facet_grid(INH_ID ~ PEP_CONC + FBP_CONC) +
      scale_y_continuous(expression(Delta ~ "Absorbance")) +
      scale_x_continuous("Time elapsed (s)")
  )
  ggsave(paste0(one_plate, "_scatter.pdf"), height = 8, width = 10)
  
  # with filtered applied, changing initial time
  filtered_plate_data <- plate_data %>% filter(!filtered) %>% mutate(A_change = OD - OD[timePoint == start],
                                                                     cumTime = cumTime - cumTime[timePoint == start])
   print(
    ggplot(filtered_plate_data, aes(x = cumTime, y = A_change, color = factor(INH_CONC), group = wellID)) + 
      geom_line(size = 1) + scatter_facet_theme + facet_grid(INH_ID ~ PEP_CONC + FBP_CONC) +
      scale_y_continuous(expression(Delta ~ "Absorbance")) +
      scale_x_continuous("Time elapsed (s)")
  )
  ggsave(paste0(one_plate, "_filtered_scatter.pdf"), height = 8, width = 10)
  
}

### Filter samples with exceptional non-linearity (relative to other assay samples) and flawed samples

all_plate_rates <- filter_bounds %>% left_join(scatter_data %>% dplyr::select(plate, wellID, INH_ID:ENZ_CONC) %>% unique(), by = "wellID") %>%
  left_join(quad_fit, by = "wellID")




ggplot(all_plate_rates, aes(x = RMSE, y = ..density..)) + geom_bar() + facet_wrap(~ plate, scale = "free")
all_plate_rates %>% filter(RMSE > 0.03)

ggplot( all_plate_data %>% filter(type == "PL"), aes(x = OD, y = ..density..)) + geom_bar()
all_plate_data %>% filter(type == "PL") %>% filter(OD< 0.6)

#all_plate_rates <- all_plate_rates %>% filter(!(wellID %in% c("PYK_regTest_D_1", "PYK_modeofInhibition_E_12", "PYK_modeofInhibition_F_12")))

### 

all_plates_PL <- all_plate_data %>% filter(type == "PL") %>% dplyr::select(plate, Row, Col, PL = OD) %>% mutate(wellID = paste(plate, Row, Col, sep = "_"))

all_plate_summaries <- left_join(all_plate_rates, all_plates_PL %>% dplyr::select(-plate), by = "wellID") %>% tbl_df() %>%
  dplyr::select(wellID, plate, INH_ID:ENZ_CONC, rate = slope, qrate = quad_fit, PL)

# subtract blank
# invert negative rates
# path-length correct

#tmp <- all_plate_summaries %>% filter(plate == "PYK_noCoupling")

all_plates_merge <- all_plate_summaries %>%
  #tmp %>%
  mutate(rate = -1 * qrate) %>%
  #mutate(rate = rate / PL) %>%
  group_by(plate) %>% mutate(relative_rate = rate - mean(rate[ENZ_CONC == 0])) %>%
  group_by(plate, INH_ID, INH_CONC, FBP_CONC, PEP_CONC, ENZ_CONC) %>% mutate(Rep = 1:n()) %>%
  group_by(plate, PEP_CONC, FBP_CONC) %>% 
  mutate(frac_max = relative_rate / mean(relative_rate[ENZ_CONC == 1 & INH_ID == "None"]))


barplot_theme <- theme(text = element_text(size = 20), title = element_text(size = 20), 
                       panel.background = element_rect(fill = "gray80"), legend.position = "right", 
                       axis.ticks = element_line(color = "black", size = 1),
                       axis.text = element_text(color = "black"), axis.text.x = element_text(size = 20, angle = 90, hjust = 1, vjust = 0.5),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1.5),
                       axis.line = element_line(color = "black", size = 1), strip.background = element_rect(fill = "cadetblue2"),
                       panel.margin = unit(1, "lines")
                       )


for(one_plate in all_plate_names){
  
  plot_data <- all_plates_merge %>% filter(plate == one_plate, ENZ_CONC == 1 & INH_CONC != 0)
  
  print(
    ggplot(plot_data, aes(x = factor(INH_CONC), y = relative_rate, group = Rep)) + 
      geom_bar(stat = "identity", position = "dodge", color = "black", fill = "coral2") + 
      barplot_theme + facet_grid(PEP_CONC ~ FBP_CONC + INH_ID) +
      scale_x_discrete("Inhibitior Concentration (mM)") + scale_y_continuous("Rate")
  )
  ggsave(paste0(one_plate, "_barplot.pdf"), height = 8, width = 10)
 
  if(any(!is.na(plot_data$frac_max))){
  print(
    ggplot(plot_data, aes(x = factor(INH_CONC), y = frac_max, group = Rep)) + 
      geom_bar(stat = "identity", position = "dodge", color = "black", fill = "coral2") + 
      barplot_theme + facet_grid(PEP_CONC ~ FBP_CONC + INH_ID) +
      scale_x_discrete("Inhibitior Concentration (mM)") + scale_y_continuous("Fraction of ``No Inhibitor``")
  )
  ggsave(paste0(one_plate, "_frac_max_barplot.pdf"), height = 8, width = 10)
   
  }
  
}

# V(S)

PEP_STD <- 0.769230769230769
ADP_STD <- 0.1538462

inhibitor_mode_data <- all_plates_merge %>% filter(plate == "PYK_modeofInhibition", ENZ_CONC == 1) %>%
  mutate(FIXED_PEP = ifelse(abs(PEP_CONC - PEP_STD) < 0.00001, T, F),
         FIXED_ADP = ifelse(abs(ADP_CONC - ADP_STD) < 0.00001, T, F))

ggplot(inhibitor_mode_data %>% filter(FIXED_ADP), aes(x = PEP_CONC, y = relative_rate, group = factor(INH_CONC), color = factor(INH_CONC))) + 
  geom_point(size = 3) + geom_smooth(method = "lm", fullrange = T, se = F) + barplot_theme +
  scale_x_continuous("[PEP]") +
  scale_color_brewer("[Citrate] (mM)", palette = "Set1")

ggplot(inhibitor_mode_data %>% filter(FIXED_PEP), aes(x = ADP_CONC, y = relative_rate, group = factor(INH_CONC), color = factor(INH_CONC))) + 
  geom_point(size = 3) + geom_smooth(method = "lm", fullrange = T, se = F) + barplot_theme +
  scale_x_continuous("[ADP]") +
  scale_color_brewer("[Citrate] (mM)", palette = "Set1")

### Lineweaver-burke ###


ggplot(inhibitor_mode_data %>% ungroup() %>% filter(FIXED_ADP) %>% filter(PEP_CONC > 0.24), aes(x = 1/PEP_CONC, y = 1/relative_rate, group = factor(INH_CONC), color = factor(INH_CONC))) + 
  geom_point(size = 3) + geom_smooth(method = "lm", fullrange = T, se = F) + barplot_theme +
  scale_x_continuous("1/[PEP]") +
  scale_color_brewer("[Citrate] (mM)", palette = "Set1") + ggtitle("Lineweaver-Burke plot testing regulation\nof pyruvate kinase") +
  expand_limits(x = 0)
ggsave("Lineweaver_pepPYK.pdf", height = 8, width = 8)




ggplot(inhibitor_mode_data %>% filter(FIXED_PEP) %>% filter(INH_CONC < 15), aes(x = 1/ADP_CONC, y = 1/relative_rate, group = factor(INH_CONC), color = factor(INH_CONC))) + 
  geom_point(size = 3) + geom_smooth(method = "lm", fullrange = T, se = F) + barplot_theme +
  scale_x_continuous("1/[ADP]") +
  scale_color_brewer("[Citrate] (mM)", palette = "Set1") + ggtitle("Lineweaver-Burke plot testing regulation\nof pyruvate kinase")
ggsave("Lineweaver_adpPYK.pdf", height = 8, width = 8)

### Hanes plot ###

ggplot(inhibitor_mode_data %>% filter(FIXED_ADP), aes(x = PEP_CONC, y = PEP_CONC/rate, group = factor(INH_CONC), color = factor(INH_CONC))) + 
  geom_point(size = 3) + geom_smooth(method = "lm", fullrange = T, se = F) + expand_limits(x = 0) + barplot_theme +
  scale_color_brewer("[Citrate] (mM)", palette = "Set1") + ggtitle("Hanes plot testing regulation\nof pyruvate kinase") +
  scale_x_continuous('[PEP] (mM)') + scale_y_continuous('[PEP] / rate (a.u.)')
ggsave("Hanes_phePDC.pdf", height = 8, width = 8)







# Estimate inhibitory constants

inh_eval_const <- function(inh_k, INH, RATE, inh_hill = 1, frac_const = 0){
  # with supplied vector of INH concnetrations and same length vector of RATE [0,1] evaluate the fit of proposed kinetic constant
  fitted_RATE <- inh_k ^ inh_hill / (inh_k ^ inh_hill + INH ^ inh_hill)
  fitted_RATE <- frac_const + (1-frac_const)*fitted_RATE
  sum((RATE - fitted_RATE)^2)
  }

# consider first just straight MM-kinetics

#inhibited_pairs <- all_plates_merge %>% filter(plate == "GsixPD_followup", INH_ID != "None") %>%
#  mutate(pair = paste(plate, INH_ID, sep = "-")) %>% dplyr::select(pair, INH_CONC, relative_rate)

inhibited_pairs <- all_plates_merge %>% filter(plate == "PDC_phenpyr", INH_ID != "None") %>%
  mutate(pair = paste(plate, INH_ID, sep = "-")) %>% dplyr::select(pair, INH_CONC, relative_rate)


test_params <- 10^(seq(-3, 3, by = 0.01))

mm_param_eval <- expand.grid(Ki = test_params, pair = unique(inhibited_pairs$pair))
mm_param_eval <- mm_param_eval %>% left_join(inhibited_pairs) %>% tbl_df() %>% group_by(Ki, pair) %>%
  summarize(SS = inh_eval_const(Ki, INH = INH_CONC, RATE = relative_rate))

MM_inhibitory_constants <- mm_param_eval %>% group_by(pair) %>% filter(SS == min(SS))

print(MM_inhibitory_constants)

# look at ultrasensitivity

test_params <- 10^(seq(-3, 3, by = 0.01))
hill_test <- seq(1, 6, by = 0.1)

hill_param_eval <- expand.grid(Ki = test_params, hill = hill_test, pair = unique(inhibited_pairs$pair))
hill_param_eval <- hill_param_eval %>% left_join(inhibited_pairs) %>% tbl_df() %>% group_by(Ki, hill, pair) %>%
  summarize(SS = inh_eval_const(Ki, INH = INH_CONC, inh_hill = hill, RATE = relative_rate))

hill_inhibitory_constants <- hill_param_eval %>% group_by(pair) %>% filter(SS == min(SS))

print(hill_inhibitory_constants)

# Determine rate at each measured concentration

inhibited_pairs <- inhibited_pairs %>% left_join(MM_inhibitory_constants %>% dplyr::select(-SS)) %>% mutate(mm_rate = Ki / (Ki + INH_CONC)) %>% dplyr::select(-Ki) %>%
  left_join(hill_inhibitory_constants %>% dplyr::select(-SS)) %>% mutate(hill_rate = Ki ^ hill / (Ki ^ hill + INH_CONC ^ hill)) %>% dplyr::select(-Ki, -hill)

inhibited_pairs_bar <- all_plates_merge %>% mutate(pair = paste(plate, INH_ID, sep = "-")) %>%
                      filter(pair %in% inhibited_pairs$pair)


ggplot() + geom_bar(data = inhibited_pairs_bar, aes(x = factor(INH_CONC), y = relative_rate, group = Rep),
                    stat = "identity", position = "dodge", fill = "cornflowerblue", color = "black") + 
  geom_line(data = inhibited_pairs, aes(x = factor(INH_CONC), y = mm_rate, group = pair), color = "coral", size = 1.3) + 
  geom_point(data = inhibited_pairs, aes(x = factor(INH_CONC), y = mm_rate), fill = "coral", shape = 21, size = 3) +
  geom_line(data = inhibited_pairs, aes(x = factor(INH_CONC), y = hill_rate, group = pair), color = "chartreuse3", size = 1.3) + 
  geom_point(data = inhibited_pairs, aes(x = factor(INH_CONC), y = hill_rate), fill = "chartreuse3", shape = 21, size = 3) +
  barplot_theme + facet_grid(INH_ID ~ plate) +
  scale_x_discrete("Inhibitior Concentration (mM)") + scale_y_continuous("Fractional Activity", label = percent_format(), breaks = c(0, 0.25, 0.5, 0.75, 1))
ggsave("PDC_inhibition.pdf", height = 7, width = 10)


# Estimate activation

act_eval_const <- function(act_k, ACT, HILL = 1, frac_const, RATE){
  # with supplied vector of INH concnetrations and same length vector of RATE [0,1] evaluate the fit of proposed kinetic constant
  # since the activity will be a fraction of Vmax, find
  # v = vmax * (c + (1-c)*a/(a + ka))
  # vo = vmax/c
  # v = (1 - c)/c * a /(a + ka)
  
  fitted_RATE <- 1 + (1-frac_const)/frac_const * ACT^HILL / (ACT^HILL + act_k^HILL)
  sum((RATE - fitted_RATE)^2)
  
}

activated_pairs <- all_plates_merge %>% filter(plate == "SixPGD_asp_glu", INH_ID != "None") %>%
  mutate(pair = paste(plate, INH_ID, sep = "-")) %>% dplyr::select(pair, INH_CONC, relative_rate)

test_params <- 10^(seq(-1, 2, by = 0.01))
const_act <- seq(0, 1, by = 0.01)

act_param_eval <- expand.grid(Ka = test_params, const = const_act, pair = unique(activated_pairs$pair))
act_param_eval <- act_param_eval %>% left_join(activated_pairs) %>% tbl_df() %>% group_by(Ka, const, pair) %>%
  summarize(SS = act_eval_const(Ka, ACT = INH_CONC, frac_const = const, RATE = relative_rate))

MM_activation_constants <- act_param_eval %>% group_by(pair) %>% filter(SS == min(SS))

# also look at ultrasensitivity

hill_test <- seq(1, 6, by = 0.1)

hill_act_param_eval <- expand.grid(Ka = test_params, hill = hill_test, const = const_act, pair = unique(activated_pairs$pair))
hill_act_param_eval <- hill_act_param_eval %>% left_join(activated_pairs) %>% tbl_df() %>% group_by(Ka, const, hill, pair) %>%
  summarize(SS = act_eval_const(Ka, ACT = INH_CONC, frac_const = const, HILL = hill, RATE = relative_rate))

hill_activation_constants <- hill_act_param_eval %>% group_by(pair) %>% filter(SS == min(SS))

# fit data with reaction parameters

activated_pairs <- activated_pairs %>% left_join(MM_activation_constants %>% dplyr::select(-SS)) %>%
  mutate(mm_rate = 1 + (1-const)/const * INH_CONC / (INH_CONC + Ka)) %>% dplyr::select(-Ka, -const) %>%
  left_join(hill_activation_constants) %>% mutate(hill_rate = 1 + (1-const)/const * INH_CONC^hill / (INH_CONC^hill + Ka^hill)) %>%
  dplyr::select(-Ka, -const, -hill)
  
activated_pairs_bar <- all_plates_merge %>% mutate(pair = paste(plate, INH_ID, sep = "-")) %>%
                      filter(pair %in% activated_pairs$pair)

ggplot() + geom_bar(data = activated_pairs_bar, aes(x = factor(INH_CONC), y = relative_rate, group = Rep),
                    stat = "identity", position = "dodge", fill = "cornflowerblue", color = "black") + 
  geom_line(data = activated_pairs, aes(x = factor(INH_CONC), y = mm_rate, group = pair), color = "coral", size = 1.3) + 
  geom_point(data = activated_pairs, aes(x = factor(INH_CONC), y = mm_rate), fill = "coral", shape = 21, size = 3) +
  geom_line(data = activated_pairs, aes(x = factor(INH_CONC), y = hill_rate, group = pair), color = "chartreuse3", size = 1.3) + 
  geom_point(data = activated_pairs, aes(x = factor(INH_CONC), y = hill_rate), fill = "chartreuse3", shape = 21, size = 3) +
  barplot_theme + facet_grid(INH_ID ~ plate) +
  scale_x_discrete("Inhibitior Concentration (mM)") + scale_y_continuous("Fractional Activity", label = percent_format(), breaks = seq(0, 2.25, by = 0.25))
ggsave("SixPGD_activation.pdf", height = 10, width = 10)





