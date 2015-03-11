setwd("~/Desktop/Rabinowitz/FBA_SRH/Validations/DAHP")

options(stringsAsFactors = F)

library(reshape2)
library(data.table)
library(scales)
library(dplyr)
library(ggplot2)

#library(devtools)
#install_github("dgrtwo/broom")
library(broom)

INHinteract <- read.delim('11.7-DAHParoAA.txt', fill = NA, header = F) 
DAHP_plate_names <- c("INHinteract")

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
DAHP_plates <- DAHP_plates %>% select(plate_ID, Row, Col, Absorbance, plate = L1)

plate_IDs <- DAHP_plates %>% group_by(plate_ID, plate) %>% summarize(type = "Kinetic", time = NA, cumTime = NA)
plate_IDs$type[grep('Pathlength', plate_IDs$plate_ID)] <- "PL"
plate_IDs$time[plate_IDs$type == "Kinetic"] <- regmatches(plate_IDs$plate_ID, regexpr('[0-9]{1}:[0-9]{2}:[0-9]{2}', plate_IDs$plate_ID))

formatted_time <- as.POSIXlt(plate_IDs$time[plate_IDs$type == "Kinetic"], format = "%H:%M:%S")

plate_IDs$cumTime[plate_IDs$type == "Kinetic"] <- formatted_time - formatted_time[1]

DAHP_plates <- tbl_df(merge(DAHP_plates, plate_IDs))

### Define experimental design of each plate ###

# INHinteract - lower concentrations of carbamoyl-phosphate

PHE <- melt(matrix(c(rep(2^(0:-7), 2), rep(c(2^(-8:-11), rep(0,4)), 2), rep(c(rep(0, 16), rep(0.5, 16)), 2)), ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = F)*(148 * 20 / 120), value.name = "PHE")
TYR <- melt(matrix(c(rep(0, 64), rep(2^(0:-7), 4)) * 1.8, ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = F), value.name = "TYR")
TRP <- melt(matrix(c(rep(0, 32), rep(2^(0:-7), 4), rep(0, 32)) * 2.4, ncol = 12, nrow = 8, dimnames = list(Row = LETTERS[1:8], Col = as.character(1:12)), byrow = F), value.name = "TRP")

All_INH_design <- merge(merge(PHE, TYR), TRP)
DAHP_plate_design <- All_INH_design

####### Plot OD(t) labelling #########

DAHP_plates_timecourse <- tbl_df(merge(DAHP_plates, DAHP_plate_design))
DAHP_plates_timecourse <- DAHP_plates_timecourse %>% filter(type == "Kinetic") %>% mutate(wellID = paste(plate, Row, Col, sep = "_")) %>%
  select(plate, wellID, Absorbance, cumTime:TRP)
DAHP_plates_timecourse <- DAHP_plates_timecourse %>% group_by(wellID) %>% mutate(A_change = Absorbance - Absorbance[cumTime == 0])

scatter_facet_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                       panel.background = element_rect(fill = "gray80"), legend.position = "top", 
                       axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "black"),
                       axis.text = element_text(color = "black"), axis.text.x = element_text(size = 10, angle = 90, hjust = 0, vjust = 0.5),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), strip.background = element_rect(color = "black", fill = "burlywood1", size = 1.2)
                       )


# INHinteract
ggplot(DAHP_plates_timecourse, aes(x = cumTime, y = A_change, color = log10(TRP), group = wellID)) + 
  geom_line(size = 1) + scatter_facet_theme +
  scale_color_continuous("[Alanine] (mM)", na.value = "black", low = "yellow", high = "blue") +
  scale_y_continuous(expression(Delta ~ "Absorbance")) +
  scale_x_continuous("Time elapsed (s)")



#### Estimate intial rate for each well ####
  
# filter outlier wells and pre-steady-state time-points
DAHP_plates_summarize <- DAHP_plates_timecourse %>% filter(
  A_change > -0.145,
  !(wellID %in% paste0("INHinteract_A_", 4:12)),
  !(plate == "INHinteract" & cumTime < 200)
  )


# take the slope by condition

DAHP_plates_rates <- DAHP_plates_summarize %>% group_by(wellID) %>% do(tidy(lm(A_change ~ cumTime, data=.))) %>%
  filter(term == "cumTime") %>% select(wellID, rate = estimate)

DAHP_plates_PL <- DAHP_plates %>% filter(type == "PL") %>% select(plate, Row, Col, PL = Absorbance) %>% mutate(wellID = paste(plate, Row, Col, sep = "_"))

DAHP_plate_summaries <- inner_join(DAHP_plates_rates, DAHP_plates_PL)

DAHP_merged <- tbl_df(merge(DAHP_plate_summaries, DAHP_plate_design)) %>% select(Row, Col, PHE:TRP, rate, PL) %>%
  mutate(corrected_rate = -1*rate / PL)

DAHP_merged <- DAHP_merged %>% mutate(relative_rate = corrected_rate / mean(corrected_rate[PHE == 0 & TYR == 0 & TRP == 0])) %>%
  filter(!(PHE == 0 & TYR == 0 & TRP == 0))

####

barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 20, face = "bold"), 
                       panel.background = element_rect(fill = "gray80"), legend.position = "right", 
                       axis.ticks = element_line(color = "black", size = 1.5),
                       axis.text = element_text(color = "black"), axis.text.x = element_text(size = 20, angle = 90, hjust = 1, vjust = 0.5),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1.5),
                       axis.line = element_line(color = "black", size = 1), strip.background = element_rect(fill = "cadetblue2"),
                       panel.margin = unit(0.6, "lines")
                       )


ggplot(DAHP_merged %>% filter(TRP == 0 & TYR == 0), aes(x = factor(round(PHE, 2)), y = relative_rate)) + geom_bar(stat = "identity", position = "dodge") + 
  barplot_theme + scale_y_continuous("V / Vmax") + scale_x_discrete("[PHE] (mM)")
ggsave("ARO3_PHE_11_7.pdf", height = 8, width = 8)

ggplot(DAHP_merged %>% filter(TRP != 0) %>% mutate(PHE = paste("[PHE] =", round(PHE, 2), "mM")) %>% group_by(PHE, TYR, TRP) %>% mutate(Rep = 1:length(relative_rate)), 
       aes(x = factor(round(TRP, 2)), y = relative_rate, group = Rep)) + geom_bar(stat = "identity", position = "dodge", color = "black", fill = "cornflowerblue") + 
  barplot_theme + facet_grid( ~ PHE, scale = "free_y") + scale_y_continuous("Fractional Activity", label = percent_format(), breaks = c(0, 0.25, 0.5, 0.75, 1)) + scale_x_discrete("[TRP] mM")
ggsave("ARO3TRP.eps", width = 9, height = 5)
ggsave("ARO3TRP.pdf", width = 9, height = 5)

ggplot(DAHP_merged %>% filter(TYR != 0) %>% mutate(PHE = paste("[PHE] =", round(PHE, 2), "mM")) %>% group_by(PHE, TYR, TRP) %>% mutate(Rep = 1:length(relative_rate)), 
       aes(x = factor(round(TYR, 2)), y = relative_rate, group = Rep)) + geom_bar(stat = "identity", position = "dodge", color = "black", fill = "cornflowerblue") + 
  barplot_theme + facet_grid( ~ PHE, scale = "free_y") + scale_y_continuous("Fractional Activity", label = percent_format(), breaks = c(0, 0.25, 0.5, 0.75, 1)) + scale_x_discrete("[TYR] mM")
ggsave("ARO3TyR.eps", width = 9, height = 5)
ggsave("ARO3TyR.pdf", width = 9, height = 5)


### Test for MM kinetics

inh_fit_data <- DAHP_merged %>% filter(TRP == 0 & TYR == 0) %>% select(INH = PHE, RATE = relative_rate)

# standard inhibition 
inh_conc <- seq(1:100)
inh_k <- 20
plot(inh_k / (inh_k + inh_conc) ~ log10(inh_conc))
plot(inh_k^3 / (inh_k^3 + inh_conc^3) ~ inh_conc)
plot(inh_k^0.6 / (inh_k^0.6 + inh_conc^0.6) ~ inh_conc)
  
inh_eval <- function(inh_k, inh_hill = 1){
  fitted_RATE <- inh_k ^ inh_hill / (inh_k ^ inh_hill + inh_fit_data$INH ^ inh_hill)
  sum((inh_fit_data$RATE - fitted_RATE)^2)
  }

k_test <- seq(0.05, 10, by = 0.01)
k_test_SS <- data.frame(k = k_test, TSS = sapply(k_test, inh_eval))
kinh_MLE <- k_test_SS$k[which.min(k_test_SS$TSS)]

inh_fit_data <- inh_fit_data %>% mutate(FITTED = kinh_MLE / (kinh_MLE + INH))

ggplot(inh_fit_data, aes(x = INH)) + geom_point(aes(y = RATE), color = "RED") + geom_line(aes(y = FITTED), color = "BLUE") +
  scale_x_log10()

# hill
hill_test <- 2^seq(-2, 2, by = 0.05)
kinetic_params <- expand.grid(K = k_test, N = hill_test)
kinetic_params_SS <- data.frame(kinetic_params, TSS = 
                                  mapply(function(x,y){
                                   inh_eval(x,y)
                                  }, x = kinetic_params$K, y = kinetic_params$N))
kinetic_params_cast <- acast(kinetic_params_SS, K ~ N, value.var = "TSS")
kinetic_params_cast_log10 <- log10(kinetic_params_cast)
heatmap.2(kinetic_params_cast_log10, Colv = F, Rowv = F, trace = "none", col = greenred(1000))

hillParams <- kinetic_params_SS[which.min(kinetic_params_SS$TSS),]

inh_fit_data <- inh_fit_data %>% mutate(HILL = hillParams$K^hillParams$N / (hillParams$K^hillParams$N + INH^hillParams$N))

ggplot(inh_fit_data, aes(x = INH)) + geom_point(aes(y = RATE), color = "RED") + geom_line(aes(y = FITTED), color = "BLUE") +
  geom_line(aes(y = HILL), color = "GREEN") + scale_x_log10()

# consitutive activity

inh_eval_const <- function(inh_k, inh_hill = 1, frac_const = 1){
  fitted_RATE <- inh_k ^ inh_hill / (inh_k ^ inh_hill + inh_fit_data$INH ^ inh_hill)
  fitted_RATE <- frac_const + (1-frac_const)*fitted_RATE
  sum((inh_fit_data$RATE - fitted_RATE)^2)
  }

kinetic_params_const <- expand.grid(K = k_test, N = 1, frac = seq(0,1,by =0.01))
kinetic_params_const_SS <- data.frame(kinetic_params_const, TSS = 
                                  mapply(function(x,y,z){
                                   inh_eval_const(x,y,z)
                                  }, x = kinetic_params$K, y = kinetic_params$N, z = kinetic_params_const$frac))
hillParams_const <- kinetic_params_const_SS[which.min(kinetic_params_const_SS$TSS),]
kinetic_params_cast <- acast(kinetic_params_const_SS, K ~ frac, value.var = "TSS")
kinetic_params_cast_log10 <- log10(kinetic_params_cast)
heatmap.2(kinetic_params_cast_log10, Colv = F, Rowv = F, trace = "none", col = greenred(1000))

inh_fit_data <- inh_fit_data %>% mutate(RATE_CONST = hillParams_const$frac + (1-hillParams_const$frac)*(hillParams_const$K/ (hillParams_const$K + INH)))

ggplot(inh_fit_data, aes(x = INH, size = 2)) + geom_point(aes(y = RATE), color = "RED") + geom_line(aes(y = FITTED), color = "BLUE") +
  geom_line(aes(y = HILL), color = "GREEN") + geom_line(aes(y = RATE_CONST), color = "ORANGE") + scale_x_log10() + 
  scale_size_identity()

kinetic_params_const_hill <- expand.grid(K = seq(0.05, 2, by = 0.02), N = seq(0.4, 2, by = 0.05), frac = seq(0,1,by =0.01))
kinetic_params_const_hill_SS <- data.frame(kinetic_params_const_hill, TSS = 
                                  mapply(function(x,y,z){
                                   inh_eval_const(x,y,z)
                                  }, x = kinetic_params_const_hill$K, y = kinetic_params_const_hill$N, z = kinetic_params_const_hill$frac))
kinetic_params_const_hill_SS <- kinetic_params_const_hill_SS[which.min(kinetic_params_const_hill_SS$TSS),]

inh_fit_data <- inh_fit_data %>% mutate(RATE_CONST_HILL = kinetic_params_const_hill_SS$frac + (1-kinetic_params_const_hill_SS$frac)*(kinetic_params_const_hill_SS$K^kinetic_params_const_hill_SS$N/ (kinetic_params_const_hill_SS$K^kinetic_params_const_hill_SS$N + INH^kinetic_params_const_hill_SS$N)))


ggplot(inh_fit_data %>% group_by(INH) %>% mutate(Rep = 1:length(INH)), aes(x = INH, size = 2)) + geom_bar(aes(y = RATE, group = Rep), stat = "identity", position = "dodge", color = "black", width = 0.2, fill = "cornflowerblue") +
  geom_line(aes(y = RATE_CONST), color = "chartreuse3") +
  geom_line(aes(y = RATE_CONST_HILL), color = "RED") +
  geom_line(aes(y = FITTED), color = "CORAL") +
  geom_point(aes(y = RATE_CONST), fill = "chartreuse3", color = "black", shape = 21, size = 4) +
  geom_point(aes(y = RATE_CONST_HILL), fill = "RED", color = "black", shape = 21, size = 4) +
  geom_point(aes(y = FITTED), fill = "CORAL", color = "black", shape = 21, size = 4) +
  scale_x_log10("[PHE] (mM)", breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)) + 
  scale_size_identity() + barplot_theme + expand_limits(y = c(0,1)) + scale_y_continuous("Fractional Activity", label = percent_format(), breaks = c(0, 0.25, 0.5, 0.75, 1))
ggsave("PHEfittedKinetics.pdf", height = 8, width = 8)



# statistical model comparison

# log-likelihood of MM-kinetics
# log-likelihood of consitutitive activity
# log-likelihood of constitutive activity + hill

nparams <- data.frame(Model = c("FITTED", "HILL", "RATE_CONST", "RATE_CONST_HILL"), npar = c(1, 2, 2, 3))

inh_model_evaluation <- inh_fit_data %>% gather("Model", "Activity", FITTED:RATE_CONST_HILL) %>% left_join(nparams)
inh_model_evaluation <- inh_model_evaluation %>% group_by(Model) %>% mutate(residualSD = sqrt(sum((RATE - Activity)^2)/(n() - npar)))

#sqrt(sum((inh_fit_data$RATE - inh_fit_data$FITTED)^2)/(nrow(inh_fit_data) - 1))

inh_model_evaluation <- inh_model_evaluation %>% mutate(normalDensity = dnorm(RATE, Activity, residualSD, log = T)) %>% summarize(logL = sum(normalDensity))

# significance of consitutive activity

# significance of constitutive activity and hill

pchisq(2 * (inh_model_evaluation$logL[inh_model_evaluation$Model == "RATE_CONST"] - inh_model_evaluation$logL[inh_model_evaluation$Model == "FITTED"]), df = 1, lower.tail = F) # < 1e-4
pchisq(2 * (inh_model_evaluation$logL[inh_model_evaluation$Model == "RATE_CONST_HILL"] - inh_model_evaluation$logL[inh_model_evaluation$Model == "RATE_CONST"]), df = 1, lower.tail = F) # < 1e-13

kinh_MLE # MLE
hillParams_const # constitutive activity
kinetic_params_const_hill_SS # constitutive + Hill

plot(pchisq(q = seq(0, 20, by = 0.1), df = 1, lower.tail = F) ~ seq(0, 20, by = 0.1))

