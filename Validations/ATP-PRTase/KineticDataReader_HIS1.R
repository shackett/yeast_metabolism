setwd("~/Desktop/Rabinowitz/FBA_SRH/Validations/ATP-PRTase")

options(stringsAsFactors = F)

library(reshape2)
library(data.table)
library(ggplot2)
library(broom)

kineticData <- read.delim('ATPprtase_10_20_14_A.txt', fill = NA, header = F) 

plate_starts <- grep('^Read', kineticData[,1])
plates <- list()

for(a_plate in 1:length(plate_starts)){
  
  if(a_plate == length(plate_starts)){
   plate_data <- kineticData[plate_starts[a_plate]:nrow(kineticData),]
  }else{
    plate_data <- kineticData[plate_starts[a_plate]:(plate_starts[a_plate + 1]-1),]
  }
  
  header_rows <- apply(is.na(plate_data), 1, sum) != 0
    
  plate_matrix <- plate_data[!header_rows,1:13]
  plate_matrix[1,1] <- "Row"
  colnames(plate_matrix) <- plate_matrix[1,]; plate_matrix <- plate_matrix[-1,]
  
  plate_summary <- melt(plate_matrix, variable.name = "Col", value.name = "OD", id.vars = "Row")
  plate_summary <- data.frame(plate_ID = last(plate_data[header_rows,1]), plate_summary)
  
  plates[[a_plate]] <- plate_summary
  
}
plates_aggregate <- do.call("rbind", plates)


plate_IDs <- sapply(plates, function(x){x$plate_ID[1]})
plate_ID_info <- data.frame(ID = plate_IDs, type = NA, time = NA)
plate_ID_info$type <- 'Kinetic'
plate_ID_info$type[grep('Pathlength', plate_ID_info$ID)] <- "PL"

plate_ID_info$time[plate_ID_info$type == "Kinetic"] <- regmatches(plate_ID_info$ID, regexpr('[0-9]{1}:[0-9]{2}:[0-9]{2}', plate_ID_info$ID))
formatted_time <- as.POSIXlt(plate_ID_info$time[plate_ID_info$type == "Kinetic"], format = "%H:%M:%S")

plate_ID_info$cumTime <- NA
plate_ID_info$cumTime[plate_ID_info$type == "Kinetic"] <- formatted_time - formatted_time[1]

plates_aggregate$secondsElapsed <- plate_ID_info$cumTime[chmatch(plates_aggregate$plate_ID, plate_ID_info$ID)]

ODmatrix <- acast(plates_aggregate[!is.na(plates_aggregate$secondsElapsed),], Row + Col ~ secondsElapsed, value.var = "OD")
OD_fit <- rep(NA, nrow(ODmatrix)); names(OD_fit) <- rownames(ODmatrix)

for(a_well in 1:nrow(ODmatrix)){
  #OD_fit[a_well] <- unname(lm(ODmatrix[a_well,] ~ as.numeric(colnames(ODmatrix)))$coef[2])
  OD_fit[a_well] <- unname(lm(ODmatrix[a_well,] ~ as.numeric(colnames(ODmatrix)) + I(as.numeric(colnames(ODmatrix))^2))$coef[2])  
}

### Inhibitor concentration ###

prpp_conc <- melt(matrix(c(rep(c(rep(1, 4), rep(0.04, 4), rep(0.2, 4), rep(0.008, 4), rep(0.4, 16)), 2), rep(0.4, 32)), nrow = 8, ncol = 12))
colnames(prpp_conc) <- c("Row", "Col", "PRPP")
atp_conc <- melt(matrix(c(rep(c(rep(0.8, 16), rep(10, 4), rep(0.4, 4), rep(2, 4), rep(0.08, 4)), 2), rep(0.8, 32)), nrow = 8, ncol = 12))
colnames(atp_conc) <- c("Row", "Col", "ATP")
his_conc <- melt(matrix(c(c(rep(c(5, 1, 0.2, 0.04), times = 8), rep(0, 32)), rep(c(5, 1, 0.2, 0.04), times = 8)), ncol = 12, nrow = 8))
colnames(his_conc) <- c("Row", "Col", "HIS")
AMP_conc <- melt(matrix(c(rep(0, 32), rep(c(1, 0.2, 0.04, 0.008), times = 8), rep(rep(c(1, 0.04, 0.2, 0.008), each = 4), times = 2)), ncol = 12, nrow = 8))
colnames(AMP_conc) <- c("Row", "Col", "AMP")

effector_conc <- merge(merge(prpp_conc, atp_conc), merge(his_conc, AMP_conc))
effector_conc$Row <- LETTERS[effector_conc$Row]

### Add pathlengths ###

pathlength_plate <- plates[plate_ID_info$type == "PL"][[1]]
colnames(pathlength_plate)[colnames(pathlength_plate) == "OD"] <- "PL"
pathlength_plate <- pathlength_plate[,colnames(pathlength_plate) %in% c("Row", "Col", "PL")]

effector_conc <- merge(effector_conc, pathlength_plate)

effector_conc$Well <- paste(effector_conc$Row, effector_conc$Col, sep = "_")
effector_conc$Activity <- OD_fit[chmatch(effector_conc$Well, names(OD_fit))]
effector_conc$n_Activity <- effector_conc$Activity/effector_conc$PL

# Effects on initial rate versus single inhibitors

ggplot(effector_conc[effector_conc$HIS != 0 & effector_conc$AMP == 0 & effector_conc$ATP != 0.8,], aes(x = 1/ATP, y = 1/n_Activity, color = factor(HIS))) + geom_point() +
  scale_x_log10() + scale_y_log10()
ggplot(effector_conc[effector_conc$HIS != 0 & effector_conc$AMP == 0 & effector_conc$PRPP != 0.4,], aes(x = 1/PRPP, y = 1/n_Activity, color = factor(HIS))) + geom_point() +
  scale_x_log10() + scale_y_log10()
ggplot(effector_conc[effector_conc$HIS == 0 & effector_conc$AMP != 0 & effector_conc$ATP != 0.8,], aes(x = 1/ATP, y = 1/n_Activity, color = factor(AMP))) + geom_point() +
  scale_x_log10() + scale_y_log10()
ggplot(effector_conc[effector_conc$HIS == 0 & effector_conc$AMP != 0 & effector_conc$PRPP != 0.4,], aes(x = 1/PRPP, y = 1/n_Activity, color = factor(AMP))) + geom_point() +
  scale_x_log10() + scale_y_log10()

# Effect with both inhibitors

ggplot(effector_conc[effector_conc$HIS != 0 & effector_conc$AMP != 0,], aes(x = factor(AMP), y = factor(HIS), fill = n_Activity)) + geom_tile()

# Looking at inhibition as product accumulates

OD_pl <- merge(data.frame(Well = rownames(ODmatrix)), effector_conc[,c('Well', 'PL')], sort = F)
OD_pl_corrected <- (ODmatrix / OD_pl$PL)

### plot time-course

OD_timecourse <- melt(OD_pl_corrected) %>% tbl_df() %>% dplyr::select(Well = Var1, Seconds = Var2, OD = value) %>%
  left_join(effector_conc) %>% dplyr::select(Well, Seconds, OD, PRPP, ATP, HIS, AMP) %>% group_by(Well) %>% mutate(ODchange = OD - OD[Seconds == 0])

ggplot(OD_timecourse, aes(x = Seconds, y = ODchange, group = Well, color = PRPP)) + geom_line() +
  facet_grid(AMP ~ HIS)


ggplot(OD_timecourse %>% filter(HIS == 0), aes(x = Seconds, y = ODchange, group = Well, color = AMP)) + geom_line() +
  facet_grid(PRPP ~ ATP)
ggplot(OD_timecourse %>% filter(AMP == 0), aes(x = Seconds, y = ODchange, group = Well, color = HIS)) + geom_line() +
  facet_grid(PRPP ~ ATP)

ggplot(OD_timecourse %>% filter(HIS != 0 & AMP != 0), aes(x = Seconds, y = ODchange, group = Well, color = factor(AMP))) + geom_line() +
  facet_grid(~ HIS) + barplot_theme + scale_y_continuous(expression(Delta ~ "OD"))
ggsave("AMPtimecourse.pdf", height = 6, width = 6)

ggplot(OD_timecourse %>% filter(HIS != 0 & AMP != 0), aes(x = Seconds, y = ODchange, group = Well, color = factor(HIS))) + geom_line() +
  facet_grid(~ AMP) + barplot_theme + scale_y_continuous(expression(Delta ~ "OD"))
ggsave("HIStimecourse.pdf", height = 6, width = 6)


ggplot(OD_timecourse %>% filter(HIS != 0 & AMP != 0), aes(x = Seconds, y = ODchange, group = Well, color = factor(HIS))) + geom_line() +
  facet_grid(~ AMP) + barplot_theme + scale_y_continuous(expression(Delta ~ "OD"))

#### Reanalysis ignoring AMP as its effect seems to be very minor ###

regulator_TC <- OD_timecourse %>% filter(PRPP == 0.4, ATP == 0.8) %>% dplyr::select(-PRPP, -ATP, -AMP) %>% mutate(HIS = factor(HIS, levels = c("0.04", "0.2", "1", "5"))) %>%
  mutate(netFlux = ODchange/3600*1e6) # umolar
# use quadratic regression to fit curve

quadratic_reg_TC <- regulator_TC %>% group_by(HIS) %>% do(tidy(lm(netFlux ~ Seconds + I(Seconds^2) + I(Seconds^3), data=.))) %>% ungroup()

fitted_TC <- lapply(unique(quadratic_reg_TC$HIS), function(x){
  seconds_test <- range(unique(regulator_TC$Seconds))
  seconds_test <- seq(seconds_test[1], seconds_test[2], by = 10)
  
  param_subset <- quadratic_reg_TC %>% filter(HIS == x)
  reg_param <- c(
    param_subset %>% filter(term == "(Intercept)") %>% dplyr::select(estimate) %>% unlist() %>% unname(),
    param_subset %>% filter(term == "Seconds") %>% dplyr::select(estimate) %>% unlist() %>% unname(),
    param_subset %>% filter(term == "I(Seconds^2)") %>% dplyr::select(estimate) %>% unlist() %>% unname(),
    param_subset %>% filter(term == "I(Seconds^3)") %>% dplyr::select(estimate) %>% unlist() %>% unname()
  )
  
  time_pts <- data.frame(HIS = x, Seconds = seconds_test) %>% tbl_df()
  time_pts <- time_pts %>% mutate(netFlux = reg_param[1] + Seconds*reg_param[2] + Seconds^2*reg_param[3] + Seconds^3*reg_param[4],
                                  Flux = reg_param[2] + 2*reg_param[3]*Seconds + 3*reg_param[4]*Seconds^2)
  
})
fitted_TC <- do.call("rbind", fitted_TC)

barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                       panel.background = element_rect(fill = "80"), legend.position = "bottom", 
                       axis.ticks = element_line(color = "black", size = 1),
                       axis.text = element_text(color = "black", size = 20),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1)
                       )

his1plots <- list()

his1plots[["cumulative"]] <- ggplot() + geom_line(data = regulator_TC, aes(x = Seconds, y = netFlux, group = Well, color = HIS)) +
  geom_line(data = fitted_TC, aes(x = Seconds, y = netFlux, group = HIS, color = HIS), color = "black", size = 2.6) +
  geom_line(data = fitted_TC, aes(x = Seconds, y = netFlux, color = HIS), size = 1.7) +
  barplot_theme + scale_y_continuous(expression(Delta ~ mu ~ "M")) +
  scale_color_brewer("Histidine [mM]", palette = "YlOrRd") +
  scale_x_continuous("Seconds Elapsed", breaks = seq(0, 3500, by = 500))
  
his1plots[["flux"]] <- ggplot() + geom_line(data = fitted_TC, aes(x = Seconds, y = Flux, color = HIS), size = 3) +
  barplot_theme + scale_y_continuous(expression("Flux ("~mu ~ "M/s)")) +
  scale_color_brewer("Histidine [mM]", palette = "YlOrRd") +
  scale_x_continuous("Seconds Elapsed", breaks = seq(0, 3500, by = 500))
  
library(gridExtra)

pdf("HIStimecourse.pdf", height = 8, width = 14)
grid.arrange(his1plots[["cumulative"]], his1plots[["flux"]], ncol = 2)
dev.off()




#PR-ATP (ε290 = 3600 / M / cm)
# A = e*l*c
# delta A = e*l* delta c

conc_pl_corrected <- OD_pl_corrected/3600*1000 # mmolar
# determine flux @ different concentrations of accumulated product
product_conc <- c(0, 0.003, 0.006, 0.009, 0.012, 0.015)
conc_fit <- matrix(NA, ncol = length(product_conc), nrow = nrow(conc_pl_corrected))
conc_fit_coefs <- matrix(NA, ncol = 3,  nrow = nrow(conc_pl_corrected))
colnames(conc_fit_coefs) <- c("int", "slope", "quad")
  
for(a_well in 1:nrow(conc_pl_corrected)){
  
  well_data <- data.frame(conc = conc_pl_corrected[a_well,], time_elapsed = as.numeric(colnames(ODmatrix))/60)
  well_qr <- lm(data = well_data, conc ~ time_elapsed + I(time_elapsed^2))
  rate_calc_concentrations <- product_conc[max(well_qr$fitted - well_qr$coef[1]) > product_conc]
  # determine time when product accumulation threshold was reached
  rate_calc_times <- sapply(rate_calc_concentrations, function(x){
    min(well_data$time_elapsed[fitted(well_qr) - coef(well_qr)[1] > x])
  })
  
  conc_fit[a_well, 1:length(rate_calc_concentrations)] <- coef(well_qr)[2] + 2*coef(well_qr)[3]*rate_calc_times
  conc_fit_coefs[a_well,] <- coef(well_qr)
  
}
colnames(conc_fit) <- product_conc
conc_fit <- cbind(conc_fit, conc_fit_coefs)
conc_fit <- as.data.frame(conc_fit)
conc_fit$Well <- rownames(conc_pl_corrected)

effector_conc_withProd <- merge(effector_conc, conc_fit)
effector_conc_predictors <- melt(effector_conc_withProd, measure.vars = c(as.character(product_conc), colnames(conc_fit_coefs)), variable.name = "PR_ATP", value.name = "flux")

effector_conc_melt <- effector_conc_predictors[effector_conc_predictors$PR_ATP %in% as.character(product_conc),]
effector_lm_melt <- effector_conc_predictors[!(effector_conc_predictors$PR_ATP %in% as.character(product_conc)),]

# Effects on initial rate versus single inhibitors

ggplot(effector_conc_melt[effector_conc_melt$HIS != 0 & effector_conc_melt$AMP == 0 & effector_conc_melt$ATP != 0.8,], aes(x = 1/ATP, y = 1/flux, color = factor(HIS))) + 
  geom_point() + facet_wrap(~ PR_ATP) + geom_smooth(method = "lm", se = F) + expand_limits(x = 0, y = 0)
ggplot(effector_conc_melt[effector_conc_melt$HIS != 0 & effector_conc_melt$AMP == 0 & effector_conc_melt$PRPP != 0.4,], aes(x = 1/PRPP, y = 1/flux, color = factor(HIS))) +
  geom_point() + facet_wrap(~ PR_ATP) + geom_smooth(method = "lm", se = F)
ggplot(effector_conc_melt[effector_conc_melt$HIS == 0 & effector_conc_melt$AMP != 0 & effector_conc_melt$ATP != 0.8,], aes(x = 1/ATP, y = 1/flux, color = factor(AMP))) +
  geom_point() + facet_wrap(~ PR_ATP) + geom_smooth(method = "lm", se = F)
ggplot(effector_conc_melt[effector_conc_melt$HIS == 0 & effector_conc_melt$AMP != 0 & effector_conc_melt$PRPP != 0.4,], aes(x = 1/PRPP, y = 1/flux, color = factor(AMP))) +
  geom_point() + facet_wrap(~ PR_ATP) + geom_smooth(method = "lm", se = F)

# Effect with both inhibitors

barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                       panel.background = element_rect(fill = "gray80"), legend.position = "right", 
                       axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "black"),
                       axis.text = element_text(color = "black"), axis.text.x = element_text(size = 10, angle = 90, hjust = 0, vjust = 0.5),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), legend.title=element_blank()
                       )

ggplot(effector_conc_melt[effector_conc_melt$HIS != 0 & effector_conc_melt$AMP != 0,], aes(x = factor(AMP), y = factor(HIS), fill = flux)) + geom_tile() + facet_wrap(~ PR_ATP) +
  barplot_theme + scale_fill_continuous(low = "blue", high = "yellow") + 
  scale_x_discrete('[AMP] (mM)') + scale_y_discrete('[HIS] (mM)') + ggtitle('Effect of inhibitors on flux at \n different product levels')
ggsave("HisAmpHM.pdf", height = 8, width = 8)


library(dplyr)
effector_lm <- tbl_df(effector_lm_melt[effector_lm_melt$HIS != 0 & effector_lm_melt$AMP != 0 & effector_lm_melt$PR_ATP == "quad",])

effector_summary <- effector_lm %>% arrange(HIS, AMP) %>% mutate(curvature = abs(flux)) %>% select(HIS, AMP, curvature)
effector_summary <- effector_summary %>% group_by(HIS, AMP) %>% mutate(name = paste0(paste0("HIS", round(HIS, 2)),":",paste0("AMP", round(AMP, 3)),":R",c(1:length(HIS))))
                                                                
ggplot(effector_summary, aes(x = name, y = curvature, fill = factor(HIS))) + geom_bar(stat = "identity", width = 0.8)

###

barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                       panel.background = element_rect(fill = "gray80"), legend.position = "right", 
                       axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "black"),
                       axis.text = element_text(color = "black"), axis.text.x = element_text(size = 10, angle = 90, hjust = 0, vjust = 0.5),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), legend.title=element_blank()
                       )


effector_fixed_PR_ATPconc <- tbl_df(effector_conc_melt[effector_conc_melt$HIS != 0 & effector_conc_melt$AMP != 0 & effector_conc_melt$PR_ATP == "0.006",])
effector_fixed_summary <- effector_fixed_PR_ATPconc %>% arrange(HIS, AMP) %>% select(HIS, AMP, flux)
effector_fixed_summary <- effector_fixed_summary %>% group_by(HIS, AMP) %>% mutate(name = paste0(paste0("HIS", round(HIS, 2)),":",paste0("AMP", round(AMP, 3)),":R",c(1:length(HIS))))

ggplot(effector_fixed_summary, aes(x = name, y = flux, fill = factor(HIS))) + geom_bar(stat = "identity", width = 0.8) + barplot_theme +
  scale_fill_brewer(palette = "Set1") + scale_x_discrete("Sample") + scale_y_continuous("Flux (mM/min)") +
  ggtitle("Flux at 6uM PR-ATP")
ggsave("InhibitorRole_6uM.pdf", height = 6, width = 8)

effector_fixed_PR_ATPconc2 <- tbl_df(effector_conc_melt[effector_conc_melt$HIS != 0 & effector_conc_melt$AMP != 0 & effector_conc_melt$PR_ATP == "0",])
effector_fixed_summary2 <- effector_fixed_PR_ATPconc2 %>% arrange(HIS, AMP) %>% select(HIS, AMP, flux)
effector_fixed_summary2 <- effector_fixed_summary2 %>% group_by(HIS, AMP) %>% mutate(name = paste0(paste0("HIS", round(HIS, 2)),":",paste0("AMP", round(AMP, 3)),":R",c(1:length(HIS))))

ggplot(effector_fixed_summary2, aes(x = name, y = flux, fill = factor(HIS))) + geom_bar(stat = "identity", width = 0.8) + barplot_theme +
  scale_fill_brewer(palette = "Set1") + scale_x_discrete("Sample") + scale_y_continuous("Flux (mM/min)") +
  ggtitle("Initial Flux")
ggsave("InhibitorRole_initial.pdf", height = 6, width = 8)

effector_fixed_combo <- rbind(effector_fixed_summary2 %>% mutate(time = "Initial rate"), effector_fixed_summary %>% mutate(time = "At 6uM PR-ATP"))
effector_fixed_combo$time <- factor(effector_fixed_combo$time, levels = unique(effector_fixed_combo$time))

ggplot(effector_fixed_combo, aes(x = name, y = flux*1000, fill = factor(HIS))) + geom_bar(stat = "identity", width = 0.8) + barplot_theme +
  scale_fill_brewer(palette = "Set1") + scale_x_discrete("Sample") + scale_y_continuous("Flux (uM/min)") +
  ggtitle("Initial Flux") + facet_grid(time ~ ., scale = "free_y")
ggsave("InhibitorRole_combo.pdf", height = 10, width = 8)








