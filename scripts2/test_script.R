
library(readxl)
library(grid)
library(gridExtra)
library(kableExtra)
library(formattable)
library(tibble)
library(dplyr)
#library(xcms)
library(Matrix)
#library(MSnbase)
#library(MsExperiment)
#library(doParallel)
library(scales)
library(ggpmisc)

## defining functions and parameters in this chunk
source("/home/analysis/scripts/test_functions.R")

print("Functions were read")

##########################################################################
# defining other user inputs
#########################################################################
plot_eics_peakscores <- Sys.getenv("INCLUDE_EICS_PEAKSCORES", "No")

Study <- Sys.getenv("STUDY", "LTRS")

path_output <- Sys.getenv("OUTPUT_PATH")

mz_rt_input <- Sys.getenv("MZ_RT_LIST", "/home/analysis/data")

metadata_input <- Sys.getenv("METADATA", "/home/analysis/data")

##########################################################################
# reading in the mz RT list
#########################################################################

mzml_path <- Sys.getenv("MZML_PATH", "/home/analysis/data") 

mz_rt_list <- read_excel(mz_rt_input)


if (all(grepl("H+", mz_rt_list$adduct))) {
  ionization <- "Pos"
  mz_rt_list <- as.data.frame(mz_rt_list)
  mz_rt_list$protonated_mass <- mz_rt_list$m_z + 1.00727
  mz_rt_list$error <- (50*mz_rt_list$protonated_mass)/1000000 # 50 ppm calculations
  mz_rt_list$mz1 <- mz_rt_list$protonated_mass - mz_rt_list$error
  mz_rt_list$mz2 <- mz_rt_list$protonated_mass + mz_rt_list$error
  
  ## calculating RT window
  mz_rt_list$min_RT <- mz_rt_list$RT - 0.5
  mz_rt_list$max_RT <- mz_rt_list$RT + 0.5
  
  
}else if (all(grepl("H-", mz_rt_list$adduct))) {
  ionization <- "Neg"
  mz_rt_list <- as.data.frame(mz_rt_list)
  mz_rt_list$deprotonated_mass <- mz_rt_list$m_z - 1.00727
  mz_rt_list$error <- (50*mz_rt_list$protonated_mass)/1000000 # 50 ppm calculations
  mz_rt_list$mz1 <- mz_rt_list$protonated_mass - mz_rt_list$error
  mz_rt_list$mz2 <- mz_rt_list$protonated_mass + mz_rt_list$error
  
  ## calculating RT window
  mz_rt_list$min_RT <- mz_rt_list$RT - 0.5
  mz_rt_list$max_RT <- mz_rt_list$RT + 0.5
  
}else{
  print("Please specify either H+ or H- in adduct column of mz_RT_list.xlsx")
  stop()
}



##########################################################################
# determining column name from user 
#########################################################################

if (all(grepl("C18", mz_rt_list$column_used, ignore.case=TRUE))) {
  
  column_used <- "C18"
  
}else if (all(grepl("Hilic", mz_rt_list$column_used, ignore.case=TRUE))) {
  
  column_used <- "Hilic"
}else{
  
  print("Please specify either C18 or Hilic under column_used of mz_RT_list.xlsx")
}


##########################################################################
# reading in data and processing it: untargeted peak detection and integration with XCMS
#########################################################################

#datafiles_list <- list.files(mzml_path, recursive = FALSE, full.names = TRUE, pattern = ".mzML")

#if(length(datafiles_list) == 0){
#print(".mzml files not found")
#stop()
#knitr::knit_exit()
#} else{
#print(".mzml files found")
#}


print(paste0("Processing data from ", Study, " study."))




#num_files <- length(datafiles_list)
#print(paste0("Number of files being processed: ", num_files))


#raw_data_processed <- perform_xcms_analysis(column_used)
data_df <- initial_processing(paste0("/home/analysis/data", "/","CIL_std_compound_quant_XCMS.txt"),
                              paste0("/home/analysis/data", "/","datafiles_list_processed_XCMS.txt"),
                              metadata_input)

print(head(data_df))

print("Generating report...")

filtered_data_df <- data_df %>%
  filter(case_when(
    compound_ID == mz_rt_list$compound_ID[1] ~ rt > mz_rt_list$min_RT[1] & rt < mz_rt_list$max_RT[1],
    compound_ID == mz_rt_list$compound_ID[2] ~ rt > mz_rt_list$min_RT[2] & rt < mz_rt_list$max_RT[2],
    compound_ID == mz_rt_list$compound_ID[3] ~ rt > mz_rt_list$min_RT[3] & rt < mz_rt_list$max_RT[3]
    
  ))

filtered_data_df <- filtered_data_df %>% group_by(compound_ID) %>% 
  mutate(median_rt = median(rt)) %>% ungroup()

filtered_data_df <- filtered_data_df %>%
  group_by(compound_ID, sample) %>%  # Group by compound and sample
  mutate(rt_diff = abs(median_rt - rt)) %>%  # Calculate the absolute difference from median rt
  filter(rt_diff == min(rt_diff)) %>%  # Filter for the closest rt value to the median
  filter(sn == max(sn)) %>%          # Keep only the rows with the highest sn value
  ungroup()

filtered_data_df <- left_join(filtered_data_df, mz_rt_list)
# check if group column is defined, and based on that generate the kable tables

if("group" %in% colnames(filtered_data_df)){
  table_xcms <- as.data.frame(table(filtered_data_df$Drug, filtered_data_df$sample_type,
                                    filtered_data_df$group), stringsAsFactors = FALSE)
  colnames(table_xcms)  <- c("Drug", "sample_type", "group", "Freq")
  all_drugs_in_data <- unique(table_xcms$Drug)
  std_3ormore <- filter(table_xcms, sample_type == "std") %>% filter(Freq > 3)
  
}else{
  table_xcms <- as.data.frame(table(filtered_data_df$Drug, filtered_data_df$sample_type), 
                              stringsAsFactors = FALSE)
  colnames(table_xcms)  <- c("Drug", "sample_type","Freq")
  all_drugs_in_data <- unique(table_xcms$Drug)
  std_3ormore <- filter(table_xcms, sample_type == "std") %>% filter(Freq > 3)
}

mz_rt_list <- filter(mz_rt_list, compound_ID %in% unique(std_3ormore$Drug))

mz_rt_list$rt1 <- (mz_rt_list$min_RT*60) - 20
mz_rt_list$rt2 <- (mz_rt_list$max_RT*60) + 20


filtered_data_df <- filtered_data_df %>% group_by(compound_ID) %>% 
  mutate(median_intb = median(intb),
         mean_intb = mean(intb),
         std_deviation = sd(intb)) %>% ungroup()


if (grepl("Pos", ionization)) {
  filtered_data_df <- filtered_data_df %>% group_by(compound_ID) %>% 
    mutate(ppm = ((protonated_mass - mz)*1000000)/protonated_mass) %>% ungroup()
}else{
  filtered_data_df <- filtered_data_df %>% group_by(compound_ID) %>% 
    mutate(ppm = ((deprotonated_mass - mz)*1000000)/deprotonated_mass) %>% ungroup()
}


filtered_data_df <- mutate(filtered_data_df, RunOrder = as.numeric(RunOrder))
filtered_data_df <- mutate(filtered_data_df, RunOrder = as.character(RunOrder))

filtered_data_df$RunOrder <- factor(filtered_data_df$RunOrder, 
                                    levels = sort(unique(data_df$RunOrder)), ordered = TRUE)

RT <- unique((filter(filtered_data_df, compound_ID == unique(std_3ormore$Drug)[1]))$median_rt)
png("/home/analysis/output/RT_drift.png", units = "in", res = 600, width = 12,
    height = 6)
filtered_data_df %>% filter(compound_ID %in% c(unique(std_3ormore$Drug)[1])) %>% 
  ggplot(., aes(x = RunOrder, y = rt)) +
  geom_point(size = 4) +
  theme_classic() +
  geom_hline(yintercept = RT[[1]], col = "blue") +
  geom_hline(yintercept = RT[[1]] + RT[[1]]*0.1) +
  geom_hline(yintercept = RT[[1]] - RT[[1]]*0.1) +
  ggtitle(paste0(Study, " ", unique(std_3ormore$Drug)[1])) +
  theme(text = element_text(size = 15),
        axis.ticks.y = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  labs(y = "Retention Time (min)") +
  scale_x_discrete(labels = function(x) {
    # Only show labels for every other tick
    x[seq_along(x) %% 2 != 1] <- ""
    return(x)
  }, drop = FALSE)
dev.off()



expected_mass <- round(unique((filtered_data_df %>% filter(compound_ID %in% c(unique(std_3ormore$Drug)[1])))$protonated_mass), 4)

adduct <- unique((filtered_data_df %>% filter(compound_ID %in% c(unique(std_3ormore$Drug)[1])))$adduct)

png("/home/analysis/output/ppm.png", units = "in", res = 600, width = 12,
    height = 6)
filtered_data_df %>% filter(compound_ID %in% c(unique(std_3ormore$Drug)[1])) %>% 
  ggplot(., aes(x = RunOrder, y = ppm)) +
  geom_point(size = 4) +
  theme_classic() +
  geom_hline(yintercept = 0, col = "blue") +
  geom_hline(yintercept = 8) +
  geom_hline(yintercept = -8) +
  ggtitle(paste0(Study, " ", unique(std_3ormore$Drug)[1], "; Expected mass - ", expected_mass, " m/z ", "(", adduct, ")")) +
  theme(text = element_text(size = 15),
        axis.ticks.y = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  scale_x_discrete(labels = function(x) {
    # Only show labels for every other tick
    x[seq_along(x) %% 2 != 1] <- ""
    return(x)
  }, drop = FALSE)
dev.off()















