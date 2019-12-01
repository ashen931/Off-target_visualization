#20191101_summarize_off-target_editing.R
# base environment
# last edited: 2019_11_01 Anne Shen

# Need to add:
# -show mismatches and bulges

library(tidyverse)
library(tidyselect)
#install.packages("extrafont")
library(extrafont)
library(gtable)
library(scales)
library(gridExtra)
library(grid)
#install.packages("effsize")
library(effsize)
#font_import()
loadfonts()
options(scipen=999) #turn off scientific notation



source("20191101_summarize_off-target_editing_functs.R")

setwd("/Users/anneshen/Documents/local_working/local_Jing_BE/1620/rhAMPSeq_1620_OT_integrated_2")

#setwd("~/Documents/local_working/off_target_crispresso")
#######################################################################################################################

#get all summary file names
list_summary_files <- list.files(pattern = "summary.csv", recursive = TRUE)
#get mock summary file names
mock_summary_files <- grep("Mock", list_summary_files, value = TRUE)

#get edited summary file names
list_1620_files <- grep("1620", list_summary_files, value = TRUE)
#list_1620_2x_files <- grep("2x", list_summary_files, value = TRUE)

#get master table linking off-target names and sequences
# The ref_seqs.csv table should have 4 columns: ot_id, amplicon_sequence, guide_sequence, and pam
master_guide_tb <- read.csv("../1620_ref_seqs_tb.csv") %>%
  filter(grepl("OT", ot_id)) %>%
  rename(off_target = ot_id) %>%
  filter(!duplicated(off_target))

master_guide_tb$off_target<- gsub("_(?=[1-9]{1}[0-9]{1}$)", "_0", 
                                     gsub("_(?=[1-9]{1}$)", "_00", master_guide_tb$off_target, perl = TRUE),
                                     perl = TRUE)

############ plot the mock off-target editing summaries ############
all_mock_tb <- get_formatted_summary(mock_summary_files, master_guide_tb, condition = "control")

all_1620_tb <- get_formatted_summary(list_1620_files, master_guide_tb, condition = "edited") 

all_samples_from_file <- rbind(all_mock_tb, all_1620_tb) %>%
  transform(sample = gsub("all_", "", sample))

#View(all_samples_from_file[grepl("1620", all_samples_from_file$off_target),])

unique_samples_idx <- unique(all_samples_from_file$sample)
unique_samples <- unique_samples_idx[which(!is.na(unique_samples_idx))]
n_control <- length(mock_summary_files)
n_edited <- length(list_1620_files)
n_samples <- n_control + n_edited

all_samples <- data.frame(off_target = rep(master_guide_tb$off_target, times = n_samples, each = 2),
                          amplicon_sequence = rep(master_guide_tb$amplicon_sequence, times = n_samples, each = 2),
                          guide_sequence = rep(master_guide_tb$guide_sequence, times = n_samples, each = 2),
                          pam = rep(master_guide_tb$pam, times = n_samples, each = 2),
                          indel = rep(c("Unedited", "Edited"), nrow(master_guide_tb) * n_samples),
                          sample = rep(c(unique_samples), times = 1, each = nrow(master_guide_tb) * 2),
                          condition = rep(c(rep("control", n_control), rep("edited", n_edited)), times = 1, each = nrow(master_guide_tb) * 2),
                          name_seq = rep(paste(paste(master_guide_tb$off_target, 
                                                 master_guide_tb$guide, sep = "  "),
                                           master_guide_tb$pam, sep = "  "),times = n_samples, each = 2 ))

all_samples_tb <- left_join(all_samples, all_samples_from_file, by = names(all_samples)) %>%
  transform(sample = gsub("merged_", "", gsub("(?<=Mock).*", "", 
                          gsub("(?<=[1,2]x).*", "", sample, perl = TRUE),
                          perl = TRUE))) %>%
  mutate(group = paste(condition, sample, sep = " "))


# all_samples_tb <- all_samples_tb[order(all_samples_tb$group, decreasing = TRUE),]
# all_samples_tb$group <- factor(all_samples_tb$group, 
#                               levels = unique(all_samples_tb$group), ordered = TRUE)

#View(all_samples_tb[grepl("1620", all_samples$off_target),])

### 1620
summary_tbs_1620 <- get_table_by_guide("1620", all_samples_tb) %>%
  filter(grepl("1620", sample) | grepl("Mock", sample))

edited_summary_1620 <- summary_tbs_1620 %>%
  filter(indel == "Edited") 

### statistical test
# Compare % edited in Mock v. Edited samples
#get list of off-targets
off_targets <- unique(edited_summary_1620$off_target)
#generate table of Edited v. Mock editing frequency results
ot_ttest_tb <- get_ttest_table(off_targets, edited_summary_1620)
sig_ots <- ot_ttest_tb$off_target[which(ot_ttest_tb$significant)] %>% droplevels()

write.csv(ot_ttest_tb, "20191130_1620_OTs_ttest.csv", row.names = FALSE, quote = FALSE)

for(n in seq(1, nrow(edited_summary_1620))){
  
  edited_summary_1620$name_seq <- as.character(edited_summary_1620$name_seq)
  
  if(edited_summary_1620$off_target[n] %in% sig_ots){
    edited_summary_1620$name_seq[n] <- paste(edited_summary_1620$name_seq[n],
                                             "*", 
                                             sep = " ")
  }else{
    edited_summary_1620$name_seq[n] <- paste(edited_summary_1620$name_seq[n],
                                             " ", 
                                             sep = " ")
  }
}

edited_summary_1620 <- edited_summary_1620 %>%
  order_alpha("name_seq", decreasing_bool = TRUE)


### generate dotplot fill and color values
fill <- list("Group", c("control BE0108-1-Mock" = "firebrick",
                        "control BE1215-1-Mock" = "royalblue",
                        "control EP1116-Mock" = "darkgoldenrod2",
                        #"" = "tomato",
                        "edited BE0108-9-1620-2x" = "coral1",
                        "edited BE1215-2-1620-1x" = "cornflowerblue",
                        "edited BE1215-3-1620-2x" = "skyblue",
                        "edited BE1116-1620-1x" = "gold",
                        "edited BE1116-1620-2x" = "lightgoldenrod"))
color <- list("Group", c("control BE0108-1-Mock" = "firebrick",
                         "control BE1215-1-Mock" = "royalblue",
                         "control EP1116-Mock" = "darkgoldenrod2",
                         #"" = "tomato",
                         "edited BE0108-9-1620-2x" = "coral1",
                         "edited BE1215-2-1620-1x" = "cornflowerblue",
                         "edited BE1215-3-1620-2x" = "skyblue",
                         "edited BE1116-1620-1x" = "gold",
                         "edited BE1116-1620-2x" = "lightgoldenrod"))
shape <- list("Group", c("control BE0108-1-Mock" = 16,
                         "control BE1215-1-Mock" = 16,
                         "control EP1116-Mock" = 16,
                         #"" = 1,
                         "edited BE0108-9-1620-2x" = 1,
                         "edited BE1215-2-1620-1x" = 1,
                         "edited BE1215-3-1620-2x" = 1,
                         "edited BE1116-1620-1x" = 1,
                         "edited BE1116-1620-2x" = 1))


#### generate editing dotplot
dotplot_1620 <- make_compiled_OT_editing_dotplot(edited_summary_1620, fill, color, shape, 6, 0.5)



#### get the total counts from "Unedited" reads 
total_counts_1620 <-get_table_by_guide("1620", summary_tbs_1620 ) %>%
  filter(indel == "Unedited" ) %>%
  mutate(total_reads = (reads *100) / frequency ) %>%
  order_alpha("group", decreasing_bool = TRUE)


#generate coverage heatmap 
heatmap_1620 <- make_compiled_OT_coverage_heatmap(total_counts_1620, guide_font_size = 6, 
                                                  group_font_size = 8, tile_font_size = 2)


### make and save composite plots
### 1620
sg1620_composite<- make_composite_grobPlot(heatmap_1620, dotplot_1620, plot_height_inch = c(7.5, 0.1), 
                        title = "sg1620 rhAmpSeq 1-60 2019_11_29", 
                        note_text = "NOTE: rhAMPSeq amplification failed for OT_001 (on-target).", 
                        note_hjust = -0.75)

save_composite_plot("20191129_1620_off_targets_1-60_Rplot", sg1620_composite, 14.5, 8.5)






############## FUNCTIONS ######################################################################
split_summary_table <- function(edited_summary_tb){
  
  #edited_summary_tb <- summary_tbs_1450
  
  #order summary table by OT number
  edited_summary_tb <- edited_summary_tb[order(edited_summary_tb$ot_num),] %>%
    put_on_target_first()
  
  #get number of samples in table
  n_samples <- length(unique(edited_summary_tb$sample))
  #get (number of guides)/2 (rounded down)
  half_guides <- round(length(unique(edited_summary_tb$ot_num)) / 2, 0)
  
  #split summary table into 2 tables, each containing half of guides
  edited_summary_1 <- edited_summary_tb[seq(1, (half_guides * n_samples * 2)),] %>%
    order_alpha("name_seq", decreasing_bool = FALSE)
  edited_summary_2 <- edited_summary_tb[seq((half_guides * n_samples * 2)+1, nrow(edited_summary_tb)),] %>%
    order_alpha("name_seq", decreasing_bool = FALSE)
  
  #return list of two tables
  return(list(edited_summary_1, edited_summary_2))
}


