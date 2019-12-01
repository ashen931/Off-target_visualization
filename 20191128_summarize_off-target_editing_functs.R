#20191102_summarize_off-target_editing_functs.R
# base environment
# last edited: 2019_11_02 Anne Shen

# Contains all the functions needed to run 20191007_summarize_off-target_editing.R

# Functions included:
# 
# get_formatted_summary()
#   get_group_summary_table()
# collapse_duplicate_OTs()
#   get_num_dups()
#   get_table_by_guide()
# make_compiled_OT_editing_dotplot()
# make_compiled_OT_editing_heatmap()
# make_compiled_OT_coverage_heatmap()
# make_composite_grobPlot()

# Useses the following packages:
# library(tidyverse)
# library(tidyselect)
# #install.packages("extrafont")
# library(extrafont)
# library(gtable)
# library(scales)
# library(gridExtra)
# library(grid)
#install.packages("effsize")
#library(effsize)
# #font_import()
# loadfonts()
# options(scipen=999) #turn off scientific notation


#Changed sample_name to 1,25
#Changed get_group_summary_table to remove "X" from beginning of OT id if it is there
# 


######################################## FUNCTIONS ########################################
# Generate composite editing summary table of all samples in list_summary_files
get_group_summary_table <- function(pool_tb, sample_name){
  
  #remove "X" from beginning of off-target names (if they begin with numbers or special symbols)
  names(pool_tb) <- gsub("^X", "", names(pool_tb))
  
  unedited <- pool_tb %>% 
    filter(indel == "Unedited") %>%
    select(vars_select(names(pool_tb), contains("_O"))) %>%
    mutate(indel = "Unedited", sample = sample_name)
  
  edited <- pool_tb %>% 
    filter(indel != "Unedited") %>%
    select(vars_select(names(pool_tb), contains("_O"))) %>%
    colSums(na.rm = TRUE)
  
  pool_summary <- rbind(unedited, edited) 
  pool_summary$indel[2] <- "Edited"
  pool_summary$sample[2] <- sample_name
  
  return(pool_summary)
}

# 1. shape the summary table for plotting
# 2. add the off-target motif sequence
get_formatted_summary <- function(list_summary_files, master_guide_tb, condition){
  
  #loop through all the summary files
  for(n in seq(1, length(list_summary_files))){
    
    #read the summary tables and get the names
    summary_tb_raw <- read.csv(list_summary_files[n])
    sample_name <- substr(list_summary_files[n], 1, 
                          min(25, str_length(list_summary_files[n])))
    
    summary_tb <- get_group_summary_table(summary_tb_raw, sample_name)
    
    #gather columns into off_target and read columns
    summary_table_reads <-summary_tb %>%
      select(vars_select(names(summary_tb), -matches("[0-9]$"))) %>%
      gather(key = "off_target", value = "reads", 
             vars_select(names(summary_tb), contains("_reads")))
    
    summary_table_reads$off_target <- gsub("_reads", "", summary_table_reads$off_target)
    
    #gather columns into off_target and read columns
    summary_table_freqs <- summary_tb %>% 
      select(vars_select(names(summary_tb), -ends_with("_reads"))) %>%
      gather(key = "off_target", value = "frequency",
             vars_select(grep("_reads", names(summary_tb), invert = TRUE, value = TRUE),
                         contains("_O")))
    
    summary_table <-full_join(summary_table_freqs, summary_table_reads, 
                              by =c("off_target", "indel", "sample"))
    #remove NA frequencies
    summary_table[is.na(summary_table)] <- 0
    
    if(n == 1){
      all_summary_table <- summary_table
    }else{
      all_summary_table <- rbind(all_summary_table, summary_table)
    }
  }
  
  all_summary_table$off_target <- gsub("_(?=[1-9]{1}[0-9]{1}$)", "_0", 
                            gsub("_(?=[1-9]{1}$)", "_00",  all_summary_table$off_target, perl = TRUE),
                            perl = TRUE)
  
  #add off-target sequence to summary table
  summary_and_seq_tb <- left_join(master_guide_tb, all_summary_table, by = "off_target")
  
  #make all 0-frequency 
  if(length(which(summary_and_seq_tb$frequency == 0)) > 0){
    
    summary_and_seq_tb[which(summary_and_seq_tb$frequency == 0),]$frequency <- 0.001
    
  }
  
  #log10 transform frequency
  summary_and_seq_tb$log_freq <- log10(summary_and_seq_tb$frequency)
  
  #merge off_target name and sequence for plotting
  summary_and_seq_tb$name_seq <- paste(paste(summary_and_seq_tb$off_target, 
                                             summary_and_seq_tb$guide, sep = "  "),
                                       summary_and_seq_tb$pam, sep = "  ")
  
  #order table by off_target
  summary_and_seq_tb <- summary_and_seq_tb[order(summary_and_seq_tb$off_target),]
  
  #add "condition" column to indicate whether this was an edited or control sample
  summary_and_seq_tb$condition <- rep(condition, nrow(summary_and_seq_tb))
  
  return(summary_and_seq_tb)
}

#collapse duplicate off-target samples by taking mean of all frequencies
collapse_duplicate_OTs <- function(summary_table){
  
  #order off_targets for collapsing
  summary_table <- summary_table[order(summary_table$off_target, decreasing = FALSE),]
  #generate is_duplicate column for collapsing
  summary_table$is_duplicate <- duplicated(summary_table$off_target)
  
  #take mean of all frequencies pertaining to one off_target sample
  n <- 1
  while(n < nrow(summary_table)){
    
    #if there is at least one duplicated off_target sample
    if(summary_table$is_duplicate[n+1]){
      
      #get the number of duplicates
      num_dups <- get_num_dups(summary_table$is_duplicate, n)
      #calculate the mean of all duplicated off-targets
      summary_table$frequency[n] <- mean(summary_table$frequency[n:(n+num_dups)])
      #skip to last duplicate row
      n <- n + num_dups
    }
    #increase counter
    n <- n + 1
  }
  
  #filter out duplicate off-targets
  summary_table_colps <- summary_table %>%
    filter(! is_duplicate) %>%
    select(-is_duplicate)
  
  return(summary_table_colps)
}


#returns the number of rows that are duplicates of row [n]
get_num_dups <- function(is_duplicate, n){
  
  #intialize the number of duplicated off-target samples
  # (must be at least 1 for this function to be called)
  num_dups <- 1 
  
  #while is_duplicate[n+num_dups+1] is TRUE
  while(is_duplicate[n+num_dups+1] & (n+num_dups+1 < length(is_duplicate))){ 
    num_dups <- num_dups + 1 #add another duplicate
  }
  
  return(num_dups)
}


# 1.filter all_data_table to get relevant off_targets for the guide
# 2.format the name_seq to the correct order (descending) for the dotplot
get_table_by_guide <- function(guide_name, all_data_table){
  
  summary_tb <- all_data_table %>%
    filter(grepl(guide_name, all_data_table$off_target))
  
  summary_tb$ot_num <- as.numeric(gsub("_", "",
                                       str_extract(summary_tb$name_seq, "_[0-9]{1,3}")))
  
  #added
  #summary_tb[is.na(summary_tb$ot_num),]$ot_num <- 0
  
  #order off-targets (row 1 should be on-target)
  summary_tb <- summary_tb[order(summary_tb$ot_num, na.last = FALSE, decreasing = FALSE),]
  
  
  summary_tb$name_seq <- factor(summary_tb$name_seq, 
                                levels = unique(summary_tb$name_seq), ordered = TRUE)
  
  # summary_tb <- summary_tb %>%
  #   select(-ot_num)
  
  return(summary_tb)
}

# Generate alphabetical levels for a specific column (for plotting purposes)
order_alpha <- function(data_tb, colname, decreasing_bool){
  
  ordered_tb <- data_tb[order(data_tb[,colname], decreasing = decreasing_bool),]
  
  ordered_tb[,colname] <- factor(ordered_tb[,colname], 
                                 levels = unique(ordered_tb[,colname]), ordered = TRUE)
  return(ordered_tb)
}

#Takes a list of unique off-target names and edited summary table and performs a t-test
# comparing the mean editing frequency between Edited and Mock samples. Returns a table
# of off-targets and their corresponding t-test p-values, as well as whether the 
# difference between Edited and Mock is significant.
get_ttest_table <- function(off_targets, edited_summary_tb){
  
  #initalize p-value vector
  p_vals <- c()
  #initalize effect size vector
  d_vals <- c()
  #initalize median vectors
  control_median <- c()
  edited_median <- c()
  
  #calculate p-values with t-test comparing Mock v. Edited samples for each off-target
  for(n in seq(1, length(off_targets))){
    
    #get rows that correspond with the off-target of interest
    ot_only_tb <- edited_summary_tb%>% 
      filter(complete.cases(edited_summary_tb)) %>%
      filter(off_target == off_targets[n])
    
    edited_freq <-ot_only_tb[ot_only_tb$condition == "edited",]$frequency
    control_freq <-ot_only_tb[ot_only_tb$condition == "control",]$frequency
    
    if(length(edited_freq) > 1 & length(control_freq) > 1){
      p_vals <- c(p_vals, 
                  t.test(edited_freq, control_freq)$p.value)
      d_vals <- c(d_vals,
                  cohen.d(edited_freq,control_freq, pooled=TRUE, paired=FALSE,
                                 na.rm=TRUE, hedges.correction=TRUE)$magnitude)
    }else{
      p_vals <- c(p_vals, NA)
      d_vals <- c(d_vals, NA)
    }
    
    #calculate median
    edited_median <- c(edited_median, median(edited_freq, na.rm = TRUE))
    control_median <- c(control_median, median(control_freq, na.rm = TRUE))
  }
  
  #generate data frame with off_targets, the t-test p-value associated with it,
  # and boolean (significant or not)
  ots_pval_tb <- data.frame(off_target = off_targets, 
                            edited_median = edited_median,
                            control_median = control_median,
                            ttest_p_value = p_vals,
                            eff_size = d_vals) %>%
    mutate(significant = ttest_p_value < 0.05)
  
  return(ots_pval_tb[order(ots_pval_tb$off_target),])
}


#Takes a summary table with all samples, all off-targets, and all editing outcomes and 
# generates editing summary dotplot. Also takes fill and color parameters for ggplot2.
# Returns a dotplot that is meant to be placed in a composite graph with a coverage heatmap.
make_compiled_OT_editing_dotplot <- function(summary_tb, fill, color, shape, 
                                             guide_font_size = 12, pointsize = 1.5){
  
  dotplot <- ggplot(data = summary_tb,
                    aes(x = name_seq,
                        y = frequency)) +
    geom_jitter(aes(fill = group,
                    color = group,
                    shape = group),
                width = 0.1,
                height = 0.1,
                size = pointsize) +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray")+
    xlab("") +
    ylab("% Editing Frequency\n") +
    scale_y_continuous(position = "right", trans="log10",
                       limits = c(0.0007, 100),
                       #breaks = scales::trans_breaks("log10", function(x) 10^x),
                       breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                       labels = c(0, 0.01, 0.1, 1.0, 10, 100)) +
    coord_flip() +
    scale_x_discrete(limits = unique(summary_tb$name_seq),
                     labels = levels(summary_tb$name_seq)) +
    scale_fill_manual(name = fill[[1]],
                      values = fill[[2]]) +
    scale_color_manual(name = color[[1]],
                       values = color[[2]]) +
    scale_shape_manual(name = shape[[1]],
                values = shape[[2]]) +
    theme_classic() +
    theme(axis.text.y = element_text(family = "Courier" , size = guide_font_size))
  
  return(dotplot)
}


#Takes a summary table with all samples, all off-targets, and all editing outcomes and 
# generates an editing outcome heatmap.
make_compiled_OT_editing_heatmap <- function(summary_tb){
  
  heatmap <- ggplot(data = summary_tb,
                    aes(x = name_seq,
                        y = sample)) +
    geom_tile(aes(fill = log_freq),
              width = 1,
              height = 1,
              color = "white") +
    xlab("") +
    ylab("Sample") +
    scale_x_discrete(limits = levels(summary_tb$name_seq),
                     expand = c(0, 0)) +
    scale_y_discrete(position = "right",
                     limits = levels(summary_tb$sample),
                     expand = c(0, 0)) +
    scale_fill_gradientn(colours = c("grey95", "skyblue", "cornflowerblue"),
                         breaks=c(-3, -2, -1, 0, 1, 2),
                         na.value = "white") +
    coord_flip() +
    theme_classic() +
    theme(axis.text.y = element_text(family = "Courier" , size = 12),
          axis.text.x = element_text(angle = 330, hjust = 0, colour = "grey50"))
  
  return(heatmap)
}


#Takes a summary table with all samples, all off-targets, and all editing outcomes and 
# generates read coverage heatmap.
# Returns a heatmap that is meant to be placed in a composite graph with an editing dotplot.
make_compiled_OT_coverage_heatmap <- function(summary_tb, guide_font_size = 12,
                                              group_font_size = 14, tile_font_size = 3){
  
  max_log10 <- round(log10(max(summary_tb$total_reads, na.rm = TRUE)), 0)
  legend_breaks <- seq(0, 10^max_log10 , 10^(max_log10 - 1))
  
  heatmap <- ggplot(data = summary_tb,
                    aes(x = name_seq,
                        y = group)) +
    geom_tile(aes(fill = total_reads)) +
    geom_text(aes(label = gsub("NA", "",
                               comma(round(total_reads, 1), trim = FALSE))),
              #position = position_nudge(y = 0.17),
              size = tile_font_size,
              family = "Courier") +
    xlab("") +
    ylab("Total Reads per Sample\n") +
    scale_x_discrete(limits = rev(unique(summary_tb$name_seq)),
                     #labels = levels(summary_tb$name_seq),
                     #breaks = levels(summary_tb$name_seq),
                     breaks = c(),
                     expand = c(0, 0)) +
    scale_y_discrete(position = "right",
                     limits = levels(summary_tb$group),
                     labels = c(gsub("[ ,-]{1}", "\n", unique(summary_tb$group))),
                     expand = c(0, 0)) +
    scale_fill_gradientn(name = "",
                         colours = c("grey95", "skyblue", "cornflowerblue"),
                         breaks = legend_breaks,
                         labels = comma(legend_breaks),
                         na.value = "white")  +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "left",
          axis.text.y = element_text(family = "Courier" , size = guide_font_size),
          axis.text.x = element_text(colour = "grey50", size = group_font_size))
  
  return(heatmap)
}

### Generates a composite plot showing coverage as a heatmap on the left and % editing as 
### a dotplot on the right.
#
# ARGUMENTS:
# heatmap = ggplot heatmap object to be plotted
# dotplot = ggplot dotplot object to be plotted
# plot_height_inch = c(plot_in, note_in) vector containing height of plot and note components in inches
# title = string of title to be displayed on top of composite graph
# note_text = string of note to be displayed in italics at bottom of composite graph
# note_hjust = hjust value for note at bottom of plot
#
# Returns composite plot as grob and displays the grob object.
make_composite_grobPlot <- function(heatmap, dotplot, plot_height_inch, title, note_text, note_hjust){
  
  g_heatmap <- ggplotGrob(heatmap)
  g_dotplot <- ggplotGrob(dotplot)
  g <- cbind(g_heatmap, g_dotplot, size = "first")
  
  grid.newpage()
  grid.arrange(g,
               ncol = 1,
               heights = unit( plot_height_inch, c("in", "in")),
               top = textGrob( title, gp=gpar(fontsize=12,font=1)),
               bottom = textGrob(note_text,
                                 gp = gpar(fontsize=10,font=3), 
                                 hjust = note_hjust))
  
  composite_grob <- arrangeGrob(g,
                                ncol = 1,
                                heights = unit( plot_height_inch, c("in", "in")),
                                top = textGrob( title, gp=gpar(fontsize=12,font=1)),
                                bottom = textGrob(note_text,
                                                  gp = gpar(fontsize=10,font=3), 
                                                  hjust = note_hjust))
  
  return(composite_grob)
}


# Takes file_name (without extension), composite plot object, and width/height of saved plot (in inches)
# and saves the plot as .pdf and png
save_composite_plot <- function(file_name, composite_plot, plot_width_in, plot_height_in){
  
  ggsave(paste(file_name, ".pdf", sep = ""), plot = composite_plot, 
         width = plot_width_in, height = plot_height_in,
         units = "in")
  
  ggsave(paste(file_name, ".png", sep = ""), plot = composite_plot, 
         width = plot_width_in, height = plot_height_in,
         units = "in")
  
}

# move the on-target to the top of the table
put_on_target_first <- function(table){
  
  on_target_idx <- which(grepl("_ON", table$off_target, ignore.case = TRUE))
  
  if(length(on_target_idx) == 0){
    print("no on-target")
    return(table) #return original table
    
  }else{
    
    ordered_table <- rbind(table[on_target_idx,], table[-(on_target_idx),])
    
    return(ordered_table) #return ordered table
  }
}
