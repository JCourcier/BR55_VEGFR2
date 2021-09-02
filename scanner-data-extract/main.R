require(stringr)
require(data.table)
require(dplyr)
require(ggplot2)
require(plotly)
library(RColorBrewer)
library(ggpubr)
library(psych)



#source_dir <- "../datas/COHORTE BR55 PRINCIPALE/"


#files <- list.files(source_dir, recursive = T, full.names = T, pattern = "Graph.txt")


# Find the position of the first string matching with a given pattern in a vector of strings
find_the_line <- function(vector, pattern){
  stringr::str_detect(vector, pattern = pattern) %>% which()
}

# 
extract_the_val <- function(vec, row_pattern, match_pattern){
  tryCatch({
    row <- find_the_line(vec, row_pattern)
    matching <- stringr::str_match_all(vec[row], match_pattern)
    return(matching[[1]][,2])
  }, error = function(e) NULL)
  
}



parse_echo_graphs <- function(x, test_cohorte = F){
  datas <- data.table::fread(x, sep = "\n", col.names = "data")
  datas_vec <- datas[,data]
  # Extract mouse id and day
  if(test_cohorte){
    mouse_id_and_day_pattern <- "Patient ID: BR55 (\\d*).* (J\\d*)"
    matching <- stringr::str_match_all(datas_vec[1], mouse_id_and_day_pattern)
    mouse_id <- matching[[1]][, 2]
    day <- matching[[1]][, 3]
    # In some files, pattern to match differs... Matching pattern must be different..
    if(is.na(as.numeric(mouse_id))){
      mouse_id_and_day_pattern <- "Patient ID: BR55 .* (\\d*) (J\\d*)"
      matching <- stringr::str_match_all(datas_vec[1], mouse_id_and_day_pattern)
      mouse_id <- matching[[1]][, 2]
      day <- matching[[1]][, 3]
    }
  }else{
    mouse_id_and_day_pattern <- "Patient ID: BR55 (\\d*).*"
    matching <- stringr::str_match_all(datas_vec[1], mouse_id_and_day_pattern)
    mouse_id <- matching[[1]][, 2]
    day <- NA
  }
 
  
  # Peak intensity
  peak <- extract_the_val(datas_vec,  "Peak intensity", "Peak intensity.*:.?(\\d*\\.?\\d*)") %>% as.numeric()
  
  # Time to peak 
  time_to_peak <- extract_the_val(datas_vec, "Time to peak", "Time to peak.*:.?(\\d*\\.?\\d*)") %>% as.numeric()
  
  # Mean transit time
  mean_transit_time <- extract_the_val(datas_vec, "Mean transit", "Mean transit.*:.?(\\d*\\.?\\d*)") %>% as.numeric()
  
  # Slope
  slope <- extract_the_val(datas_vec, "Slope", "Slope.*:.?(\\d*\\.?\\d*)") %>% as.numeric()
  
  # Area
  area <- extract_the_val(datas_vec, "Area\\(", "Area.*:.?(\\d*\\.?\\d*)") %>% as.numeric()
  
  # Area wash in
  area_wash_in <- extract_the_val(datas_vec, "Area wash in", "Area wash in.*:.?(\\d*\\.?\\d*)") %>% as.numeric()
  
  # Area wash out
  area_wash_out <- extract_the_val(datas_vec, "Area wash out", "Area wash out.*:.?(\\d*\\.?\\d*)") %>% as.numeric()
  
  # Detect headers position
  table_start <- find_the_line(datas_vec, "\\[msec\\]")
  table_headers <- datas_vec[table_start]
  
  # Get colums names from detected string from headers
  splitted <- stringr::str_split(table_headers, "\t")
  splitted <- unlist(splitted)[which(unlist(lapply(splitted, nchar)) != 0)] # 
  splitted <- stringr::str_replace_all(splitted, "\\[", "")
  splitted <- stringr::str_replace_all(splitted, "\\]", "")
  
  # Read the table
  table_datas <- data.table::fread(text = datas_vec[(table_start):(length(datas_vec))],
                             dec = ".",
                             sep = " ",
                             encoding = "UTF-8",
                             header = T,
                             col.names = splitted) 

  # Add additionnal parsed data from files
  table_datas <- cbind(table_datas,
                       mouse_id = mouse_id,
                       day = day,
                       peak_intensity = peak,
                       time_to_peak = time_to_peak,
                       mean_transit_time = mean_transit_time,
                       slope = slope,
                       area = area, 
                       area_wash_in = area_wash_in,
                       area_wash_out = area_wash_out,
                       source_file = x)
}

adjusted_mean <- function(x, min_destructive_pulse = 614000, max_destructive_pulse = 617000){
  mean_1 <- x[msec < min_destructive_pulse, mean(`T`)]
  mean_2 <- x[msec > max_destructive_pulse, mean(`T`)]
  return(c(mean_1 - mean_2, mean_1, mean_2))
}





plot_variable_per_day_and_groups <- function(df,
                                             var_to_plot,
                                             yLab = NULL,
                                             saving_path = NULL,
                                             custom_ylims = NULL,
                                             fill_colors = c("#4373c5", "#ffc002", "#ed7d31")){
  # Build the main plot
  g <- ggplot(df, aes_string(y= var_to_plot, fill = "Groupe", x = "Groupe"))+
    geom_boxplot()+
    geom_jitter(color="black", size=0.5, alpha=0.5, width = 0.15) +
    facet_wrap(day ~ .)+
    scale_fill_manual(values=fill_colors)+
    theme_minimal()+
    ylab(yLab)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  if(!is.null(custom_ylims)){
    stopifnot(is.numeric(custom_ylims))
    stopifnot(is.vector(custom_ylims))
    stopifnot(length(custom_ylims) == 2)
    message("Custom y limits applied...")
    g <- g + ylim(custom_ylims)
  }
  
  if(!is.null(saving_path)){
    ggsave(saving_path, g)
    message(paste0("Plot save: ", saving_path))
  }
  
  g
}

variables_comparison_plots <- function(datas, y_var, yLab){
  my_comparisons <- list( c("Placebo", "axi 7.5"),c("axi 7.5", "axi 15"),  c("Placebo", "axi 15") )
  p <- ggboxplot(datas,
                 x = "Groupe",
                 y = y_var,
                 color = "Groupe",
                 fill = "#d8dde6",
                 palette = c("#4373c5", "#ffc002", "#ed7d31"),
                 add = "jitter",
                 shape = "Groupe",
                 facet.by = "day", 
                 short.panel.labs = F,
                 panel.labs.background = list(fill = "steelblue", color = "steelblue"),
                 ylab = yLab)
  
  p <- p + stat_compare_means(label.x.npc = 0.8,label.y.npc = 1) +
    stat_compare_means(comparisons = my_comparisons, method = "t.test") # Add pairwise comparisons
  return(p)
}


plot_graph_signals <- function(tardive, sel_day, saving_location, colors_scales = c("#ed7d31", "#ffc002", "#4373c5")){
  temp <- tardive %>% filter(day == sel_day) %>% dplyr::arrange(enc_groups, mouse_id)
  
  temp$Ordered <- factor(temp$mouse_id,      # Reordering group factor levels
                         levels = temp$mouse_id %>% unique())
  
  temp$Groupe <- factor(temp$Groupe, levels=c("Placebo", "axi 7.5", "axi 15"), labels=c("Placebo", "axi 7.5", "axi 15"))
  
  g <- ggplot(temp, aes(Time, `T`))+
    geom_line(aes(col = Groupe), lwd = 0.8)+
    facet_wrap(Ordered ~ . , ncol = 4, drop = T)+
    theme_minimal()+
    theme(strip.text.x = element_blank())+
    scale_color_manual(values = colors_scales)+
    # scale_colour_viridis_d(alpha = 1, begin = 0.2, end = 0.85)+
    xlab("Time (s)") + ylab("Signal intensity (a.u)") + ggtitle(sel_day)
  
  ggsave(file = saving_location, plot=g, width=10, height=8)
  
  return(g)
}






