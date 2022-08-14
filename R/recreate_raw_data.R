# function to recreate raw data from the frequency data
recreate_raw <- function(d) {
  
  recreated_raw_data <- data.frame("category"=character(0), 
                                   "count"=double(0), 
                                   "gene.count"=double(0)) # create empty data frame to bind to
  #fuck me, do this in the worst way possible
  for(i in 1:nrow(d)) {
    this_row <- d[i, ]
    for(j in 1:this_row$count) {
      if(this_row$count > 0) {
        recreated_raw_data <- rbind(recreated_raw_data, this_row)
      }
    }
  }
  recreated_raw_data 
}

# Create plots of raw data
raw_plot <- function(df) {
  p <- ggplot(df, aes(x=category, y=gene.count)) + 
    geom_boxplot() + 
    geom_point(position=position_jitter(height = 0.3), alpha = 0.5) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle=-45, hjust=0))
  p
}

# fit poisson data 
fit_poisson <- function(df) {
  dist_obj <- MASS::fitdistr(df$gene.count, 
                 densfun=dpois, 
                 start=list(lambda = 1)) # Gotta figure out how to handle that warning
  dist_obj
}

get_lambda <- function(fit_obj) {
  fit_obj$estimate
}
