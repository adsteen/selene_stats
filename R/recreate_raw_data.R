read_gene_freq_data <- function(filename) {
  d <- read_csv(filename)
  
  # Remove the "total" row if there is one
  total.row <- str_detect(str_to_lower(d$`gene number`), "total")
  d <- d %>% filter(!total.row)
  
  # Make the data set long, change the `gene number` column name to gene.count, and calulate frequency of gene count within a category
  d <- d %>%
    pivot_longer(cols = -1, names_to = "category", values_to = "count") %>%
    mutate(gene.count = as.numeric(`gene number`)) %>%
    dplyr::select(-`gene number`) %>%
    group_by(category) %>%
    mutate(freq = count / sum(count, na.rm = TRUE)) %>%
    arrange(category, count)
  d
}



# function to recreate raw data from the frequency data
recreate_raw <- function(d) {
  
  # Initialize a zero-row 
  recreated_raw_data <- data.frame("category"=character(0), 
                                   "count"=double(0), 
                                   "gene.count"=double(0)) # create empty data frame to bind to
  #browser()
  # This is a really inelegant way to recreate the raw count data from the summary counts 
  # "count" is the number of genomes that has a given number of relevant genes
  # gene.count is the number of relevant genes that genome has
  # So this loop runs through each row of the original data frame
  # and appends a number of rows equal to the "count" to the new df
  # Thus giving us a "recreated" data frame with a row for every genome observed
  
  # Note how I'm growing the df in the loop, very cool stuff
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

# fit gene counts to poisson distribution
fit_poisson <- function(df) {
  dist_obj <- MASS::fitdistr(df$gene.count, 
                 densfun=dpois, 
                 start=list(lambda = 1)) # Gotta figure out how to handle that warning
  dist_obj
}

# function to get lambda from poisson fit 
get_lambda <- function(fit_obj) {
  fit_obj$estimate
}

barplot_2 <- function(df) {
  max.num.levels <- max(df$gene.count, na.rm=TRUE)
  df <- df %>% 
    mutate(gene.count = factor(gene.count, levels=0:max.num.levels))
  
  p <- ggplot(df, aes(x=category, y=freq)) +
    geom_bar(aes(fill=gene.count), 
             position="dodge", stat="identity", colour="black", linewidth=0.1) + 
    scale_fill_brewer(name="gene number", type="seq") + 
    scale_y_continuous(name="frequency", labels = scales::percent) + 
    theme(axis.text.x = element_text(angle=-45, hjust=0))
  p
}

shuf_calc_f <- function(df) {
  nrow.df <- nrow(df)
  df <- df %>% 
    ungroup() %>%
    mutate(shuf.cat = sample(category, size = nrow.df, replace=FALSE))
  df
  #browser()
  #f <- tryCatch({
    m <- aov(gene.count ~ shuf.cat, data = df)
    f <- summary(m)[[1]][1,4]
  #}, 
  #  error=function(e) NA
  #)
  #m <- aov(gene.count ~ shuf.cat, data = df)
  
  f
}

# Summarise raw data for dot-standard deviation plots
summ_for_dotplot <- function(df) {
  df %>%
    ungroup() %>%
    group_by(category) %>%
    summarise(mean.gene.count = mean(gene.count),
              sd.gene.count = sd(gene.count, 0.25)) %>%
    mutate(cat.ordered = fct_reorder(category, mean.gene.count, mean, .desc=TRUE))
}

dotplot <- function(df) {
  ggplot(region_summ, aes(x=cat.ordered, y=mean.gene.count)) + 
    geom_pointrange(aes(ymin=mean.gene.count-sd.gene.count,
                        ymax=mean.gene.count+sd.gene.count)) 
}