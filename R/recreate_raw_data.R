# function to recreate raw data from the frequency data
recreate_raw <- function(d) {
  
  recreated_raw_data <- data.frame("category"=character(0), 
                                   "frequency"=double(0), 
                                   "gene.count"=double(0)) # create empty data frame to bind to
  #fuck me, do this in the worst way possible
  for(i in 1:nrow(d)) {
    this_row <- d[i, ]
    for(j in 1:this_row$frequency) {
      if(this_row$frequency > 0) {
        recreated_raw_data <- rbind(recreated_raw_data, this_row)
      }
    }
  }
  recreated_raw_data
}
