#anomaly detection
#period = the level at which the seasonality should be extracted, for.e.g 52 weeks of data
# you need to have at least 2 periods i.e. (52*2)+1 observations of data
#alpha = confidence interval default - 0.05 i.e. 95%
#max_anoms = number of observations you think will be anomalous in percentage
#if your data is longer than 12 weeks then longterm is suggested[ long tems splits your data into smaller windows]
#window_size = number of weeks by how you want to split your data
#threshold = med_max or p95 or p99 ; this allows you to choose anomalies higher than 95 percentile or 99 percentile

raw_data = read.csv(file.choose())
x = raw_data
anomaly_detection<-function(x, period = 52 , alpha = 0.05 , max_anoms = 0.49, longterm = TRUE, window_size = 26, threshold = "med_max"){
library('lubridate')
num_obs<-length(x[[2]])
#max_anoms = 1/num_obs
#max_anoms = 0.49
piecewise_median_period_weeks = 2
#period = 52
#alpha = 0.05

data = x

data_decomp <- stl(ts(data[[2L]], frequency = period),s.window = "periodic", robust = TRUE)
# Remove the seasonal component, and the median of the data to create the univariate remainder
data <- data.frame(timestamp = data[[1L]], count = (data[[2L]]-data_decomp$time.series[,"seasonal"]))
# Store the smoothed seasonal component, plus the trend component for use in determining the "expected values" option
data_decomp <- data.frame(timestamp=data[[1L]], count=(as.numeric(trunc(data_decomp$time.series[,"trend"]+data_decomp$time.series[,"seasonal"]))))

####################
# If longterm is enabled, break the data into subset data frames and store in all_data
if(longterm){
  # Pre-allocate list with size equal to the number of piecewise_median_period_weeks chunks in x + any left over chunk
  # handle edge cases for daily and single column data period lengths
  
  
  # Store last date in time series
  last_date <- data[[1]][num_obs]
  num_days_in_period = window_size #weeks
  num_obs_in_period = window_size
  all_data <- vector(mode="list", length=ceiling(length(data[[1]])/(num_obs_in_period))) 
  # Subset x into piecewise_median_period_weeks chunks
  for(j in seq(1,length(data[[1]]), by=num_obs_in_period)){
    start_date <- ymd(data[[1]][j]) 
    end_date <- min(ymd(start_date) + weeks(num_days_in_period), ymd(data[[1]][length(data[[1]])]))
    # if there is at least 14 days left, subset it, otherwise subset last_date - 14days
    if(as.integer(difftime(end_date, start_date, units = "weeks")) == as.integer(as.difftime(num_days_in_period, units="weeks"))){
      all_data[[ceiling(j/(num_obs_in_period))]] <- subset(data, as.Date(data[[1]]) >= as.Date(start_date) & as.Date(data[[1]]) < as.Date(end_date))
    }else{
      all_data[[ceiling(j/(num_obs_in_period))]] <- subset(data, as.Date(data[[1]]) > (as.Date(ymd(last_date))-weeks(num_days_in_period)) & as.Date(data[[1]]) <= as.Date(last_date))
    }
  }
}else{
  # If longterm is not enabled, then just overwrite all_data list with x as the only item
  all_data <- list(data)
}

# Create empty data frames to store all anoms and seasonal+trend component from decomposition
all_anoms <- data.frame(timestamp=numeric(0), count=numeric(0))
seasonal_plus_trend <- data.frame(timestamp=numeric(0), count=numeric(0))

# Detect anomalies on all data (either entire data in one-pass, or in 2 week blocks if longterm=TRUE)
for(k in 1:length(all_data)) {
  
  # detect_anoms actually performs the anomaly detection and returns the results in a list containing the anomalies
  # as well as the decomposed components of the time series for further analysis.
  data = all_data[[k]]
  # Maximum number of outliers that S-H-ESD can detect (e.g. 49% of data)
  num_obs = nrow(data)
  max_outliers <- trunc(num_obs*0.49)
  
  func_ma <- match.fun(median)
  func_sigma <- match.fun(mad)
  
  ## Define values and vectors.
  n <- length(data[[2L]])
  R_idx <- 1L:max_outliers
  
  num_anoms <- 0L
  verbose = FALSE
  # Compute test statistic until r=max_outliers values have been
  # removed from the sample.
  R_idxp<-vector("list", max_outliers) 
  for (i in 1L:max_outliers){
    if(verbose) print(paste(i,"/", max_outliers,"completed"))
    ares = abs(data[[2L]] - func_ma(data[[2L]]))
    ares <- ares/func_sigma(data[[2L]])
    R <- max(ares)
    temp_max_idx <- which(ares == R)[1L]
    R_idxp[i] <- as.character(as.Date(data[[1L]][temp_max_idx]))
    data <- data[-which(data[[1L]] == (R_idxp[i])), ]
    
    ## Compute critical value.
    p <- 1 - alpha/(2*(n-i+1))
    t <- qt(p,(n-i-1L))
    lam <- t*(n-i) / sqrt((n-i-1+t**2)*(n-i+1))
    
    if(R > lam)
      num_anoms <- i
  }
  
  if(num_anoms > 0) {
    R_idxp <- R_idxp[1L:num_anoms]
  } else {
    R_idxp = NULL
  }
  
  s_h_esd_timestamps = list(anoms = R_idxp,stl = data_decomp)
  
  # store decomposed components in local variable and overwrite s_h_esd_timestamps to contain only the anom timestamps
  data_decomp <- s_h_esd_timestamps$stl
  s_h_esd_timestamps <- s_h_esd_timestamps$anoms
  if(!is.null(s_h_esd_timestamps)){
    anoms <- subset(all_data[[k]], (all_data[[k]][[1]] %in% s_h_esd_timestamps))
  } else {
    anoms <- data.frame(timestamp=numeric(0), count=numeric(0))
  }
  
  # Filter the anomalies using one of the thresholding functions if applicable
  if(threshold != "None"){
    # Calculate daily max values
    periodic_maxs <- tapply(x[[2]],as.Date(x[[1]]),FUN=max)
    
    # Calculate the threshold set by the user
    if(threshold == 'med_max'){
      thresh <- median(periodic_maxs)
    }else if (threshold == 'p95'){
      thresh <- quantile(periodic_maxs, .95)
    }else if (threshold == 'p99'){
      thresh <- quantile(periodic_maxs, .99)
    }
    # Remove any anoms below the threshold
    anoms <- subset(anoms, anoms[[2]] >= thresh)
  }
  
  all_anoms <- rbind(all_anoms, anoms)
  seasonal_plus_trend <- rbind(seasonal_plus_trend, data_decomp)
}

# Cleanup potential duplicates
all_anoms <- all_anoms[!duplicated(all_anoms[[1]]), ]
seasonal_plus_trend <- seasonal_plus_trend[!duplicated(seasonal_plus_trend[[1]]), ]


# Calculate number of anomalies as a percentage
anom_pct <- (length(all_anoms[[2]]) / nrow(x)) * 100

# If there are no anoms, then let's exit
if(anom_pct == 0){
  print("No anomalies detected.")
  return (list("anoms"=NULL, "plot"=NULL))
}
plot = TRUE
title = "Anomaly Detection"
#plot
require(ggplot2)
library(scales)
library(ggplot2)
x$new_date = as.Date(x$timestamp)
all_anoms$new_date = as.Date(all_anoms$timestamp)
xplot = ggplot(x, aes(new_date, count)) + geom_line()+
  scale_x_date(labels= date_format( "%Y-%m-%d"), breaks = "1 week") + xlab("Date") + ylab("Bookings") 

all_anoms$actual_booking = x[which(x$new_date %in% all_anoms$new_date),2]
colnames(all_anoms) <- c("timestamp","count","new_date","actual_booking")
all_anoms <- subset(all_anoms, select = c(new_date,actual_booking))
xplot <- xplot + ggplot2::geom_point(data=all_anoms, ggplot2::aes(new_date, actual_booking), size = 3, shape = 1)
result <-list(anomalous_data=all_anoms,plot=xplot)
return (result)
}

output = anomaly_detection(x, period = 52 , alpha = 0.05 , max_anoms = 0.49, longterm = TRUE, window_size = 26, threshold = "med_max")
output$anomalous_data
output$plot

