# AnomalyDetection
This code uses generalized ESD technique to detect anomalies
This code also takes into account any seasonality in the data

All you need to pass into the code is a raw data with 2 columns
<column1 - date> and <column2 - metric value>

Define your needs:
confidence interval --> alpha [ 0.05 (95%) , 0.01(99%) ]

Constraints : 
We only consider weekly data for now 
You need to at least have more then (2*56)+1 observations of data for this code to work 

if your data has more than 150 observations and you want to very closely remove outliers then please set the longterm flag as TRUE and also decide on what the window size will be.

