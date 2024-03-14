

setwd("C:/Users/zachr/Desktop/6970_Asn_3")

library(readxl)

data <- read_xlsx("Immunologic profiles of patients with COVID-19.xlsx")

#Data Exploration

mean(data$AGE)
sd(data$AGE)
(data$BTLA)
