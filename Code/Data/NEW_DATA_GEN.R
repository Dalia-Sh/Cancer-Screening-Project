# New patients values

SS
RR
Platelets
Smoking
Age

# Age
Age_HCC <- trunc(rnorm(82,53,7))
Age_Non_HCC <- trunc(rnorm(246,50,7))
Age = c(Age_HCC, Age_Non_HCC)

# Platelets
Platelets_HCC <- trunc(rnorm(82, 126, 51))
Platelets_Non_HCC <- trunc(rnorm(246, 169, 65))
Platelets = c(Platelets_HCC, Platelets_Non_HCC)

# Smoked
Smoked_HCC <- c(rep(1, trunc(41/100*82)), rep(0, (82-trunc(41/100*82))))
Smoked_Non_HCC <- c(rep(1, trunc(24/100*246)), rep(0, (246-trunc(24/100*246))))
Smoked = c(Smoked_HCC,Smoked_Non_HCC)

# Pid
Pid <- 1:(82+246)

# HCC
C_N <- c(rep("C", 82), rep("N", 246))

# tumor sbar
sbar_HCC <- trunc(runif(82,1,5))
sbar_Non_HCC <- rep(0,246)
sbar = c(sbar_HCC, sbar_Non_HCC)

# t_sbar
tbar_HCC <- trunc(runif(82,1,10))
tbar_Non_HCC <- rep(0,246)
tbar = c(tbar_HCC, tbar_Non_HCC)

# Static dataframe
static_df <- data.frame(Pid=Pid,Age=Age, Platelets=Platelets, Smoked=Smoked, C_N = C_N)

# Vector of static c1Bi ouptut
c1Bi_func <- function(Age, Platelets, Smoked){
  return(1/(1+(0.99*Platelets+3.06*Smoked+1.05*Age)))
} 
  
c1Bi <- c1Bi_func(static_df$Age, static_df$Platelets, static_df$Smoked)
static_df$c1Bi <- c1Bi

# Patient dataframe (STATIC at time 0)
patient_df <- static_df[,c(1,6,5)]
patient_df$SS <- rep(0, (82+246))
patient_df$RR <- rep(0, (82+246))
patient_df$sbar <- sbar
patient_df$t_sbar <- tbar
patient_df$t <- rep(0, (82+246))
#write.csv(patient_df, "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv")

# SS
SS_HCC <- rnorm(82*10, 51, 86)
SS_Non_HCC <- rnorm(246*10, 9, 19)
SS = c(SS_HCC, SS_Non_HCC)
C_N <- c(rep("C", 82*10), rep("N",246*10))
t = rep(1:10, (82+246))
SS_df <- data.frame(SS, t, C_N)

# RR
RR_HCC <- rnorm(82*10, 5, 11)
RR_Non_HCC <- rnorm(246*10, 0.11, 2.1)
RR = c(RR_HCC, RR_Non_HCC)
C_N <- c(rep("C", 82*10), rep("N",246*10))
RR_df <- data.frame(RR, C_N)

# Dynamic Data
dynamic_data <- cbind(RR=RR_df[,-2],SS_df)
#write.csv(dynamic_data, "/Users/daliashanshal/Desktop/CancerScreening/Code/dynamic_data.csv")