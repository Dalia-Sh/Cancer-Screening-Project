# First run the all Functions_Modules_Code2

## Greedy Simulation

Greedy_Function <- function(Patient_list, N_Patient_list, T_max, k, epsilon){
  # Policy Module
  K_Patients <<- policy_mod_greedy(N_Patient_list, k, epsilon) 
  X_increment(K_Patients)
  # Imaging Module
  Cancer_Free <<- list() 
  Replace_Patient <<- list()
  imaging_mod(K_Patients, mysim_time=simulation_time)
  if(length(Cancer_Free)!=0){
    # AFP Reading and Update Module if Cancer-Free
    N_Patient_list <<- AFP_Update_mod2(N_Patient_list, Cancer_Free)
  }
  if(length(Replace_Patient)!=0){
    # New Patient Module if Early or Late (and E and L incremented)
    New_Patient_mod(Patient_list, Replace_Patient, N_Patient_list)
    Remaining_Patients <- Update_remaining_patients(N_Patient_list, Patient_list)
  }
  # Advance to time t+1
  simulation_time <<- simulation_time + 1
  if(simulation_time<T_max){
    # Patient Exit Module
    T_max_list <<- list()
    Exit_Set_func(N_Patient_list, T_max = T_max)
    if(length(T_max_list)!=0){
      New_Patient_mod(Patient_list, Replace_Patient =  T_max_list, N_Patient_list)
      Remaining_Patients <- Update_remaining_patients(N_Patient_list, Patient_list)
      L_increment(T_max_list)
    }
  }
  Canc <- list()
  Nocanc <- list()
  C_and_N(N_Patient_list)
}


Greedy_Sim <- function(T_max, N_N, k_prop, Patients_Data, epsilon){
  # Reading Patient Data
  Patients_Data <<- "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv" # Path to Data
  df <- read.csv(Patients_Data)
  Patient_list <<- setNames(split(df, seq(nrow(df))), rownames(df))
  k <<- k_prop*N_N 
  # Creating AFP_values
  AFP_values <<- data.frame(
    RR = c(rnorm(82*10, 5, 11), rnorm(246*10, 0.11, 2.1)),
    SS = c(rnorm(82*10, 51, 86), rnorm(246*10, 9, 19)),
    C_N = c(rep("C", 82*10), rep("N",246*10)))
  # Initialization
  X <<- 0
  E <<- 0
  L <<- 0 
  simulation_time <<- 0 # Simulation start time
  # Saving Patient ids in a vector (needed for the: update_remaining_patients function)
  All_Pids <<- c()
  for (i in 1:length(Patient_list)){
    All_Pids[i] <<- Patient_list[[i]][1,"Pid"]
  }
  # Initial Panel Module
  N_Patient_list <<- initial_panel(Patient_list, N=N_N)
  Remaining_Patients <<- setdiff(Patient_list, N_Patient_list)
  # Creation of C and N
  Canc <<- list()
  Nocanc <<- list()
  C_and_N(N_Patient_list)
  while(simulation_time < T_max){
    Greedy_Function(Patient_list, N_Patient_list, T_max, k, epsilon)
  }
}

Greedy_many_iteration <- function(iteration, T_max, N_N, k_prop, Patients_Data, epsilon){
  E_vec <<- c()
  L_vec <<- c()
  X_vec <<- c()
  for(i in 1:iteration){
    Greedy_Sim(T_max, N_N, k_prop, Patients_Data, epsilon)
    E_vec[i] <<- E
    L_vec[i] <<- L
    X_vec[i] <<- X
  }
}


## RESULTS (N_N = 100, iteration = 10, T_max =15)

Greedy_many_iteration(
  iteration = 10,
  T_max = 15, 
  N_N = 100, 
  k_prop = 0.1, 
  epsilon = 0.025,
  Patients_Data = "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv")

# 10% Prop: mean(E_vec/(E_vec+L_vec)) = 0.3767379
# 10% Prop: mean(X_vec/((k_prop*N_N)*T_max)) = 0.338

Greedy_many_iteration(
  iteration = 10,
  T_max = 15, 
  N_N = 100, 
  k_prop = 0.2, 
  epsilon = 0.05,
  Patients_Data = "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv")

# 20% Prop: mean(E_vec/(E_vec+L_vec)) = 0.5284845
# 20% Prop: mean(X_vec/((k_prop*N_N)*T_max)) = 0.2663333

Greedy_many_iteration(
  iteration = 10,
  T_max = 15, 
  N_N = 100, 
  k_prop = 0.3, 
  epsilon = 0.10,
  Patients_Data = "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv")

# 30% Prop: mean(E_vec/(E_vec+L_vec)) = 0.5771196
# 30% Prop: mean(X_vec/((k_prop*N_N)*T_max)) = 0.2473333

Greedy_many_iteration(
  iteration = 10,
  T_max = 15, 
  N_N = 100, 
  k_prop = 0.4, 
  epsilon = 0.10,
  Patients_Data = "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv")

# 40% Prop: mean(E_vec/(E_vec+L_vec)) = 0.5828866
# 40% Prop: mean(X_vec/((k_prop*N_N)*T_max)) = 0.2138333

Greedy_many_iteration(
  iteration = 10,
  T_max = 15, 
  N_N = 100, 
  k_prop = 0.5, 
  epsilon = 0.5, # This parameter was changed (not the optimal one as per paper) and performed better here0.
  Patients_Data = "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv")

# 50% Prop: mean(E_vec/(E_vec+L_vec)) = 0.6601321
# 50% Prop: mean(X_vec/((k_prop*N_N)*T_max)) = 0.1538667

## RESULT GRAPH COMPARISON

library(ggplot2)
library(tidyr)

test_data <-
  data.frame(
    Paper_Greedy = c(0.19, 0.36, 0.47, 0.53, 0.62),
    My_Greedy =  c(0.37, 0.52, 0.57, 0.58, 0.66),
    Screening_Capacity = c(10, 20, 30, 40, 50)
  )

test_data %>%
  gather(key,value, Paper_Greedy, My_Greedy) %>%
  ggplot(aes(x=Screening_Capacity, y=value, colour=key)) +
  geom_line() +
  ggtitle("Greedy Policy Early Detection")

test_data2 <-
  data.frame(
    Paper_Greedy = c(0.165, 0.15, 0.14, 0.13, 0.11),
    My_Greedy =  c(0.33, 0.26, 0.24, 0.22, 0.15),
    Screening_Capacity = c(10, 20, 30, 40, 50)
  )

test_data2 %>%
  gather(key,value, Paper_Greedy, My_Greedy) %>%
  ggplot(aes(x=Screening_Capacity, y=value, colour=key)) +
  geom_line() +
  ggtitle("Greedy Policy Screening Rate")

