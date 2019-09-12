# First run the all Functions_Modules_Code2

## Boltzmann Simulation Function

Boltzmann_Function <- function(Patient_list, N_Patient_list, T_max, k, boltz_t){
  # Policy Module
  K_Patients <<- policy_mod_boltzmann(N_Patient_list, k, boltz_t)
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

Boltzmann_Sim <- function(T_max, N_N, k_prop, Patients_Data, boltz_t){
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
    Boltzmann_Function(Patient_list, N_Patient_list, T_max, k, boltz_t)
  }
}

Boltzmann_many_iteration <- function(iteration, T_max, N_N, k_prop, Patients_Data, boltz_t){
  E_vec <<- c()
  L_vec <<- c()
  X_vec <<- c()
  for(i in 1:iteration){
    Boltzmann_Sim(T_max, N_N, k_prop, Patients_Data, boltz_t)
    E_vec[i] <<- E
    L_vec[i] <<- L
    X_vec[i] <<- X
  }
}

Boltzmann_many_iteration(
  iteration = 10, # number of simulation iteration
  T_max = 10, # Time each simulation iteration will run
  N_N = 100, # N (Populaiton)
  k_prop = 0.5, # k Proportion of N
  Patients_Data = "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv", # Static Patient Value csv path
  boltz_t = 0.25 # boltz_t parameter
  )

## RESULTS (N_N = 100, iteration = 10, T_max =15)

Boltzmann_many_iteration(
  iteration = 10, # number of simulation iteration
  T_max = 15, # Time each simulation iteration will run
  N_N = 100, # N (Populaiton)
  k_prop = 0.1, # k Proportion of N
  Patients_Data = "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv", # Static Patient Value csv path
  boltz_t = 0.2 # boltz_t parameter
)

# 10% Prop: mean(E_vec/(E_vec+L_vec)) = 0.7648148
# 10% Prop: mean(X_vec/((k_prop*N_N)*T_max)) = 0.07066667

Boltzmann_many_iteration(
  iteration = 10, # number of simulation iteration
  T_max = 15, # Time each simulation iteration will run
  N_N = 100, # N (Populaiton)
  k_prop = 0.2, # k Proportion of N
  Patients_Data = "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv", # Static Patient Value csv path
  boltz_t = 0.25 # boltz_t parameter
)

# 20% Prop: mean(E_vec/(E_vec+L_vec)) = 0.7037879
# 20% Prop: mean(X_vec/((k_prop*N_N)*T_max)) = 0.06166667

Boltzmann_many_iteration(
  iteration = 10, # number of simulation iteration
  T_max = 15, # Time each simulation iteration will run
  N_N = 100, # N (Populaiton)
  k_prop = 0.3, # k Proportion of N
  Patients_Data = "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv", # Static Patient Value csv path
  boltz_t = 0.3 # boltz_t parameter
)

# 30% Prop: mean(E_vec/(E_vec+L_vec)) = 0.7471429
# 30% Prop: mean(X_vec/((k_prop*N_N)*T_max)) = 0.08022222

Boltzmann_many_iteration(
  iteration = 10, # number of simulation iteration
  T_max = 15, # Time each simulation iteration will run
  N_N = 100, # N (Populaiton)
  k_prop = 0.4, # k Proportion of N
  Patients_Data = "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv", # Static Patient Value csv path
  boltz_t = 0.325 # boltz_t parameter
)

# 40% Prop: mean(E_vec/(E_vec+L_vec)) = 0.7473112
# 40% Prop: mean(X_vec/((k_prop*N_N)*T_max)) = 0.0855

Boltzmann_many_iteration(
  iteration = 10, # number of simulation iteration
  T_max = 15, # Time each simulation iteration will run
  N_N = 100, # N (Populaiton)
  k_prop = 0.5, # k Proportion of N
  Patients_Data = "/Users/daliashanshal/Desktop/CancerScreening/Code/patient_static_data.csv", # Static Patient Value csv path
  boltz_t = 0.4 # boltz_t parameter
)

# 50% Prop: mean(E_vec/(E_vec+L_vec)) = 0.7979181
# 50% Prop: mean(X_vec/((k_prop*N_N)*T_max)) = 0.08666667


## RESULT GRAPH

library(ggplot2)
library(tidyr)

test_data <-
  data.frame(
    Paper_Boltzmann = c(0.25, 0.4, 0.5, 0.55, 0.63),
    My_Boltzmann =  c(0.76, 0.70, 0.74, 0.74, 0.8),
    Screening_Capacity = c(10, 20, 30, 40, 50)
  )

test_data %>%
  gather(key,value, Paper_Boltzmann, My_Boltzmann) %>%
  ggplot(aes(x=Screening_Capacity, y=value, colour=key)) +
  geom_line() +
  ggtitle("Boltzmann Policy Early Detection")


test_data2 <-
  data.frame(
    Paper_Boltzmann = c(0.15, 0.14, 0.13, 0.12, 0.11),
    My_Boltzmann =  c(0.07, 0.062, 0.08, 0.08, 0.08),
    Screening_Capacity = c(10, 20, 30, 40, 50)
  )

test_data2 %>%
  gather(key,value, Paper_Boltzmann, My_Boltzmann) %>%
  ggplot(aes(x=Screening_Capacity, y=value, colour=key)) +
  geom_line() +
  ggtitle("Boltzmann Policy Screening Rate")




