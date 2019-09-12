# Simulation Function

#### SIM_DATA2
Run_Simulation2 <- function(Patient_list, N_Patient_list, T_max, k){
  # Policy Module
  K_Patients <<- policy_mod(N_Patient_list, k) 
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
  C <- list()
  N <- list()
  C_and_N(N_Patient_list)
}


Run_Simulation_greedy2 <- function(Patient_list, N_Patient_list, T_max, k, epsilon){
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
  C <- list()
  N <- list()
  C_and_N(N_Patient_list)
}

Run_Simulation_interv2<- function(Patient_list, N_Patient_list, T_max, k, z){
  # Policy Module
  K_Patients <<- policy_mod_interv(N_Patient_list, k, z) 
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
  C <- list()
  N <- list()
  C_and_N(N_Patient_list)
}


#### DRAFT

Run_Simulation <- function(Patient_list, N_Patient_list, T_max, k){
  # Policy Module
  K_Patients <<- policy_mod(N_Patient_list, k) 
  X_increment(K_Patients)
  # Imaging Module
  Cancer_Free <<- list() 
  Replace_Patient <<- list()
  imaging_mod(K_Patients, mysim_time=simulation_time)
  if(length(Cancer_Free)!=0){
    # AFP Reading and Update Module if Cancer-Free
    N_Patient_list <<- AFP_Update_mod(N_Patient_list, Cancer_Free)
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
  C <- list()
  N <- list()
  C_and_N(N_Patient_list)
}

Run_Simulation_greedy <- function(Patient_list, N_Patient_list, T_max, k, epsilon){
  # Policy Module
  K_Patients <<- policy_mod_greedy(N_Patient_list, k, epsilon) 
  X_increment(K_Patients)
  # Imaging Module
  Cancer_Free <<- list() 
  Replace_Patient <<- list()
  imaging_mod(K_Patients, mysim_time=simulation_time)
  if(length(Cancer_Free)!=0){
    # AFP Reading and Update Module if Cancer-Free
    N_Patient_list <<- AFP_Update_mod(N_Patient_list, Cancer_Free)
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
  C <- list()
  N <- list()
  C_and_N(N_Patient_list)
}
