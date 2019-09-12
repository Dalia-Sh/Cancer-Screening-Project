## Functions and Modules Codes 

### All functions

initial_panel <- function(Patient_list, N){
  N_idx <- sample(1:length(Patient_list), N, replace=F)
  N_Patient_list <- Patient_list[N_idx]
  return(N_Patient_list)
}

C_and_N <- function(N_Patient_list){
  for (i in 1:length(N_Patient_list)){
    if(N_Patient_list[[i]][1,"C_N"]=="C"){
      C[[i]] <- N_Patient_list[[i]] 
    }
    else if (N_Patient_list[[i]][1,"C_N"]=="N"){
      N[[i]] <- N_Patient_list[[i]]
    }
  }
  N <<- N[!sapply(N, is.null)]
  C <<- C[!sapply(C, is.null)]
}

policy_modHCC <- function(N_Patient_list, k){
  PHCC <- c()
  for (i in 1:length(N_Patient_list)){
    c1Bi <- N_Patient_list[[i]][(nrow(N_Patient_list[[i]])), "c1Bi"]
    c2SD <- 1.02*N_Patient_list[[i]][(nrow(N_Patient_list[[i]])), "SS"]
    c3RR <- 1.14*N_Patient_list[[i]][(nrow(N_Patient_list[[i]])), "RR"]
    PHCC[i] <- (1+exp(-c1Bi-c2SD-c3RR))^(-1)
  }
  names(PHCC) <- as.character(1:length(N_Patient_list)) # name is the patient order in the list of patients
  PHCC_ordered <- sort(PHCC, decreasing = T)
  PHCC_ordered <- PHCC_ordered[1:k]
  return(PHCC_ordered)
}

policy_mod <- function(N_Patient_list, k){
  PHCC <- c()
  for (i in 1:length(N_Patient_list)){
    c1Bi <- N_Patient_list[[i]][(nrow(N_Patient_list[[i]])), "c1Bi"]
    c2SD <- 1.02*N_Patient_list[[i]][(nrow(N_Patient_list[[i]])), "SS"]
    c3RR <- 1.14*N_Patient_list[[i]][(nrow(N_Patient_list[[i]])), "RR"]
    PHCC[i] <- (1+exp(-c1Bi-c2SD-c3RR))^(-1)
  }
  names(PHCC) <- as.character(1:length(N_Patient_list)) # name is the patient order in the list of patients
  PHCC_ordered <- sort(PHCC, decreasing = T)
  PHCC_ordered <- PHCC_ordered[1:k]
  K_Patients <- N_Patient_list[as.numeric(names(PHCC_ordered))]
  return(K_Patients)
}

policy_mod_greedy <- function(N_Patient_list, k, epsilon){
  PHCC <- c()
  for (i in 1:length(N_Patient_list)){
    c1Bi <- N_Patient_list[[i]][(nrow(N_Patient_list[[i]])), "c1Bi"]
    c2SD <- 1.02*N_Patient_list[[i]][(nrow(N_Patient_list[[i]])), "SS"]
    c3RR <- 1.14*N_Patient_list[[i]][(nrow(N_Patient_list[[i]])), "RR"]
    PHCC[i] <- (1+exp(-c1Bi-c2SD-c3RR))^(-1)
  }
  names(PHCC) <- as.character(1:length(N_Patient_list)) # name is the patient order in the list of patients
  top = (1-epsilon)*k
  bottom = epsilon*k
  PHCC_ordered <- sort(PHCC, decreasing = T)
  PHCC_ordered_top <- PHCC_ordered[1:top]
  PHCC_rest <- PHCC_ordered[(top+1):length(PHCC_ordered)]
  rdm_idx <- sample(length(PHCC_rest), k)
  PHCC_ordered_random <- PHCC_ordered[rdm_idx]
  PHCC_K <- c(PHCC_ordered_top, PHCC_ordered_random)
  K_Patients <- N_Patient_list[as.numeric(names(PHCC_K))]
  return(K_Patients)
}

X_increment <- function(K_Patients){
  for(i in 1:length(K_Patients)){
    if(K_Patients[[i]][1,"C_N"]=="C"){
      X <<- X + 1
    }
  }
}

estimated_tumor <- function(time, t_sbar, sbar){
  return((2^((time-t_sbar)/0.5))*sbar)
}

imaging_mod <- function(K_Patients, mysim_time){
  for (i in 1:length(K_Patients)){
    if(mysim_time==0){
      Cancer_Free[[i]] <<- K_Patients[[i]]
    }
    else{
      sim_time = mysim_time
      t_sbar = K_Patients[[i]][nrow(K_Patients[[i]]), "t_sbar"]
      sbar = K_Patients[[i]][nrow(K_Patients[[i]]), "sbar"]
      s <- estimated_tumor(time=sim_time, t_sbar, sbar) 
      if (sim_time>=t_sbar & (s!=0 & s<=5)){
        E <<- E + 1
        Replace_Patient[[i]] <<- K_Patients[[i]]
      }
      else if (sim_time>=t_sbar & (s!=0 & s>5)){
        L <<- L+1
        Replace_Patient[[i]] <<- K_Patients[[i]]
      }
      else {
        Cancer_Free[[i]] <<- K_Patients[[i]]
      }
    }
  }
  Cancer_Free <<- deleteNULLs(Cancer_Free)
  Replace_Patient <<- deleteNULLs(Replace_Patient)
}

AFP_Update_mod <- function(N_Patient_list, Cancer_Free){
  xrow <<- data.frame(Pid = 0, AFP = 0, SS = 0, RR=0, c1Bi=0, sbar=0, t_sbar=0, t = 0, C_N="X")
  for(i in 1:length(N_Patient_list)){
    for(j in 1:length(Cancer_Free)){
      if (N_Patient_list[[i]][1,"Pid"]==Cancer_Free[[j]][1,"Pid"]){
        N_Patient_list[[i]] <- rbind(N_Patient_list[[i]], xrow)
        N_Patient_list[[i]][nrow(N_Patient_list[[i]]),] <- N_Patient_list[[i]][nrow(N_Patient_list[[i]])-1,]
        N_Patient_list[[i]][nrow(N_Patient_list[[i]]),"t"] <- N_Patient_list[[i]][nrow(N_Patient_list[[i]]), "t"]+1
        N_Patient_list[[i]][nrow(N_Patient_list[[i]]),"AFP"] <- AFP_values[which(AFP_values$Pid==N_Patient_list[[i]][1,"Pid"]&AFP_values$t==N_Patient_list[[i]][nrow(N_Patient_list[[i]]),"t"]),"AFP"]
        N_Patient_list[[i]][nrow(N_Patient_list[[i]]),"SS"] <- sd(N_Patient_list[[i]][,"AFP"])
        N_Patient_list[[i]][nrow(N_Patient_list[[i]]),"RR"] <- N_Patient_list[[i]][nrow(N_Patient_list[[i]]),"AFP"] - N_Patient_list[[i]][nrow(N_Patient_list[[i]])-1,"AFP"]
      }
    }
  }
  return(N_Patient_list)
}

AFP_Update_mod2 <- function(N_Patient_list, Cancer_Free){
  xrow <<- data.frame(Pid = 0, c1Bi=0, C_N="X", SS = 0, RR=0, sbar=0, t_sbar=0, t = 0)
  for(i in 1:length(N_Patient_list)){
    for(j in 1:length(Cancer_Free)){
      if (N_Patient_list[[i]][1,"Pid"]==Cancer_Free[[j]][1,"Pid"]){
        N_Patient_list[[i]] <- rbind(N_Patient_list[[i]], xrow)
        N_Patient_list[[i]][nrow(N_Patient_list[[i]]),] <- N_Patient_list[[i]][nrow(N_Patient_list[[i]])-1,]
        N_Patient_list[[i]][nrow(N_Patient_list[[i]]),"t"] <- N_Patient_list[[i]][nrow(N_Patient_list[[i]]), "t"]+1
        if(N_Patient_list[[i]][1,"C_N"]=="C"){
          N_Patient_list[[i]][nrow(N_Patient_list[[i]]),"SS"] <- sample(AFP_values$SS[AFP_values$C_N=="C"], 1)
          N_Patient_list[[i]][nrow(N_Patient_list[[i]]),"RR"] <- sample(AFP_values$RR[AFP_values$C_N=="C"], 1)
        }
        else{
          N_Patient_list[[i]][nrow(N_Patient_list[[i]]),"SS"] <- sample(AFP_values$SS[AFP_values$C_N=="N"], 1)
          N_Patient_list[[i]][nrow(N_Patient_list[[i]]),"RR"] <- sample(AFP_values$RR[AFP_values$C_N=="N"], 1)
        }
      }
    }
  }
  return(N_Patient_list)
}

New_Patient_mod <- function(Remaining_Patients, Replace_Patient, N_Patient_list){
  smple <-sample(1:length(Remaining_Patients), length(Replace_Patient), replace = F)
  new_patient <- Remaining_Patients[smple] # I select a random new patient from the remaining_patient_list
  idxx <- c()
  for(i in 1:length(N_Patient_list)){
    for(j in 1:length(Replace_Patient)){
      if (N_Patient_list[[i]][1,"Pid"]==Replace_Patient[[j]][1,"Pid"]){
        idxx[i] <- i # This is the index of the Replace_Patient in the N_Patient_list
      }
    }
  }
  idxx<-idxx[!(is.na(idxx))] # Removing any NAs
  N_Patient_list[idxx] <<- new_patient # Replacing the 'Replace_Patient' in the N_Patient_list with new_patient
}

Update_remaining_patients <- function(N_Patient_list, Patient_list){
  New_Pids <- c()
  for (i in 1:length(N_Patient_list)){
    New_Pids[i] <- N_Patient_list[[i]][1,"Pid"]
  }
  Rem_Pids <- setdiff(All_Pids, New_Pids) # These are the Patient ids that are remaining after replacement (not selected in the N patients)
  Rem_idx <- c()
  for(i in 1:length(Patient_list)){
    if(Patient_list[[i]]$Pid %in% Rem_Pids){
      Rem_idx[[i]] <- i
    }
    else{Rem_idx[[i]] <- F}
  }
  Remaining_Patients <- c()
  return(Remaining_Patients <- Patient_list[Rem_idx])
}

Exit_Set_func <- function(N_Patient_list, T_max){
  idxxx <- c()
  for(i in 1:length(N_Patient_list)){
    if (N_Patient_list[[i]][nrow(N_Patient_list[[i]]),"t"]>=T_max){
      idxxx[i] <- i
    }
  }
  idxxx<<-idxxx[!(is.na(idxxx))]
  T_max_list <<-N_Patient_list[idxxx]
  T_max_list <<- deleteNULLs(T_max_list)
}

L_increment <- function(T_max_list){
  D_C <- list()
  for(i in 1:length(T_max_list)){
    if(T_max_list[[i]][nrow(T_max_list[[i]]),"C_N"]=="C"){
      D_C[[i]] <- T_max_list[[i]]
    }
  }
  D_C <<- deleteNULLs(D_C)
  L <<- L + length(D_C)
}

policy_mod_interv <- function(N_Patient_list, k, z){
  PHCC <- c()
  for (i in 1:length(N_Patient_list)){
    c1Bi <- N_Patient_list[[i]][(nrow(N_Patient_list[[i]])), "c1Bi"]
    c2 <- 1.02
    SD <- N_Patient_list[[i]][(nrow(N_Patient_list[[i]])), "SS"]
    c3 <- 1.14
    RR <- N_Patient_list[[i]][(nrow(N_Patient_list[[i]])), "RR"]
    w <- var(N_Patient_list[[i]]$SS)
    if(is.na(w)==T){
      w <- 1
    }
    v <- var(N_Patient_list[[i]]$RR)
    if(is.na(v)==T){
      v <- 1
    }
    PHCC[i] <- (1+exp(-c1Bi-(c2*(SD+z*sqrt(v)))-(c3*(RR+z*sqrt(w)))))^(-1)
  }
  names(PHCC) <- as.character(1:length(N_Patient_list)) # name is the patient order in the list of patients
  PHCC_ordered <- sort(PHCC, decreasing = T)
  PHCC_ordered <- PHCC_ordered[1:k]
  K_Patients <- N_Patient_list[as.numeric(names(PHCC_ordered))]
  return(K_Patients)
}



