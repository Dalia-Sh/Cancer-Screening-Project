# Evaluation of Learning-Based Screening Policies for Early Detection of Hepatocellular Carcinoma (HCC)

## About This Project
This project is a replication of the following paper: Lee, Elliot, et al. "Applying reinforcement learning techniques to detect hepatocellular carcinoma under limited screening capacity." Health care management science 18.3 (2015): 363-375.
I replicate the simulation flow presented in the paper using a toy dataset that I created. All code is written in R, and a presentation report is available under the presentation folder. 

## Introduction
Cancer, when detected early, increase the chance of survivals significantly. 50% of patients diagnosed with early-stage cancer survive beyond 5 years, as opposed to 10% of patients diagnosed with late-stage cancer. 
Hepatocellular Carcinoma (HCC) is commonly known as liver cancer, its risk has been increasing over the past decades. In Canada, HCC frequency tripled in men and doubles in women since 1970.

## Current Challenges

#### Screening Process
The standard screening process consists of: (1) Ultrasound image of patient liver; (2) Blood-test measure for alphafeto-protein level (AFP). This standard (type of) screening does not consider the changes in blood tests done at each screening visit. However, it has been shown that AFP past patterns significantly improve estimation of patient's risk of developing HCC. 

#### Screening Cost
Screening is costly, and the number of patients at risk is often greater than the number of screening resources available.

## Proposed Solution
To overcome the above challenges, Lee et al. [1] proposes a framework to test different strategies for selecting individuals to be screened for HCC in an optimal manner, by assigning a probability of cancer risk based on each individual's demographic information and blood-tests, while taking into account the AFP blood-test fluctuations.

### Simulation input
- n: List of n patients
- k: number of k patients to subset (i.e. from those n patients, how many individuals do we want to select for screening?)
- T: iteration time

### Simulation Output
- E: number of early stage cancers detected
- L: number of late stage cancer detected
- X: number of screenings spent on patients who would eventually develop HCC
 
 ### Other Variables
 - C: set of patients who will develop cancer
 - N: set of patients who will not develop cancer
 - D: Departing Patients (Patients are eliminated and replaced)
 
 ### Performance Measures
 - E/(E + L): Proportion of cancers detected in early stage
 - X/(k * T): Proportion of resources spent on patients who will eventually develop cancer
 
 ### Simulation Flow Graph (Refer to Lee, Elliot, et al. paper [1])
 
 ![Simulation Flow](https://github.com/Dalia-Sh/Cancer-Screening/blob/master/Other/simulationflow.png)
 
 ## Related References
[1] Lee, Elliot, et al. "Applying reinforcement learning techniques to detect
hepatocellular carcinoma under limited screening capacity." Health care management
science 18.3 (2015): 363-375.

[2] Lee, Elliot, et al. "Improving screening for hepatocellular carcinoma by
incorporating data on levels of α-fetoprotein, over time." Clinical Gastroenterology
and Hepatology 11.4 (2013): 437-440.

[3] O'Connor, S., et al. "Hepatocellular carcinoma-United States, 2001-
2006." Morbidity and Mortality Weekly Report 59.17 (2010): 517-520.

[4] Okada, S., et al. "Follow-up examination schedule of postoperative HCC patients
based on tumor volume doubling time." Hepato-gastroenterology 40.4 (1993): 311-
315.

## Some Thoughts
It appears that in Ontario, cervical cancer is predominant in women who experience imprisonement, and a lack of cancer screening is contributing to this phenomena as shown in the following paper: 

Kouyoumdjian, Fiona G., et al. "Cervical Cancer Screening Access for Women Who Experience Imprisonment in Ontario, Canada." JAMA network open 1.8 (2018): e185637-e185637.

Implementing a similar type of strategy to assess risk of cancer within imprisoned women in Ontario, and choosing the most at risk subset of women for screening could be a potential solution to the scarce screening resources available. 
