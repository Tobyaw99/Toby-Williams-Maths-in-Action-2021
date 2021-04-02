import numpy as np

#########Constants############

betaA = 0.079      #prob of contact causing infection for strain A
betaB = 1.7*betaA           #strain B beta

kappa = 1/7.7     #rate of preinfectious becoming infectious
rho = [0.016,0.054,0.124,0.208,0.318]       #prob of preinfectious going to ISC1
rhoprime = [0.0348,0.0544,0.129,0.249,0.322]     #prob of ISC1 going into IC1
gammaC = 1/10     #rate of going from ISC1 to IC1
gammaR = 1/10    #rate of recovery from SC2
gammaRC = 1/10  #rate of recovery from IC2
mu = 1/6           #rate of leaving the ICU
tau = 0.9       #vaccine effectiveness
mvu = [0.7,0.66,0.68,0.75,0.85]     #max vaccine uptake proportion for each age
phi = 0.316     #Proportion of ICU patients who die

#labels are used for plotting
labels = np.array([["0 - 19 Susceptible", "0 - 19 Preinfectious", "0 - 19 Infectious","0 - 19 Hospitalized (all)", "0 - 19 Hospitalised (ICU)", "0 - 19 Resistant", "0 - 19 Deaths"],          
          ["20 - 39 Susceptible", "20 - 39 Preinfectious", "20 - 39 Infectious","20 - 39 Hospitalized (all)", "20 - 39 Hospitalized (ICU)", "20 - 39 Resistant", "20 - 39 Deaths"],
          ["40 - 59 Susceptible", "40 - 59 Preinfectious", "40 - 59 Infectious","40 - 59 Hospitalized (all)", "40 - 59 Hospitalized (ICU)", "40 - 59 Resistant", "40 - 59 Deaths"],
          ["60 - 79 Susceptible", "60 - 79 Preinfectious", "60 - 79 Infectious","40 - 59 Hospitalized (all)", "40 - 59 Hospitalized (ICU)", "60 - 79 Resistant", "60 - 79 Deaths"],
          ["80+ Susceptible", "80+ Preinfectious", "80+ Infectious", "80+ Hospitalized (all)", "80+ Hospitalized (ICU)", "80+ Resistant", "80+ Deaths"],
          ["Total Susceptible", "Total Preinfectious", "Total Infectious", "Total Hospitalized (all)","Total Hospitalized (ICU)", "Total Resistant", "Total Deaths"]])
