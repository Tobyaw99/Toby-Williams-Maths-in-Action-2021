##File contains the simulation for simulationg covid with the basic SEIR model



import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def model(z,t):

    dzdt = [0,0,0,0]

    dzdt[0] = -beta*z[0]*z[2]/sum(z)
    dzdt[1] = beta*z[0]*z[2]/sum(z) - sigma*z[1]
    dzdt[2] = sigma*z[1] - gamma*z[2]
    dzdt[3] = gamma*z[2]
    return dzdt





z0 = [1000000,0,1,0]
beta =  0.215
gamma = 1/7
sigma = 1/3
t = np.linspace(0,400)
results = odeint(model,z0,t)

plt.plot(t,results[:,0],'b',label="Susceptible")
plt.plot(t,results[:,1],'y',label="Exposed")
plt.plot(t,results[:,2],'r',label="Infected")
plt.plot(t,results[:,3],'g',label="Resistant")

plt.ylabel('Number of People')
plt.xlabel('Time (days)')
plt.legend(loc='best')
plt.show()
