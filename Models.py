#File contains code for the models used in ther project

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
import ContactMatrices as cm
import Constants as c

############Standard Model Class##################

class Model():

    def __init__(self, name):
        self.name = name
        self.set_C()
        self.fig = 1
        #Create arrays for Re(t)
        self.RearrayA =  []
        self.RearrayB = []

    def set_C(self):
        self.C = cm.School + cm.Home + cm.Work + cm.Other #contact matrix

    def set_t(self,duration):                   #sets timespan for model
        self.t = np.array(range(0,duration+1))

    def set_conditions(self, conditions):       #sets initial conditions of model
        self.z0 = conditions
        self.N = [0,0,0,0,0]
        for i in range(5):
            self.N[i] = sum(self.z0[11*i:11*i+11])

    def equations(self,z,t):

        dzdt = [0]*55
        
        for i in range(5):
            j = 11*i
            sumSCA = 0
            sumSCB = 0
            
            for k in range(5):
                sumSCA += self.C[i,k] * (z[11*k+2] + z[11*k+3])/self.N[k]
                sumSCB += self.C[i,k] * (z[11*k+5] + z[11*k+6])/self.N[k]

  
            dzdt[j]   = -c.betaA * z[j] * sumSCA - c.betaB * z[j] * sumSCB                              #0 susceptible

            #Strain A
            dzdt[j+1] = c.betaA * z[j] * sumSCA - c.kappa * z[j+1]                                      #1 pre infectious
            dzdt[j+2] = c.rho[i] * c.kappa * z[j+1] - c.gammaC * z[j+2]                                 #2 infected subclinical 1
            dzdt[j+3] = (1 - c.rho[i])*c.kappa*z[j+1] - c.gammaR*z[j+3]                                 #3 infected subclinical 2

            #Strain B
            dzdt[j+4] = c.betaB * z[j] * sumSCB - c.kappa * z[j+4]                                      #4 pre infectious
            dzdt[j+5] = c.rho[i] * c.kappa * z[j+4] - c.gammaC * z[j+5]                                 #5 infected subclinical 1
            dzdt[j+6] = (1 - c.rho[i])*c.kappa*z[j+4] - c.gammaR*z[j+6]                                 #6 infected subclinical 2
            
            dzdt[j+7] = c.rhoprime[i] * c.gammaC *(z[j+2] + z[j+5]) - c.mu * z[j+7]                     #7 infected clinical 1
            dzdt[j+8] = (1-c.rhoprime[i])*c.gammaC*(z[j+2]+z[j+5])-c.gammaRC*z[j+8]+(1-c.phi)*c.mu*z[j+7]    #8 infected clinical 2
            dzdt[j+9] = c.gammaRC * z[j+8] + c.gammaR * (z[j+3] + z[j+6])                               #9 resistant
            dzdt[j+10] = c.phi*c.mu * z[j+7]                                                            #10 deceased

        return dzdt

    #########Formatting##############
    #Takes the values from z and formats them into a more intuitive array with rows of the form
    #[susceptibe, exposed, infected, hospitalized (all), hospitalized (ICU), recovered, deceased] for each age group
    #and the final row will be for the total over all ages. [time,age,field]

    def format(self,z):

        new_array = np.array([[[0] * 7] *6] * len(z))

        #t = time, i = age group
        for t in range(len(z)):
            for i in range(5):
                k=11*i
                new_array[t,i,0] = z[t,k]                            
                new_array[t,i,1] = z[t,k+1] + z[t,k+4]               
                new_array[t,i,2] = z[t,k+2] + z[t,k+3] + z[t,k+5] + z[t,k+6] + z[t,k+7] +z[t,k+8]
                new_array[t,i,3] = z[t,k+7] + z[t,k+8]
                new_array[t,i,4] = z[t,k+7]
                new_array[t,i,5] = z[t,k+9]
                new_array[t,i,6] = z[t,k+10]

                #Totals
                new_array[t,5,0] += z[t,k]
                new_array[t,5,1] += z[t,k+1] + z[t,k+4]               
                new_array[t,5,2] += z[t,k+2] + z[t,k+3] + z[t,k+5] + z[t,k+6] + z[t,k+7] +z[t,k+8]
                new_array[t,5,3] += z[t,k+7] + z[t,k+8]
                new_array[t,5,4] += z[t,k+7]
                new_array[t,5,5] += z[t,k+9]
                new_array[t,5,6] += z[t,k+10]
    

        return new_array
    


    def run_simulation(self):                #runs simluation: solves equations
        self.raw_results = odeint(self.equations,self.z0,self.t)
        self.results = self.format(self.raw_results)

        #Sort our Rearrays
        if self.RearrayA == []:      #First compute R0
            self.RearrayA.append(self.get_ReA(0))
            self.RearrayB.append(self.get_ReB(0))
            
        for t in self.t[1:]:    #Don't do this for the first t or you'll get repeated values
            self.RearrayA.append(self.get_ReA(t))
            self.RearrayB.append(self.get_ReB(t))


    #Move one case from Si into ISC1
    def add_variantB(self):
        for i in range(5):
            self.raw_results[-1,11*i] = self.raw_results[-1,11*i] - 10
            self.raw_results[-1,11*i+5] = self.raw_results[-1,11*i+5] + 10


    ###################Caclulating Metrics#####################
    #R0 values
    def get_R0A(self):
        evalues = [0,0,0,0,0]

        for i in range(5):
            evalues[i] = (c.rho[i]*c.betaA/c.gammaC + (1-c.rho[i])*c.betaA/c.gammaR)*self.N[i]*(self.C[i,1]/self.N[1]+self.C[i,2]/self.N[2]+self.C[i,3]/self.N[3]+self.C[i,4]/self.N[4]+self.C[i,0]/self.N[0])
        return max(evalues)

    def get_R0B(self):
        evalues = [0,0,0,0,0]

        for i in range(5):
            evalues[i] = (c.rho[i]*c.betaB/c.gammaC + (1-c.rho[i])*c.betaB/c.gammaR)*self.N[i]*(self.C[i,1]/self.N[1]+self.C[i,2]/self.N[2]+self.C[i,3]/self.N[3]+self.C[i,4]/self.N[4]+self.C[i,0]/self.N[0])
    
        return max(evalues)

        
    #Re(t) values
    def get_ReA(self,t):
        evalues = [0,0,0,0,0]

        for i in range(5):
             evalues[i] = (c.rho[i]*c.betaA/c.gammaC + (1-c.rho[i])*c.betaA/c.gammaR)*self.results[t,i,0]*(self.C[i,1]/self.N[1]+self.C[i,2]/self.N[2]+self.C[i,3]/self.N[3]+self.C[i,4]/self.N[4]+self.C[i,0]/self.N[0])
        return max(evalues)

    def get_ReB(self,t):
        evalues = [0,0,0,0,0]

        for i in range(5):
             evalues[i] = (c.rho[i]*c.betaB/c.gammaC + (1-c.rho[i])*c.betaB/c.gammaR)*self.results[t,i,0]*(self.C[i,1]/self.N[1]+self.C[i,2]/self.N[2]+self.C[i,3]/self.N[3]+self.C[i,4]/self.N[4]+self.C[i,0]/self.N[0])
        return max(evalues)

    #Calculate final size
    def get_final_size(self):
        return (self.results[-1,5,5]+ self.results[-1,5,6]) - (self.results[0,5,5]+ self.results[0,5,6])

    #Peak infections
    def get_peak_infections(self):
        return [max(self.results[:,5,2]), np.where(self.results[:,5,2] == max(self.results[:,5,2]))[0][0]]

    #Peak hospitalisations
    def get_peak_hospitalisations(self):
        return [max(self.results[:,5,3]), np.where(self.results[:,5,3] == max(self.results[:,5,3]))[0][0]]

    #Duration
    def get_duration(self):
        return np.where(self.results[:,5,2] > 1)[-1][-1] + 1

    #######################Plotting#########################
    #Plot specific fields 
    def plot_specific(self,nums):

        plt.figure(self.fig)
        for i in nums:
            try:
                plt.plot(self.t,self.results[:,i[0],i[1]],label = c.labels[i[0],i[1]])
            except:
                print("Error, invalid index.")
                return

        plt.ylabel('People')
        plt.xlabel('Time')
        plt.legend(loc='best')
        plt.show()
        self.fig = self.fig + 1


    #Plot all fields for a specified age range or the totals
    def plot_age(self,age):

        plt.figure(self.fig)
        for i in range(7):
            plt.plot(self.t,self.results[:,age,i],label = c.labels[age,i])

        plt.ylabel('People')
        plt.xlabel('Time')
        plt.legend(loc='best')
        plt.show()
        self.fig = self.fig + 1

    #Plot all ages for a specified field not including the total
    def plot_field(self,field):

        plt.figure(self.fig)
        for i in range(5):
            plt.plot(self.t,self.results[:,i,field],label = c.labels[i,field])

        plt.ylabel('People')
        plt.xlabel('Time')
        plt.legend(loc='best')
        plt.show()
        self.fig = self.fig + 1

    #Plot the SEIR classes
    def plot_SEIR(self,age):

        plt.figure(self.fig)
        plt.plot(self.t,self.results[:,age,0],label = c.labels[age,0])
        plt.plot(self.t,self.results[:,age,1],label = c.labels[age,1])
        plt.plot(self.t,self.results[:,age,2],label = c.labels[age,2])
        plt.plot(self.t,self.results[:,age,5],label = c.labels[age,5])

        plt.ylabel('People')
        plt.xlabel('Time')
        plt.legend(loc='best')
        plt.show()
        self.fig = self.fig + 1

    #Plot the Re value over time
    def plot_Re(self):

        plt.figure(self.fig)

        plt.plot(self.t,self.RearrayA,label = "Re(t) variant A",color='r')
        plt.plot(self.t,self.RearrayB,label = "Re(t) variant B",color='b')

        plt.ylabel('Reproduction Number')
        plt.xlabel('Time')
        plt.legend(loc='best')
        plt.show()
        self.fig = self.fig + 1


##############Tiered Model###########################

class TierModel(Model):

    def __init__(self, name):
        self.name = name
        self.tier = 0       #default tier, tier 0 is no restrictions, tier 5 is full lockdown
        self.fig = 1
        #Create arrays for Re(t)
        self.RearrayA =  []
        self.RearrayB = []


    def set_C(self):
        newC = cm.SchoolRestrictions[self.tier] *cm.School
        newC += cm.HomeRestrictions[self.tier] * cm.Home
        newC += cm.WorkRestrictions[self.tier] * cm.Work
        newC += cm.OtherRestrictions[self.tier] * cm.Other
        self.C = newC
        

    def set_tier(self,new_tier):
        if new_tier >= 0 and new_tier <= 5:
            self.tier = new_tier
            self.set_C()
        else:
            print("Invalid tier, must be a number between 0 and 5.")

    def continue_simulation(self,duration,new_tier):        #runs the simulation again with new tier and time

        #Sorts out time for the plotting
        list_old_t = self.t.tolist()        
        list_new_t = range(list_old_t[-1],list_old_t[-1] + duration+1)
        del list_old_t[-1]
        list_old_t.extend(list_new_t)
        arr_new_t = np.asarray(list_old_t)
        
        #Setting the new simulation parameters
        self.z0 = self.raw_results[-1]
        old_results = self.raw_results[0:-1]
        self.set_t(duration)
        self.set_tier(new_tier)

        self.run_simulation()

        #Combines the old results with the new
        old = old_results.tolist()
        new = self.raw_results.tolist()
        old.extend(new)
        
        self.raw_results = np.asarray(old)
        self.results = self.format(self.raw_results)

        self.t = arr_new_t



###############Vaccine Model####################

class VaccineModel(TierModel):

    def __init__(self, name):
        self.name = name
        self.numvaxxed = 0   #used for finding number of people vaccinated for section 8

        self.tier = 0 
        self.order = [4,3,2,1,0]     #default order
        self.day = 0        #keeps track of what day of the simulation we are on
        self.temp_t = np.array([0,1])
        self.t = [0,1]
        self.fig = 1
        #Create arrays for Re(t)
        self.RearrayA =  []
        self.RearrayB = []
        
        

    def set_t(self,duration):
        self.final_t = np.array([0]*(duration+1))
        for i in range(duration + 1):
            self.final_t[i] = i

    def set_temp_t(self):
        self.temp_t = np.array([self.day,self.day+1])

    def set_conditions(self, conditions):       #now also sets  the initial unvaccinated vaules
        self.z0 = conditions
        self.N = [0,0,0,0,0]
        for i in range(5):
            self.N[i] = sum(self.z0[11*i:11*i+11])
        
        

    def set_vaccine_order(self,new_order):
        if len(new_order) != 5:
            print("The array parameter must have length 5. Default order = [4,3,2,1,0]")
            return
        for i in new_order:
            if i > 4 or i < 0:
                print("All elements of the array must be between 0 and 4. Default order = [4,3,2,1,0]")
                return
        self.order = new_order

        
    def get_vaxxed(self,t):

        v = [0] * 5
        to_be_vaxxed = self.rollout(t)   #keeps track of remaining vaccinations left to distribute
        if to_be_vaxxed == 0:
            return v            #ends early if no vaccines to be distributed

        for j in range(5):        #vaccinates specified age groups first people first

            k = self.order[j]
        
            if self.unvaxxed[k] != 0:           #moves into the next age group if an age group is all vaccinated
                if self.unvaxxed[k] < to_be_vaxxed:
                    v[k] = self.unvaxxed[k]
                    to_be_vaxxed -= self.unvaxxed[k]
                    self.unvaxxed[k] = 0

                else:                       #oteherwise only members of one group need to be vaccinated
                    v[k] = to_be_vaxxed
                    self.unvaxxed[k] -= to_be_vaxxed
                    to_be_vaxxed = 0
        return v

    def rollout(self,t):        #Number of vaccines rolled out every day
        return 8300

    def vax_simulation(self,duration):

        #Sets up the initial values for unvaxxed, the maximum number of people in each
        #age group who can be vaccinated
        self.unvaxxed = np.array([0]*5)
        for i in range(5):
            self.unvaxxed[i] = self.z0[11*i] * c.mvu[i]

        v = self.get_vaxxed(0)

        for k in range(5):
            self.z0[11*k] -= v[k]*c.tau *(self.z0[11*k]/self.N[k])
            self.z0[11*k+9] += v[k]*c.tau *(self.z0[11*k]/self.N[k])
            self.numvaxxed += v[k]*c.tau

        self.run_simulation()           #run the simulation once

        for day in range (duration):        #run the rest of the days
            self.day = day+1     
            self.continue_simulation(self.tier)

        self.t = self.final_t

        new_raw_results = np.array(self.raw_results[0:-1])
        new_results = np.array(self.results[0:-1])
        
        self.raw_results = new_raw_results
        self.results = new_results
        

    def continue_simulation(self,new_tier):        #runs the simulation again with new tier and time
    
        #Setting the new simulation parameters
        self.z0 = self.raw_results[-1]
        old_results = self.raw_results[0:-1]

        N = [0] * 5
        v = self.get_vaxxed(self.day) 
        for k in range(5):
            N[k] = np.sum(self.z0[11*k:11*k+10])
            self.z0[11*k] -= v[k]*c.tau *(self.z0[11*k]/N[k])
            self.z0[11*k+9] += v[k]*c.tau *(self.z0[11*k]/N[k])
            self.numvaxxed += v[k]*c.tau
        
        self.set_tier(new_tier)

        self.run_simulation()


        #Combines the old results with the new
        old = old_results.tolist()
        new = self.raw_results.tolist()
        old.extend(new)
        
        self.raw_results = np.asarray(old)
        self.results = self.format(self.raw_results)


class BadVaccineModel(VaccineModel):

    def __init__(self, name):
        self.name = name
        
        self.numvaxxed = 0
        self.tier = 0 
        self.order = [4,3,2,1,0]     #default order
        self.day = 0        #keeps track of what day of the simulation we are on
        self.temp_t = np.array([0,1])
        self.t = [0,1]
        self.fig = 1
        #Create arrays for Re(t)
        self.RearrayA =  []
        self.RearrayB = []


    def get_vaxxed(self,t):

        n=0     #number of groups to divide by

        v = [0] * 5
        redist = 0   #Vaccines to be redistributed
        to_be_vaxxed = self.rollout(t)   #keeps track of remaining vaccinations left to distribute
        if to_be_vaxxed == 0:
            return v            #ends early if no vaccines to be distributed

        for i in range(5):
            if t != 0:
                v[i] = to_be_vaxxed*(self.raw_results[t,11*i]/sum(self.z0))

            else:
                v[i] = to_be_vaxxed*(self.z0[11*i]/sum(self.z0))

            if v[i] > self.unvaxxed[i]:
                redist = redist + v[i]
                v[i] = self.unvaxxed[i]
                self.unvaxxed[i] = 0

            else:
                n=n+1

        for i in range(5):
            if self.unvaxxed[i] != 0:
                v[i] = v[i] + redist/n
        return v

    
        

#########Testing##############
#initial conditions
"""      
z0 = [23500,0,100,0,0,100,0,0,0,0,0,   #0-19
      26500,0,100,0,0,0,0,0,0,0,0,      #20-39
      26500,0,100,0,0,0,0,0,0,0,0,      #40-59
      19100,0,100,0,0,0,0,0,0,0,0,      #60-79
      4400,0,100,0,0,0,0,0,0,0,0]       #80+

"""




