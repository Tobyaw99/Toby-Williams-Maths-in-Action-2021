#File contains code for the lockdown case 2 simulation from section 7

import Models as m
import matplotlib.pyplot as plt

#initial conditions

zB = [239790,0,10,0,0,0,0,0,0,0,0,      #0-19
      269790,0,10,0,0,0,0,0,0,0,0,      #20-39
      267090,0,10,0,0,0,0,0,0,0,0,      #40-59
      177190,0,10,0,0,0,0,0,0,0,0,      #60-79
      45890,0,10,0,0,0,0,0,0,0,0]       #80+


D = m.TierModel("UK Govt Response")
D.set_t(30)
D.set_conditions(zB)
D.set_C()
D.run_simulation()

D.continue_simulation(75,5)
D.continue_simulation(95,1)

D.add_variantB()
D.continue_simulation(1000,4)
print("Final Size =" , D.get_final_size())
print("Final Death Count =" , D.results[-1,5,6])

print("Peak Infections =" , D.get_peak_infections())
print("Peak Hospitalisations =" , D.get_peak_hospitalisations())
print("Duration =" , D.get_duration())

D.plot_Re()
D.plot_SEIR(5)


plt.figure(1)

plt.plot(range(0,30),D.RearrayA[0:30],label = "Re(t) variant A",color='r')
plt.plot(range(31,105),D.RearrayA[31:105],color='r')
plt.plot(range(106,200),D.RearrayA[106:200],color='r')
plt.plot(range(201,1200),D.RearrayA[201:1200],color='r')
plt.plot(range(0,30),D.RearrayB[0:30],label = "Re(t) variant B",color='b')
plt.plot(range(31,105),D.RearrayB[31:105],color='b')
plt.plot(range(106,200),D.RearrayB[106:200],color='b')
plt.plot(range(201,1200),D.RearrayB[201:1200],color='b')
plt.ylabel('Reproduction Number')
plt.xlabel('Time')
plt.legend(loc='best')
plt.show()

