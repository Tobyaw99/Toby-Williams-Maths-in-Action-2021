#File contains code for the lockdown case 1 simulation from section 7

import Models as m
import matplotlib.pyplot as plt

#initial conditions

zB = [239790,0,10,0,0,0,0,0,0,0,0,      #0-19
      269790,0,10,0,0,0,0,0,0,0,0,      #20-39
      267090,0,10,0,0,0,0,0,0,0,0,      #40-59
      177190,0,10,0,0,0,0,0,0,0,0,      #60-79
      45890,0,10,0,0,0,0,0,0,0,0]       #80+


C = m.TierModel("National Lockdown")
C.set_t(30)
C.set_conditions(zB)
C.set_C()
C.run_simulation()

C.continue_simulation(170,5)

C.add_variantB()
C.continue_simulation(1000,5)
print("Final Size =" , C.get_final_size())
print("Final Death Count =" , C.results[-1,5,6])

print("Peak Infections =" , C.get_peak_infections())
print("Peak Hospitalisations =" , C.get_peak_hospitalisations())
print("Duration =" , C.get_duration())

C.plot_Re()
C.plot_SEIR(5)
C.plot_specific([[5,3]])
C.plot_field(6)

plt.figure(1)

plt.plot(range(0,30),C.RearrayA[0:30],label = "Re(t) variant A",color='r')
plt.plot(range(31,200),C.RearrayA[31:200],color='r')
plt.plot(range(201,1200),C.RearrayA[201:1200],color='r')
plt.plot(range(0,30),C.RearrayB[0:30],label = "Re(t) variant B",color='b')
plt.plot(range(31,200),C.RearrayB[31:200],color='b')
plt.plot(range(201,1200),C.RearrayB[201:1200],color='b')
plt.ylabel('Reproduction Number')
plt.xlabel('Time')
plt.legend(loc='best')
plt.show()
