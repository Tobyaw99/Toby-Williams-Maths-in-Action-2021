#File contains code for the control case 1 simulation from section 7
import Models as m

#initial conditions
zA = [239790,0,10,0,0,0,0,0,0,0,0,      #0-19
      269790,0,10,0,0,0,0,0,0,0,0,      #20-39
      267090,0,10,0,0,0,0,0,0,0,0,      #40-59
      177190,0,10,0,0,0,0,0,0,0,0,      #60-79
      45890,0,10,0,0,0,0,0,0,0,0]       #80+


A = m.Model("No restricitons")
A.set_t(500)
A.set_conditions(zA)
A.set_C()
A.run_simulation()

print("Final Size =" , A.get_final_size())
print("Final Death Count =" , A.results[-1,5,6])

print("Peak Infections =" , A.get_peak_infections())
print("Peak Hospitalisations =" , A.get_peak_hospitalisations())
print("Duration =" , A.get_duration())

A.plot_Re()
A.plot_SEIR(5)
