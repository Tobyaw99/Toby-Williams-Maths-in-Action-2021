#File contains code for the control case 2 simulation from section 7
import Models as m

#initial conditions

zB = [239790,0,10,0,0,0,0,0,0,0,0,      #0-19
      269790,0,10,0,0,0,0,0,0,0,0,      #20-39
      267090,0,10,0,0,0,0,0,0,0,0,      #40-59
      177190,0,10,0,0,0,0,0,0,0,0,      #60-79
      45890,0,10,0,0,0,0,0,0,0,0]       #80+


B = m.TierModel("No restricitons: Vairant B introduced")
B.set_t(200)
B.set_conditions(zB)
B.set_C()
B.run_simulation()

B.add_variantB()
B.continue_simulation(300,0)
print("Final Size =" , B.get_final_size())
print("Final Death Count =" , B.results[-1,5,6])

print("Peak Infections =" , B.get_peak_infections())
print("Peak Hospitalisations =" , B.get_peak_hospitalisations())
print("Duration =" , B.get_duration())

B.plot_Re()
B.plot_SEIR(5)
B.plot_specific([[5,3]])
