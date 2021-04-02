#File contains the simulations used to generate the vaccine graph used in the poster

import Models as m
import matplotlib.pyplot as plt


z0=[2.35446336e+05, 1.06129535e+02, 2.44740713e+00, 1.50515538e+02,
       5.23539520e+01, 8.95966732e-01, 5.50740326e+01, 7.08978203e-02,
       3.41912381e+00, 3.98186172e+03, 8.95383799e-01, 2.65857994e+05,
       8.96902524e+01, 6.98121757e+00, 1.22300589e+02, 4.42592070e+01,
       2.55568222e+00, 4.47638127e+01, 3.16149244e-01, 9.69457591e+00,
       3.61784855e+03, 3.59634526e+00, 2.63658904e+05, 7.87713285e+01,
       1.40770415e+01, 9.94474871e+01, 3.89099330e+01, 5.16153734e+00,
       3.64605565e+01, 1.51222229e+00, 1.91070109e+01, 3.13139765e+03,
       1.62510649e+01, 1.76236705e+05, 2.06835842e+01, 6.20126140e+00,
       2.36124953e+01, 1.05889808e+01, 2.40913854e+00, 9.17152956e+00,
       1.31005870e+00, 8.28240539e+00, 8.65627166e+02, 1.54086118e+01,
       4.56897227e+04, 3.83176015e+00, 1.75760464e+00, 3.76945397e+00,
       2.09722596e+00, 7.59060223e-01, 1.62694792e+00, 4.97996712e-01,
       2.39552649e+00, 1.85955152e+02, 7.58658111e+00]
z0s = [z0] * 4

V0 = m.Model("No Vaccines")
V0.set_t(300)
V0.set_conditions(z0s[0])
V0.run_simulation()

V1 = m.VaccineModel("V1")
V1.set_t(300)
V1.set_conditions(z0s[1])
V1.set_C()
V1.set_vaccine_order([4,3,2,1,0])
V1.vax_simulation(300)

V2 = m.VaccineModel("V2")
V2.set_t(300)
V2.set_conditions(z0s[2])
V2.set_C()
V2.set_vaccine_order([0,1,2,3,4])
V2.vax_simulation(300)

V3 = m.BadVaccineModel("V3")
V3.set_t(300)
V3.set_conditions(z0s[3])
V3.set_C()
V3.set_vaccine_order([0,1,2,3,4])
V3.vax_simulation(300)

print(V3.results[-1,5,5]-V3.results[-1,5,6])
print(V3.results[-1,5,6])


V3.plot_field(6)

plt.figure(1)

plt.plot(V0.t,V0.results[:,5,6],label = 'No Vaccines',color='r')
plt.plot(range(0,301),V2.results[:,5,6],label = 'Alternative Vaccination Plan',color='g')

plt.plot(range(0,301),V1.results[:,5,6],label = 'Optimal Vaccination Plan',color='b')
plt.ylabel('Deaths')
plt.xlabel('Time')
plt.legend(loc='best')
plt.show()


