#File contains code for finding the values of R0


import Constants as c
import ContactMatrices as cm

def get_C(tier):
        newC = cm.SchoolRestrictions[tier] *cm.School
        newC += cm.HomeRestrictions[tier] * cm.Home
        newC += cm.WorkRestrictions[tier] * cm.Work
        newC += cm.OtherRestrictions[tier] * cm.Other
        return newC

N = [239800,269800,267100,177200,45900]

def get_R0(C):

    evalues = [0,0,0,0,0]

    for i in range(5):
        evalues[i] = (c.rho[i]*c.betaA/c.gammaC + (1-c.rho[i])*c.betaA/c.gammaR)*N[i]*(C[i,1]/N[1]+C[i,2]/N[2]+C[i,3]/N[3]+C[i,4]/N[4]+C[i,0]/N[0])

    return max(evalues)

def get_R0B(C):

    evalues = [0,0,0,0,0]

    for i in range(5):
        evalues[i] = (c.rho[i]*c.betaB/c.gammaC + (1-c.rho[i])*c.betaB/c.gammaR)*N[i]*(C[i,1]/N[1]+C[i,2]/N[2]+C[i,3]/N[3]+C[i,4]/N[4]+C[i,0]/N[0])

    return max(evalues)

R0 = [0,0,0]
R0B = [0,0,0]


C = get_C(0)
print(get_R0(C))
print(get_R0B(C))
print('\n')

C = get_C(1)
print(get_R0(C))
print(get_R0B(C))
print('\n')

C = get_C(4)
print(get_R0(C))
print(get_R0B(C))
print('\n')
