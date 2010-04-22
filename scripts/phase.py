#!/usr/bin/python

from math import *
from numpy import *
from scipy.integrate import quad
from scipy import interpolate
from scipy.optimize import fsolve
from scipy.optimize import fmin
from scipy.optimize import fmin_powell
from scipy.optimize import fmin_cg
from scipy.optimize import fmin_ncg
from scipy.optimize import fmin_bfgs


################################################################################
################################################################################
def betaA_res(rho):
    return interpolate.splint(0.0, rho, Zrho_fit)


def intU_func(beta, rho):
    return interpolate.splev(rho, intU_fit[beta])


def Z_func(beta, rho):
    return interpolate.splev(rho, Z_fit[beta])


def Zrho_func(rho):
    return interpolate.splev(rho, Zrho_fit)


def betamu_func(beta, rho):
    betamu_ig  = log(rho)
    betamu_res = betaA_res(rho) + intU_func(beta, rho) + Z_func(beta, rho)-1.0
    return betamu_ig + betamu_res


def f_res(x, beta):
    rho_g,  rho_l = x
    Z_g = Z_func(beta, rho_g)
    Z_l = Z_func(beta, rho_l)
    p_g = rho_g*Z_g
    p_l = rho_l*Z_l
    I1 = betaA_res(rho_l) - betaA_res(rho_g)
    I2 = intU_func(beta, rho_l) - intU_func(beta, rho_g)
    I3 = log(rho_l/rho_g)
    residual = []
    residual.append(p_g-p_l)
    residual.append(I1+I2+I3+Z_l-Z_g)
    print "-----------"
    print x[0]*16.0*pi/6.0, x[1]*16.0*pi/6.0
    print residual
    print "I1=", I1, "I2=", I2
    return residual


################################################################################
################################################################################
# import simulation data
fit_data = []
infile = open('data.dat', 'r')
for line in infile:
    line = line.strip()
    data = line.split()
    if (len(data)==5):
        beta = float(data[0]) 
        rho  = float(data[1])/16.0
        y    = 16.0*rho*pi/6.0
        p    = float(data[3])
        U    = float(data[4])/27.0
        Z = beta*p/rho
        Zrho = (Z-1.0)/rho
        fit_data.append([beta, rho, Z, U])
#        print T, rho, Z, U
infile.close()
#print fit_data


# sort data into isotherms and isochores
rho_isotherm = {}
isotherm = {}
isochore = {}
Zrho_data = []
for data in fit_data:
    beta = data[0]
    rho  = data[1]
    Z    = data[2]
    U    = data[3]
    
    if (beta not in isotherm):
        isotherm[beta] = [[], [], []]
    isotherm[beta][0].append(rho)
    isotherm[beta][1].append(Z)
    isotherm[beta][2].append(U)

    if (rho not in isochore):
        isochore[rho] = [[], [], []]
    isochore[rho][0].append(beta)
    isochore[rho][1].append(Z)
    isochore[rho][2].append(U)

    if (beta == 0.001):
        Zrho = (Z-1.0)/rho
        Zrho_data.append(Zrho)

#
Zrho_fit = interpolate.splrep(isotherm[0.001][0], Zrho_data, s=0)


# fit U along isochores
U_fit = {}
for rho in isochore.keys():
    print "rho=", rho
    print len(isochore[rho][0]), isochore[rho][0]
    print len(isochore[rho][2]), isochore[rho][2]
    x = isochore[rho][0]
    y = isochore[rho][2]
    U_fit[rho] = interpolate.splrep(x, y, s=0)
#    U_fit[rho] = interpolate.splrep(isochore[rho][0], isochore[rho][2], s=0)


# fit Z along isotherms
Z_fit = {}
intU_fit = {}
for beta in isotherm.keys():
    Z_fit[beta] = interpolate.splrep(isotherm[beta][0], isotherm[beta][1], s=0)
    intU_data = []
    for rho in isotherm[beta][0]:
        intU = interpolate.splint(0.0, beta, U_fit[rho])
        intU_data.append(intU)
    intU_fit[beta] = interpolate.splrep(isotherm[beta][0], intU_data, s=0)


beta = 0.45
rho_g = 0.025/16.*6.0/pi
rho_l = 0.225/16.*6.0/pi
x = [rho_g, rho_l]
f = fsolve(f_res, x, args=(beta,))
print f


# create a plot of mu vs p to make first initial guess
for rho in isotherm[beta][0]:
    p = rho*Z_func(beta, rho)/beta
    print rho, p, betamu_func(beta, rho) 


