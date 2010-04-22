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
#    print "-----------"
#    print x[0]*16.0*pi/6.0, x[1]*16.0*pi/6.0
#    print residual
#    print "I1=", I1, "I2=", I2
    return residual


################################################################################
################################################################################
import sys
from optparse import OptionParser


# parse command line arguments
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("--beta", type="float", dest="beta", default=0.454545)
parser.add_option("--y_g", type="float", dest="y_g", default=0.025)
parser.add_option("--y_l", type="float", dest="y_l", default=0.225)

(options, args) = parser.parse_args()



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


# determine minimum value of beta
beta_min = min(isotherm.keys())
#print beta_min, len(isotherm[beta_min][0])
for i in arange( len(isotherm[beta_min][0]) ):
    rho = isotherm[beta_min][0][i]
    Z   = isotherm[beta_min][1][i]
    Zrho = (Z-1.0)/rho
    Zrho_data.append(Zrho)
x = isotherm[beta_min][0]
y = Zrho_data
Zrho_fit = interpolate.splrep(x, y, s=0)


# fit U along isochores
U_fit = {}
for rho in isochore.keys():
    x = isochore[rho][0]
    y = isochore[rho][2]
    U_fit[rho] = interpolate.splrep(x, y, s=0)
#    U_fit[rho] = interpolate.splrep(isochore[rho][0], isochore[rho][2], s=0)


# fit Z along isotherms
Z_fit = {}
intU_fit = {}
for beta in isotherm.keys():

    intU_data = []
    for rho in isotherm[beta][0]:
        intU = interpolate.splint(beta_min, beta, U_fit[rho])
        intU_data.append(intU)
    intU_fit[beta] = interpolate.splrep(isotherm[beta][0], intU_data, s=0)

    x = isotherm[beta][0]
    y = isotherm[beta][1]
    x.insert(0, 0.0)
    y.insert(0, 1.0)
    Z_fit[beta] = interpolate.splrep(x, y, s=0)


beta = options.beta
#rho_g = options.y_g/16.*6.0/pi
#rho_l = options.y_l/16.*6.0/pi
#x = [rho_g, rho_l]
#x_soln = fsolve(f_res, x, args=(beta,))
#print 1.0/beta, x_soln[0], x_soln[1]
#print 1.0/beta, x_soln[0]*16.0*pi/6.0, x_soln[1]*16.0*pi/6.0

beta_list = []
for item in isotherm.keys():
    if (item >= beta):
        beta_list.append(item)
beta_list.sort()
#print beta_list

rho_g = options.y_g/16.*6.0/pi
rho_l = options.y_l/16.*6.0/pi
for beta in beta_list:
    x = [rho_g, rho_l]
    x_soln, infodict, ier, mesg = fsolve(f_res, x, args=(beta,), full_output=True)
    rho_g = x_soln[0]
    rho_l = x_soln[1]
#    print 1.0/beta, rho_g, rho_l
    if (ier == 1):
        print 1.0/beta, rho_g*16.0*pi/6.0, rho_l*16.0*pi/6.0
    else:
#        print mesg
        sys.exit(2)


## create a plot of mu vs p to make first initial guess
#for rho in isotherm[beta][0]:
#    p = rho*Z_func(beta, rho)/beta
#    print rho, p, betamu_func(beta, rho) 


