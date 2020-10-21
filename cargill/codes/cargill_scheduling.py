# Import various libraries

import sys
import numpy as np
from numpy import linalg as la
import gurobipy as gp
from gurobipy import *

# Set parameters
NumLoads = 10
NumDrivers = 3
NumPeriods = 30
NumLines = 2




# Define problem parameters
GL = {
     1 : 1,
     2 : 2,
     3 : 3,
     4 : 4,
     5 : 3,
     6 : 3,
     7 : 2,
     8 : 1,
     9 : 5,
     10 : 4
    }

GU = {
     1 : 25,
     2 : 27,
     3 : 27,
     4 : 25,
     5 : 23,
     6 : 23,
     7 : 28,
     8 : 30,
     9 : 30,
     10 : 30
    }

SL = {
    1 : 1,
    2 : 5,
    3 : 10
    }

SU = {
    1 : 22,
    2 : 25,
    3 : 30
    }

# Set up necessary sets
I = set(range(1, NumLoads + 1)) # set of loads
J = set(range(1, NumDrivers + 1))  # set of drivers
T = list(range(1, NumPeriods + 1)) # set of periods
L = set(range(1, NumLines + 1)) # set of lines

Im = { 1, 2, 3 }
Ir = { 4, 5, 6, 7 }
In = I.difference(Im.union(Ir))

ld = {i : 1 for i in I}
dt = {
     1 : 1,
     2 : 2,
     3 : 1,
     4 : 3,
     5 : 4,
     6 : 1,
     7 : 1,
     8 : 1,
     9 : 3,
     10 : 2
    }


ul = {i : 1 for i in I}

TR = { key : ld[key] + 2*dt[key] + ul[key] for key in dt.keys() }
HTR = { key : ld[key] + dt[key] for key in dt.keys() }


key_min = min(TR.keys(), key=(lambda k: TR[k]))

HTR_key_min = min(HTR.keys(), key=(lambda k: HTR[k]))
    
P = {
     (t,l) : 1 for t in T for l in L
    }

# Create the model

Cargill_model = gp.Model("Cargill Scheduling Model")

# Define variables
x = Cargill_model.addVars( I, J, T, L, vtype = GRB.BINARY, name = "Driver start")
y = Cargill_model.addVars( I, J, T, L, vtype = GRB.BINARY, name = "Driver finish")
w = Cargill_model.addVars( J, T, vtype = GRB.BINARY, name = "Driver busy?")


constr_1 = Cargill_model.addConstrs( x.sum('*', j, t, '*') <= 1 for j in J for t in T )
constr_2 = Cargill_model.addConstrs( y.sum('*', j, t, l) <= 1 for j in J for t in T for l in L )
constr_3 = Cargill_model.addConstrs( y[i, j, t + TR[i], l] >= x[i, j, t, l]
                                     for i in I for j in J for t in T for l in L
                                     if t + TR[i] <= T[-1])
constr_4 = Cargill_model.addConstrs( x.sum('*', '*', t, l) <= P[t,l] for t in T for l in L )
constr_5 = Cargill_model.addConstrs( x[i, j, t, l]*GL[i] <= (t + HTR[i])*x[i, j, t, l]
                                     for i in I for j in J for t in T for l in L )
constr_6 = Cargill_model.addConstrs( x[i, j, t, l]*GU[i] >= (t + HTR[i])*x[i, j, t, l]
                                     for i in I for j in J for t in T for l in L )
constr_7 = Cargill_model.addConstrs( SL[j]*w[j, t] <= t*w[j, t] for j in J for t in T )
constr_8 = Cargill_model.addConstrs( SU[j]*w[j, t] >= t*w[j, t] for j in J for t in T )
constr_9 = Cargill_model.addConstrs( quicksum(x[i, j, t - ld[i] - dt[i], l] for i in I if t >= 1 + ld[i] + dt[i] ) <= 1
                                     for j in J for l in L for t in T )
constr_10 = Cargill_model.addConstrs( quicksum(w[j, tp] for tp in T if tp >= t and tp <= t + dt[i] - 1)
                                      >= dt[i]*x[i, j, t, l] for i in I for j in J for t in T for l in L )


obj = w.sum('*', '*') + quicksum( TR[i]*x[i, j, t, l] for i in I for j in J for t in T for l in L)

Cargill_model.setObjective(obj, GRB.MAXIMIZE)

def printSolution():
    if Cargill_model.status == GRB.OPTIMAL:
        print('\nCost: %g' % Cargill_model.objVal)
        print('\nStart:')
        xx = Cargill_model.getAttr('x', x)
        for i in I:
            for j in J:
                for t in T:
                    for l in L:
                        if x[i, j, t, l].x > 0.0001:
                            print('%s %s %s %s %g' % (i, j, t, l, xx[i, j, t, l]))
    else:
        print('No solution')
        
#Chemo_model.modelSense = GRB.MINIMIZE
Cargill_model.optimize()
printSolution()



