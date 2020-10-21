# Import various libraries

import sys
import numpy as np
import pandas as pd

from numpy import linalg as la
import gurobipy as gp
from gurobipy import *

# Set parameters

SlotLen = 15 # minutes (unit)
NumLines = 2

df_customers = pd.read_csv("data/processed/customers.csv")
df_drivers = pd.read_csv("data/processed/drivers.csv")
df_transport_costs = pd.read_csv("data/processed/transport_costs.csv")

NumDrivers = len(df_drivers['driver__c'])

startTime = df_drivers['starttime']
endTime = df_drivers['endtime']

startTime = pd.to_datetime(startTime)
endTime = pd.to_datetime(endTime)

startTime_hour = list(startTime.dt.hour)
startTime_minute = list(startTime.dt.minute)
endTime_hour = list(endTime.dt.hour)
endTime_minute = list(endTime.dt.minute)


start_minute = [ int(startTime_minute[i]) + int(startTime_hour[i])*60 for i in range(len(startTime_hour)) ]
end_minute = [ int(endTime_minute[i]) + int(endTime_hour[i])*60 for i in range(len(endTime_hour)) ]


startTimeSlot = [ int(item / SlotLen) for item in start_minute ]
endTimeSlot = [ int(item / SlotLen) for item in end_minute ]


col_customerid = df_customers['unique_customer_id']
col_medicated = list(df_customers['non_medicated_yard__c'])
col_gl = list(df_customers['begin__c'])
col_gu = list(df_customers['end__c'])



transports = df_transport_costs[ df_transport_costs['unique_customer_id'].isin(col_customerid) ]
col_miles = transports['miles_round_trip__c']
col_mph = transports['avg_mph__c']
col_miles = list(col_miles)
col_mph = list(col_mph)

distance = [ int(np.ceil((col_miles[i] / col_mph[i])*30/SlotLen)) for i in range(len(col_mph)) ]

col_loads = df_customers['loads']
loads = list(col_loads)


I = set(range(1,sum(loads)+1))
Ir = { 1, 2, 3 }
In = { 94, 95, 96, 97, 98, 99, 100, 101, 115, 116, 117, 118, 119, 120 }
Im = I.difference(In.union(Ir))


T = list(range(min(startTimeSlot), max(endTimeSlot)))
J = list(df_drivers['driver__c'])
L = set(range(1, NumLines + 1)) # set of lines

for i in range(len(col_gl)):
    if col_gl[i] == ' ':
        col_gl[i] = min(T)
    else:
        col_gl[i] = int(((float(col_gl[i])/100)*60)/SlotLen)
        
for i in range(len(col_gu)):
    if col_gu[i] == ' ':
        col_gu[i] = max(T)
    else:
        col_gu[i] = int(((float(col_gu[i])/100)*60)/SlotLen)


dt = {}
GL = {}
GU = {}

counter = 1
for i in range(len(loads)):
    for j in range(loads[i]):
        dt[counter] = distance[i]
        GL[counter] = col_gl[i]
        GU[counter] = col_gu[i]
        counter += 1
        
ld = {i : 1 for i in I}
ul = {i : 1 for i in I}

SL = {
      J[i] : startTimeSlot[i] for i in range(len(J))
      }

SU = {
      J[i] : endTimeSlot[i] for i in range(len(J))
      }


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



# transports.to_csv('data/processed/tranports.csv', encoding='utf-8', index=False)











