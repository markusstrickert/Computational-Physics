#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 13:26:04 2019

@author: Markus Strickert
"""

from scipy import *
from matplotlib import pyplot as plt 
import numpy as np
import math
from pylab import figure, axes, pie, title, show

alpha = 1/2 #Constant recovery factor

#Defining function for our parameters
def SIS(infection,N,beta,h):
    time = 60 #Years on x axis. 
    m = round(time/h) #Number of iterations for the given time interval
    
    I = np.zeros(m) #Infection list.
    S = np.zeros(m) #Susceptible list.
    H = np.zeros(m) #Time step list.
    I_fraction = np.zeros(m) #Fraction infected population list. 

#introducing the initial values and infection fraction
    I[0] = infection  
    S[0] = N - I[0]
    I_fraction[0] = I[0]/N
    
#Looping the infected values over m iterations with given conditions
    for i in range(1,m): 
        I[i] = I[i-1] + h*beta*((S[i-1]*I[i-1])/N) - h*alpha*I[i-1] #Recursion formula
        if I[i] > N: 
            I = N
        if I[i] < 0:
            I = 0
#Dividing every element with number of people
        I_fraction[i] = I[i]/N
#Looping the susceptible values over m iterations with given conitions
        S[i] = S[i-1] - h*beta*((S[i-1]*I[i-1])/N) + h*alpha*I[i-1]
        if S[i] > N:
            S = N
        if S[i] < 0:
            S = 0
#Time steps        
        H[i] = h*i
#Plotting        
    plt.plot(H,I, 'r', label = 'Infected population')
    plt.plot(H,S, 'b', label = 'Susceptible population')
    plt.xlabel(r'Time [$t$]')
    plt.ylabel(r'Population [$N$]')
    plt.legend()
    
    plt.figure()
    plt.plot(H,I_fraction, 'b', label = 'Fraction infected people')
    plt.axhline(y=1-alpha/beta, color = 'r', label = r'$1 - \alpha/ \beta$' )
    plt.xlabel(r'Time [$t$]')
    plt.ylabel(r'Infection fraction [$I/N$]')
    plt.legend()
    plt.show()
    
    
    return " "
#Changebale input values
infection = 100
stepsize = 0.01
population = 40000
infectionrate = 1.35
recovery_time = 1/2

#X = SIS(infection, population, infectionrate, stepsize)
#print (X)


time = 40 #Years  

def SIR(w1,w2,w3,w4,infection,N,beta,gamma,h):
    
    n = len(N) #Number of input cities.
    m = round(time/h) 
#Creating mxn matrices. Rows = city and columns = time steps
    I = np.zeros((n,m))
    S = np.zeros((n,m))
    R = np.zeros((n,m))
    I_c = np.zeros((n)) #Chronically ill with no recovery. 
    H = np.zeros(m) #Time step matrix. 
    
    w = np.row_stack((w1,w2,w3,w4)) #Travel probabilty matrix. 
    
#initial conditions
    for i in range(n):  
        I[(i,0)] = infection[i]
        S[(i,0)] = N[i] - I[(i,0)] 

    for j in range(1,m): #Updates columns. 
        for i in range(n): #Updates rows. 
            T = np.zeros(2) #Number of travelers to one given city. First element is travelers which are infected, second are susceptible.
            for k in range(n):                 
                if i !=k:
                    T[0] += w[i,k]*I[k,j-1] - w[k,i]*I[i,j-1] #Infected people leaving/enering city. 
                    T[1] += w[i,k]*S[k,j-1] - w[k,i]*S[i,j-1] #Susceptible people leaving/enering city. 
#Conditions and recursion formulas     
            I[(i,j)] = I[(i,j-1)] + h*beta[i] * (S[(i,j-1)]*I[(i,j-1)])/(N[i]) - h*alpha*I[(i,j-1)] + h*T[0] + I_c[i] 
            if I[(i,j)] < 0: 
                I[(i,j)] = 0
            if I[(i,j)] > N[i]:
                I[(i,j)] = N[i]
                
            dI = (I[(i,j)]-I[(i,j-1)])
            if dI < 0:
                dI = 0
            I_c[i] += 0.01*dI #1% of infected become chronically ill. 

            S[(i,j)] = S[(i,j-1)] - h*beta[i] * (S[(i,j-1)]*I[(i,j-1)])/(N[i]) - h*gamma[i]*S[(i,j-1)] + h*T[1]
            if S[(i,j)] < 0:
                S[(i,j)] = 0
            if S[(i,j)] > N[i]:
                S[(i,j)] = N[i]
            
            R[(i,j)] = N[i] - I[(i,j)] - S[(i,j)]
            if R[(i,j)] < 0:
                R[(i,j)] = 0
            if R[(i,j)] > N[i]:
                R[(i,j)] = N[i]
            
            H[j] = j*h
            
#Code for plotting in the SIR model
    '''
    cities = ('Kalmar', 'Tanta', 'Las Vegas', 'Kabul') #City names
    for i in range(len(Population_N)):
        plt.plot(H,S[i],'b', label = 'Susceptible population') #Plots each row (=city) from the Susceptible matrix to time axis. 
        plt.plot(H,I[i],'r', label = 'Infected population')
        plt.plot(H,R[i],'g', label = 'Recovered population') #Comment out when running Error function.
        
        plt.title(cities[i])
        plt.xlabel(r'Time [$t$]')
        plt.ylabel(r'Infection fraction [$I/N$]')
        plt.legend()
        plt.show()'''
        
    return I,S,R  #Uncomment when running Error function. 
      

def Error(w1,w2,w3,w4,infection,N,beta,gamma,h,city):    
    A = SIR(w1,w2,w3,w4,infection,N,beta,gamma,h) #Tuple containing matrices I,S,R when running SIR with step size h. 
    B = SIR(w1,w2,w3,w4,infection,N,beta,gamma,2*h) #Tuple containing matrices I,S,R when running SIR with step size 2h. 
    
    m = round(time/(2*h))
    n = len(N)


    InfectionMatrix_h = A[0] #Takes the infection matrix from A. 
    InfectionMatrix_2h = B[0] 
    
    #II = np.array(InfectionMatrix_h) #Used to call certian elements in tuples. 
    #InfectionList_h = II[0]
    
    SusceptibilityMatrix_h = A[1]
    SusceptibilityMatrix_2h = B[1] #Takes the susceptible matrix. 
    
    RecoveryMatrix_h = A[2]
    RecoveryMatrix_2h = B[2] #Takes the recovery part. 

    H = np.zeros(m)
    L = linspace(0,time,2*m)
    
    E = np.zeros((n,m))
#Magnitude of errors    
    for j in range(m): 
        E[0,j] = abs(InfectionMatrix_2h[city,j] - InfectionMatrix_h[city,2*j])
        E[1,j] = abs(SusceptibilityMatrix_2h[city,j] - SusceptibilityMatrix_h[city,2*j])
        E[2,j] = abs(RecoveryMatrix_2h[city,j] - RecoveryMatrix_h[city,2*j])
        H[j] = h*2*j
    
    I = InfectionMatrix_h[city]
    S = SusceptibilityMatrix_h[city]
    R = RecoveryMatrix_h[city]
    
    plt.title('Kalmar')
    plt.plot(H,E[0],'grey', label = 'Magnitude of error') #Plots error between infection graphs for step sized h and 2h, for a given city. 
    plt.plot(L,InfectionMatrix_h[city],'r', label = 'Infected inhabitants') #plots infection graph for step size h, for a given city. 
    plt.ylabel(r' Population $[N]$')
    plt.xlabel(r' time $[t]$ ')
    plt.legend()
    plt.show()
    #plt.plot(H,InfectionMatrix_2h[city], 'r') #plots infection graph for step size h, for a given city.
    
    plt.title('Kalmar')
    plt.plot(H,E[1],'grey', label = 'Magnitude of error')
    plt.plot(L,SusceptibilityMatrix_h[city],'b', label = 'Susceptible inhabitants') 
    plt.ylabel(r' Population $[N]$')
    plt.xlabel(r' time $[t]$' )
    plt.legend()
    plt.show()
    
    #plt.plot(H,SusceptibilityMatrix_2h[city], 'r')
    plt.title('Kalmar')
    plt.plot(H,E[2],'grey', label = 'Magnitude of error')
    plt.plot(L,RecoveryMatrix_h[city],'g', label = 'Recovered inhabitants')
    plt.ylabel(r' Population $[N]$')
    plt.xlabel(r' time $[t]$ ')
    plt.legend()
    plt.show()
    #plt.plot(H,RecoveryMatrix_2h[city], 'r')
    
    plt.title('Kalmar')
    plt.plot(H,E[0]/(I[::2]),'r', label = 'Relative error for infected')
    plt.plot(H,E[1]/(S[::2]),'b', label = 'Relative error for susceptible')
    plt.plot(H,E[2]/(R[::2]),'g', label = 'Relative error for recovered')
    plt.ylabel(r'Relative global error $[\%]$')
    plt.xlabel(r' time $[t]$')
    plt.legend()
    plt.show()
    


    
    return #I[::2], I #" "
    
'Input for SIR and Error'
H = 0.2 #Timestep. 
InfectionPopulation = (100,0,0,0) #Number of infected people in each city at t = 0.
Population_N = (40000,400000,4000000,4000000)
Beta = (1.35,2,2,3)
Gamma = (0.001,0.001,0.01,0.0)

W1 = (0,0.003,0.0004,0.0004) #Probability of people travelling to first city from other cities. First element set to zero.  
W2 = (0.01,0,0.001,0.001) 
W3 = (0.05,0.009,0,0.005) 
W4 = (0.05,0.009,0.005,0)
Y = SIR(W1,W2,W3,W4,InfectionPopulation, Population_N, Beta, Gamma,H)
#print(Y)


City = 0 #Error estimation of city 1. 
Z = Error(W1,W2,W3,W4,InfectionPopulation, Population_N, Beta, Gamma,H, City)
print(Z)

