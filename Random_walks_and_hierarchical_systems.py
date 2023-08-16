#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 20 21:37:15 2019

@author: Markus Strickert
"""


import numpy as np
import pylab 
import random 
import matplotlib.pyplot as plt

'''Code is divided into different functions.

OneDdist, twoDdist
and threeDdist returns distribution of first passage time, in
different dimensions.

FoodConc1 return equillibrium population, given a 
food concentration/food amount (including reproduction)
PlotFoodConc simply plots the returned values from FoodConc1, 
manually.
Foodconc2 returns population survival time, given a food 
concentration/ food amount (excluding reproduction).
PlotFoodConc2 plots values.

NaturalSelection1 simulates two herbivore species,
with different traits. Returnes time evolution of populations.
NaturalSelection2 includes three species, two herbivore, and
one carnivore.

The FoodConc functions and NaturalSelection functions have similiar
structure. Read the first ones of each, most comments are there.
The functions are simply slightly more developed than the ones before,
including different things, so the code is quite long.
Have fun and Good luck ;)

To run functions, uncoment code underneath them. 

'''


'''First passage, 1D'''
'''Function runs one dimensional random walk. When walker
 on target site, breaks. Simulation run m number of times. 
 Returns distribution.'''
def oneDdist(m):
    distlist = []
    for j in range(m):
        walk = 0
        for i in range(5000):
            if walk == 10:
                break
        
            step = random.choice(['W','E'])
            if step == 'W':
                walk -= 1
            if step == 'E':
                walk += 1
        if i != 4999:
            distlist.append(i)
    plt.hist(distlist, bins = 30)

#oneD = oneDdist(1000)
#print(oneD)

'''First passage, 2D'''
def twoDdist(m):
    distlist=[]
    l = 10
    
    for j in range(m):
        xwalk = 0
        ywalk = 0
        for i in range(700):
            if xwalk and ywalk == 7:
                break
            xstep = random.choice(['W','E'])
            ystep = random.choice(['N','S'])
            if xstep == 'W':
                if xwalk - 1 < -l:
                    xwalk += 1
                else: 
                    xwalk -= 1
            if xstep == 'E':
                if xwalk + 1 > l:
                    xwalk -= 1
                else: 
                    xwalk += 1    
            if ystep == 'S':
                if ywalk - 1 < -l:
                    ywalk += 1
                else: 
                    ywalk -= 1
            if ystep == 'N':
                if ywalk + 1 > l:
                    ywalk -= 1
                else: 
                    ywalk += 1
        if i != 699:
            distlist.append(i)
            
    plt.show()   
    plt.hist(distlist, bins = 60)

twoD = twoDdist(10000)
print(twoD)

'''First passage, 3D'''
def threeDdist(m):
    distlist = []
    l = 5
    
    for j in range(m):
        xwalk = 0
        ywalk = 0
        zwalk = 0
        for i in range(400):
            if xwalk and ywalk and zwalk == 4:
                break
            xstep = random.choice(['W','E'])
            ystep = random.choice(['N','S'])
            zstep = random.choice(['I','O'])
            if xstep == 'W':
                if xwalk - 1 < -l:
                    xwalk += 1
                else: 
                    xwalk -= 1
            if xstep == 'E':
                if xwalk + 1 > l:
                    xwalk -= 1
                else: 
                    xwalk += 1    
            if ystep == 'S':
                if ywalk - 1 < -l:
                    ywalk += 1
                else: 
                    ywalk -= 1
            if ystep == 'N':
                if ywalk + 1 > l:
                    ywalk -= 1
                else: 
                    ywalk += 1
            if zstep == 'O':
                if zwalk - 1 < -l:
                    zwalk += 1
                else: 
                    zwalk -= 1
            if zstep == 'I':
                if zwalk + 1 > l:
                    zwalk -= 1
                else: 
                    zwalk += 1
        if i != 399:
            distlist.append(i)
            
    plt.show()   
    plt.hist(distlist, bins = 60)

#threeD = threeDdist(50000)
#print(threeD)

'''Function gives a mean equillibrium population, given a food concentration.
Choose appropriate number of step iterations, such that the population can stabilize.
And apporopriate simulationIterations, such that simulation is run sufficient number
of times such that the mean converges '''
def FoodConc1(stepIterations, simulationIterations, FoodConcentration):
    meanlist = [] #will append equillibrium populations in this list
    walkers = 25 #number of walkers
    for k in range(simulationIterations):
        N = walkers 
        n = stepIterations
        totalwalkers = N #dead or alive
        alivewalkers = N #alive
        alivewalkerslist = [N]
        playerlist = np.ones(N).tolist() #Consists of either 1 or 0. 1 
        #are alive walkers. 0 are dead walkers. 
        
        l = 100 # lxl 2D plane
        
        foodproducerate = 1 #food given after every -- iterations. 
        foodamount = FoodConcentration
        
        xpos = np.zeros((N,n)) #x position of players.
        #Each row corresponds to one player. Every columnt
        #to a iteration.
        ypos = np.zeros((N,n)) # y position
        
        #Initial positions of walkers, spread over plane.
        sqrtN = int(np.sqrt(N))
        xposinitial = np.linspace(-l,l,sqrtN).tolist() * sqrtN
        yposinitial = []
        for y in range(sqrtN):
            y = -l + (y+0.5)*(2*l)/(sqrtN)
            yposinitial.extend([y for i in range(sqrtN)])
        ypos[:,0] = yposinitial
        xpos[:,0] = xposinitial
        
        hungermatrix = np.ones((N,n)) #hungermatrix conisists of 0,1 or 2.
        #0 if walker is dead. 1 if walker alive but has not eaten in
        #that iteration. 2 if walker  recently ate. 
        foodlistx = np.zeros(foodamount).tolist()
        foodlisty = np.zeros(foodamount).tolist()
        
        
        m = l/10 #speed of walker
        walkerradar = 8 #walkers food detection radar. 
        survivaltime = 25 
        
        ilist = [0]
        emptylist = np.zeros(n) #just some lists
        onelist = np.ones(n)
        
        #iterating over iterations
        for i in range(1,n):
            ilist.append(i)
            
            '''Producing food'''
            if i % foodproducerate == 0:
                for k in range(foodamount):
                    food = (random.randint(-l, l), random.randint(-l, l))
                    foodlistx[k] = food[0]
                    foodlisty[k] = food[1]
                
            for j in range(totalwalkers):
     
                '''deathrow''' #if walker has not eaten
                #in a period equal to its survival time, it dies. 
                deathrow = 0
                if playerlist[j] != 0:
                    for q in range(survivaltime):
                        if q <= i:
                            deathrow += hungermatrix[j,i-q]
                            
                    if int(deathrow) == min(q, survivaltime):
                        alivewalkers -= 1
                        playerlist[j] = 0
                        hungermatrix[j] = 0
                        xpos[j] = 0
                        ypos[j] = 0            
                
                '''simulating random walks'''
                if playerlist[j] != 0:
                    xval = random.randint(0,100*m) 
                    #Player takes somewhere between 0 and m steps
                    xval = xval/100
                    yval = np.sqrt(m**2 - xval**2)
                    xsign = random.randint(1, 2) 
                    #Randomly choosing in which direction
                    #player takes step.
                    ysign = random.randint(1, 2)
                    '''direction of walk'''
                    if xsign == 1:
                        xval = -xval
                    if ysign  == 1:
                        yval = -yval
                    '''reflective boundaries'''
                    if abs(xpos[j,i-1] + xval) > l:
                        xpos[j,i] = xpos[j,i-1] - xval  
                    else:
                        xpos[j,i] = xpos[j,i-1] + xval 
                    if abs(ypos[j,i-1] + yval) > l:
                        ypos[j,i] = ypos[j,i-1] - yval  
                    else:
                        ypos[j,i] = ypos[j,i-1] + yval 
                        
                
    
                '''food and multiplication of walkers'''
                for k in range(foodamount):
                    if (foodlistx[k] and foodlisty[k] and playerlist[j] != 0):
                        distance = np.sqrt((xpos[j,i]-foodlistx[k])**2 + (ypos[j,i]-foodlisty[k])**2)
                        if distance < walkerradar:
                            totalwalkers += 1
                            alivewalkers += 1
                            playerlist.append(1)
    
                            hungermatrix = np.row_stack((hungermatrix,emptylist))
                            hungermatrix[j,i] = 2
                            foodlistx[k] = 0
                            foodlisty[k] = 0
                            
                            xpos = np.row_stack((xpos,emptylist))
                            ypos = np.row_stack((ypos,emptylist))
                            dim = len(xpos)
                            xpos[dim-1,i] = xpos[j,i]
                            ypos[dim-1,i] = ypos[j,i]                        
                        else:
                            hungermatrix[j,i] = 1
                
        meanlist.append(alivewalkers)
    mean = sum(meanlist)/(len(meanlist))

    return mean

#Y = FoodConc(400, 100, 8)
#print(Y)

def PlotFoodConc():
    foodlist = [6,7,8,9,10,11]
    walkerlist = [0.34,0.79,1.27,6.65,38.34,81.83]
    plt.plot(foodlist,walkerlist, 'yo')
    plt.ylabel('Equillibrium population [N]')
    plt.xlabel('Food amount')
    plt.legend()

#Y = PlotFoodConc()
#print(Y)

def FoodConc2(stepIterations, simulationIterations, FoodConcentration):
    meanlist=[]
    walkers = 25
    for k in range(simulationIterations):
        N = walkers
        totalwalkers = N #dead or alive
        alivewalkers = N
        alivewalkerslist = [N]
        playerlist = np.ones(N).tolist()
        
        l = 33
        n = stepIterations
        
        foodproducerate = 1 #food given after every -- iterations. 
        foodamount = FoodConcentration #ammount of food given 
        m = l/10
        walkerradar = 10
        survivaltime = 10
        xpos = np.zeros((N,n)) #x position of players.
        #Each row corresponds to one player. 
        ypos = np.zeros((N,n)) # y position
        
        sqrtN = int(np.sqrt(N))
        xposinitial = np.linspace(-l,l,sqrtN).tolist() * sqrtN
        yposinitial = []
        for y in range(sqrtN):
            y = -l + (y+0.5)*(2*l)/(sqrtN)
            yposinitial.extend([y for i in range(sqrtN)])
        ypos[:,0] = yposinitial
        xpos[:,0] = xposinitial
        
        hungermatrix = np.ones((N,n))
        foodlistx = np.zeros(foodamount).tolist()
        foodlisty = np.zeros(foodamount).tolist()
        
        ilist = [0]
        emptylist = np.zeros(n) #just some lists
        onelist = np.ones(n)
        
        for i in range(1,n):
            ilist.append(i)
            
            if alivewalkers == 0:
                meanlist.append(i)
                break
            
            '''Producing food'''
            if i % foodproducerate == 0:
                for k in range(foodamount):
                    food = (random.randint(-l, l), random.randint(-l, l))
                    foodlistx[k] = food[0]
                    foodlisty[k] = food[1]
                    
            for j in range(totalwalkers):
                '''deathrow'''
                deathrow = 0
                if playerlist[j] != 0:
                    if i > survivaltime:
                        for q in range(survivaltime):
                            deathrow += hungermatrix[j,i-q]
                            
                        if int(deathrow) == survivaltime:
                            alivewalkers -= 1
                            playerlist[j] = 0
                            hungermatrix[j] = 0
                            xpos[j] = 0
                            ypos[j] = 0            
                    
                '''simulating random walks'''
                if playerlist[j] != 0:
                    xval = random.randint(0,100*m) 
                    #Player takes somewhere between 0 and m steps
                    xval = xval/100
                    yval = np.sqrt(m**2 - xval**2)
                    xsign = random.randint(1, 2) 
                    #Randomly choosing in which direction player takes step.
                    ysign = random.randint(1, 2)
                    '''direction of walk'''
                    if xsign == 1:
                        xval = -xval
                    if ysign  == 1:
                        yval = -yval
                    '''reflective boundaries'''
                    if abs(xpos[j,i-1] + xval) > l:
                        xpos[j,i] = xpos[j,i-1] - xval  
                    else:
                        xpos[j,i] = xpos[j,i-1] + xval 
                    if abs(ypos[j,i-1] + yval) > l:
                        ypos[j,i] = ypos[j,i-1] - yval  
                    else:
                        ypos[j,i] = ypos[j,i-1] + yval 
                            
                    
                if playerlist[j] != 0:
                    '''food and multiplication of walkers'''
                    a = 0 
                    for k in range(foodamount):
                        if (foodlistx[k] and foodlisty[k] != 0):
                            distance = np.sqrt((xpos[j,i]-foodlistx[k])**2 + (ypos[j,i]-foodlisty[k])**2)
                            if distance < walkerradar:
                                a += 1
                                hungermatrix[j,i] = 2
                                foodlistx[k] = 0
                                foodlisty[k] = 0   
    
                    if a == 0: 
                        hungermatrix[j,i] = 1
                    
    populationSurvival = sum(meanlist)/(len(meanlist))
 
    return populationSurvival, len(meanlist)

#Foodconc = FoodConc2(4000,300,11)
#print(Foodconc)

def plotFoodConc2():
    foodconc = [1,2,3,4,5,6,7,8,9,10]
    PopulationSurvivalTime = [29,46.9,69.4,105.4,155.5,223.9,334,526.1,737.7,1119]
    plt.xlabel('Food amount')
    plt.ylabel('Population survival time')
    plt.plot(foodconc,PopulationSurvivalTime, 'go')
    plt.legend()
    
#Foodconc2 = plotFoodConc2()
#print(Foodconc2)


'''One species - 2D plane - reflective boundaries. 1 species. Walker survives, 
but does not multiplicate when food is found. '''
'''Not mentioned in the report. Runs simulation for FoodConc 2.
NaturalSelection 1 and 2 part of the report. '''
def NaturalSelection(walkers, planedim, iterations):

    N = walkers
    totalwalkers = N #numver of dead or alive
    alivewalkers = N #numver of alive walkers
    alivewalkerslist = [N] #Appends values of alivewalkers to this list.
    playerlist = np.ones(N).tolist()
    
    l = planedim #dimension of plane
    n = iterations
    
    foodproducerate = 1 #food given after every -- iterations. 
    foodamount = 6 #ammount of food given 
    m = l/10
    walkerradar = 10
    survivaltime = 10
    xpos = np.zeros((N,n)) #x position of players. Each row corresponds to one player. 
    ypos = np.zeros((N,n)) # y position
    
    #Initial positions of walkers.
    sqrtN = int(np.sqrt(N))
    xposinitial = np.linspace(-l,l,sqrtN).tolist() * sqrtN
    yposinitial = []
    for y in range(sqrtN):
        y = -l + (y+0.5)*(2*l)/(sqrtN)
        yposinitial.extend([y for i in range(sqrtN)])
    ypos[:,0] = yposinitial
    xpos[:,0] = xposinitial
    
    hungermatrix = np.ones((N,n)) 
    foodlistx = np.zeros(foodamount).tolist() #xpos of food.
    foodlisty = np.zeros(foodamount).tolist() #ypos of food.
    
    ilist = [0]
    emptylist = np.zeros(n) #just some lists
    onelist = np.ones(n)
    
    for i in range(1,n):
        ilist.append(i)
        
        '''Producing food'''
        if i % foodproducerate == 0:
            for k in range(foodamount):
                food = (random.randint(-l, l), random.randint(-l, l))
                foodlistx[k] = food[0]
                foodlisty[k] = food[1]
                
        for j in range(totalwalkers):
            '''deathrow'''
            deathrow = 0
            if playerlist[j] != 0:
                if i > survivaltime:
                    for q in range(survivaltime):
                        deathrow += hungermatrix[j,i-q]
                        
                    if int(deathrow) == survivaltime:
                        alivewalkers -= 1
                        playerlist[j] = 0
                        hungermatrix[j] = 0
                        xpos[j] = 0
                        ypos[j] = 0            
                
            '''simulating random walks'''
            if playerlist[j] != 0:
                xval = random.randint(0,100*m) #Player takes somewhere between 0 and m steps
                xval = xval/100
                yval = np.sqrt(m**2 - xval**2)
                xsign = random.randint(1, 2) #Randomly choosing in which direction player takes step.
                ysign = random.randint(1, 2)
                '''direction of walk'''
                if xsign == 1:
                    xval = -xval
                if ysign  == 1:
                    yval = -yval
                '''reflective boundaries'''
                if abs(xpos[j,i-1] + xval) > l:
                    xpos[j,i] = xpos[j,i-1] - xval  
                else:
                    xpos[j,i] = xpos[j,i-1] + xval 
                if abs(ypos[j,i-1] + yval) > l:
                    ypos[j,i] = ypos[j,i-1] - yval  
                else:
                    ypos[j,i] = ypos[j,i-1] + yval 
                        
                
            if playerlist[j] != 0:
                '''food and multiplication of walkers'''
                a = 0
                for k in range(foodamount):
                    if (foodlistx[k] and foodlisty[k] != 0):
                        distance = np.sqrt((xpos[j,i]-foodlistx[k])**2 + (ypos[j,i]-foodlisty[k])**2)
                        if distance < walkerradar:
                            a += 1
                            hungermatrix[j,i] = 2
                            foodlistx[k] = 0
                            foodlisty[k] = 0   

                if a == 0: 
                    hungermatrix[j,i] = 1
                
        alivewalkerslist.append(alivewalkers)
        
    plt.plot(ilist,alivewalkerslist, label = 'Total walkers')
    plt.ylim(0)
    plt.legend()
    #plt.plot(xposinitial,yposinitial, 'bo')
    plt.show()
 
    return alivewalkers
    
#Naturalselection = NaturalSelection(25,33,500)
#print(Naturalselection)
    
    
def NaturalSelection1(walkers, planedim, iterations):
    N = walkers #alive walkers
    totalwalkers = N #dead or alive
    alivewalkers = N #Total numver of walkers
    alivespecies1 = N-25 #number of alive species 1
    alivespecies2 = 25
    alivewalkerslist = [alivewalkers] #Number of alive walkers will be
    #appended here.
    alivespecies1list = [alivespecies1] #Number of alive walkers of species 1
    #will be appended here. 
    alivespecies2list = [alivespecies2] 
    
    species1list = [1] * (alivespecies1)
    species2list = [2] * (alivespecies2)     #Not important
    specieslist = species1list + species2list #Equal to 'playerlist' used
    #before. Consists of either 0,1 or 2. 0 dead walkers. 1 alive species 1,
    #2 alive species 2. 
    
    l = planedim #dimensions of plane on which walkers live on
    n = iterations
    
    emptylist = np.zeros(n) #just some lists
    onelist = np.ones(n)

    
    foodproducerate = 1 #food given after every -- iterations. 
    foodamount = 6 #ammount of food given 
    
    
    #survivaltime = 5 #number of iterations walkers survive without food.
    xpos = np.zeros((N,n)) 
    #x position of players. Each row corresponds to one player. 
    ypos = np.zeros((N,n)) # y position
    
    sqrtN = int(np.sqrt(N))
    xposinitial = np.linspace(-l,l,sqrtN).tolist() * sqrtN
    yposinitial = []
    for y in range(sqrtN):
        y = -l + (y+0.5)*(2*l)/(sqrtN)
        yposinitial.extend([y for i in range(sqrtN)])
    ypos[:,0] = yposinitial
    xpos[:,0] = xposinitial
    
    hungermatrix = np.ones((N,n))
    foodlistx = np.zeros(foodamount).tolist()
    foodlisty = np.zeros(foodamount).tolist()
    
    ilist = [0]

    for i in range(1,n):
        ilist.append(i)
        #if alivewalkers > 52:
            #foodamount = 0
            #print('Its ZERO!')
        #print(foodamount)
        '''Producing food'''
        if i % foodproducerate == 0:
            for k in range(foodamount):
                food = (random.randint(-l, l), random.randint(-l, l))
                foodlistx[k] = food[0]
                foodlisty[k] = food[1]

        for j in range(totalwalkers):
            '''species'''
            if specieslist[j] != 0:
                if specieslist[j] == 1:
                    m = l/10
                    walkerradar = 5
                    species = 1
                    survivaltime = 21
                if specieslist[j] == 2:
                    m = (l)/10
                    walkerradar = 11
                    species = 2
                    survivaltime = 10
            
            '''deathrow'''
            deathrow = 0
            if specieslist[j] !=0 :
                if i > survivaltime:
                    for q in range(survivaltime):
                        deathrow += hungermatrix[j,i-q]
                        
                    if int(deathrow) == (survivaltime):
                        #if i > 800:
                            #print(int(deathrow))
                            #print(survivaltime)
                            #print(species)
                        alivewalkers -= 1
                        if species == 1:
                            alivespecies1 -= 1
                        if species == 2:
                            alivespecies2 -= 1
                        specieslist[j] = 0
                        hungermatrix[j] = 0
                        xpos[j] = 0
                        ypos[j] = 0
                        
            
            
            '''simulating random walks'''
            if specieslist[j] != 0:
                xval = random.randint(0,100*m) #Player takes somewhere between 0 and m steps
                xval = xval/100
                yval = np.sqrt(m**2 - xval**2)
                xsign = random.randint(1, 2) #Randomly choosing in which direction player takes step.
                ysign = random.randint(1, 2)
                '''direction of walk'''
                if xsign == 1:
                    xval = -xval
                if ysign  == 1:
                    yval = -yval
                '''reflective boundaries'''
                if abs(xpos[j,i-1] + xval) > l:
                    xpos[j,i] = xpos[j,i-1] - xval  
                else:
                    xpos[j,i] = xpos[j,i-1] + xval 
                if abs(ypos[j,i-1] + yval) > l:
                    ypos[j,i] = ypos[j,i-1] - yval  
                else:
                    ypos[j,i] = ypos[j,i-1] + yval 
                    
            

            '''food and multiplication of walkers'''
            if specieslist[j] != 0:
                a = 0
                for k in range(foodamount):
                    if foodlistx[k] and foodlisty[k] != 0:
                        distance = np.sqrt((xpos[j,i]-foodlistx[k])**2 + (ypos[j,i]-foodlisty[k])**2)
                        if distance < walkerradar:
                            a += 1
                            totalwalkers += 1
                            alivewalkers += 1
                            if species == 1:
                                alivespecies1 += 1
                                specieslist.append(1)
                            if species == 2:
                                alivespecies2 += 1
                                specieslist.append(2)
                            hungermatrix = np.row_stack((hungermatrix,emptylist))
                            hungermatrix[j,i] = 2
                            foodlistx[k] = 0
                            foodlisty[k] = 0
                            
                            xpos = np.row_stack((xpos,emptylist))
                            ypos = np.row_stack((ypos,emptylist))
                            dim = len(xpos)
                            xpos[dim-1,i] = xpos[j,i]
                            ypos[dim-1,i] = ypos[j,i]    
                            hungermatrix[dim-1,i] = 2
                    else:
                        hungermatrix[j,i] = 1 #Remove
                    
                if a == 0:
                    hungermatrix[j,i] = 1
            

        alivewalkerslist.append(alivewalkers)
        alivespecies1list.append(alivespecies1)
        alivespecies2list.append(alivespecies2)
        
        fraction = 0.1 * n
        if i % fraction == 0:
            print('{} procent done, there are {} alive walkers'.format(((i/n)*100),alivewalkers))
    
    #plt.plot(ilist,alivewalkerslist, label = 'Total walkers')
    plt.plot(ilist,alivespecies1list,'r', label = 'species1')
    plt.plot(ilist,alivespecies2list,'seagreen', label = 'species2')
    plt.ylim(0)
    plt.legend()
    plt.xlabel(r'time [days] ')
    plt.ylabel(r'alive walkers')

    #plt.plot(xposinitial,yposinitial, 'bo')
    plt.show()
    
    #uncomment part below, if you want to see walkers themselves.
# =============================================================================
#     plt.show()
#     pylab.title("Random Walk ($n = " + str(n) + "$ steps)") #plotting everythin.
#     for i in range (totalwalkers): #plottar alla walkers
#         if xpos[i,n-1] and ypos[i,n-1] != 0: #plottar inte punker om walkern,
#             #i sista iteration hamnar på punkten (0,0), då jag satt positionen av
#             #alla döda walkers till (0,0) vid den tidpunkten
#             if specieslist[i] == 1: #plottar första arted med blåa punkter
#                 pylab.plot(xpos[i,n-10:], ypos[i, n-10:])  #plottar sista 10 positionerna av walkern.
#                 plt.plot(xpos[i,n-1],ypos[i,n-1], 'bo') #plottar sista positionen med cirkel,
#                 #så man ser vart walkern hamnar sist. 
#             if specieslist[i] == 2: #plottar andra arted med röda
#                 pylab.plot(xpos[i,n-10:], ypos[i, n-10:]) 
#                 plt.plot(xpos[i,n-1],ypos[i,n-1], 'ro')
#     plt.plot(foodlistx,foodlisty, 'gx') #plottar mat
#     
#     plt.xlim((-l-0.1*l,l+0.1*l))
#     plt.ylim((-l-0.1*l,l+0.1*l))
#     pylab.savefig("rand_walk"+str(n)+".png",bbox_inches="tight",dpi=600) 
#     pylab.show() 
# =============================================================================
    
    return alivewalkers, alivespecies1, alivespecies2
    

#X = NaturalSelection1(49,35,400) #do 20 or 40 planedim. Blue in 20
#print(X)
    

def NaturalSelection2(walkers, planedim, iterations):
    N = walkers #alive walkers
    totalwalkers = N #dead or alive
    alivewalkers = N
    alivespecies1 = 400
    alivespecies2 = 400
    alivespecies3 = 100
    alivewalkerslist = [alivewalkers]
    alivespecies1list = [alivespecies1] #See code above
    alivespecies2list = [alivespecies2]
    alivespecies3list = [alivespecies3]
    
    species1list = [1] * (alivespecies1)
    species2list = [2] * (alivespecies2)
    species3list = [3] * (alivespecies3)    
    specieslist = species1list + species2list + species3list #See code above
    random.shuffle(specieslist)
    
    l = planedim #dimensions of plane on which walkers live on
    n = iterations

    emptylist = np.zeros(n) #just some lists
    onelist = np.ones(n)
    
    foodproducerate = 1 #food given after every -- iterations. 
    foodamount = 52 #ammount of food given 
    #survivaltime = 5 #number of iterations walkers survive without food.
    xpos = np.zeros((N,n)) #x position of players. Each row corresponds to one player. 
    ypos = np.zeros((N,n)) # y position
    
    #Initial positions
    sqrtN = int(np.sqrt(N))
    xposinitial = np.linspace(-l,l,sqrtN).tolist() * sqrtN
    yposinitial = []
    for y in range(sqrtN):
        y = -l + (y+0.5)*(2*l)/(sqrtN)
        yposinitial.extend([y for i in range(sqrtN)])
    ypos[:,0] = yposinitial
    xpos[:,0] = xposinitial
    
    hungermatrix = np.ones((N,n))
    foodlistx = np.zeros(foodamount).tolist()
    foodlisty = np.zeros(foodamount).tolist()
    
    ilist = [0]
    
    for i in range(1,n):
        ilist.append(i)
        
        '''Producing food'''
        if i % foodproducerate == 0:
            for k in range(foodamount):
                food = (random.randint(-l, l), random.randint(-l, l))
                foodlistx[k] = food[0]
                foodlisty[k] = food[1]
            
        for j in range(totalwalkers):
            '''species'''
            if specieslist[j] != 0:
                if specieslist[j] == 1:
                    m = l/10
                    walkerradar = 3
                    species = 1
                    survivaltime = 15
                if specieslist[j] == 2:
                    m = (5*l)/10
                    walkerradar = 8
                    species = 2
                    survivaltime = 7
                if specieslist[j] == 3:
                    m = (3*l)/10
                    walkerradar = 3.7
                    species = 3
                    survivaltime = 35
                    attackradar = 2.5
            
                '''deathrow'''
                deathrow = 0
                if specieslist[j] != 0:
                    if i > survivaltime:
                        for q in range(survivaltime):
                            deathrow += hungermatrix[j,i-q]
                                    
                        if int(deathrow) == survivaltime:
                            alivewalkers -= 1
                            if species == 1:
                                alivespecies1 -= 1
                            if species == 2:
                                alivespecies2 -= 1
                            if species == 3:
                                alivespecies3 -= 1
                                #print('Hunger')
                            specieslist[j] = 0
                            hungermatrix[j] = 0
                            xpos[j] = 0
                            ypos[j] = 0            
                
                
                '''simulating random walks'''
                xval = random.randint(0,100*m) #Player takes somewhere between 0 and m steps
                xval = xval/100
                yval = np.sqrt(m**2 - xval**2)
                xsign = random.randint(1, 2) #Randomly choosing in which direction player takes step.
                ysign = random.randint(1, 2)
                '''direction of walk'''
                if xsign == 1:
                    xval = -xval
                if ysign  == 1:
                    yval = -yval
                '''reflective boundaries'''
                if abs(xpos[j,i-1] + xval) > l:
                    xpos[j,i] = xpos[j,i-1] - xval  
                else:
                    xpos[j,i] = xpos[j,i-1] + xval 
                if abs(ypos[j,i-1] + yval) > l:
                    ypos[j,i] = ypos[j,i-1] - yval  
                else:
                    ypos[j,i] = ypos[j,i-1] + yval 
                        
                
    
                '''food and multiplication of walkers'''
                '''if species 1 or 2 find food'''
                if species == 1 or species == 2:
                    a = 0
                    for k in range(foodamount):
                        if foodlistx[k] and foodlisty[k] != 0:
                            distance = np.sqrt((xpos[j,i]-foodlistx[k])**2 + (ypos[j,i]-foodlisty[k])**2)
                            if distance < walkerradar:
                                a += 1
                                totalwalkers += 1
                                alivewalkers += 1
                                if species == 1:
                                    alivespecies1 += 1
                                    specieslist.append(1)
                                if species == 2:
                                    alivespecies2 += 1
                                    specieslist.append(2)
                                hungermatrix = np.row_stack((hungermatrix,emptylist)) #why empty
                                hungermatrix[j,i] = 2
                                foodlistx[k] = 0
                                foodlisty[k] = 0
                                
                                xpos = np.row_stack((xpos,emptylist))
                                ypos = np.row_stack((ypos,emptylist))
                                dim = len(xpos)
                                xpos[dim-1,i] = xpos[j,i]
                                ypos[dim-1,i] = ypos[j,i]    
                                hungermatrix[dim-1,i] = 2
                            else:
                                hungermatrix[j,i] = 1
                    if a == 0:
                        hungermatrix[j,i] = 1
                
                '''if species 3 finds species 1 or 2 as food'''
                if i > 5:
                    if species == 3:
                        a = 0
                        b = 0
                        for k in range(totalwalkers):
                            if specieslist[k] == 1 or specieslist[k] == 2:
                                distance = np.sqrt((xpos[j,i]-xpos[k,i])**2 + (ypos[j,i]-ypos[k,i])**2)
                                if distance < walkerradar:
                                    a += 1
                                    hungermatrix[j,i] += 1
                                    
                                    if hungermatrix[j,i] == 4:
                                        specieslist.append(3)
                                        alivespecies3 += 1
                                        alivewalkers += 1
                                        totalwalkers += 1
                                        hungermatrix = np.row_stack((hungermatrix,emptylist)) #why empty
                                        hungermatrix[j,i] = 2
                                        
                                        xpos = np.row_stack((xpos,emptylist))
                                        ypos = np.row_stack((ypos,emptylist))
                                        dim = len(xpos)
                                        xpos[dim-1,i] = xpos[j,i]
                                        ypos[dim-1,i] = ypos[j,i]
                                        hungermatrix[dim-1,i] = 2
                                    
                                    if specieslist[k] == 1:
                                        specieslist[k] = 0
                                        alivespecies1 -= 1
                                        alivewalkers -= 1
        
                                        hungermatrix[k] = 0
                                        xpos[k] = 0
                                        ypos[k] = 0   
                                    
                                    
                                    if specieslist[k] == 2:
                                        specieslist[k] = 0
                                        alivespecies2 -= 1
                                        alivewalkers -= 1
        
                                        hungermatrix[k] = 0
                                        xpos[k] = 0
                                        ypos[k] = 0  
                                    
                        
                        if a == 0:
                            hungermatrix[j,i] = 1
                        
                        '''if two walkers of species 3 collide, one dies'''
                        for k in range(totalwalkers):
                            if specieslist[k] == 3:
                                if j > k:
                                    distance = np.sqrt((xpos[j,i]-xpos[k,i])**2 + (ypos[j,i]-ypos[k,i])**2)
                                    if distance < attackradar*walkerradar:
                                        b += 1
                                        specieslist[k] = 0
                                        alivespecies3 -= 1
                                        alivewalkers -= 1
                                        xpos[k] = 0
                                        ypos[k] = 0
                                        hungermatrix[k] = 0
                                        #print('True')
                        #Uncomment if you want both species 3 walkers to die.                
                        #if b > 0:
                            #specieslist[j] = 0
                            #alivewalkers -= 1
                            #alivespecies3 -= 1
                            #xpos[j] = 0
                            #ypos[j] = 0
                            #hungermatrix[j] = 0
                            #print('Fight!!')
                            
                        
                        if b == 0:
                            hungermatrix[j,i] = 1
                
    
        alivewalkerslist.append(alivewalkers)
        alivespecies1list.append(alivespecies1)
        alivespecies2list.append(alivespecies2)
        alivespecies3list.append(alivespecies3)
        
        
    
    
    #plt.plot(ilist,alivewalkerslist, label = 'Total walkers')
    plt.plot(ilist,alivespecies1list, label = 'species1')
    plt.plot(ilist,alivespecies2list, label = 'species2')
    '''1000 lines!! Hurray!! Sorry :/ '''
    plt.plot(ilist,alivespecies3list, label = 'species3')
    plt.ylabel('alive walkers')
    plt.xlabel('time [days]')
    plt.ylim(0)
    plt.legend()
    #plt.plot(xposinitial,yposinitial, 'bo')
    plt.show()

    #Uncomment if you want to see walkers themselves.
# =============================================================================
#     plt.show()
#     pylab.title("Random Walk ($n = " + str(n) + "$ steps)") #plotting everythin.
#     for i in range (totalwalkers): #plottar alla walkers
#         if xpos[i,n-1] and ypos[i,n-1] != 0: #plottar inte punker om walkern,
#             #i sista iteration hamnar på punkten (0,0), då jag satt positionen av
#             #alla döda walkers till (0,0) vid den tidpunkten
#             if specieslist[i] == 1: #plottar första arted med blåa punkter
#                 pylab.plot(xpos[i,n-10:], ypos[i, n-10:])  #plottar sista 10 positionerna av walkern.
#                 plt.plot(xpos[i,n-1],ypos[i,n-1], 'bo') #plottar sista positionen med cirkel,
#                 #så man ser vart walkern hamnar sist. 
#             if specieslist[i] == 2: #plottar andra arted med röda
#                 pylab.plot(xpos[i,n-10:], ypos[i, n-10:]) 
#                 plt.plot(xpos[i,n-1],ypos[i,n-1], 'ro')
#     plt.plot(foodlistx,foodlisty, 'gx') #plottar mat
#     
#     plt.xlim((-l-0.1*l,l+0.1*l))
#     plt.ylim((-l-0.1*l,l+0.1*l))
#     pylab.savefig("rand_walk"+str(n)+".png",bbox_inches="tight",dpi=600) 
#     pylab.show() 
# =============================================================================
    
    
    
    return alivewalkers, alivespecies1, alivespecies2, alivespecies3

#Z = NaturalSelection2(900,48,100)
#print(Z)
