0# -*- coding: utf-8 -*-
"""
Created on Wed Nov 09 14:08:28 2016

@author: Markus
"""

from matplotlib import cm
from itertools import cycle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#functions

def velo(rho):
    global rhomax, vmax
    if rho<rhomax:
        return vmax*(1-rho/rhomax)
    return 0
    
def f(rho):
    global vmax, rhomax
    return rho*velo(rho)
    
def init(x,rho,dx):
    """defines the inital values for given n-dim vector x and n-1-dim vector rho"""
    if len(x)<2 or len(rho)!=len(x)-1:
        print "check dimensions of rho and x"
        return None
    xstart = x[0]
    xend = x[-1]
    X = np.arange(xstart,xend,dx)
    RHO = []
    ind=0
    for i in X:
        if i>=x[ind] and i<x[ind+1]:
            RHO.append(rho[ind])
        else:
            RHO.append(rho[ind+1])
            ind=ind+1
    return X, RHO   


        
def godunovFlux(f,rhol,rhor):
    global sigma
    if rhol<= rhor:
        return min([f(rhol),f(rhor)])
    elif rhol<sigma:
        return f(rhol)
    elif rhor<=sigma and rhol>=sigma:
        return f(sigma)
    else:
        return f(rhor)
        
def godunovFluxes(f,rho):
    global sigma
    roadfluxes = []
    for cell in range(1,len(rho)):   
        roadfluxes.append(godunovFlux(f,rho[cell-1],rho[cell]))
    return roadfluxes
        
def getMaxInflux(f,rho):
    global sigma
    if rho[-1]<sigma:
        return f(rho[-1])
    else:
        return f(sigma)
        
def getMaxOutflux(f,rho):
    global sigma
    if rho[0]<sigma:
        return f(sigma)
    else:
        return f(rho[0])
        
def getIncomingFluxesSB(maxIn,junction):
    #incoming fluxes using a single buffer zone
    incomingFluxes = []
    _c = junction["_c"]
    _M = junction["_M"]
    _Buffer = junction["Buffer"]
    for i in range(len(maxIn)):
        _temp = min(maxIn[i],_c[i]*(_M-sum(_Buffer[-1])))
        incomingFluxes.append(_temp)
    return incomingFluxes

def getIncomingFluxesMB(maxIn,junction):
    #incoming fluxes using separate buffers for every outgoing road
    incomingFluxes = []
    _c = junction["_c"]
    _M = junction["_M"]
    _TM = junction["matrix"]
    _Buffer = junction["Buffer"]
    for i in range(len(maxIn)):
        bufferFluxes = []
        for j in range(len(_TM)):
            bufferFluxes.append(_c[i]*(_M[j]-_Buffer[-1][j])/_TM[j][i])
        _temp = min(maxIn[i],min(bufferFluxes))
        incomingFluxes.append(_temp)
    #print "incoming fluxes are {}".format(incomingFluxes)
    return incomingFluxes
        
def getOutgoingFluxes(maxOut,incomingFluxes,junction):
    #outgoing fluxes: same when using single buffers and multiple buffers
    TM = junction["matrix"]
    Buffer = junction["Buffer"]
    outgoingFluxes = []
    for j in range(len(maxOut)):
        if Buffer[-1][j]>0:
            outgoingFluxes.append(maxOut[j])
        elif Buffer[-1][j]==0:
            temp=np.dot(incomingFluxes,TM[j,:])
            outgoingFluxes.append(min(maxOut[j],temp))
            #print "jup buffer ist leer"
        else:
            print "problem detected: your buffer at road {} is negative".format(j)
    return outgoingFluxes

    
    
def updateBuffer(junction,incomingFluxes,outgoingFluxes):
    #todo: for multiple junctions there are multiple bufferArrays
    global dt
    TM = junction["matrix"]
    Buffer = junction["Buffer"]
    newBuffer = []
    #print len(Buffer)
    for j in range(len(Buffer[-1])):
        temp = np.dot(incomingFluxes,TM[j,:])
        newBuffer.append(Buffer[-1][j]+dt*(temp-outgoingFluxes[j]))
    junction["Buffer"].append(newBuffer)

def godunovStep(fluxes,rho):
    global dx, dt    
 #   print "data = {}".format(data)
    temp = []
    for i in range(0,len(rho)):
        discstep = rho[i]-dt/dx*(fluxes[i+1]-fluxes[i])
        temp.append(discstep)
    return temp

    
def test22(tend):
    global vmax,rhomax,sigma,TM,initIn,initOut,dt,dx,c,M
    #1 initialization
    roadsin = len(TM[0])
    roadsout = len(TM) 
    roads = roadsin+roadsout
    
    maxinflow = [0 for s in range(roadsin)]
    maxoutflow = [0 for s in range(roadsout)]
    
    Densities = [0 for s in range(roads)] 
    X = [0 for s in range(roads)]
    
    X[0],Densities[0] = init([0,100],[0.5],dx)
    X[1],Densities[1] = init([0,100],[0.5],dx)
    X[2],Densities[2] = init([100,160],[0.5],dx)
    X[3],Densities[3] = init([100,160],[0.5],dx)

    R = [[s] for s in Densities]
    print "init: {}".format(R[1])
    flows = [0 for s in range(roads)]
    
    t=dt
    while t<tend:
    #2 compute fluxes    
    #2.1 compute road fluxes
        
        for i in range(roads):
            flows[i]=godunovFluxes(f,R[i][-1])
        for i in range(roadsin):
            maxinflow[i]=getMaxInflux(f,R[i][-1])
        for i in range(roadsout):
            maxoutflow[i]=getMaxOutflux(f,R[i+2][-1])
    
    #2.2 compute boundary fluxes    
    
        incomingFluxes = getIncomingFluxesSB(maxinflow)
        outgoingFluxes = getOutgoingFluxes(maxoutflow,incomingFluxes,TM)
        
        for i in range(roadsin):
            flows[i].append(incomingFluxes[i])
            flows[i].insert(0,f(initIn))
        for i in range(roadsout):
            flows[i+2].append(f(initOut))
            flows[i+2].insert(0,outgoingFluxes[i])
        
    #3 compute Godunov step
        for i in range(roads):
            R[i].append(godunovStep(flows[i],R[i][-1]))
    
    #4 update buffer
        updateBuffer(TM,incomingFluxes,outgoingFluxes)
    
        t+=dt
    plot2D(X,R)

# functions for arbitrary networks

def network2junctions(networkMatrix,c,M,mode="SB"):
    """ networkMatrix: Transition matrix representing the network
        c: vector of length #junctions where each element is a array consisting of the preference values for each outgoing road
        M: vector of length #junctions where each element is a array consisting of the maximum buffer sizes for each outgoing road"""
        
    roads = len(networkMatrix)
    # list to be updated: roads already assigned to a junction are deleted
    roadList = range(roads)
    junctions = {}
    i = 0
    while len(roadList)>0:
        #print roadList
        nonzero_out_tupel = np.nonzero(networkMatrix[:,roadList[0]])
        nonzero_out=nonzero_out_tupel[0]
        if len(nonzero_out)>0:
            junctions[i]= {}
            junctions[i]["out"]=nonzero_out
            nonzero_in_tupel = np.nonzero(networkMatrix[nonzero_out[0],:])
            nonzero_in = nonzero_in_tupel[0]
            #print nonzero_in
            junctions[i]["in"] = nonzero_in
            for el in nonzero_in:
                roadList.remove(el)
            junctions[i]["matrix"] = networkMatrix[np.ix_(nonzero_out,nonzero_in)]
            if (mode=="SB"):
                junctions[i]["_M"]=sum(M[i])
            else:
                junctions[i]["_M"]=M[i]
            junctions[i]["_c"]=c[i]
            junctions[i]["eta"] = [0.95,0.05]
            junctions[i]["Buffer"]=[[0]*len(nonzero_out)]
            i+=1 
        else:
            roadList.remove(roadList[0])
    return junctions
            
def solveArbitraryNetworks(tend,network,_c,_M,initX,initRho,ext_inflow, mode="SB"):
    #TODO: somehow i have to include different buffers and c_i for the different junctions    
    #parameter: #tend: final time
                # network: (n+m)x(n+m) distribution matrix
                # initX: (n+m)-dim array consisting of the space discretizations of the roads dependent on dx 
                # initRho: (n+m)-dim array consisting of the density at discretization points initX
                # ext_inflow: (n+m)-dim vector consisting of the external density entering the roads (0 if no influx, e.g. road is exiting from a junction)
                # mode: "SB" for single buffer junctions, "MB" for multiple buffer junctions    
    global dt,dx,vmax,rhomax,sigma,c,M,RhoPlot, XPlot
    _t = np.arange(0,tend+dt,dt)
    roads = len(network)
    junctions = network2junctions(network,_c,_M,mode)
    #print "junctions in network: {}".format(junctions)
    #init    
    flows = [0 for s in range(roads)]
    R = [[s] for s in initRho]
    t=0
    while t<tend:
        #1 road internal fluxes
        for i in range(roads):
            flows[i]=godunovFluxes(f,R[i][-1])
        #2 fluxes at junction
        for junc in junctions:
            indices_in = junctions[junc]["in"]
            rho_in = []
            for i in indices_in:
                rho_in.append(R[i][-1])
            indices_out = junctions[junc]["out"]
            rho_out = []
            for i in indices_out:
                rho_out.append(R[i][-1])
            jfluxes_in,jfluxes_out = getJunctionFluxes(rho_in,rho_out,junctions[junc],mode) 
            updateBuffer(junctions[junc],jfluxes_in,jfluxes_out)        
            j=0    
            k=0
            for i in junctions[junc]["in"]:
                flows[i].append(jfluxes_in[j])
                j+=1
            for i in junctions[junc]["out"]:
                flows[i].insert(0,jfluxes_out[k])
                k+=1
        #3 external fluxes
        rowSum = np.sum(network,axis=1)
        columnSum = np.sum(network,axis=0)
        for i in range(roads):
            if rowSum[i]==0: #has external influx
                flows[i].insert(0,f(ext_inflow[i]))
            elif columnSum[i]==0: #has external outflow
                #flows[i].append(f(ext_outflow[i]))
                flows[i].append(flows[i][-1])
        #4 godunov step
        for i in range(roads):
            R[i].append(godunovStep(flows[i],R[i][-1]))
        t+=dt
    #update global var for animation
    RhoPlot = R
    XPlot = initX
    
    doAnimation("test")
    return _t,R,initX,junctions

def solveNetworkWithLifetime(tend,network,_c,_M,initX,initRho,ext_inflow,lifetime, mode="SB"):
    global dt,dx,vmax,tend_input,rhomax,tend_input,sigma,c,M,RhoPlot, XPlot, totalflows,ext_in2
    _t = np.arange(0,tend+dt,dt)
    roads = len(network)
    junctions = network2junctions(network,_c,_M,mode)
    #print "junctions in network: {}".format(junctions)
    #init    
    totalflows = []
    flows = [0 for i in range(roads)]
    R = [[s] for s in initRho]
    t=0
    controls = [[],[]]
    while t<tend:
        #1 road internal fluxes
        for i in range(roads):
            flows[i]=godunovFluxes(f,R[i][-1])
        #2 fluxes at junction
        for junc in junctions:
            indices_in = junctions[junc]["in"]
            rho_in = []
            for i in indices_in:
                rho_in.append(R[i][-1])
            indices_out = junctions[junc]["out"]
            rho_out = []
            for i in indices_out:
                rho_out.append(R[i][-1])
            jfluxes_in,jfluxes_out,c= getJunctionFluxesWithLifetime(rho_in,rho_out,junctions[junc],lifetime,t,mode) 
            controls[junc].append(c)
            updateBuffer(junctions[junc],jfluxes_in,jfluxes_out)        
            j=0    
            k=0
            for i in junctions[junc]["in"]:
                flows[i].append(jfluxes_in[j])
                j+=1
            for i in junctions[junc]["out"]:
                flows[i].insert(0,jfluxes_out[k])
                k+=1
        #3 external fluxes
        rowSum = np.sum(network,axis=1)
        columnSum = np.sum(network,axis=0)
        for i in range(roads):
            if rowSum[i]==0: #has external influx
                if t<tend_input:
                    if t%120<=60:
                        flows[i].insert(0,f(ext_inflow[i]))
                    else:
                        flows[i].insert(0,f(ext_in2[i]))
                else:
                    flows[i].insert(0,0)
            elif columnSum[i]==0: #has external outflow
                #flows[i].append(f(ext_outflow[i]))
                flows[i].append(flows[i][-1])
        #4 godunov step
		totalflows.append(flows[:])
        for i in range(roads):
            R[i].append(godunovStep(flows[i],R[i][-1]))
        t+=dt
    #update global var for animation
    RhoPlot = R
    XPlot = initX
    #traveltime = computeTravelTime(0,2,tend,R)
    return _t,R,initX,junctions, totalflows,controls

def getJunctionFluxesWithLifetime(rho_in,rho_out,junction,lifetime,t,mode="SB"):
    global switchingTimes, switchingValues,switchingPeriod
    #TODO: adjust computation of TM for multiple outgoing roads
    TM = junction["matrix"]
    #for i in range(len(TM)):
      #  TM[i][0]=switchingValues[findIndex(switchingTimes,t-tao)]
       # TM[i][1]=(1-switchingValues[findIndex(switchingTimes,t-tao)])
    #eta = [switchingValues[findIndex(switchingTimes,t-tao)],1-switchingValues[findIndex(switchingTimes,t-tao)]]
    eta = [1.,0.0]
    roadsin = len(TM[0])
    roadsout = len(TM)
    maxinflow = [0 for s in range(roadsin)]
    maxoutflow = [0 for s in range(roadsout)]
    
    for i in range(roadsin):
        maxinflow[i]=getMaxInflux(f,rho_in[i])
    for j in range(roadsout):
        maxoutflow[j]=getMaxOutflux(f,rho_out[j])
    if(mode=="SB"):
        incomingFluxes = getIncomingFluxesSB(maxinflow,junction)
        if t%(2*lifetime)>lifetime:
            eta=list(reversed(eta))
        for i in range(len(incomingFluxes)):
            incomingFluxes[i]=incomingFluxes[i]*eta[i]
    elif(mode=="MB"):
        incomingFluxes = getIncomingFluxesMB(maxinflow,junction)
        for i in range(len(incomingFluxes)):
            incomingFluxes[i]=incomingFluxes[i]*eta[i]
    outgoingFluxes = getOutgoingFluxes(maxoutflow,incomingFluxes,junction)
    
    return incomingFluxes, outgoingFluxes,eta[0]

def getJunctionFluxes(rho_in,rho_out,junction,mode="SB"):
    
    TM = junction["matrix"]
    roadsin = len(TM[0])
    roadsout = len(TM)
    maxinflow = [0 for s in range(roadsin)]
    maxoutflow = [0 for s in range(roadsout)]
    
    for i in range(roadsin):
        maxinflow[i]=getMaxInflux(f,rho_in[i])
    for j in range(roadsout):
        maxoutflow[j]=getMaxOutflux(f,rho_out[j])
    if(mode=="SB"):
        incomingFluxes = getIncomingFluxesSB(maxinflow,junction)
    elif(mode=="MB"):
        incomingFluxes = getIncomingFluxesMB(maxinflow,junction)
    outgoingFluxes = getOutgoingFluxes(maxoutflow,incomingFluxes,junction)
    
    return incomingFluxes, outgoingFluxes

def computeTravelTime(roadIndexStart,roadIndexEnd,tend,fluxes):
    global dt
    n=0
    totalInflux = 0
    travelTime = 0
    while n*dt<tend:
        totalInflux += dt*fluxes[n][roadIndexStart][0]
        n+=1
    m=0
    while m*dt<tend:
        travelTime +=dt*m*dt*fluxes[m][roadIndexEnd][-2]
        m+=1
    nettoTravelTime =1./totalInflux*travelTime
    return nettoTravelTime

def findIndex(switchingTimes,currentTime):
    index = 0
    while currentTime > switchingTimes[index]  and index<len(switchingTimes)-1:
        index+=1
    return index

def getTotalFluxesOnNetwork(fluxes):
    totalFluxes = []
    for i in range(len(fluxes)):
        totalFluxes.append(getRoadFluxes(fluxes[i]))
    return totalFluxes
            
def getRoadFluxes(fluxes):
    return [sum(i) for i in fluxes]     

def getOverallFlux(fluxes):
    fluxes_on_roads = getTotalFluxesOnNetwork(fluxes)
    totalFlux = sum([sum(i) for i in fluxes_on_roads])
    return totalFlux
def avgTravelTimes(fluxes,tend):
    times = []
    for index in range(len(fluxes[0])):
        times.append(computeTravelTime(index,index,tend,fluxes))
    return times
    
# plot functions
 
def plot2D(x,rho):

    f,ax = plt.subplots(len(rho))
    plt.xlabel('x')
    f.set_size_inches(10, 8)
    plt.ylabel('\rho')
    #plt.ylim([0,0.8])
    for i in range(len(rho)):
        ax[i].plot(x[i][:-1],np.array(rho[i][-1][:-1]))  
        ax[i].set_ylabel(r'$\rho$')
        ax[i].set_yticks([0.0,0.25,0.5,0.75])
        
def plotBuffers(t,junctions):
    fig = plt.figure()
    plt.xlabel("t")
    plt.ylabel("total buffers")
    legend = []
    i=0
    for junc in junctions:
        Buffer = junctions[junc]["Buffer"]
        #print junctions[junc]["matrix"]
        bufferSum = [sum(buf) for buf in Buffer]
        #print bufferSum
        ax = fig.add_subplot(111)
        ax.plot(np.array(t),np.array(bufferSum))
        legend.append("junction {}".format(i))
        i=i+1
    plt.legend(legend)
   
def plotTotalDensity(t,rho):
    totalRhoRoad = []
    totalRho = []
    for j in range(len(rho)):
        totalRhoRoad.append([dx*sum(rho[j][i]) for i in range(len(rho[0]))])
    for s in range(len(totalRhoRoad[0])):
        total =0
        for i in range(len(totalRhoRoad)):
            total+=totalRhoRoad[i][s]
        totalRho.append(total)
    fig = plt.figure()
    plt.xlabel("t")
    plt.ylabel("total density on the network")
    ax = fig.add_subplot(111)
    ax.plot(t,totalRho)
    #plt.ylim([100,200])

def animate(i):
    global RhoPlot,XPlot,line,fig,ax
    ax.clear()    
    ax.plot(XPlot[0],RhoPlot[0][i])
    ax.set_ylim([0,1]) 
    return line,

def initAnim():
    global line
    line.set_data([],[])
    return line,


def doAnimation(filename="test_animation"):
    global RhoPlot,XPlot,fig,ax,line,tend,dt
    frames = int(tend/dt)
    filenamefull = filename+".mp4"
    print filenamefull
    anim = animation.FuncAnimation(fig, animate, frames=frames,
                              blit=True, init_func=initAnim)
    anim.save(filenamefull, fps=30, extra_args=['-vcodec', 'libx264'])
    plt.show()

def doAnimationNew(filename="test_animation"):
    global RhoPlot,XPlot,tend,dt,linenew, axnew
    f,axnew = plt.subplots(len(RhoPlot))
    #f.figsize=(5,4*len(RhoPlot)) 
    linenew = []

    def initAnimNew():
        global linenew, RhoPlot,axnew
        return linenew
    def animateNew(i):
        global RhoPlot,axnew,linenew
        for j in range(len(RhoPlot)):
            axnew[j].clear()
            axnew[j].plot(XPlot[j][:-2],RhoPlot[j][i][:-2])
            axnew[j].set_ylim([0,1])
        return linenew
    frames = int(tend/dt)
    filenamefull = filename+".mp4"   
    anim = animation.FuncAnimation(f, animateNew, frames=frames,
                              blit=True, init_func=initAnimNew)
    anim.save(filenamefull, fps=30, extra_args=['-vcodec', 'libx264'])
    plt.show()
    
        
        
def simu(network,_c,_M,_vmax,_rhomax,_sigma,_ext_in,_dt,_dx,_tend,_X,_Rho,mode="SB"):
    global switchingPeriod,RhoPlot, XPlot
    vmax = _vmax
    rhomax = _rhomax
    sigma = _sigma
    dt = _dt
    tend = _tend
    M = _M
    c = _c
    ext_in = _ext_in
    
    t = np.arange(0,tend,dt)
    
    Densities = [0]*len(_Rho) 
    X = [0]*len(_X)
    
    for i in range(len(X)):
        X[i],Densities[i] = init(_X[i],_Rho[i],dx)
    
    #lifetime = np.arange(1,150.,1)  
    lifetime = [34]
    #tauValues = np.arange(2*switchingPeriod) 
    #tauValues = [34]
    traveltimes = []
    for lt in lifetime:
        time,R,XX,junctions, totalflows,controls = solveNetworkWithLifetime(tend,network,c,M,X,Densities,ext_in,lt,mode)
        traveltimes.append((lt,getOverallFlux(totalflows)))       
    #print traveltimes
    plot2D(XX,R)
    #ff = plt.figure()
    #plt.plot(t,controls[0],label="$\eta_C(t)$")
    #plt.plot(t,controls[1],label="$\eta_C(t-34s)$")
    plt.legend()
    plotBuffers(time,junctions)
    plotTotalDensity(time,R)
    XPlot = XX
    RhoPlot =R
    return traveltimes
    
#global variables   
fig = plt.figure()
ax = fig.add_subplot(111)     
axnew = fig.add_subplot(111) 
linenew, = axnew.plot([],[])
line, = ax.plot([],[])
ax.set_ylim([0,1]) 
RhoPlot = []
XPlot = []
vmax=1
rhomax =1
sigma=.5
initIn = .5
initOut = 0.5
dt=.5
totalflows = []
tend = 500
tend_input = +200
switchingPeriod = 60
t = np.arange(0,tend,dt)
dx=.5    
#M=[0.5]
M=[[0.5,0.5],[0.5,0.5]]
c=[[1.,1.],[1.,1.]]
control = [1,0]
switchingTimes = np.arange(0,tend,switchingPeriod)
myIterator = cycle([0.95,0.05])
switchingValues = [0.05]+[myIterator.next() for i in range(len(switchingTimes))]

#R=test22(100)

#maybe include following in global init function
_X = [[0,50],[0,50],[50,100]]
_Rho = [[0.2],[0.2],[0.2]]
network = np.array([[0,0,0],[0,0,0],[1.,1.,0]])
ext_in = [0.4,0.2,0.]
ext_in2 = [0.4,0.2,0]
#_X = [[0,50],[0,50],[50,150],[100,150],[150,200]]
#_Rho = [[0.2],[0.2],[0.2],[0.2],[0.2]]
#X[4],Densities[4] = init([160,220],[0.5],dx)
#
#network = np.array([[0,0,0,0,0],[0,0,0,0,0],[1.,1.,0,0,0],[0,0,0,0,0],[0,0,1.,1.,0]]) 
#
#ext_in = [0.4,0.2,0,0.2,0]
#ext_in2 = [0.4,0.2,0,0.2,0]

#ext_in5 = [0.5,0.5,0,0,0]
#ext_out5 = [0,0,0.5,0,0.5]
fluxes=simu(network,c,M,vmax,rhomax,sigma,ext_in,dt,dx,tend,_X,_Rho,"SB")
doAnimationNew("34")