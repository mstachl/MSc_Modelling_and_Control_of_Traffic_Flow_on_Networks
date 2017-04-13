# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 13:39:51 2017

@author: Markus
"""

from matplotlib import cm
from itertools import cycle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import minimize
import copy
import time
import pickle

#functions

def velo(rho):
    global rhomax, vmax
    if rho<rhomax:
        return vmax*(1-rho/rhomax)
    return 0
    
def f(rho):
    global vmax, rhomax
    flux = rho*velo(rho)
    return flux
    
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

def getIncomingFluxesSB_MPC(Buffer,maxIn,junction):
    #incoming fluxes using a single buffer zone
    incomingFluxes = []
    _c = junction["_c"]
    _M = junction["_M"]
    #_Buffer = junction["Buffer"]
    _Buffer = Buffer
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

def getIncomingFluxesMB_MPC(Buffer,maxIn,junction):
    #incoming fluxes using separate buffers for every outgoing road
    incomingFluxes = []
    _c = junction["_c"]
    _M = junction["_M"]
    _TM = junction["matrix"]
    #_Buffer = junction["Buffer"]
    _Buffer = Buffer
    for i in range(len(maxIn)):
        bufferFluxes = []
        for j in range(len(_TM)):
            bufferFluxes.append(_c[i]*(_M[j]-_Buffer[-1][j])/_TM[j][i])
        _temp = min(maxIn[i],min(bufferFluxes))
        incomingFluxes.append(_temp)
    #print "incoming fluxes are {}".format(incomingFluxes)
    return incomingFluxes
        
def getOutgoingFluxes_MPC(Buffer,maxOut,incomingFluxes,junction):
    #outgoing fluxes: same when using single buffers and multiple buffers
    TM = junction["matrix"]
    #Buffer = junction["Buffer"]
    Buffer = Buffer
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
   # print incomingFluxes,outgoingFluxes,Buffer[-1]
    return outgoingFluxes

    
    
def updateBuffer(junction,incomingFluxes,outgoingFluxes):
    #todo: for multiple junctions there are multiple bufferArrays
    global dt
    TM = junction["matrix"]
    Buffer = junction["Buffer"]
    newBuffer = []
    for j in range(len(Buffer[-1])):
        temp = np.dot(incomingFluxes,TM[j,:])
        if Buffer[-1][j]+dt*(temp-outgoingFluxes[j])>0:
            newBuffer.append(Buffer[-1][j]+dt*(temp-outgoingFluxes[j]))
        else:
            newBuffer.append(0)
    junction["Buffer"].append(newBuffer)
    #print junction["Buffer"]
    
    
def computeBuffer_MPC(Buffer,TM,incomingFluxes,outgoingFluxes):
    #todo: for multiple junctions there are multiple bufferArrays
    global dt
    newBuffer = []
    for j in range(len(outgoingFluxes)):
        temp = np.dot(incomingFluxes,TM[j,:])
        if Buffer[-1][j]+dt*(temp-outgoingFluxes[j])>0:
            newBuffer.append(Buffer[-1][j]+dt*(temp-outgoingFluxes[j]))
        else:
            newBuffer.append(0)
    return newBuffer
    #print junction["Buffer"]

def updateBuffer_MPC(junctions):
    #todo: for multiple junctions there are multiple bufferArrays
    global MPC_buffers
    for junc in junctions:
        Buffer = junctions[junc]["Buffer"]
        Buffer.extend(MPC_buffers[junc])

def godunovStep(fluxes,rho):
    global dx, dt    
 #   print "data = {}".format(data)
    temp = []
    for i in range(0,len(rho)):
        discstep = rho[i]-dt/dx*(fluxes[i+1]-fluxes[i])
        if discstep>1:
            print "!!!!!!density >1!!!!!"
        temp.append(discstep)
    return temp
# functions for arbitrary networks

def network2junctions(network,c,M,mode="SB"):
    """ networkMatrix: Transition matrix representing the network
        c: vector of length #junctions where each element is a array consisting of the preference values for each outgoing road
        M: vector of length #junctions where each element is a array consisting of the maximum buffer sizes for each outgoing road"""
    networkMatrix = np.array(network)
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
            junctions[i]["Buffer"]=[[0 for s in range(len(nonzero_out))]]
            i+=1 
        else:
            roadList.remove(roadList[0])
    return junctions

def network2Adjescency(network):
    #removes 0 rows and columns of the network matrix 
    # returns matrix of dim (n_in x n_out) where n_in number of all outgoing 
    # roads of all junctions in the network
    temp = np.array(network)
    for i in range(len(network)):
        for j in range(len(network)):
            if network[i][j]>0:
                temp[i,j]=1
    rowSum = np.sum(temp,axis=1)
    colSum = np.sum(temp,axis=0)
    rows_to_slice = [i for i,j in enumerate(rowSum) if j==0]
    cols_to_slice = [i for i,j in enumerate(colSum) if j==0]
    temp=np.delete(temp,rows_to_slice,axis=0)
    temp = np.delete(temp,cols_to_slice,axis=1)
    return temp

            

def getJunctionFluxesMPC(Buffer,rho_in,rho_out,u,junction,mode="SB"):
    ## u is control vector
    #u=np.array(u)
    TM = junction["matrix"]
    inroads = junction["in"]
    #outroads = junction["out"]
    roadsin = len(TM[0])
    roadsout = len(TM)
    maxinflow = [0 for s in range(roadsin)]
    maxoutflow = [0 for s in range(roadsout)]
    for i in range(roadsin):
        maxinflow[i]=getMaxInflux(f,rho_in[i])
    for j in range(roadsout):
        maxoutflow[j]=getMaxOutflux(f,rho_out[j])
    if(mode=="SB"):
        incomingFluxes = getIncomingFluxesSB_MPC(Buffer,maxinflow,junction)
        for i in range(len(incomingFluxes)):
            incomingFluxes[i]=incomingFluxes[i]*u[inroads[i]]
    elif(mode=="MB"):
        incomingFluxes = getIncomingFluxesMB_MPC(Buffer,maxinflow,junction)
        for i in range(len(incomingFluxes)):
            incomingFluxes[i]=incomingFluxes[i]*u[inroads[i]]
    outgoingFluxes = getOutgoingFluxes_MPC(Buffer,maxoutflow,incomingFluxes,junction)
    #print "junction fluxes: ",incomingFluxes,outgoingFluxes
    return incomingFluxes, outgoingFluxes

def getTotalFluxesOnNetwork(fluxes):
    totalFluxes = []
    for i in range(len(fluxes)):
        totalFluxes.append(getRoadFluxes(fluxes[i]))
    return totalFluxes
            
def getRoadFluxes(fluxes):
    return [sum(i) for i in fluxes]     

def getOverallFlux(fluxes):
    fluxes_on_roads = getTotalFluxesOnNetwork(fluxes)
    totalFlux = sum([sum(i)*dt*dx for i in fluxes_on_roads])
    return totalFlux

def dynamicEvolution(rho_init,control,t0,n_pred):
    global junctions,n_p,n_c,vmax,tend_input,tend,ext_in,rhomax,sigma,c,M,RhoPlot, XPlot, Rho_MPCstep
    roads = len(rho_init)    
    flows = [0 for i in range(roads)]
    #u=np.array(u)
    #print u
    tm = [junctions[junc]["matrix"] for junc in junctions]
    totalflows = []
    Current_Buffers = [junctions[junc]["Buffer"][-1] for junc in junctions]
    #print rho_init
    R = [[s] for s in rho_init]
    New_Buffers=[[Current_Buffers[i]] for i in range(len(Current_Buffers))]
    t=t0
    sum_of_incoming_roads = 0
    for junc in junctions:
        sum_of_incoming_roads+=len(junctions[junc]["in"])
    u_laststep = control[-sum_of_incoming_roads:]
    unpacked_controls = []
    for i in range(n_p):
        temp = np.asarray([control[i*sum_of_incoming_roads:(i+1)*sum_of_incoming_roads]])
        #print temp
        temp = temp.repeat(signal_horizon,axis=0)
        #print temp
        temp_list_multidim = temp.tolist()
        temp_list_1d = sum(temp_list_multidim,[])
        unpacked_controls.append(temp_list_1d)
    u = sum(unpacked_controls,[])
    #print "len u in DE: {}".format(len(u))
    norm_u_dot = 0
    w_eta = 0
    timestep_in_MPC=0
    #print signal_horizon,n_pred,dt,t0
    #print np.arange(t0,t0+n_pred*signal_horizon*dt,dt)
    for timestep_in_MPC in range(n_pred*signal_horizon):
    #for t in np.arange(t0,t0+n_pred*signal_horizon*dt,dt):
    #while t<t0+n_pred*signal_horizon*dt:
        #1 road internal fluxes
        #print timestep_in_MPC,t
        for i in range(roads):
            flows[i]=godunovFluxes(f,R[i][-1])
            #print len(flows[i])
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
            #print column(u,timestep_in_MPC)
            jfluxes_in,jfluxes_out = getJunctionFluxesMPC(New_Buffers[junc],rho_in,rho_out,u[sum_of_incoming_roads*timestep_in_MPC:sum_of_incoming_roads*(timestep_in_MPC+1)],junctions[junc]) 
            #updateBuffer(junctions[junc],jfluxes_in,jfluxes_out)        
            New_Buffers[junc].append(computeBuffer_MPC(New_Buffers[junc],tm[junc],jfluxes_in,jfluxes_out))       
            j=0    
            k=0
            for i in junctions[junc]["in"]:
                flows[i].append(jfluxes_in[j])
                j+=1
            for i in junctions[junc]["out"]:
                flows[i].insert(0,jfluxes_out[k])
                k+=1
            ### compute W(\eta)
            for road_index in range(len(indices_in)):
                temp = 1.
                temp*=(u[timestep_in_MPC*sum_of_incoming_roads+indices_in[road_index]]-1.)
                for i in range(len(indices_in)-1):
                    temp*=(u[timestep_in_MPC*sum_of_incoming_roads+indices_in[road_index]]-0.)                
                w_eta +=pow(temp,2)
                #print u[timestep_in_MPC*sum_of_incoming_roads+indices_in[road_index]]
        ## compute ||\dot(\eta)||
        u_with_last_control = np.append(u_laststep,u)       
        for i in range(sum_of_incoming_roads):
            norm_u_dot+=pow(u_with_last_control[sum_of_incoming_roads*(timestep_in_MPC+1)+i]-u_with_last_control[sum_of_incoming_roads*timestep_in_MPC+i],2)               
                
        #timestep_in_MPC+=1
        #3 external fluxes
       
        rowSum = np.sum(network,axis=1)
        columnSum = np.sum(network,axis=0)
        for i in range(roads):
            if rowSum[i]==0: #has external influx
                #print "influx at {}".format(i)
                if t<tend_input:
                    flows[i].insert(0,f(ext_in[i]))
                else:
                    flows[i].insert(0,0)
            elif columnSum[i]==0: #has external outflow
                #flows[i].append(f(ext_outflow[i]))
                flows[i].append(flows[i][-1])
        ##4 godunov step
        #print "flows at road 0:",flows[0]
        #print "flows: ", flows
        #if np.sum(flows)>0.5:
        #    finalTravelTime=t
        totalflows.append(flows[:])
        #print totalflows
        for i in range(roads):
            R[i].append(godunovStep(flows[i],R[i][-1]))
        #print t0,t
        t+=dt
    R_new = [R[i][1:] for i in range(len(R))]
    #print len(R_new)
    Rho_MPCstep = copy.deepcopy(R_new)
    New_Buffers = [New_Buffers[i][1:] for i in range(len(New_Buffers))]
    #print Rho_MPCstep
    return R_new,totalflows, New_Buffers,norm_u_dot,w_eta
    


def solveNetworkWithMPC(_X,_Rho,network):
    global XPlot,RhoPlot,tend,dx,dt,c,M,signal_horizon, Rho_init,Rho_MPCstep,junctions
    junctions = network2junctions(network,c,M)
    sum_of_incoming_roads = 0
    
    for junc in junctions:
        sum_of_incoming_roads+=len(junctions[junc]["in"])
    u_init=sum([[0 for i in range(sum_of_incoming_roads)] for i in range(n_p)],[])
    t0=0.
    Densities = [0 for i in range(len(_Rho))] 
    X = [0 for i in range(len(_X))]
    controls_sol = []
    for i in range(len(X)):
        X[i],Densities[i] = init(_X[i],_Rho[i],dx)
    Rho_init = Densities
    Rho_sol = [[Densities[i]] for i in range(len(Densities))]
    XPlot = X
    res_fluxes = []

    t=[]
    ftt = 0
    while t0<tend:
        time_vector = np.arange(t0,t0+n_c*signal_horizon*dt,dt).tolist()
        t.extend(time_vector)
        R,flows,buffers,norm_u_dot,w_eta = dynamicEvolution(Rho_init,u_init,t0,n_p)
        controls=u_init
        res_fluxes.append(flows)
        #print "len rho_sol before: {}".format(len(Rho_sol[0]))
        for i in range(len(R)):
            Rho_sol[i].extend(R[i])
        #print len(Rho_sol[0])
        Rho_init = [R[i][-1] for i in range(len(R))]
        t0+=n_c*signal_horizon*dt
        for i in range(n_c):
            temp = np.asarray([controls[i*sum_of_incoming_roads:(i+1)*sum_of_incoming_roads]])
            #print temp
            temp_v = temp.repeat(signal_horizon,axis=0)
            #print temp
            temp_list_multidim = temp_v.tolist()
            temp_list_1d = sum(temp_list_multidim,[])
            controls_sol.extend(temp_list_1d)
        #print controls_sol, len(controls_sol)
        #for i in range(len(Rho_MPCstep)):
            #Rho_sol[i].extend(Rho_MPCstep[i][:n_c*signal_horizon])
        #for junc in junctions:
           # print junctions[junc]["Buffer"]
        RhoPlot = Rho_sol
    control_vectors=computeControlVectors(controls_sol,network)  
    RhoPlot = Rho_sol[:]
    return junctions,t,XPlot, RhoPlot, control_vectors, res_fluxes


def computeControlVectors(all_controls,network):
    global c,M,t
    junctions = network2junctions(network,c,M)
    roads_incoming = []
    for junc in junctions:
        roads_incoming.extend(junctions[junc]["in"])
    number_incoming_roads=len(roads_incoming)
    control_vector = [[] for i in range(number_incoming_roads)]
    for control in range(len(all_controls)):
        index = control%number_incoming_roads
        control_vector[index].append(all_controls[control])
    return control_vector


#### plot functions
def plot2D(x,rho):
    f,ax = plt.subplots(len(rho))
    plt.xlabel('x')
    #plt.tight_layout()
    f.set_size_inches(10, 20)
    #plt.yticks([0,0.5])
    #plt.ylim([0,0.8])
    for i in range(len(rho)):
        ax[i].plot(x[i][:-1],np.array(rho[i][-1][:-1]))  
        ax[i].set_ylabel(r'$\rho$')
        ax[i].set_yticks([0.0,0.25,0.5,0.75])
        ax[i].set_xlim([x[i][0],x[i][-2]])
        
def plot2DatTime(x,rho,s):
    global dt
    f,ax = plt.subplots(len(rho))
    plt.xlabel('x')
    f.set_size_inches(6, 24)

    if isinstance(s,int):
        s=[s]
    #plt.tight_layout()
    
    #plt.yticks([0,0.5])
    #plt.ylim([0,0.8])
    for t in s:
        index = int(t/dt)
        for i in range(len(rho)):
            ax[i].plot(x[i][:-1],np.array(rho[i][index][:-1]))  
            ax[i].set_ylabel(r'$\rho$')
            ax[i].set_yticks([0.0,0.25,0.5,0.75])
            ax[i].set_xlim([x[i][0],x[i][-2]])
   # plt.savefig("plot_at_{}.eps".format(t),format="eps")
        
def plotControls(control_vector,t,index):
    print len(t),len(control_vector[index])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('time t')
    ax.set_ylabel(r'$\eta$')
    ax.plot(t,control_vector[index][:len(t)])
    fig.savefig("output_fdp_{}_{}_{}_control_{}.eps".format(n_p,n_c,gamma,index),format='eps')
    
def plotBuffers(t,junctions):
    fig = plt.figure()
    plt.xlabel("t")
    plt.ylabel("total buffers")
    legend = []
    i=0
    for junc in junctions:
        Buffer = junctions[junc]["Buffer"]
        bufferSum = [sum(buf) for buf in Buffer]
        ax = fig.add_subplot(111)
        ax.plot(np.array(t),np.array(bufferSum[:len(t)]))
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
    print len(t),len(totalRho)
    ax.plot(t,totalRho[:len(t)])
    #plt.ylim([100,200])
    
    
def doAnimationNew(x,rho,filename="test_animation"):
    global RhoPlot,XPlot,tend,dt,linenew, axnew
    f,axnew = plt.subplots(len(rho))
    #f.figsize=(5,4*len(RhoPlot)) 
    linenew = []
    def initAnimNew():
        global linenew, RhoPlot,axnew
        return linenew
    def animateNew(i):
        global RhoPlot,axnew,linenew
        for j in range(len(rho)):
            axnew[j].clear()
            axnew[j].plot(x[j][:-2],rho[j][i][:-2])
            axnew[j].set_ylim([0,1])
        return linenew
    frames = int(tend/dt)
    filenamefull = filename+".mp4"   
    anim = animation.FuncAnimation(f, animateNew, frames=frames,
                              blit=True, init_func=initAnimNew)
    anim.save(filenamefull, fps=30, extra_args=['-vcodec', 'libx264'])
    plt.show()
junctions = []

RhoPlot = []
XPlot = []
vmax=50./3.6
rhomax =1
sigma=.5
dt=0.5
tend = 50
tend_input = 200.
t = np.arange(0,tend,dt)
dx=10    

n_junctions = 3
M=[[2,2] for i in range(n_junctions)]
c=[[1,1] for i in range(n_junctions)]

MPC_flows_laststep=[]
MPC_flows=[]

MPC_buffers=[]


#_X=[[-400,0,400],[-200,0,400],[0,350],[-700,-500,350],[0,400],[0,400]]
#_Rho = [[1,0.0],[0.8,0.0],[0.0],[0.8, 0.0],[0.0],[0.0]]
ext_in=[0,0,0,0,0,0]

#_X=[[0,600],[0,400],[0,350],[0,350],[0,400],[0,400]]
#_Rho = [[0.0],[0.0],[0.0],[0.0],[0.0],[0.0]]
#ext_in = [0.4,0.2,0,0.2,0,0]

network = np.array([[0,0,0,0,0,0],[0,0,0,0,0,0],[1,1,0,0,0,0],[0,0,0,0,0,0],
                    [0,0,0.8,0.6,0,0],[0,0,0.2,0.4,0,0]]) 



_X=[[-1000,0,400],[-1000,0,500],[0,500],[-800,0,860],[0,400],[0,400]]

#_X = [[-400,200],[-400,200],[200,400],[200,700],[500,700],[700,900],[700,1550],[1350,1550],[1550,1750],[1550,1700]]
_Rho = [[1,0.0],[1,0.0],[0.0],[1, 0.0],[0.0],[0.0]]
ext_in = [0,0,0,0,0,0]

### 2. two in one out
#_X=[[-400,0,400],[-200,0,500],[0,350]]
#_Rho=[[1,0.0],[1,0.0],[1.0]]
#ext_in = [0.,0.,0]
#network=np.array([[0,0,0],[0,0,0],[1,1,0]])


### 3. one in one out
#_X=[[-300,0,100],[0,100]]
#_Rho=[[1,0],[1]]
#network=np.array([[0,0],[1,0]])
#ext_in=[0.,0]
Rho_init = []
Rho_MPCstep = []


gamma=10
eps=5
signal_horizon = 50
n_p=2
n_c=2
ftt = 0
starttime = time.time()
junctions,t,X,R,u,fluxes = solveNetworkWithMPC(_X,_Rho,network)
totaltime = time.time()-starttime
overallFlux = np.array(fluxes)
while isinstance(overallFlux,(list,np.ndarray)):
    overallFlux = np.sum(overallFlux, axis=-1)
print overallFlux
#data = {"t": t,"x": X, "rho": R, "fluxes": fluxes, "total flux": overallFlux, "total time": totaltime, "controls": u, "T": ftt, "feval": feval,"junctions": junctions}
#with open('output_buffer_dp_{}_{}_{}_{}.txt'.format(signal_horizon,n_p,n_c,gamma), 'wb') as handle:
#    pickle.dump(data, handle)
plot2D(X,R)
plotControls(u,t,0)
plotBuffers(t,junctions)