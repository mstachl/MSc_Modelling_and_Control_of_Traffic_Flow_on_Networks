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
    for j in range(len(Buffer[-1])):
        temp = np.dot(incomingFluxes,TM[j,:])
        if Buffer[-1][j]+dt*(temp-outgoingFluxes[j])>0:
            newBuffer.append(Buffer[-1][j]+dt*(temp-outgoingFluxes[j]))
        else:
            newBuffer.append(0)
    junction["Buffer"].append(newBuffer)
    #print junction["Buffer"]
    
    
def computeBuffer_MPC(junction,incomingFluxes,outgoingFluxes):
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
            junctions[i]["Buffer"]=[[0]*len(nonzero_out)]
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

            

def getJunctionFluxesMPC(rho_in,rho_out,u,junction,mode="SB"):
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
        incomingFluxes = getIncomingFluxesSB(maxinflow,junction)
        for i in range(len(incomingFluxes)):
            incomingFluxes[i]=incomingFluxes[i]*u[inroads[i]]
    elif(mode=="MB"):
        incomingFluxes = getIncomingFluxesMB(maxinflow,junction)
        for i in range(len(incomingFluxes)):
            incomingFluxes[i]=incomingFluxes[i]*u[inroads[i]]
    outgoingFluxes = getOutgoingFluxes(maxoutflow,incomingFluxes,junction)
    #print "junction fluxes: ",incomingFluxes,outgoingFluxes
    return incomingFluxes, outgoingFluxes

def dynamicEvolution(rho_init,u,t0,n_pred):
    global junctions,n_p,n_c,vmax,tend_input,tend,ext_in,rhomax,sigma,c,M,RhoPlot, XPlot, Rho_MPCstep
    roads = len(rho_init)    
    flows = [0 for i in range(roads)]
    #u=np.array(u)
    #print u
    totalflows = []
    
    Buffer = [[] for i in range(len(junctions))]
    #print rho_init
    R = [[s] for s in rho_init]
    t=t0
    sum_of_incoming_roads = 0
    for junc in junctions:
        sum_of_incoming_roads+=len(junctions[junc]["in"])
    u_laststep = u[:sum_of_incoming_roads]
#    def column(matrix, i):
#        return [row[i] for row in matrix]
    norm_u_dot = 0
    w_eta = 0
    timestep_in_MPC=0
    while t<=t0+n_pred*dt:
        #1 road internal fluxes
        #print timestep_in_MPC,t,n_p
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
            #print column(u,timestep_in_MPC)
            jfluxes_in,jfluxes_out = getJunctionFluxesMPC(rho_in,rho_out,u[sum_of_incoming_roads*timestep_in_MPC:sum_of_incoming_roads*(timestep_in_MPC+1)],junctions[junc]) 
            #updateBuffer(junctions[junc],jfluxes_in,jfluxes_out)        
            Buffer[junc].append(computeBuffer_MPC(junctions[junc],jfluxes_in,jfluxes_out))       
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
        ### compute ||\dot(\eta)||
        u_with_last_control = np.append(u_laststep,u)
        #print u_with_last_control        
        for i in range(sum_of_incoming_roads):
            norm_u_dot+=pow(u_with_last_control[sum_of_incoming_roads*(timestep_in_MPC+1)+i]-u_with_last_control[sum_of_incoming_roads*timestep_in_MPC+i],2)               
                
        timestep_in_MPC+=1
        #3 external fluxes
       
        rowSum = np.sum(network,axis=1)
        columnSum = np.sum(network,axis=0)
        for i in range(roads):
            if rowSum[i]==0: #has external influx
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
        totalflows.append(flows[:])
        #print totalflows
        for i in range(roads):
            R[i].append(godunovStep(flows[i],R[i][-1]))
        t+=dt
    
    Rho_MPCstep = copy.deepcopy(R)
    #print Rho_MPCstep
    return R,totalflows, Buffer,norm_u_dot,w_eta

def MPC_functional(u):
    global gamma, eps,junctions, n_p, n_c, MPC_flows, MPC_buffers, t0, Rho_init,dt,dx,vmax,tend_input, rhomax,sigma,c,M,RhoPlot, XPlot, totalflows
    Rho,flows,buffers_in_n_p,norm_u_dot,w_eta=dynamicEvolution(Rho_init,u,t0,n_p)
    #print w_eta
    val = -getOverallFlux(flows)
    #gamma = 0.
    #eps = 0.
    fval = val+gamma*norm_u_dot+eps*w_eta
    MPC_flows = flows[:n_c+1]
    MPC_buffers=[buffers_in_n_p[i][:n_c+1] for i in range(len(junctions))]   
    return fval

def MPC_main(_X,_Rho,network):
    global MPC_buffers, MPC_flows, Rho_MPCstep, t0,tend,n_p,n_c,junctions,Rho_init,dt, dx, vmax, tend_input, rhomax, simga, c,M,RhoPlot, XPlot, totalflows
    _t = np.arange(0,tend+dt,dt) #time for plots
    junctions = network2junctions(network,c,M)
    sum_of_incoming_roads = 0
    
    for junc in junctions:
        sum_of_incoming_roads+=len(junctions[junc]["in"])
    u_init=sum([[0 for i in range(sum_of_incoming_roads)] for i in range(n_p+1)],[])
    t0=0.
    Densities = [0 for i in range(len(_Rho))] 
    X = [0 for i in range(len(_X))]
    controls_sol = []
    for i in range(len(X)):
        X[i],Densities[i] = init(_X[i],_Rho[i],dx)
    Rho_init = Densities
    XPlot = X
    
    Adm = network2Adjescency(network)
    cnts = adjescency2constraint(Adm,n_p)
    bnds = sum_of_incoming_roads*(n_p+1)*[(0,1)]
    bnds = tuple(bnds)
    totalflux=0
    res_fluxes=[]
    feval = 0
    Rho_sol = [[] for i in range(len(X))] #initialize solution array of length #roads    
    while t0<tend:
        res = minimize(MPC_functional,u_init,bounds = bnds,constraints = cnts,method="SLSQP")
        if t0%50==0:
            print "current time: {}".format(t0)
        print res
        res_fluxes.extend(MPC_flows) 
        feval+=res.nfev
        updateBuffer_MPC(junctions)
        totalflux+=getOverallFlux(MPC_flows)                
        controls_sol.append(res.x[:sum_of_incoming_roads*(n_c+1)])
        Rho_init=[]
        for i in range(len(Rho_MPCstep)): #update initial rho for next MPC step
            Rho_init.append(Rho_MPCstep[i][n_c+1])
        t0+=n_c*dt+dt
        #update u_init to the a n_p-extension of the controls at time n_c
        u_init=[]
        u_init.extend(sum([controls_sol[-1][-sum_of_incoming_roads:].tolist() for i in range(n_p+1)],[]))
        print u_init
        for i in range(len(Rho_MPCstep)):
            Rho_sol[i].extend(Rho_MPCstep[i][:n_c+1])
        RhoPlot = Rho_sol
    #print "Buffer: {}".format(junctions[0]["Buffer"])
    control_vectors=computeControlVectors(controls_sol,network)
    return _t,X,Rho_sol,control_vectors,res_fluxes,feval
    #print Rho,X

def MPC_main_new(_X,_Rho,network):
    global MPC_buffers, MPC_flows, Rho_MPCstep, t0,tend,n_p,n_c,junctions,Rho_init,dt, dx, vmax, tend_input, rhomax, simga, c,M,RhoPlot, XPlot, totalflows
    _t = np.arange(0,tend+dt,dt) #time for plots
    junctions = network2junctions(network,c,M)
    sum_of_incoming_roads = 0
    
    for junc in junctions:
        sum_of_incoming_roads+=len(junctions[junc]["in"])
    u_init=sum([[0 for i in range(sum_of_incoming_roads)] for i in range(n_p+1)],[])
    t0=0.
    Densities = [0 for i in range(len(_Rho))] 
    X = [0 for i in range(len(_X))]
    controls_s = []
    for i in range(len(X)):
        X[i],Densities[i] = init(_X[i],_Rho[i],dx)
    Rho_init = Densities
    XPlot = X
    
    Adm = network2Adjescency(network)
    cnts = adjescency2constraint(Adm,n_p)
    bnds = sum_of_incoming_roads*(n_p+1)*[(0,1)]
    bnds = tuple(bnds)
    totalflux=0
    res_fluxes=[]
    feval = 0
    Rho_sol = [[] for i in range(len(X))] #initialize solution array of length #roads    
    while t0<tend:
        #fluxes,rho,controls=doOptimizationStep(u_init,bnds,cnts,feval,sum_of_incoming_roads,junctions)
        fluxes,rho,controls,feval=doOptimizationStep_updated(u_init,bnds,cnts,feval,sum_of_incoming_roads,junctions,t0)
        res_fluxes.extend(fluxes) 
        totalflux+=getOverallFlux(fluxes)  
        controls_s.append(controls)
        Rho_init=[]
        for i in range(len(rho)): #update initial rho for next MPC step
            Rho_init.append(rho[i][n_c+1])
        #update u_init to the a n_p-extension of the controls at time n_c
        u_init=[]
        u_init.extend(sum([controls_s[-1][-sum_of_incoming_roads:] for i in range(n_p+1)],[]))
        print u_init
        t0+=n_c*dt+dt
        for i in range(len(rho)):
            Rho_sol[i].extend(rho[i])
    RhoPlot = Rho_sol
    #print "Buffer: {}".format(junctions[0]["Buffer"])
    control_vectors=computeControlVectors(controls_s,network)
    return _t,X,Rho_sol,control_vectors,res_fluxes,feval
    #print Rho,X

def doOptimizationStep(u_init,bnds,cnts,feval,sum_of_incoming_roads,junctions):
    global MPC_flows,dt,n_c,RhoPlot,Rho_MPCstep,Rho_init
    res = minimize(MPC_functional,u_init,bounds = bnds,constraints = cnts,method="SLSQP")
    if t0%50==0:
        print "current time: {}".format(t0)
    print res
    feval+=res.nfev
    updateBuffer_MPC(junctions)                
    controls_sol=res.x[:sum_of_incoming_roads*(n_c+1)]
    controls_sol = controls_sol.tolist()
    Rho_init=[]
    for i in range(len(rho)): #update initial rho for next MPC step
        Rho_init.append(rho[i][n_c])

    #update u_init to the a n_p-extension of the controls at time n_c
    temp = np.array(Rho_MPCstep)
    Rho_s= temp[:,:n_c+2]
    Rho_sol = Rho_s.tolist()
    res_fluxes=MPC_flows

    return res_fluxes,Rho_sol,controls_sol

def doOptimizationStep_updated(u_init,bnds,cnts,feval,sum_of_incoming_roads,junctions,t0):
    global MPC_buffers,MPC_flows,dt,n_c,RhoPlot,Rho_MPCstep,Rho_init
    res = minimize(MPC_functional,u_init,bounds = bnds,constraints = cnts,method="SLSQP")
    if t0%50==0:
        print "current time: {}".format(t0)
    print res
    feval+=res.nfev
    updateBuffer_MPC(junctions)                
    controls_sol=res.x[:sum_of_incoming_roads*(n_c+1)]
    controls_sol = controls_sol.tolist()
    #controls_bin = [round(x,0) for x in controls_sol]  
    R,totalflows, Buffer,norm_u_dot,w_eta=dynamicEvolution(Rho_init,controls_sol,t0,n_c)    
    Rho_init=[]
    for i in range(len(R)): #update initial rho for next MPC step
        Rho_init.append(R[i][n_c])

    #update u_init to the a n_p-extension of the controls at time n_c
    temp = np.array(R)
    Rho_s= temp[:,:n_c+2]
    Rho_sol = Rho_s.tolist()
    res_fluxes=totalflows
    MPC_buffers=[Buffer[i][:n_c+1] for i in range(len(junctions))] 
    return res_fluxes,Rho_sol,controls_sol,feval

def computeFluxesFromControls(u,_X,_Rho,network):
    global MPC_buffers, MPC_flows, Rho_MPCstep, t0,tend,n_p,n_c,junctions,Rho_init,dt, dx, vmax, tend_input, rhomax, simga, c,M,RhoPlot, XPlot, totalflows
    _t = np.arange(0,tend+dt,dt) #time for plots
    junctions = network2junctions(network,c,M)
    sum_of_incoming_roads = 0
    
    for junc in junctions:
        sum_of_incoming_roads+=len(junctions[junc]["in"])
    Densities = [0 for i in range(len(_Rho))] 
    X = [0 for i in range(len(_X))]
    for i in range(len(X)):
        X[i],Densities[i] = init(_X[i],_Rho[i],dx)
    R = [[s] for s in Densities]
    totalflows = []
    roads = len(Rho_init) 
    XPlot = X
    flows = [0 for i in range(roads)]
    Buffer = [[] for i in range(len(junctions))]
    #Rho_sol = [[] for i in range(len(X))] #initialize solution array of length #roads    
    for timestep in range(len(_t)):
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
            #print column(u,timestep_in_MPC)
                
            def column(matrix, i):
                return [row[i] for row in matrix]
            jfluxes_in,jfluxes_out = getJunctionFluxesMPC(rho_in,rho_out,column(u,timestep),junctions[junc]) 
            #updateBuffer(junctions[junc],jfluxes_in,jfluxes_out)        
            Buffer[junc].append(computeBuffer_MPC(junctions[junc],jfluxes_in,jfluxes_out))       
            j=0    
            k=0
            for i in junctions[junc]["in"]:
                flows[i].append(jfluxes_in[j])
                j+=1
            for i in junctions[junc]["out"]:
                flows[i].insert(0,jfluxes_out[k])
                k+=1
        rowSum = np.sum(network,axis=1)
        columnSum = np.sum(network,axis=0)
        for i in range(roads):
            if rowSum[i]==0: #has external influx
                if _t[timestep]<tend_input:
                    flows[i].insert(0,f(ext_in[i]))
                else:
                    flows[i].insert(0,0)
            elif columnSum[i]==0: #has external outflow
                #flows[i].append(f(ext_outflow[i]))
                flows[i].append(flows[i][-1])
        ##4 godunov step
        #print "flows at road 0:",flows[0]
        #print "flows: ", flows
        totalflows.append(flows[:])
        MPC_buffers=Buffer
        updateBuffer_MPC(junctions)
        #print totalflows
        for i in range(roads):
            R[i].append(godunovStep(flows[i],R[i][-1]))
    
    #Rho_MPCstep = copy.deepcopy(R)
    #print Rho_MPCstep
    return R,totalflows

def roundControls(u):
    u_rounded = []
    for control in u:
        eta_rounded = [int(round(x,0)) for x in control]
        u_rounded.append(eta_rounded)
    return u_rounded
    
    
def adjescency2constraint(Adm,n_p):
    # returns necessary constrains for optimization of form \sum_{i in \E^{in}} u_i <=1
    #Adm adjescency matrix, n_p predictive horizon
    cnts = {}
    M = copy.deepcopy(Adm)
    if type(M).__module__ == np.__name__:
        M=M.tolist()
    temp = copy.deepcopy(M)
    step=1
    while step-1<n_p:   
        for i in range(len(M)):
            M[i].extend(len(temp[0])*[0])
        for i in range(len(temp)):
            M.append(step*len(temp[0])*[0]+temp[i])
        step+=1
    M=np.array(M)
    cnts=[]
    for i in range(len(M)):
        def fun(u,i=i):
            return 1-np.dot(M[i,:],u)
        cnts.append({"type": "ineq", "fun": fun})
    print cnts
    return tuple(cnts)
    
    
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
    

# plot functions
 
def plot2D(x,rho):
    global _X
    f,ax = plt.subplots(len(rho))
    plt.xlabel('x')
    #plt.tight_layout()
    f.set_size_inches(10, 8)
    #plt.yticks([0,0.5])
    #plt.ylim([0,0.8])
    for i in range(len(rho)):
        ax[i].plot(x[i][:-1],np.array(rho[i][-1][:-1]))  
        ax[i].set_ylabel(r'$\rho$')
        ax[i].set_yticks([0.0,0.25,0.5,0.75])
        ax[i].set_xlim(_X[i])
        
def plotBuffers(t,junctions):
    fig = plt.figure()
    plt.xlabel("t")
    print len(t)
    plt.ylabel("total buffers")
    legend = []
    i=0
    for junc in junctions:
        Buffer = junctions[junc]["Buffer"]
        #print junctions[junc]["matrix"]
        bufferSum = [sum(buf) for buf in Buffer]
        #print bufferSum
        print len(bufferSum)
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

def plotControls(control_vector,t,index):
    global gamma, eps,n_p, n_c
    print len(t),len(control_vector[index])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('time t')
    ax.set_ylabel(r'$\eta$')
    ax.plot(t,control_vector[index][:len(t)])
    fig.savefig("control_mpc_{}_{}_{}.eps".format(n_p,n_c,gamma),format="eps")
    
def computeControlVectors(all_controls,network):
    global c,M,t
    junctions = network2junctions(network,c,M)
    roads_incoming = []
    for junc in junctions:
        roads_incoming.extend(junctions[junc]["in"])
    number_incoming_roads=len(roads_incoming)
    control_vector = [[] for i in range(number_incoming_roads)]
    for i in all_controls:
        for control in range(len(i)):
            index = control%number_incoming_roads
            control_vector[index].append(i[control])
    return control_vector

def animate(i):
    global RhoPlot,XPlot,line,fig,ax
    ax.clear()    
    ax.plot(XPlot[1][:-2],RhoPlot[1][i][:-2])
    ax.set_ylim([0,1]) 
    return line,

def initAnim():
    global line
    line.set_data([],[])
    return line,

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
##obsolte



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
dt=0.5

totalflows = []
junctions = []
t0=0.
tend=500
tend_input= 200
t = np.arange(0,tend,dt)
dx=.5   
#M=[0.5]
M=[[0.5,0.5],[0.5,0.5]]
c=[[1.,1.],[1.,1.]]

MPC_flows_laststep=[]
MPC_flows=[]

MPC_buffers=[]



######## test cases

### test 1: very very small "network"

#_X = [[0,50],[50,150]]
#_Rho = [[0.],[0.0]]
#Rho_init = _Rho
#Rho_MPCstep = []
#network = [[0.,0.],[1.,0]]
#ext_in = [0.5,0]

### test 2: single junction

_X = [[0,50],[0,50],[50,100]]
_Rho = [[0.2],[0.2],[0.2]]
Rho_init = _Rho
Rho_MPCstep = []
network = [[0,0,0],[0,0,0],[1.,1.,0]]
ext_in = [0.4,0.2,0]
#

### test 3: two junctions

#_X = [[0,50],[0,50],[50,150],[100,150],[150,200]]
#_Rho = [[0.2],[0.2],[0.2],[0.2],[0.2]]
#ext_in = [0.5,0.2,0,0.2,0]
#Rho_init = _Rho
#Rho_MPCstep = []
#network = np.array([[0,0,0,0,0],[0,0,0,0,0],[1.,1.,0,0,0],[0,0,0,0,0],[0,0,1.,1.,0]])


# voariables for MPC

n_p =20
gamma = 10.
eps = 5.
n_c= 20
start = time.time()
_t,x,rho,u_vector,totalflux,feval=MPC_main_new(_X,_Rho,network)
elapsedTime = time.time()-start
print elapsedTime
#plot2D(x,rho)

#rr,ff=computeFluxesFromControls(u_vector,_X,_Rho,network)
#
##binary controls
#u_bin=roundControls(u_vector)
#rr2,ff2=computeFluxesFromControls(u_bin,_X,_Rho,network)
#flux_from_binary_controls = getOverallFlux(ff2)
#
plotControls(u_vector,t,0)
#plotControls(u_bin,0)
flux = getOverallFlux(totalflux)
#
## write data files
input_data = open("input_{}_{}.txt".format(n_p,n_c),"w")
output_data = open("output_{}_{}_{}.txt".format(n_p,n_c,gamma),"w")
input_data.write("np: {}\n n_c: {}\n X: {}\n Rho: {} gamma: {}\n eps: {}\n network: {}\n external inflow: {}\n input end: {}\n tend:".format(n_p,n_c,_X,_Rho,gamma, eps,network,ext_in,tend_input,tend))
#output_data.write("elapsed time: {}\n t: {}\n x: {}\n rho: {}\n fluxes: {}\n controls: {}\n binary controls: {}\n binary fluxes: {}\n feval: {}\n total flux: {}\n binary flux {}".format(elapsedTime,t,x,rho,totalflux,u_vector,u_bin,ff2,feval,flux,flux_from_binary_controls))
output_data.write("elapsed time: {}\n t: {}\n x: {}\n rho: {}\n fluxes: {}\n controls: {}\n feval: {}\n total flux: {}".format(elapsedTime,t,x,rho,totalflux,u_vector,feval,flux))

input_data.close()
output_data.close()