from matplotlib import cm
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import pickle


rho = [[[0,1],[0.2,0.5]],[[1]]]
vmax_2 =50
vmax = 120
rhomax = 1
psi_0 = 2
dt = 0.5
sigma=0.5
dx = 10.
alpha=1
buffer_file = "buffer_test_v2.01.txt"
delay_file = 'buffer_delay_v2.01.txt'
mpc_file = 'output_buffer_dp_50_2_2_10.txt'

def vel(rho):
    global vmax_2
    return vmax_2*(1-rho)

def getFuelConsumptionByDensity(rho):
    global psi_0,vmax
    v= vel(rho)
#    emission = 25*math.exp(4.5-10*v)*v+100*pow(v,.5)+psi_0
    emission = 1.2*((3.5*math.exp(5-15*v/vmax)+8*v/vmax)*v/vmax+2)
    return emission

def getPollutionOnNetwork(rho_vector):
    global alpha,vmax
    pollution = []
    for road in range(len(rho_vector)):
        pollutionOnRoad = []
        for timestep in range(len(rho_vector[road])):
            pollutionOnRoad.append([alpha*0.01*(getFuelConsumptionByDensity(rho)*vel(rho)+psi_0)*rho for rho in rho_vector[road][timestep]])
        pollution.append(pollutionOnRoad)
    return pollution

def getBufferPollution(junctions):
    global dt,alpha
    totalBufferPollution = []
    junctionPollution =[]
    for junc in junctions:
        buffers = junctions[junc]["Buffer"]
        b = [sum(buf) for buf in buffers]
        junctionPollution=[alpha*0.01*(getFuelConsumptionByDensity(sigma)*vel(sigma)*bufferLevel) for bufferLevel in b]
        print junctionPollution
        totalBufferPollution.append(junctionPollution)
    return totalBufferPollution
 
def computePollutionOnNetwork(buffer_file,delay_file,mpc_file):
    data = {}
    with open(buffer_file,"rb") as handle:
        data = pickle.loads(handle.read())
    RhoPlot = data["rho"]
    t = data["t"]
    junctions_buffer = data["junctions"]
    
    data_var_lt = {}
    with open(delay_file, 'rb') as handle:
        data_var_lt = pickle.loads(handle.read())
    Rho_var_lt = data_var_lt["rho"]
    t_var_lt=data_var_lt["t"]
    junctions_lt = data_var_lt["junctions"]
    #
    data_mpc = {}
    with open(mpc_file, 'rb') as handle:
        data_mpc = pickle.loads(handle.read())
    RhoMPC = data_mpc["rho"]
    t_mpc=data_mpc["t"]
    junctions_mpc = data_mpc["junctions"]
    
    pollution=getPollutionOnNetwork(RhoPlot)
    pollutionOnRoad=[]
    bufferPollution = getBufferPollution(junctions_buffer)
    for j in range(len(pollution)):
        pollutionOnRoad.append([np.sum(pollution[j][i]) for i in range(len(pollution[j]))])
    for j in range(len(pollutionOnRoad)):
        plt.plot(t,pollutionOnRoad[j][:-1],label="road {}".format(j))
    plt.legend()
    plt.xlabel("time in s")
    plt.ylabel("pollution in mg/s")
    plt.title("No optimization")
    plt.savefig("pollution_noopt.eps", format="eps")
    totalBufferPollution = np.sum([np.sum(bufferPollution[i]) for i in range(len(bufferPollution))])
    total1 = np.sum([np.sum(pollutionOnRoad[i]) for i in range(len(pollutionOnRoad))])
    total1PlusBuffer=total1+totalBufferPollution
        
    plt.figure()
    pollution2 = getPollutionOnNetwork(Rho_var_lt)
    pollutionOnRoad2 = []
    bufferPollution_lt = getBufferPollution(junctions_lt)
    for j in range(len(pollution2)):
        pollutionOnRoad2.append([np.sum(pollution2[j][i]) for i in range(len(pollution2[j]))])
    for j in range(len(pollutionOnRoad2)):
        plt.plot(t_var_lt,pollutionOnRoad2[j][:len(t_var_lt)],label="road {}".format(j))
    plt.legend()
    plt.xlabel("time in s")
    plt.ylabel("pollution in mg/s")
    plt.title("Delay model")
    plt.savefig("pollution_delay.eps", format="eps")
    totalBufferPollution_lt = np.sum([np.sum(bufferPollution_lt[i]) for i in range(len(bufferPollution_lt))])
    total2 = np.sum([np.sum(pollutionOnRoad2[i]) for i in range(len(pollutionOnRoad2))])
    total2PlusBuffer = totalBufferPollution_lt+total2
    #    
    plt.figure()
    pollution3 = getPollutionOnNetwork(RhoMPC)
    pollutionOnRoad3 = []
    bufferPollution_mpc = getBufferPollution(junctions_mpc)
    for j in range(len(pollution3)):
        pollutionOnRoad3.append([np.sum(pollution3[j][i]) for i in range(len(pollution3[j]))])
    for j in range(len(pollutionOnRoad3)):
        plt.plot(t_mpc,pollutionOnRoad3[j][:len(t_mpc)],label="road {}".format(j))
    plt.legend()
    plt.xlabel("time in s")
    plt.ylabel("pollution in mg/s")
    plt.title("MPC")
    plt.savefig("pollution_mpc.eps", format="eps")
    total3 = np.sum([np.sum(pollutionOnRoad3[i]) for i in range(len(pollutionOnRoad3))]) 
    totalBufferPollution_mpc = np.sum([np.sum(bufferPollution_mpc[i]) for i in range(len(bufferPollution_mpc))])
    total3PlusBuffer = total3+totalBufferPollution_mpc
    return pollutionOnRoad, pollutionOnRoad2, pollutionOnRoad3, total1PlusBuffer,total2PlusBuffer,total3PlusBuffer

pollBuffer,pollDelay,pollMPC,totalBuffer,totalDelay,totalMPC=computePollutionOnNetwork(buffer_file,delay_file,mpc_file)
