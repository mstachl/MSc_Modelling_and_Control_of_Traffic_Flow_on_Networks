README file, written 2017-04-12
*******************************************************
How to use the traffic flow modelling toolbox?
##########################################################

The toolbox permits traffic flow simulation on arbitrary networks using three 
different approaches:
	- buffer_model.py: 			buffer model by Bressan without optimization 
	- buffer_mode_delay.py: 		optimisation via fixed strategy for traffic lights
	- MPC_fixed_decision_points: 	optimisation using Model Predictive Control and 
	fixed minimum green phase

The toolbox also permits computation of generated pollution on the network 
(without temporal development) via post-processing of the traffic data:
	- pollution.py

Link: www.github.com/mstachl/Traffic_flow_optimization/Experiments

##########################################################
1. Traffic Flow Modelling:

1.1 Parameters

	In all functions, several parameters have to be adjusted:
		- dt,dx: step sizes in time and space (in s, m) of type double
		- vmax: maximum allowed velocity on the road (in m/s) of type double

		Important: check that CFL condition is fulfilled

		- rhomax: maximum allowed density on the road of type double
		- tend: final simulation time (in s) of type Integer
		- network: (n x n)-traffic distribution matrix, n number of roads, ndarray
		- _X: 2d vector consisting of starting and end point of the roads and points of 
density discontinuities of initial data
		- _Rho: vector of lists with initial densities on the road, 
dim(_Rho[i])=dim(_X[i])-1 for all i=1…n
		-ext_in: n-dim vector of Doubles with incoming densities into the roads
		- M: n_J-dim array of lists of buffer capacities for all junctions, one list 
consists of buffer capacities for all outgoing roads separately
		- c: n_J dim array of priority values for incoming roads at junctions
		
		only applicable for MPC:
		- gamma: cost of control changes
		- eps: weight of quantisation by double-well
		- n_p: Integer number of predictive signal intervals
		- n_c: Integer number of control signal intervals
		- signal_horizon: number of time steps in one signal interval


1.2 Output

	The pickle-module is used for storing output. The resulting dictionaries are s
tored in non-human-readable txt files and need to be imported first 
for further processing.
	The output dictionary consists of the following elements:
	- “t”: time steps
	- “x”: space discretization
	- “rho”: densities at time steps at every road
	- “fluxes”: resulting Godunov fluxes at time steps at every road
	- “total flux”: sum of all fluxes
	- “total time”: total computational time needed
	- “junctions”: dictionary of numbered junctions with following elements:
			# “in”: list of indexes of incoming roads
			# “out”: list of indexes of outgoing roads
			# “matrix”: traffic distribution matrix of this junction, (n_in x n_out)-dim
			# “_M”: buffer capacity
			# “_c”: priority vector for incoming roads
			# “Buffer”: list of filling rates of the buffer at time steps
			only for MPC and fixed strategy optimisation:
			# “eta”: initial control
	Only for MPC and fixed strategy:
	- “controls”: control vector
	Only for MPC
	- “feval”: number of function evaluations

1.2.1 Buffer model:

	Output file is named “buffer_data.txt”

1.2.2 Fixed strategy:

	Output file is named “buffer_delay.txt”

1.2.3 MPC:

	Output file is named 'output_buffer_dp_{}_{}_{}_{}.txt'
     .format(signal_horizon,n_p,n_c,gamma)

1.3 Import

	To import a generated text file of name filename use pickle module in the following:

	with open(filename, ‘rb’) as handle:
    		data = pickle.load(handle)

	The data can be accessed by dict commands, e.g. X= data[“x”]

1.4 Functions

1.4.1 General Functions

	- velo(rho): 							computes density-dependent velocity
	- f(rho): 								fundamental diagram
	- init(x,rho,dx): 							generates initial setup on the network 
	- godunovFlux(f,rhol,rhor):					computes flux between two cells
	- godunovFluxes(f, rho):					computes Godunov fluxes on a road
	- getMaxInflux(f, rho): 					computes maximum influx into junction based on 
equation 2.3.5a
	- getMaxOutflux(f, rho):					computes maximum outflux of junction based on 
equation 2.3.5b
	- getIncomingFluxesSB(Buffer, maxIn, junction):	computes junction influxes by 
current buffer level using SBJ
	- getIncomingFluxesMB(Buffer, maxIn,junction): computes junction influxes by 
current buffer level using MBJ
	- getOutgoingFluxes(Buffer,maxOut,incomingFluxes): computes junction outfluxes
	- updateBuffer(junction,incomingFluxes,outgoingFluxes): computes changes in the 
buffer level using junction fluxes
	- godunovStep(fluxes,rho): computes one Godunov step based on equations 3.1.12
	- network2junctions(network, c,M,mode): builds junction dictionary based on 
transition matrix
	- network2Adjescency(network): builds adjacency matrix based on transition matrix
	- getJunctionFluxes(Buffer,rho_in,rho_out,mode): computes incoming and outgoing 
junction fluxes
	- getOverallFlux(fluxes): computes total flux based on flux array
- plot2D(x,rho): plots densities at time tend on all roads
- plot2DatTime(x, rho,s): plots densities at times t in s, s list of time points
- plotControls(control,t,index): plots controls of junction index 
- plotBuffers(t, junctions): plots buffer levels over time for all junctions
- plotTotalDensity(t,rho): plots total density on the network over time
- doAnimationNew(x, rho,filename): produces .mp4 file of the dynamics	

1.4.2 Buffer specific functions

- solveArbitraryNetwork(tend, network, c,M, initX,initRho,ext_in): main function

1.4.3 Fixed strategy specific functions

- getJunctionFluxesWithDelay(rho_in,rho_out,junction,tau,t,mode): computes junction
fluxes based on a signal delay of tau
- solveNetworkWithDelay(tend, network,_c,_M,initX,initRho,ext_inflow): main function

1.4.4 MPC specific functions

	- dynamicEvolution(rho_init,control,t0,n_p):	simulated the dynamics of the system 
over the predictive horizon
	- MPC_functional(u,t0):					functional used for minimization
	- adjescency2constraint: computes inequality constraints for controls based
on feasibility condition 4.0.1
- solveNetworkWithMPC(_X,_Rho,network): main function


2. Pollution modeling

2.1 Parameters

- buffer_file: path to output file of the buffer model
- delay_file: path to output file of fixed strategy model
- mpc_file: path to output file of mpc model
- vmax: maximum velocity for pollution plot
- vmax_2: maximum velocity allowed on network
- psi_0: fuel consumption in idle mode (in l/h)

2.2 Functions

- vel(rho): computes velocity from density
- getFuelConsumptionByDensity(rho): computes fuel consumption given the density at 
some point
- getPollutionOnNetwork(rho_vector): computes pollution given the density on the network
- getBufferPollution(junctions): computes pollution generated at the junctions

2.3 Output

- pollution distribution on roads
- total pollution


		