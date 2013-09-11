def likelihoodfunction(xdata,ydata,parameters):
	a = parameters[0]
	b = parameters[1]
	dim = len(xdata)
	chi2 = 0.0
	for i in range(dim):
		y = a*xdata[i] + b
		chi2 = chi2 + ((y-ydata[i])/(ydata[i]*0.05))**2
	return numpy.exp(-chi2/2.0)

#===================================================
import random
import numpy
import pylab

outputfile=open('chain.txt','w')	#output file
#====================================================
#reading input file
xdata = []
ydata = []
datafile = open('fort.11111','r')
while True:
	l = datafile.readline()
	if not l:
		break
	xdata.append(float(l.split()[0]))
	ydata.append(float(l.split()[1]))

aa=[]
bb=[]
#print xdata, ydata
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#input parameters
#parametersvalues = [[Mean, minimum, maximum, std.dev.],....]
parametersvalues = [[5.0,0.0,10.0,0.1],[10.0,0.0,20.0,0.1]]#,[10.0,0.0,20.0,0.1],[10.0,0.0,20.0,0.1]]
numberofparameters = 2		#number of parameters
numberofsteps = 100000		#number of steps
random.seed(216.0)		#random seed to start the chain, if this line is commented, random seen is taken from system clock
stepmultiplier = 7.0		#how big the steps should be (for each parameter)
burn_in=0			#how many initial accepted points to remove
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#generating initial point
parametersi = []
parametersf = []
for i in range(numberofparameters):
	parameterini = parametersvalues[i][1] + random.random() * (parametersvalues[i][2] - parametersvalues[i][1])
	parametersi.append(parameterini)	#initial point for each parameter
	parametersf.append(parameterini)	#new point for each parameter, just to make the array, not using here


# Calculate likelihood at parametersi as likelihoodold
likelihoodold = likelihoodfunction(xdata,ydata,parametersi)#0.5 #initial likelihood choosen randomly, need to replace

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
multiplier = 0			#how many times the chain come back to the point
acceptedpoints = 0		#how many points are accepted
acceptedratio =0.0		#Acceptance ratio of the chain
for i in range(numberofsteps):
	multiplier += 1
	for j in range(numberofparameters):	#generating new point in parameter space
		parametervalue = parametersi[j] + random.gauss(0.0,1.0) * parametersvalues[j][3] * stepmultiplier
		parametersf[j]=(parametervalue)

#Calculate likelihood at parametersf as likelihood new
	likelihoodnew = likelihoodfunction(xdata,ydata,parametersf)

	posteriorratio = likelihoodnew/likelihoodold	#ratio of new likelihood to the old
	randomnumber01 = random.random()		#random number from uniform distribution (0,1)

	if posteriorratio > randomnumber01:		#If the ratio>random number, the point is accepted
		acceptedpoints += 1
		acceptedratio = float(acceptedpoints)/(i+1)

		for k in range(numberofparameters):
			parametersi[k] = parametersf[k]
		likelihoodold = likelihoodnew
		
		if acceptedpoints > burn_in:
			print >>outputfile, likelihoodold,'\t',multiplier,'\t',parametersi[0],'\t',parametersi[1]#,parametersi[2],'\t',parametersi[3]	#print to file
			aa.append(parametersi[0])
			bb.append(parametersi[1])
	
		multiplier = 0
#Chain finished
print "Acceptance ratio is: ",acceptedratio

pylab.plot(aa,bb,'r')
pylab.xlabel('a [slope]')
pylab.ylabel('b [intercept]')
#pylab.show()
pylab.savefig('chain.pdf')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
