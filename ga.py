#fitness function need to be change for particular problem, right now just fitting straight line [imp: lower fitness is better]
def fitnessfunction(xdata,ydata,parameters):
	a = parameters[0]
	b = parameters[1]
	dim = len(xdata)
	chi2 = 0.0
	for i in range(dim):
		y = a*xdata[i] + b
		chi2 = chi2 + ((y-ydata[i])/(ydata[i]*0.05))**2
#	print a,b,chi2
	return chi2/11.0

#==============================================
#selection of parents
def selection(fitness,numberofpopulations):		#Here higher fitness values are accpted so taking inverse of individual fitness values
	fitnessinverse = []				#making a list for inverse of fitness	
	for i in range(numberofpopulations):		#doing inverse of individual fitness
		fitnessinverse.append(1./fitness[i])	
	sumfitnessinverse = numpy.sum(fitnessinverse)	#Taking sum of all fitness (inverse)
	randomsum = random.random()*sumfitnessinverse	#Generating random number betwee 0 and sum of inverse fitness
	
	fitnesscompare = 0.				
	for i in range(numberofpopulations):		#Comparing summing over individual fitnesses
		fitnesscompare = fitnesscompare + fitnessinverse[i]
		if fitnesscompare > randomsum:
			break				#when summing is greater than random number, that index is returned and chain is broken
#	print randomsum,fitnesscompare,i
	return i
#==============================================
#crossover only depends on parents
def crossover(parent1,parent2,numberofparameters,Pc):
	crossat = int(numberofparameters*Pc)		#Calculating at which index crossover to be done
	newgeneration = []				
	for i in range(numberofparameters):
		if (i+1) <= crossat:
			newgeneration.append(parent1[i])#Making new generation with the two parents
		else:
			newgeneration.append(parent2[i])


	return newgeneration
#==============================================
import numpy, random, time

numberofpopulations = 100	#number of chromosomes in each generation
random.seed(21506.0)		#random seed; if this line is commented, random seed is taken from system clock
numberofparameters = 2
parametersvalues = [[5.0,0.0,10.0],[10.0,0.0,20.0]]#,[10.0,0.0,20.0],[10.0,0.0,20.0]] 	#[mean, min., max.] for each parameter
fitness_limit = 0.001		#The best fitness to reach during the process, lower is strict.
outputfile = open('output_ga.txt','w')
Pc = 0.5			#Crossover probability
Pm = 0.05			#Mutation probability
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
#generating initial point
parametersi = numpy.zeros((numberofpopulations,numberofparameters))	#making empty lists
parametersf = numpy.zeros((numberofpopulations,numberofparameters))

parent1 = numpy.ones((numberofparameters))		#making empty lists for parents
parent2 = numpy.ones((numberofparameters))

fitness = numpy.ones((numberofpopulations))		#making empty lists for fitness of each population
for j in range(numberofpopulations):			#making first random points of parameters/chromosomes/solution and corresponding fitness
	for i in range(numberofparameters):
		parametersi[j][i] = parametersvalues[i][1] + random.random() * (parametersvalues[i][2] - parametersvalues[i][1])
		parametersf[j][i] = parametersvalues[i][1] + random.random() * (parametersvalues[i][2] - parametersvalues[i][1])
	fitness[j] = fitness[j] * random.random() * 1e20

meanfitness = numpy.mean(fitness)			#mean of fitness population
minfitness = numpy.amin(fitness)			#minimum of fitness population

#print meanfitness,minfitness
step = 0						#initiating generation number
minfitnesswrite = 1e20

print "Sampling started with generation 0"
while minfitness > fitness_limit:			#To do while minimum fitness matches the limit required
	step += 1
	for j in range(numberofpopulations):		#selecting parents from the whole population of the generation
		index1 = selection(fitness,numberofpopulations)
		index2 = selection(fitness,numberofpopulations)
		parent1 = parametersi[index1,:]
		parent2 = parametersi[index2,:]
	
		newgeneration = crossover(parent1,parent2,numberofparameters,Pc)	#Crossover
		for i in range(numberofparameters):	#mutation
			newgeneration[i] = newgeneration[i] * (1.0+Pm*random.gauss(0.0,1.0))
		
		parametersf[j,:] = newgeneration
		fitness[j] = fitnessfunction(xdata,ydata,newgeneration)	#calculatin fitness of new generation
	meanfitness = numpy.mean(fitness)
	minfitness = numpy.amin(fitness)
#	ind = fitness.index(minfitness)
	if minfitness<minfitnesswrite:
		print >> outputfile, 'Current best in generation:',step,'\t with best fitness', minfitness
		minfitnesswrite = minfitness
	print meanfitness,'\t',minfitness,'\t','Generation number: \t',step

for index in range(len(fitness)):		#finding the chromosome with minimum fitness out of the final generation population
	if fitness[index] == minfitness :
		print 'Best fitness is: ', minfitness,'with parameter values [a,b] = ', parametersf[index,:], 'found in generation' , step
		break
	
