import numpy,pylab
from scipy.stats import kde
import matplotlib.pyplot as plt

burn_in = 500
data = numpy.genfromtxt('chain.txt',skip_header=burn_in)

p1 = 2
p2 = 3
dim = 20

p1min = min(data[:,p1])
p1max = max(data[:,p1])
p2min = min(data[:,p2])
p2max = max(data[:,p2])

p1mean = numpy.mean(data[:,p1])
p2mean = numpy.mean(data[:,p2])

p1std = numpy.std(data[:,p1])
p2std = numpy.std(data[:,p2])

print "Mean, St.dev., Minimum, Maximum"
print p1mean,p1std,p1min,p1max
print p2mean,p2std,p2min,p2max

p1array = numpy.linspace(p1min,p1max,dim)
p2array = numpy.linspace(p2min,p2max,dim)



counts = numpy.zeros((dim,dim))

for i in range(dim-1):
	for j in range(dim-1):
		for k in range(len(data)):
			if data[k,p1]>p1array[i] and data[k,p1]<p1array[i+1] and \
				 data[k,p2]>p2array[j] and data[k,p2]<p2array[j+1]:
				 counts[i,j] += 1

pylab.contourf(p1array,p2array,counts,(150,200))
pylab.colorbar()
pylab.show()




