import numpy
import pylab

x = numpy.genfromtxt('fort.11')
y = numpy.genfromtxt('fort.12')
z = numpy.genfromtxt('fort.13')

pylab.semilogx(x[:,0],x[:,1]/x[:,1])
pylab.semilogx(y[:,0],y[:,1]/x[:,1])
pylab.semilogx(z[:,0],z[:,1]/x[:,1])

pylab.ylim(0.5,1.5)
pylab.show()