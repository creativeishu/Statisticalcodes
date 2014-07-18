import numpy
import pylab

x = numpy.genfromtxt('fort.11')
x1 = numpy.genfromtxt('fort.12')
x2 = numpy.genfromtxt('fort.13')
x3 = numpy.genfromtxt('fort.14')
x4 = numpy.genfromtxt('fort.15')
x5 = numpy.genfromtxt('fort.16')
x6 = numpy.genfromtxt('fort.17')
x7 = numpy.genfromtxt('fort.18')
x8 = numpy.genfromtxt('fort.19')

pylab.loglog(x[:,0],x[:,1])
pylab.loglog(x1[:,0],x1[:,1])
pylab.loglog(x2[:,0],x2[:,1])
pylab.loglog(x3[:,0],x3[:,1])
pylab.loglog(x4[:,0],x4[:,1])
pylab.loglog(x5[:,0],x5[:,1])
pylab.loglog(x6[:,0],x6[:,1])
pylab.loglog(x7[:,0],x7[:,1])
#pylab.semilogx(x8[:,0],x8[:,1]/x[:,1])

#pylab.ylim(0.5,1.5)
pylab.show()


'''
x = numpy.genfromtxt('test.output')
pylab.loglog(x[:,0],x[:,1],'r')
pylab.loglog(x[:,0],x[:,2],'g')
pylab.loglog(x[:,0],x[:,3],'b')
pylab.loglog(x[:,0],x[:,4],'k')

pylab.ylim(ymin=min(x[:,2]))
pylab.show()
'''