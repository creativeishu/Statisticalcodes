import numpy
from scipy import linalg
import pylab
import subprocess

pylab.ion()


lmax = 20000

filename2 = 'data/cl_lmax_%i_Mcrit_1d13_bar.txt'%lmax
filename1 = filename2.replace('bar','bar_unnoise')
invname = filename2.replace('cl','invcov')
print filename1,filename2,invname

cl = numpy.genfromtxt(filename1)
inverse = numpy.genfromtxt(invname)

spec = ['r','b','g','m','c','k']
label = [11,12,13,22,23,33]

for i in range(6):
	pylab.loglog(cl[:,0],cl[:,1+i],spec[i],lw=2,label=str(label[i]))

for i in range(len(cl)):
	covariance = linalg.inv(inverse[i*6:(i+1)*6,:])
#	print linalg.eigh(covariance*1e25,eigvals_only=True)
	delta = numpy.random.multivariate_normal(cl[i,1:7],covariance)
	cl[i,1:7] = delta
#	print i,cl[i,1:7]
numpy.savetxt(filename2,cl)

#c1 = 'cp %s %s'%(filename1,filename1+'temp')
#c2 = 'cp %s %s'%(filename2,filename2+'temp')
#c3 = 'mv %s %s'%(filename2+'temp',filename1)
#c4 = 'mv %s %s'%(filename1+'temp',filename2)

#subprocess.call(c1,shell=True)
#subprocess.call(c2,shell=True)
#subprocess.call(c3,shell=True)
#subprocess.call(c4,shell=True)

for i in range(6):
	pylab.loglog(cl[:,0],cl[:,1+i],'--'+spec[i],lw=2)

pylab.legend(loc=3,fontsize=18)
pylab.xlabel('$\mathtt{\ell}$',fontsize=22)
pylab.ylabel('$\mathtt{C_{\ell}}$',fontsize=22)
pylab.xlim(min(cl[:,0]),max(cl[:,0]))
pylab.ylim(ymax = 2e-7)
pylab.savefig(filename2.replace('txt','eps'))
pylab.show()
