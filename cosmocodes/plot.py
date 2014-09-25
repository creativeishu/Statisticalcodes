import numpy
import pylab
'''
spec = ['r','b','g','m','c','y','k']
bins = ['1-1','1-2','1-3','2-2','2-3','3-3']
x = numpy.genfromtxt('plot_data/cl_bins_dmo.txt')
y = numpy.genfromtxt('plot_data/cl_bins_Mcrit_1e13_AC.txt')
z = numpy.genfromtxt('plot_data/cl_bins_Mcrit_1e13_noAC.txt')

pylab.semilogx(x[:,0],x[:,1]/x[:,1]-1,'k',lw=2)
for j in range(6):
		pylab.semilogx(x[:,0],y[:,1+j]/x[:,1+j]-1,spec[j],lw=2,\
				label='$\mathtt{Bin:\ %s}$'%bins[j])
		pylab.semilogx(x[:,0],z[:,1+j]/x[:,1+j]-1,'--'+spec[j],lw=2)
pylab.axvline(x=1000,color='k',ls='--')
pylab.axvline(x=5000,color='k',ls='--')
pylab.axvline(x=10000,color='k',ls='--')
#pylab.axhline(y=-0.0462/0.279,color='k',ls=':',lw=3,label='$\mathtt{\Omega_b/\Omega_m}$')
pylab.legend(loc=3,fontsize=16)
pylab.ylim(-0.1,0.05)
pylab.xlim(min(x[:,0]),max(x[:,0]))
pylab.xlabel('$\mathtt{\ell}$',fontsize=22)
pylab.ylabel('$\mathtt{C_{\ell}/C_{\ell}^{DMO}-1\ [M_{crit}=1e13 M_{sun}/h]}$',fontsize=22)
pylab.savefig('plots/cl_ratio_bins_Mcrit_1e13.eps')
pylab.show()

'''

'''
spec = ['r','b','g','m','c','k','y']
x = numpy.genfromtxt('fort.50',skip_footer=1)
j=0
pylab.semilogx(x[:,0],x[:,1+j]/x[:,1+j]-1,'k',lw=2)

for i in range(7):
	y = numpy.genfromtxt('fort.%i'%(51+i),skip_footer=1)
	pylab.semilogx(x[:,0],y[:,1+j]/x[:,1+j]-1,spec[i],lw=2,\
			label='$\mathtt{M_{crit}=%1.1e\ M_{sun}/h}$'%10**(i+9))
	z = numpy.genfromtxt('fort.%i'%(61+i),skip_footer=1)
	pylab.semilogx(x[:,0],z[:,1+j]/x[:,1+j]-1,'--'+spec[i],lw=2)
#pylab.axvline(x=1000,color='k',ls='--')
#pylab.axvline(x=5000,color='k',ls='--')
#pylab.axvline(x=10000,color='k',ls='--')
pylab.axhline(y=-0.0462/0.279,color='k',ls=':',lw=3,label='$\mathtt{\Omega_b/\Omega_m}$')
pylab.legend(loc=2,fontsize=16)
pylab.ylim(-0.2,0.5)
pylab.xlim(min(x[:,0]),max(x[:,0]))
pylab.xlabel('$\mathtt{k \ [h/Mpc]}$',fontsize=22)
pylab.ylabel('$\mathtt{P(k)/P(k)[DMO]-1\ [z=4.0]}$',fontsize=22)
#pylab.savefig('plots/pk_ratio_z4.0.eps')
pylab.show()
'''
'''
spec = ['r','b','g','m','c','y','k']
rvir = [3.512e-2,7.5658e-2,0.163,0.351173,0.75679,1.63,3.5117298]
for i in range(7):
	x = numpy.genfromtxt('fort.%i'%(51+i))
	pylab.loglog(x[:,0]/rvir[i],x[:,1],spec[i],lw=2,label='$\mathtt{M_{halo} = 10^{%i} M_{sun}/h}$'%(9+i))
	pylab.loglog(x[:,0]/rvir[i],x[:,2],'--'+spec[i],lw=2)
pylab.legend(loc=1,fontsize=12)
pylab.xlim(0.01,1)
pylab.ylim(ymax=2e17)
pylab.xlabel('$\mathtt{Radius/Rvir}$',fontsize=22)
pylab.ylabel('$\mathtt{\\rho\ [h^2 M_{sun}/Mpc^3]\ [z=0.0]}$',fontsize=22)
pylab.savefig('rho_nfw_gas_z0.eps')
pylab.show()
'''
'''
redshift = [0.0,0.4,0.8,1.2,1.6,2.0]
spec = ['r','b','g','m','c','k']
for i in range(6):
	x = numpy.genfromtxt('fort.%i'%(61+i))
	pylab.loglog(x[:,0],x[:,1]/x[:,0],spec[i],lw=2,label='$\mathtt{z=%1.1f}$'%redshift[i])
pylab.legend(loc=1,fontsize=16)
pylab.xlabel('$\mathtt{M_{halo}\ [M_{sun}/h]}$',fontsize=22)
pylab.ylabel('$\mathtt{M_{CentralGalaxy}/M_{halo}}$',fontsize=22)
#pylab.ylim(ymin=1e-4)
pylab.tight_layout()
pylab.savefig('plots/bcg_mhalo.eps')
pylab.show()
'''
'''
spec = ['r','b','g','m','c','y','k']

for i in range(7):
	x = numpy.genfromtxt('fort.%i'%(61+i))
	pylab.loglog(x[:,0],x[:,1],spec[i],lw=2,label='$\mathtt{M_{crit} = 10^{%i} M_{sun}/h}$'%(10+i))
pylab.axhline(y=0.0462/0.279,color='k',ls='--',lw=2,label='$\mathtt{\Omega_b/\Omega_m}$')
pylab.legend(loc=4,fontsize=14)
pylab.ylim(1e-4,0.5)
pylab.xlim(1e9,1e17)
pylab.xlabel('$\mathtt{M_{halo}\ [M_{sun}/h]}$',fontsize=22)
pylab.ylabel('$\mathtt{f_{gas}}$',fontsize=22)
pylab.tight_layout()
pylab.savefig('plots/fgas_mhalo.eps')
pylab.show()
'''

redshift = [0.0,0.4,0.8,1.2,1.6,2.0]
spec = ['r','b','g','m','c','k']
for i in range(6):
	x = numpy.genfromtxt('fort.%i'%(61+i))
	pylab.loglog(x[:,0],x[:,1],spec[i],lw=2,label='$\mathtt{z=%1.1f}$'%redshift[i])
pylab.axhline(y=4.0,color='k',ls='--',lw=2,label='$\mathtt{Conc=4}$')
pylab.legend(loc=1,fontsize=16)
pylab.xlabel('$\mathtt{M_{halo}\ [M_{sun}/h]}$',fontsize=22)
pylab.ylabel('$\mathtt{Concentration}$',fontsize=22)
pylab.ylim(3,20)
pylab.xlim(1e10,1e16)
pylab.tight_layout()
pylab.savefig('plots/conc_mhalo.eps')
pylab.show()

