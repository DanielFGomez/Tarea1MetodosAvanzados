import numpy as np
import matplotlib.pyplot as plt
import os

datosLF = np.loadtxt('LaxFriedrichs_finalstate.dat')
datosUG = np.loadtxt('UpwindGodunov_finalstate.dat')
datosAN = np.loadtxt('Analytic.dat')

os.system('mkdir Plots/')

results=['Density','Velocity','Energy','Pressure']
i=1

for label in results:

	fig = plt.figure()
	ax = plt.subplot(111)

	ax.plot(datosLF[:,0],datosLF[:,i],label='Lax Fridrich method')
	ax.plot(datosUG[:,0],datosUG[:,i],label='Upwind Godunov method')
	ax.plot(datosAN[:,0],datosAN[:,i],label='Analytic solution')
	ax.set_ylim(0,np.amax(datosAN[:,i])*1.15)

	lgd = ax.legend(loc=3, bbox_to_anchor=(1.3, 0.5))
	ax.set_title(label+' profile after one unit of time')
	ax.set_xlabel('x')
	ax.set_ylabel(label) 
	# Put a legend below current axis

	plt.savefig(("Plots/"+label+".pdf"),bbox_extra_artists=(lgd,), bbox_inches='tight')
	i+=1

