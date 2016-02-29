import numpy as np
import matplotlib.pyplot as plt

datos0 = np.loadtxt('LaxFriedrichs_step_0.dat')
datos1 = np.loadtxt('LaxFriedrichs_step_1.dat')
datos2 = np.loadtxt('LaxFriedrichs_step_2.dat')
datos3 = np.loadtxt('LaxFriedrichs_step_3.dat')
datos4 = np.loadtxt('LaxFriedrichs_step_4.dat')
datos5 = np.loadtxt('LaxFriedrichs_step_5.dat')

results=['Density','Velocity','Energy','Pressure']
i=1

for label in results:

	fig = plt.figure()
	ax = plt.subplot(111)

	ax.plot(datos0[:,0],datos0[:,i],label='0')
	ax.plot(datos1[:,0],datos1[:,i],label='1')
	ax.plot(datos2[:,0],datos2[:,i],label='2')
	ax.plot(datos3[:,0],datos3[:,i],label='3')
	ax.plot(datos4[:,0],datos4[:,i],label='4')
	ax.plot(datos5[:,0],datos5[:,i],label='5')

	lgd = ax.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))
	ax.set_title(label+' profile')
	ax.set_xlabel('x')
	ax.set_ylabel(label) 
	# Put a legend below current axis

	plt.savefig((label+".pdf"),bbox_extra_artists=(lgd,), bbox_inches='tight')
	i+=1

