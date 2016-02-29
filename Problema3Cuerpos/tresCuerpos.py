import numpy as np
import matplotlib.pylab as plt
from matplotlib.backends.backend_pdf import PdfPages

data=np.genfromtxt('tresCuerpos.dat');
dataRK=np.genfromtxt('tresCuerposRK.dat');
dataA=np.genfromtxt('tresCuerposA.dat');
E=np.genfromtxt('energia.dat');
ERK=np.genfromtxt('energiaRK.dat');

with PdfPages('tresCuerpos.pdf') as pdf:
    ############
    plt.figure()
    plt.scatter(data[:,0],data[:,1],s=0.01)
    plt.ylim([-2.5,2.5])
    plt.xlim([-3,3])
    plt.xlabel('q3')
    plt.ylabel('p3')
    plt.title('Integrador Simplectico')
    pdf.savefig()
    plt.close()

    plt.figure()
    plt.scatter(dataRK[:,0],dataRK[:,1],s=0.01)
    plt.ylim([-2.5,2.5])
    plt.xlim([-3,3])
    plt.xlabel('q3')
    plt.ylabel('p3')
    plt.title('Integrador Runge Kutta')
    pdf.savefig()
    plt.close()
    
    ############
    plt.figure()
    plt.scatter(data[:,0],data[:,1],s=0.1)
    plt.ylim([-0.6,0.6])
    plt.xlim([-2.4,-1.0])
    plt.title('Integrador Simplectico')
    plt.xlabel('q3')
    plt.ylabel('p3')
    pdf.savefig()
    plt.close()
    
    plt.figure()
    plt.scatter(dataRK[:,0],dataRK[:,1],s=0.1)
    plt.ylim([-0.6,0.6])
    plt.xlim([-2.4,-1.0])
    plt.xlabel('q3')
    plt.ylabel('p3')
    plt.title('Integrador Runge Kutta')
    pdf.savefig()
    plt.close()


    ############
    plt.figure()
    plt.scatter(data[:,0],data[:,1],s=0.1)
    plt.ylim([-0.3,0.3])
    plt.xlim([-0.2,0.7])
    plt.title('Integrador Simplectico')
    plt.xlabel('q3')
    plt.ylabel('p3')
    pdf.savefig()
    plt.close()    
   
    plt.figure()
    plt.scatter(dataRK[:,0],dataRK[:,1],s=0.1)
    plt.ylim([-0.3,0.3])
    plt.xlim([-0.2,0.7])
    plt.xlabel('q3')
    plt.ylabel('p3')
    plt.title('Integrador Runge Kutta')
    pdf.savefig()
    plt.close()
    
    ############
    plt.figure()
    plt.scatter(dataA[:,0],dataA[:,1],s=0.1)
    plt.ylim([-2.5,2.5])
    plt.xlim([-2,2])
    plt.xlabel('q3')
    plt.ylabel('p3')
    plt.title('Integrador Simplectico')
    pdf.savefig()
    plt.close()    
   
    ############
    
    plt.figure()
    plt.scatter(dataA[:,0],dataA[:,1],s=0.1)
    plt.ylim([-0.6,0.6])
    plt.xlim([0.2,0.7])
    plt.xlabel('q3')
    plt.ylabel('p3')
    plt.title('Integrador Simplectico')
    pdf.savefig()
    plt.close()
    
    ############
    
    plt.figure()
    plt.scatter(dataA[:,0],dataA[:,1],s=0.1)
    plt.ylim([-0.4,0.4])
    plt.xlim([-0.4,0.4])
    plt.xlabel('q3')
    plt.ylabel('p3') 
    plt.title('Integrador Simplectico') 
    pdf.savefig()
    plt.close() 
 
    ###########
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(E[::5000,0],E[::5000,1])
    plt.xlabel('t')
    plt.ylabel('E')   
    plt.title('Integrador Simplectico')
   
    plt.subplot(1,2,2)
    plt.plot(ERK[::5000,0],ERK[::5000,1])
    plt.xlabel('t')
    plt.ylabel('E')
    plt.title('Integrador Runge Kutta')
    pdf.savefig()
    plt.close()
