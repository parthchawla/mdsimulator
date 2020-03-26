import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
################################################################################
## Get data ####################################################################
################################################################################
file1 = 'U.txt'
file2 = 'K.txt'
file3 = 'E.txt'
file4 = 'P.txt'
file5 = 'sE.txt'
file6 = 'sP.txt'
with open(file1) as f:
    U = np.loadtxt(f)
with open(file2) as f:
    K = np.loadtxt(f)
with open(file3) as f:
    E = np.loadtxt(f)
with open(file4) as f:
    P = np.loadtxt(f)
with open(file5) as f:
    sE = np.loadtxt(f)
with open(file6) as f:
    sP = np.loadtxt(f) 
################################################################################
## Plot energy against time ####################################################
################################################################################
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Potential and Kinetic Energy Vs Time")
ax1.set_xlabel('Timesteps')
ax1.set_ylabel('Magnitude')
line1, = ax1.plot(U, c='r', label='U')
line2, = ax1.plot(K, c='b', label='K')
plt.gca().add_artist(ax1.legend(handles=[line2], loc=1, prop={'size':12}))
ax1.legend(handles=[line1], loc=4, prop={'size':12})
fig1 = plt.figure()
fig1.subplots_adjust(hspace=0.5)
ax2 = fig1.add_subplot(211)
ax2.set_title("Total Energy Vs Time")
ax2.set_xlabel('Timesteps')
ax2.set_ylabel('Magnitude')
ax2.plot(E, c = 'purple')
################################################################################
## Plot pressure ###############################################################
################################################################################
ax3 = fig1.add_subplot(212)
ax3.set_title("Pressure Vs Time")
ax3.set_xlabel('Timesteps')
ax3.set_ylabel('Magnitude')
ax3.plot(P, c = 'green')
################################################################################
## Energy analysis #############################################################
################################################################################
meanu = np.mean(U)
meank = np.mean(K)
variancek = np.var(K)
sigmak = np.sqrt(variancek)
mean1 = np.mean(E)
variance1 = np.var(E)
sigma1 = np.sqrt(variance1)
print "\nPotential Energy Mean = ", meanu
print "Kinetic Energy Mean = ", meank
print "Total Energy Mean = ", mean1
print "Total Energy Variance = ", variance1
print "Total Energy Standard Deviation = ", sigma1
################################################################################
## Plot energy histogram #######################################################
################################################################################
fig2 = plt.figure()
fig2.subplots_adjust(hspace=0.5)
ax4 = fig2.add_subplot(211)
ax4.hist(E, bins=100, normed="True", facecolor='purple')
ax4.set_title("Total Energy Histogram")    
ax4.set_xlabel('Magnitude')
ax4.set_ylabel('Frequency')
plt.gca().set_xlim([min(E), max(E)])
x = np.linspace(min(E), max(E), 100) # returns 100 evenly spaced numbers
ax4.plot(x, mlab.normpdf(x, mean1, sigma1), c = 'yellow', label = 'Gaussian')
ax4.legend(prop={'size':12})
fig4 = plt.figure()
ax8 = fig4.add_subplot(111)
ax8.hist(K, bins=30, normed="True", facecolor='blue')
ax8.set_title("Temperature Histogram")    
ax8.set_xlabel('Magnitude')
ax8.set_ylabel('Frequency')
plt.gca().set_xlim([min(K), max(K)])
x = np.linspace(min(K), max(K), 30) # returns 50 evenly spaced numbers
ax8.plot(x, mlab.normpdf(x, meank, sigmak), c = 'yellow', label = 'Gaussian')
ax8.legend(prop={'size':12})
################################################################################
## Pressure analysis ###########################################################
################################################################################
mean2 = np.mean(P)
variance2 = np.var(P)
sigma2 = np.sqrt(variance2)
print "\nPressure Mean = ", mean2
print "Pressure Variance = ", variance2
print "Pressure Standard Deviation = ", sigma2
################################################################################
## Plot pressure histogram #####################################################
################################################################################
ax5 = fig2.add_subplot(212)
ax5.hist(P, bins=100, normed="True", facecolor='green')
ax5.set_title("Pressure Histogram")
ax5.set_xlabel('Magnitude')
ax5.set_ylabel('Frequency')
plt.gca().set_xlim([min(P), max(P)])
x = np.linspace(min(P), max(P), 100) # returns 100 evenly spaced numbers
ax5.plot(x, mlab.normpdf(x, mean2, sigma2), c = 'yellow', label = 'Gaussian')
ax5.legend(prop={'size':12})
################################################################################
## Plot s-values ###############################################################
################################################################################
x = [1,2,5,10,20,50,100,200,500,625,1000,1250,2000,3125,5000,6250,10000,12500,
20000,25000]
x = np.sqrt(x)
fig3 = plt.figure()
fig3.subplots_adjust(hspace=0.5)
ax6 = fig3.add_subplot(211)
ax6.set_title("s values for Energy")
ax6.set_xlabel('$t_b^{1/2}$')
ax6.set_ylabel('$s$')
ax6.scatter(x, sE)
ax7 = fig3.add_subplot(212)
ax7.set_title("s values for Pressure") 
ax7.set_xlabel('$t_b^{1/2}$')
ax7.set_ylabel('$s$')
ax7.scatter(x, sP)
plt.show()