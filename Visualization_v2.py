import matplotlib.pyplot as plt
import numpy as np
import astropy.coordinates as ap
projection = 'mollweide'
fig = plt.figure(figsize=(25,12.5))
ax = fig.add_subplot(111, projection=projection)
ax.grid(True)
a = np.loadtxt("o.dat",usecols=(0,1))
ax.plot(np.radians(a[:,0]),np.radians(a[:,1]),'*',label = "Galactic Cut (+-30 Degrees)")
ax.plot(np.radians(262.92103505720036-180),np.radians(47.67440904728711),'k*',label="True Dipole")
plt.title("Dipole Distribution")
ax.set_yticklabels([])
ax.set_xticklabels([])
plt.legend()
plt.savefig("DipoleDistribution.png")
plt.show()

a = np.loadtxt("output.dat",usecols=(4,5))
y,x = a[:,0],a[:,1]
x_min = np.min(x)
x_max = np.max(x)
y_min = np.min(y)
y_max = np.max(y)
x_bins = np.linspace(x_min, x_max, 30)
y_bins = np.linspace(y_min, y_max, 30)
fig = plt.subplots(figsize =(10, 7))
# Creating plot
plt.hist2d(x, y,bins =[x_bins, y_bins])
plt.title("Magnitude vs Offset")
plt.savefig("Histogram.png")
# show plot
plt.show()