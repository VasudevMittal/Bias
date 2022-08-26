#visualizing the simulated dipoles.
import matplotlib.pyplot as plt
import numpy as np
import astropy.coordinates as ap
projection = 'mollweide'
fig = plt.figure(figsize=(25,12.5))
ax = fig.add_subplot(111, projection=projection)
ax.grid(True)
a = np.loadtxt("o.dat",usecols=(0,1))
ax.plot(np.radians(a[:,0]),np.radians(a[:,1]),'*',label = "Galactic Cut (+-30 Degrees)")
ax.plot(np.radians(262.92103505720036-360),np.radians(47.67440904728711),'k*',label="True Dipole")
plt.title("Dipole Distribution")
ax.set_yticklabels([])
ax.set_xticklabels([])
plt.legend()
plt.savefig("DipoleDistribution.png")
plt.show()