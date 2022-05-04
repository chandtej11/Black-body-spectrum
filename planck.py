import matplotlib.pyplot as plt
import numpy as np
import math

T = 10000 # IN KELVIN
h = 6.63*10**-27 # in cgs unit
c = 3*10**10 # in cm units
K = 1.38*10**-16
A = 3*10**15 # in erg/(cm^3 s Hz)
l = np.arange(1e-5, 7e-5, 0.1e-6) #l for wavelength

S = (np.exp(h*c/(l*K*T))-1)**-1*2*h*c**2/(l**5)
plt.yscale('log')
plt.xscale('log')
plt.xlabel("$\lambda$ (cm)")
plt.ylabel('Intensity(I)')

plt.xlim([1e-5, 7e-5])
plt.ylim([7.0e13, 3.4e16])

t01 = 10**-2
t02 = 10**-1
t03 = 1
t04 =10
t05 = 100
t06 = 300
C = (l/(5*10**-5))**-1.5
t1 = (t01*(l/(5*10**-5))**3)*(1-np.exp(- h*c/(K*l*T)))
t2 = (t02*(l/(5*10**-5))**3)*(1-np.exp(- h*c/(K*l*T)))
t3 = (t03*(l/(5*10**-5))**3)*(1-np.exp(- h*c/(K*l*T)))
t4 = (t04*(l/(5*10**-5))**3)*(1-np.exp(- h*c/(K*l*T)))
t5 = (t05*(l/(5*10**-5))**3)*(1-np.exp(- h*c/(K*l*T)))
t6 = (t06*(l/(5*10**-5))**3)*(1-np.exp(- h*c/(K*l*T)))
I1 = A*C*np.exp(-t1) + S*(1-np.exp(-t1))
I2 = A*C*np.exp(-t2) + S*(1-np.exp(-t2))
I3 = A*C*np.exp(-t3) + S*(1-np.exp(-t3))
I4 = A*C*np.exp(-t4) + S*(1-np.exp(-t4))
I5 = A*C*np.exp(-t5) + S*(1-np.exp(-t5))
I6 = A*C*np.exp(-t6) + S*(1-np.exp(-t6))


plt.plot(l,I1, color="brown", linestyle="--", lw=1.5, label="Optical depth= $10^{-2}$")
plt.plot(l,I2, color="b", linestyle="--", lw=1.5, label="Optical depth=$10^{-1}$")
plt.plot(l,I3, color="r", linestyle="--", lw=1.5, label="Optical depth=1")
plt.plot(l,I4, color="c", linestyle="--", lw=1.5, label="Optical depth=10")
plt.plot(l,I5, color="g", linestyle="--", lw=1.5, label="Optical depth=100")
plt.plot(l,I6, color="y", linestyle="--", lw=1.5, label="Optical depth=300")
plt.plot(l,S, color="K", linestyle="-", lw=2.0, label = "Blackbody spectrum")
plt.savefig('test.pdf')
plt.title("Radiation transfer equation compared with blackbody spectrum")
plt.legend()
plt.savefig("planck_radiation.pdf")
plt.show()
