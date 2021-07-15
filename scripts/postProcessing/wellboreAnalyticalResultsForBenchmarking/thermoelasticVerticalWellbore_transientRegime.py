from scipy.special import kn
from math import factorial, floor
import numpy as np
import matplotlib.pyplot as plt


def FFunction(s, R):
    P = kn(0,R * (s ** 0.5)) / (s * kn(0,s ** 0.5))
    U = -(R * kn(1,R * (s ** 0.5)) - kn(1,s ** 0.5)) / (s * R * (s ** 0.5) * kn(0,s ** 0.5))
    Srr = -(-R * kn(1,R * (s ** 0.5)) + kn(1,s ** 0.5)) / (R ** 2. * (s ** 1.5) * kn(0,s ** 0.5))
    Stt = -P - Srr
   
    return [P, U, Srr, Stt]


def Vfunction(i, N):
    
	sum1 = 0.
	kmin = int(floor((i + 1.) / 2.))

	kmax = min(i, N)
    
	for k in range(kmin, kmax+1):
		sum1 = sum1 + (1.*(k**N) * factorial(2 * k) / (factorial(N - k) * factorial(k) * factorial(k - 1) * factorial(i - k) * factorial(2 * k - i)))


	return ((-1.)**(N + i)) * sum1

def StehfestTransform(t, R):
	N = 5
	sum1 = 0.
	sum2 = 0.
	sum3 = 0.
	sum4 = 0.
	for j in range(1,2*N+1):
		Lresult = FFunction(j * np.log(2.) / t, R)
		sum1 = sum1 + Vfunction(j, N) * Lresult[0]
		sum2 = sum2 + Vfunction(j, N) * Lresult[1]
		sum3 = sum3 + Vfunction(j, N) * Lresult[2]
		sum4 = sum4 + Vfunction(j, N) * Lresult[3]

	return [sum1* np.log(2.) / t, sum2* np.log(2.) / t, sum3* np.log(2.) / t, sum4* np.log(2.) / t] 

t = 10.
a = 0.1
T0 = 100.

thermalStressCoeff = 1e4 #Pa/C
E = 20.6e9
nu = 0.189

K = E/3./(1.-2.*nu)
G = E/2./(1.+nu)

coefc = 0.1

tDimensionless = t*coefc/(a**2.)

listR = np.arange(1.,10.,0.01)
listP = []
listU = []
listSrr = []
listStt = []
for R in listR:

	result = StehfestTransform(tDimensionless, R)
	listP.append(result[0])
	listU.append(result[1])
	listSrr.append(result[2])
	listStt.append(result[3])

r_anal = [R*a for R in listR]
T_anal = [P*T0 for P in listP]

sigrr_anal = [1e-6*Srr*thermalStressCoeff*(1.-2.*nu)/(1.-nu)*T0 for Srr in listSrr]
sigtt_anal = [1e-6*Stt*thermalStressCoeff*(1.-2.*nu)/(1.-nu)*T0 for Stt in listStt]

fig = plt.figure(figsize=[13,10])

# Temperature
r_geosx, T_geosx = [], []
for line in open('temperature.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		rVal = values[0]
		TVal = values[1]

		r_geosx.append( rVal )
		T_geosx.append( TVal )

plt.subplot(221)
plt.plot(r_geosx, T_geosx, 'ko', label='GEOSX result')
plt.plot(r_anal, T_anal,  'r', linewidth=2, label='Analytic')
plt.ylabel('Temperature (C)')
plt.xlabel('Radial coordinate (m)')
plt.xlim(a,10*a)
plt.legend()

# Radial stress
sigrr_geosx = []
for line in open('radialStress.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		sigrrVal = values[1]
		sigrr_geosx.append( sigrrVal*1e-6 ) 

plt.subplot(223)
plt.plot(r_geosx, [val1-1e4*val2*1e-6 for val1, val2 in zip(sigrr_geosx, T_geosx) ], 'ko', label='GEOSX result')
plt.plot(r_anal, sigrr_anal,  'r', linewidth=2, label='Analytic')
plt.ylabel('Radial stress (MPa)')
plt.xlabel('Radial coordinate (m)')
plt.xlim(a,10*a)


# Tangent stress
sigtt_geosx = []
for line in open('tangentStress.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		sigttVal = values[1]
		sigtt_geosx.append( sigttVal*1e-6 ) 

plt.subplot(224)
plt.plot(r_geosx, [val1-1e4*val2*1e-6 for val1, val2 in zip(sigtt_geosx, T_geosx) ], 'ko', label='GEOSX result')
plt.plot(r_anal, sigtt_anal,  'r', linewidth=2, label='Analytic')
plt.ylabel('Tangent stress (MPa)')
plt.xlabel('Radial coordinate (m)')
plt.xlim(a,10*a)

plt.show()
