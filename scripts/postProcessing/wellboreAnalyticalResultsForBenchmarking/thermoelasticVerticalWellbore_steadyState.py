import numpy as np
import matplotlib.pyplot as plt

def temperature( Thole, Tinf, r, rhole, rinf ):
	return Thole + (Tinf - Thole) * ( np.log(r) - np.log(rhole) ) / ( np.log(rinf) - np.log(rhole) )

def intTemperature( Thole, Tinf, r, rhole, rinf ):
	val1 = Thole - (Tinf - Thole) * np.log(rhole) / ( np.log(rinf) - np.log(rhole) )
	val2 = (Tinf - Thole) / ( np.log(rinf) - np.log(rhole) )
	return val1 * r*r/2.0 + val2 * ( r*r/2.0*np.log(r) - r*r/4.0 )

def coeffC( Thole, Tinf, r, rhole, rinf ):
	return ( intTemperature( Thole, Tinf, r, rhole, rinf ) - intTemperature( Thole, Tinf, rhole, rhole, rinf ) )/r/r

rhole = 0.1
rinf = 2.0
Thole = 100.
Tinf = 0.

K = 5.5556e9
G = 4.16667e9
TEC = 1e4/K/3.0
E = 9.0*K*G / ( 3.0*K + G )
nu = ( 3.0*K - 2.0*G ) / ( 6.0*K + 2.0*G )
lam = K - 2.0/3.0*G

coeffChole = coeffC( Thole, Tinf, rhole, rhole, rinf )
coeffCinf = coeffC( Thole, Tinf, rinf, rhole, rinf )

rigidity = np.array([ [rinf, 1.0/rinf],
                      [2.0*(lam+G), -2.0*G/rhole/rhole] ])

force = np.array([ -(1.0+nu)/(1.0-nu)*TEC*rinf*coeffCinf, E/(1.0-nu)*TEC*coeffChole ])
vectorCoefficientAB = np.dot( np.linalg.inv( rigidity ), force )

coeffA = vectorCoefficientAB[0]
coeffB = vectorCoefficientAB[1] #TEC*(1.0+nu) / (1.0-nu) * coeffChole * rhole * rhole

r_anal = np.arange(rhole, 20.*rhole, 0.01*rhole)
T_anal = temperature(Thole, Tinf, r_anal, rhole, rinf)
sigzz_anal = 2.0*lam*coeffA - E*TEC/(1.0-nu) * T_anal
sigrr_anal = 2.0*(lam+G)*coeffA - 2.0*G*coeffB /r_anal /r_anal  - E*TEC/(1.0-nu) * coeffC( Thole, Tinf, r_anal, rhole, rinf )
sigtt_anal = 2.0*(lam+G)*coeffA + 2.0*G*coeffB /r_anal /r_anal  + E*TEC/(1.0-nu) * ( coeffC( Thole, Tinf, r_anal, rhole, rinf ) - T_anal )

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
#plt.xlim(rhole,10*rhole)

# Vertical stress
sigzz_geosx = []
for line in open('verticalStress.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]

		sigzzVal = values[1]
		sigzz_geosx.append( sigzzVal*1e-6 ) 

plt.subplot(222)
plt.plot(r_geosx, [val1-1e4*val2*1e-6 for val1, val2 in zip(sigzz_geosx, T_geosx) ], 'ko', label='GEOSX result')
plt.plot(r_anal, sigzz_anal*1e-6,  'r', linewidth=2, label='Analytic')
plt.ylabel('Vertical stress (MPa)')
plt.xlabel('Radial coordinate (m)')
#plt.xlim(rhole,10*rhole)

# Radial stress
sigrr_geosx = []
for line in open('radialStress.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		sigrrVal = values[1]
		sigrr_geosx.append( sigrrVal*1e-6 ) 

plt.subplot(223)
plt.plot(r_geosx, [val1-1e4*val2*1e-6 for val1, val2 in zip(sigrr_geosx, T_geosx) ], 'ko', label='GEOSX result')
plt.plot(r_anal, sigrr_anal*1e-6,  'r', linewidth=2, label='Analytic')
plt.ylabel('Radial stress (MPa)')
plt.xlabel('Radial coordinate (m)')
#plt.xlim(rhole,10*rhole)


# Tangent stress
sigtt_geosx = []
for line in open('tangentStress.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		sigttVal = values[1]
		sigtt_geosx.append( sigttVal*1e-6 ) 

plt.subplot(224)
plt.plot(r_geosx, [val1-1e4*val2*1e-6 for val1, val2 in zip(sigtt_geosx, T_geosx) ], 'ko', label='GEOSX result')
plt.plot(r_anal, sigtt_anal*1e-6,  'r', linewidth=2, label='Analytic')
plt.ylabel('Tangent stress (MPa)')
plt.xlabel('Radial coordinate (m)')
#plt.xlim(rhole,10*rhole)

plt.show()
