import numpy as np
import matplotlib.pyplot as plt

# Rotate stress from global coordinates system to the local coordinates of an inclined borehole
# See the description in fig.1 in Abousleiman and Cui 1998
def stressRotation(stress,phi_x,phi_z):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x),0.],[-np.sin(phi_x), np.cos(phi_x),0.],[0.,0.,1.]])
    rotz = np.array([[np.cos(phi_z),0., np.sin(phi_z)],[0.,1.,0.],[-np.sin(phi_z),0., np.cos(phi_z)]])

    return np.dot(np.dot(np.transpose(rotz),np.dot(np.dot(np.transpose(rotx),stress),rotx)),rotz)


a = 0.1
p0 = -10.0

K = 5.5556e9
G = 4.16667e9

#Deviation angles
phi_x = 0.
phi_z = 45./180.*3.1416

fig = plt.figure(figsize=[13,10])

# Get radial coordinate and compute analytical results
r, sig_rr_anal, sig_tt_anal = [], [], []
for line in open('stress_11.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		rval = values[0]
		sig_rr_val = p0*a*a/rval/rval
		sig_tt_val = -p0*a*a/rval/rval

		r.append( rval )
		sig_rr_anal.append( sig_rr_val )
		sig_tt_anal.append( sig_tt_val )

# Get stress_ij
stress_11, stress_12, stress_13, stress_22, stress_23, stress_33 = [], [], [], [], [], []
for line in open('stress_11.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		sigVal = values[1]*1e-6 # convert to MPa
		stress_11.append( sigVal )

for line in open('stress_12.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		sigVal = values[1]*1e-6 # convert to MPa
		stress_12.append( sigVal )

for line in open('stress_13.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		sigVal = values[1]*1e-6 # convert to MPa
		stress_13.append( sigVal )

for line in open('stress_22.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		sigVal = values[1]*1e-6 # convert to MPa
		stress_22.append( sigVal )

for line in open('stress_23.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		sigVal = values[1]*1e-6 # convert to MPa
		stress_23.append( sigVal )

for line in open('stress_33.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		sigVal = values[1]*1e-6 # convert to MPa
		stress_33.append( sigVal )

#Compute sig_rr, sig_tt
sig_rr, sig_tt = [], []
for i in range(len(stress_11)):
	stress = np.array([[stress_11[i],stress_12[i],stress_13[i]],\
	                   [stress_12[i],stress_22[i],stress_23[i]],\
	                   [stress_13[i],stress_23[i],stress_33[i]]])

	stressLocal = stressRotation(stress,phi_x,phi_z)
	sig_rr.append(stressLocal[0][0])
	sig_tt.append(stressLocal[1][1])


plt.subplot(221)
plt.plot(r, sig_rr, 'ko', label='GEOSX result')
plt.plot(r, sig_rr_anal,  'k', linewidth=2, label='Analytic')
plt.ylabel('Radial stress (MPa)')
plt.xlabel('r (m)')
plt.xlim(a,3*a)

plt.subplot(222)
plt.plot(r, sig_tt, 'ko', label='GEOSX result')
plt.plot(r, sig_tt_anal,  'k', linewidth=2, label='Analytic')
plt.ylabel('Tangent stress (MPa)')
plt.xlabel('r (m)')
plt.xlim(a,3*a)
plt.legend()
plt.show()
