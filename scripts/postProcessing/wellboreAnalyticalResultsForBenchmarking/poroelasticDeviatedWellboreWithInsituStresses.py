import numpy as np
import matplotlib.pyplot as plt
import poroelasticDeviatedWellbore_analytic as analytic

# Rotate stress from local coordinates of an inclined borehole to the global coordinates system
# See the description in fig.1 in Abousleiman and Cui 1998
def stressRotation(stress,phi_x,phi_z):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x),0.],[-np.sin(phi_x), np.cos(phi_x),0.],[0.,0.,1.]])
    rotz = np.array([[np.cos(phi_z),0., np.sin(phi_z)],[0.,1.,0.],[-np.sin(phi_z),0., np.cos(phi_z)]])

    return np.dot(np.dot(np.transpose(rotz),np.dot(np.dot(np.transpose(rotx),stress),rotx)),rotz)

# Rotate stress from global coordinates system to the local coordinates of an inclined borehole
# See the description in fig.1 in Abousleiman and Cui 1998
def stressRotationInv(stress,phi_x,phi_z):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x),0.],[-np.sin(phi_x), np.cos(phi_x),0.],[0.,0.,1.]])
    rotz = np.array([[np.cos(phi_z),0., np.sin(phi_z)],[0.,1.,0.],[-np.sin(phi_z),0., np.cos(phi_z)]])

    return np.dot(np.dot(np.transpose(rotx),np.dot(np.dot(np.transpose(rotz),stress),rotz)),rotx)


ri = 0.1 #borehole radius
theta = -90. * np.pi / 180. # Tangent angle to x-axis where data are extracted, converted degree to radial
phi_x = 0.* np.pi / 180. #Azimuth angle
phi_z = 45.* np.pi / 180. #Inclination angle

# Five poroelastic properties are required for an isotropic problem
E = 20.6e9 # in-plane Young's modulus
nu = 0.189 # in-plane Poisson's ratio
M = 15.8e9 # Biot's modulus
Ks = 48.2e9 # Bulk modulus of the solid phase
kappa = 1e-14 # (m2/Pa/s) ratio between the intrinsic permeability and the dynamic viscosity of fluid

# Anisotropy: three additional parameters are needed for a transversly isotropic problem
# Out-plane shear modulus is not required because only stresses are calculated
nE =1.0 # =1 for the isotropic case
nnu = 1.0 # =1 for the isotropic case
nkappa = 1.0 # ratio between the in-plane and out-plane permeability, =1 for the isotropic case

bBiot = 1.0 - E/3.0/(1.0-2.0*nu) / Ks

# Loading: in-situ stress, in-situ pore pressure, borehole pore pressure and mud pressure
pi = 0.0e6 # borehole pore pressure
pw = 0.0 # mud pressure
p0 = 10e6 # in-situ pore pressure

# For inclined borehole, the in-situ stress must be rotated to the local coordinates of the borehole
# The solutions of Abousleiman and Cui 1998 are restricted to the case where the borehole is oriented in the direction of the material anisotropy
S = analytic.stressRotation(29e6,20e6,25e6,phi_x,phi_z)

def analyticalResults():

    Sx  = S[0][0]
    Sxy = S[0][1]
    Sxz = S[0][2]
    Sy  = S[1][1]
    Syz = S[1][2]
    Sz  = S[2][2]

    if (Sx != Sy):
        theta_r = 0.5*np.arctan(2.*Sxy/(Sx-Sy))
    else:
        theta_r = 0.

    PP0 = (Sx+Sy)/2.
    S0  = -( ( (Sx-Sy)/2. )**2. + Sxy**2. )**0.5

    r = np.arange(ri,10.*ri,0.005*ri)

    # Elastic stiffnesses: see Eq.2 in Abousleiman and Cui 1998
    E_p = E/nE
    nu_p = nu/nnu
    kappa_p = kappa/nkappa

    M11 = E*(E_p - E*nu_p**2.)/(1.+nu)/(E_p-E_p*nu-2.*E*nu_p**2.)
    M12 = E*(E_p*nu + E*nu_p**2.)/(1.+nu)/(E_p-E_p*nu-2.*E*nu_p**2.)
    M13 = E*E_p*nu_p/(E_p-E_p*nu-2.*E*nu_p**2.)
    M33 = E_p**2.*(1.-nu)/(E_p-E_p*nu-2.*E*nu_p**2.)
    M44 = E/2./(1.+nu)
    #M55 = G_p
    G = M44

    # Anisotropic Biot's coefficients
    alpha = 1.-(M11+M12+M13)/(3.*Ks)
    alpha_p = 1.-(2.*M13+M33)/(3.*Ks)

    # Fluid diffusion coefficient
    c = kappa*M*M11/(M11+alpha**2.*M) 


    t = 78.0
    p, sig_rr, sig_tt, sig_rt = analytic.inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)
    #sig_zz,tau_rz,tau_tz = analytic.inTime_outPlane(r, ri,theta,p0,Sx,Sy,Sz,Sxz,Syz,sig_rr,sig_tt,p,nu_p,alpha,alpha_p)
    return [r, sig_rr/1e6, sig_tt/1e6, p/1e6]

r_anal, sig_rr_anal, sig_tt_anal, pPore_anal = analyticalResults()

fig = plt.figure(figsize=[13,10])

# Get radial coordinate and compute analytical results
r = []
for line in open('stress_11.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		rval = values[0]
		r.append( rval )

# Get stress_ij and pore pressure
stress_11, stress_12, stress_13, stress_22, stress_23, stress_33, pPore = [], [], [], [], [], [], []
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

for line in open('pressure.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		pPoreVal = values[1]*1e-6 # convert to MPa
		pPore.append( pPoreVal )

#Compute sig_rr, sig_tt
sig_rr, sig_tt = [], []
for i in range(len(stress_11)):
	stress = np.array([[stress_11[i],stress_12[i],stress_13[i]],\
		               [stress_12[i],stress_22[i],stress_23[i]],\
		               [stress_13[i],stress_23[i],stress_33[i]]])

	stressLocal = stressRotationInv(stress,theta+phi_x,phi_z)
	sig_rr.append(stressLocal[0][0])
	sig_tt.append(stressLocal[1][1])
	

plt.subplot(221)
plt.plot(r, pPore, 'ko', label='GEOSX result')
plt.plot(r_anal, pPore_anal,  'k', linewidth=2, label='Analytic')
plt.ylabel('Pore pressure (MPa)')
plt.xlabel('r (m)')
plt.xlim(ri,10*ri)

plt.subplot(223)
plt.plot(r, sig_rr, 'ko', label='GEOSX result')
plt.plot(r_anal, sig_rr_anal+bBiot*pPore_anal,  'k', linewidth=2, label='Analytic')
plt.ylabel('Effective radial stress (MPa)')
plt.xlabel('r (m)')
plt.xlim(ri,10*ri)
plt.legend()

plt.subplot(224)
plt.plot(r, sig_tt, 'ko', label='GEOSX result')
plt.plot(r_anal, sig_tt_anal+bBiot*pPore_anal,  'k', linewidth=2, label='Analytic')
plt.ylabel('Effective tangent stress (MPa)')
plt.xlabel('r (m)')
plt.xlim(ri,10*ri)

plt.show()
