# Objective: analytical modeling transient temperature and fluid diffusion around an inclined borehole 
#            & computation of pore pressure and stress  
#            in a transversly isotropic rock 
#            subjected to THM loading on the borehole wall 
#            and in-situ stress and pressure
# Method: analytical solutions are obtained in Laplace space and results in time space are obtained by the Stehfest inversion method
# This is the combination of the analytical results derived by Detournay and Cheng 1988, Wang and Papamichos 1994 and Abousleiman and Cui 1998 
# Input data are: in-situ stress and pressure, borehole orientation, mud pressure and temperature
#                 THM properties of rock
# Limitations: the borehole is oriented in the axis of the material anisotropy
# The code is extensively validated against results presented in the most relevent references

import warnings
warnings.filterwarnings("ignore")

from scipy.special import kv, iv, kn, erfc
from math import factorial, floor
import numpy as np

import matplotlib.pyplot as plt
font = {'size'   : 14}
plt.rc('font', **font)


# Rotate the in-situ stress to the local coordinates of an inclined borehole
# See the description in fig.1 in Abousleiman and Cui 1998
def stressRotation(Sx_p,Sy_p,Sz_p,phi_x,phi_z):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x),0.],[-np.sin(phi_x), np.cos(phi_x),0.],[0.,0.,1.]])
    rotz = np.array([[np.cos(phi_z),0., np.sin(phi_z)],[0.,1.,0.],[-np.sin(phi_z),0., np.cos(phi_z)]])
    
    S_p = np.array([[Sx_p,0.,0.],[0.,Sy_p,0.],[0.,0.,Sz_p]])
    return np.dot(np.dot(np.transpose(rotz),np.dot(np.dot(np.transpose(rotx),S_p),rotx)),rotz)


# Stehfest method for inverse Laplace transformation
def Vfunction(i, N):
    sum1 = 0.
    kmin = int(floor((i + 1.) / 2.))
    kmax = min(i, N)
    for k in range(kmin, kmax+1):
        sum1 = sum1 + (1.*(k**N) * factorial(2 * k) / (factorial(N - k) * factorial(k) * factorial(k - 1) * factorial(i - k) * factorial(2 * k - i)))
    return ((-1.)**(N + i)) * sum1


# Mode 1 loading: pure elastic problem
# Without mud pressure this loading mode is described by Eq.14 in Detournay and Cheng 1988 
# With mud pressure this loading mode is described by Eq.54 in Abousleiman and Cui 1998
# Material properties are not needed for the calculation of the stresses
def inTime_mode1(r, ri, PP0, pw):
    sig_rr_1 = (PP0 - pw)*ri**2./r**2.
    sig_tt_1 = -(PP0 - pw)*ri**2./r**2.
    return [sig_rr_1, sig_tt_1]


# Mode 2 loading: axisymmetric fluid diffusion problem
# Without initial pore pressure, this mode is described by Eq.15 in Detournay and Cheng 1988
# With initial pore pressure (pi), this mode is described by Eq.55 in Abousleiman and Cui 1998
# Analytical solutions are expressed in Laplace space
# The material property that is required to compute pore pressure is the hydraulic diffusivity
# For transversely isotropic material, the diffusivity is defined by Eq.11 in Abousleiman and Cui 1998
#
#                c = kappa*M*M11/(M11+alpha**2.*M)
#
# For isotropic material, the diffusivity is defined by Eq.9 in Detournay and Cheng 1988
# Additional properties that are required to compute the stresses are the Biot's coefficient and the elastic moduli: alpha, G and M11
# In isotropic case, M11 = K + 4*G/3 and the term G*alpha/M11 is noted by eta (see Eq.B.11 in Cheng 2016)
# s is the Laplace variable

def inLaplace_mode2(s, r, ri, p0, pi, c, alpha, G, M11):
    xi = r*(s/c)**0.5
    beta = ri*(s/c)**0.5

    p_2 = -(p0-pi)/s * kv(0.,xi) / kv(0.,beta)
    sig_rr_2 = -(p0-pi)/s * (2.*G*alpha)/M11 * ( ri/r*kv(1.,xi) - ri**2./r**2*kv(1.,beta) ) / (beta * kv(0.,beta))
    sig_tt_2 = (p0-pi)/s * (2.*G*alpha)/M11 * ( ( ri/r*kv(1.,xi) - ri**2./r**2*kv(1.,beta) ) / (beta * kv(0.,beta) ) + kv(0.,xi) / kv(0.,beta) )
    return [p_2, sig_rr_2, sig_tt_2]

def inTime_mode2(t, r, ri, p0, pi, c, alpha, G, M11):    
    N = 5
    sum1 = 0.
    sum2 = 0.
    sum3 = 0.
    for j in range(1,2*N+1):
        Lresult = inLaplace_mode2(j * np.log(2.) / t, r, ri, p0, pi, c, alpha, G, M11)

        sum1 += Vfunction(j, N) * Lresult[0]
        sum2 += Vfunction(j, N) * Lresult[1]
        sum3 += Vfunction(j, N) * Lresult[2]
    return [sum1* np.log(2.) / t, sum2* np.log(2.) / t, sum3* np.log(2.) / t] 


# Mode 3 loading: pure shearing
# This loading mode is described by Eq.16 in Detournay and Cheng 1988 for the case without in-plane in-situ shear stress
# An extension for the case with in-plane shear stress sig_xy is described by Eqs.22b, 22c and 56 in Abousleiman and Cui 1998
# Analytical solutions are expressed in Laplace space
# The material property that is required to compute pore pressure is the hydraulic diffusivity
# For transversely isotropic material, the diffusivity is defined by Eq.11 in Abousleiman and Cui 1998
#
#                c = kappa*M*M11/(M11+alpha**2.*M)
#
# For isotropic material, the diffusivity is defined by Eq.9 in Detournay and Cheng 1988
# Additional properties that are required to compute the stresses are the Biot's coefficient and modulus and the elastic moduli: alpha, M, G, M11 and M12
# In isotropic case, M11 = K + 4*G/3 and M11 = K - 2*G/3
# s is the Laplace variable
# The loading is defined by the shear stress S0 and angle theta_r that are defined by Eqs.22b, 22c in Abousleiman and Cui 1998
# For the case sig_xy = 0 that was considered by Detournay and Cheng: theta_r = 0 and S0 = (sig_h_max-sig_h_min)/2
# The results depend on the angle theta, see the description of theta in Fig.1 in Detournay and Cheng 1988 or Fig.1b in Abousleiman and Cui 1998 

def inLaplace_mode3(s, r, ri, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa):
    xi = r*(s/c)**0.5
    beta = ri*(s/c)**0.5

    A1 = alpha*M/(M11+alpha**2.*M)
    A2 = (M11+M12+ 2.*alpha**2.*M)/(M11+alpha**2.*M)
    
    M1 = M11/(2.*G*alpha)*kv(2.,beta)
    M2 = kv(1.,beta)/beta + 6.*kv(2.,beta)/beta**2.
    M3 = 2.*( kv(1.,beta)/beta + 3.*kv(2.,beta)/beta**2. )
    
    C1 = 4./(2.*A1*(M3-M2) - A2*M1)
    C2 = -4.*M1/(2.*A1*(M3-M2) - A2*M1)
    C3 = (2.*A1*(M2+M3) + 3.*A2*M1)/3./(2.*A1*(M3-M2) - A2*M1)

    p_3 = np.cos(2.*(theta-theta_r)) *S0/s *( c/(2.*G*kappa)*C1*kv(2.,xi) + A1*C2* ri**2./r**2. ) 
    sig_rr_3 = np.cos(2.*(theta-theta_r)) *S0/s *( A1*C1*( kv(1.,xi)/xi + 6.*kv(2.,xi)/xi**2. ) - A2*C2* ri**2./r**2. - 3.*C3* ri**4./r**4.  ) 
    sig_tt_3 = np.cos(2.*(theta-theta_r)) *S0/s *( -A1*C1*( kv(1.,xi)/xi + kv(2.,xi) + 6.*kv(2.,xi)/xi**2. ) + 3.*C3* ri**4./r**4.  ) 
    tau_rt_3 = np.sin(2.*(theta-theta_r)) *S0/s *( 2.*A1*C1*( kv(1.,xi)/xi + 3.*kv(2.,xi)/xi**2. ) - 0.5*A2*C2* ri**2./r**2. - 3.*C3* ri**4./r**4.  ) 
    
    return [p_3,sig_rr_3,sig_tt_3,tau_rt_3]

def inTime_mode3(t, r, ri, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa):    
    N = 5
    sum1 = 0.
    sum2 = 0.
    sum3 = 0.
    sum4 = 0.
    for j in range(1,2*N+1):
        Lresult = inLaplace_mode3(j * np.log(2.) / t, r, ri, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)

        sum1 += Vfunction(j, N) * Lresult[0]
        sum2 += Vfunction(j, N) * Lresult[1]
        sum3 += Vfunction(j, N) * Lresult[2]
        sum4 += Vfunction(j, N) * Lresult[3]
    return [sum1* np.log(2.) / t, sum2* np.log(2.) / t, sum3* np.log(2.) / t, sum4* np.log(2.) / t] 


# Plane-strain problem: sum-up modes 1, 2 and 3
def inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa):
    sig_rr_1, sig_tt_1 = inTime_mode1(r, ri, PP0, pw)
    p_2,sig_rr_2,sig_tt_2 = inTime_mode2(t, r, ri, p0, pi, c, alpha, G, M11)
    p_3, sig_rr_3, sig_tt_3, sig_rt_3 = inTime_mode3(t, r, ri, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)
    
    sig_rr = -PP0 + S0*np.cos( 2.*(theta - theta_r) ) + sig_rr_1 + sig_rr_2 + sig_rr_3
    sig_tt = -PP0 - S0*np.cos( 2.*(theta - theta_r) ) + sig_tt_1 + sig_tt_2 + sig_tt_3
    p = p0 + p_2 + p_3
    sig_rt = -S0*np.sin(2.*(theta-theta_r))  + sig_rt_3
    return [p, sig_rr, sig_tt, sig_rt]


# Reference: Abousleiman and Cui 1998, Eqs. 21c, 23a, 24a, 24b
# The out-plane Poisson ratio and Biot's coefficient are: nu_p and alpha_p
# sig_rr, sig_tt and p are the stress and pore pressure obtained by the plane-strain problem (the problem considered by Detournay and Cheng 1988)
# Sx, Sy, Sz, Sxz and Syz are the in-situ stresses rotated to the local coordinate of the borehole
# p0 is pore pressure
# angle theta is defined in Fig.1 in Detournay and Cheng 1988 or Fig.1b in Abousleiman and Cui 1998

def inTime_outPlane(r, ri,theta,p0,Sx,Sy,Sz,Sxz,Syz,sig_rr,sig_tt,p,nu_p,alpha,alpha_p):
	sig_zz = nu_p*(sig_rr+sig_tt) - (alpha_p - 2.*nu_p*alpha)*p - Sz + nu_p*(Sx+Sy) + (alpha_p - 2.*nu_p*alpha)*p0
	tau_rz = -(Sxz*np.cos(theta) + Syz*np.sin(theta))*(1 - ri**2./r**2.)
	tau_tz = (Sxz*np.sin(theta) - Syz*np.cos(theta))*(1 + ri**2./r**2.)
	return [sig_zz,tau_rz,tau_tz]


# Thermal loading: thermal and fluid diffusion by a tempareture raise on the borehole
# This problem was developed by Wang and Papamichos 1994
# It is also reformulated in Cheng 2016
# The code is programmed based on the results given in Cheng 2016
# The applied temperature is noted by T1 following Eq.11.332 in p.654
# Temperature diffusion is uncoupled from fluid diffusion


def inLaplace_thermal(s, r, ri, Ti, G, nu, alpha, alpha_d, c, kappaT, alpha_e):
    T1 = Ti # The applied temperature Ti is noted by T1 in Cheng 2016
    eta = alpha * (1-2*nu)/2/(1-nu)
    eta_d = alpha_d*(1.-2.*nu)/2./(1.-nu)

    # Temperature, see Eq.11.344 p.656
    T = T1/s* kn(0, r*(s/kappaT)** 0.5) / kn(0, ri*(s/kappaT)** 0.5)
    
    # Pore pressure, see Eq.11.346 p.656
    p = alpha_e*T1/(1-c/kappaT)/s * ( kn(0, r*(s/kappaT)** 0.5) / kn(0, ri*(s/kappaT)** 0.5) - kn(0, r*(s/c)** 0.5) / kn(0, ri*(s/c)** 0.5) )
    
    # The integral terms needed to determine the stresses, Eqs. 11.347 to 11.349
    int_rT = -( T1*r*kn(1, r*(s/kappaT)**0.5) ) / ( s*(s/kappaT)**0.5 * kn(0, ri*(s/kappaT)**0.5) )
    int_rp = -alpha_e/(1-c/kappaT) * T1*r/s * ( kn(1, r*(s/kappaT)**0.5) / ((s/kappaT)**0.5 * kn(0, ri*(s/kappaT)**0.5)) - kn(1, r*(s/c)**0.5) / ((s/c)**0.5 * kn(0, ri*(s/c)**0.5)) )
    A2 = T1*ri/(1-c/kappaT)/G/s * ( ( eta*alpha_e + eta_d*(1-c/kappaT) ) * kn(1, ri*(s/kappaT)**0.5) / ((s/kappaT)**0.5 * kn(0, ri*(s/kappaT)**0.5)) - eta*alpha_e * kn(1, ri*(s/c)**0.5) / ((s/c)**0.5 * kn(0, ri*(s/c)**0.5)) )
    
    # Stress and displacement, Eqs. 11.326 to 11.329
    ur = eta/G/r*int_rp + eta_d/G/r*int_rT + A2/r
    sig_rr = -2.*eta/r**2. * int_rp - 2.*eta_d/r**2. * int_rT - 2.*G/r**2.*A2
    sig_tt = 2.*eta/r**2. * int_rp + 2.*eta_d/r**2. * int_rT + 2.*G/r**2.*A2 - 2.*eta*p - 2.*eta_d*T
    sig_zz = -2.*eta*p - 2.*eta_d*T
    return [T,p,ur,sig_rr,sig_tt, sig_zz]

def inTime_thermal(t, r, ri, T1, G, nu, alpha, alpha_d, c, kappaT, alpha_e):    
    N = 5
    sum1 = 0.
    sum2 = 0.
    sum3 = 0.
    sum4 = 0.
    sum5 = 0.
    sum6 = 0.
    for j in range(1,2*N+1):
        Lresult = inLaplace_thermal(j * np.log(2.) / t, r, ri, T1, G, nu, alpha, alpha_d, c, kappaT, alpha_e)

        sum1 += Vfunction(j, N) * Lresult[0]
        sum2 += Vfunction(j, N) * Lresult[1]
        sum3 += Vfunction(j, N) * Lresult[2]
        sum4 += Vfunction(j, N) * Lresult[3]
        sum5 += Vfunction(j, N) * Lresult[4]
        sum6 += Vfunction(j, N) * Lresult[5]
    return [sum1* np.log(2.) / t, sum2* np.log(2.) / t, sum3* np.log(2.) / t, sum4* np.log(2.) / t, sum5* np.log(2.) / t, sum6* np.log(2.) / t] 


# 3D problem: Inclined borehole in transversely isotropic rock
# sum-up modes 1, 2, 3, thermal and out-plane loadings
# Limitations: rock anisotropy coincide with the axe of the hole

def inTime(t, r, theta, ri, Ti, pi, pw, p0, Sx_p, Sy_p, Sz_p, phi_x, phi_z, E, nu, M, Ks, kappa, nE, nnu, nkappa, alpha_d, kappaT, alpha_e):
    # For inclined borehole, the in-situ stress must be rotated to the local coordinates of the borehole
    # The solutions of Abousleiman and Cui 1998 are restricted to the case where the borehole is oriented in the direction of the material anisotropy
    S = stressRotation(Sx_p,Sy_p,Sz_p,phi_x,phi_z)
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


    # Out-plane properties
    E_p = E/nE
    nu_p = nu/nnu
    kappa_p = kappa/nkappa

    # Elastic stiffness, see Abousleiman and Cui 1998    
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
    
    p, sig_rr, sig_tt, sig_rt = inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)
    
    # Out-plane loading
    sig_zz, tau_rz, tau_tz = inTime_outPlane(r, ri,theta,p0,Sx,Sy,Sz,Sxz,Syz,sig_rr,sig_tt,p,nu_p,alpha,alpha_p)
    
    # Thermal loading, axisymmetric thermal and fluid diffusion
    T_thermal,p_thermal,ur_thermal,sig_rr_thermal,sig_tt_thermal, sig_zz_thermal = inTime_thermal(t, r, ri, Ti, G, nu, alpha, alpha_d, c, kappaT, alpha_e)
    
    p += p_thermal
    sig_rr +=sig_rr_thermal
    sig_tt +=sig_tt_thermal
    sig_zz +=sig_zz_thermal
    
    
    return [p, sig_rr, sig_tt, sig_rt, sig_zz, tau_rz, tau_tz]




################## Tests##########################
def testStressRotation():
    print('Validation vs data given by Cui et al. 1997, p.958')
    print(stressRotation(29.,20.,25.,0.,70.*np.pi/180.))

def testMode1():
    PP0 = 1e6
    pw = 0
    ri = 0.1
    r = np.arange(ri,3.*ri,0.1*ri)
    results = inTime_mode1(r,ri, PP0, pw)
    
    plt.plot(r/ri,results[0],label='sig_rr')
    plt.plot(r/ri,results[1],label='sig_tt')
    plt.xlabel('r/Rhole')
    plt.ylabel('Stress (Pa)')
    plt.legend()
    plt.title('Loading mode 1: pure elastic axisymmetric\n')
    plt.show()

    
def vsDetournay1988_mode2():
    #Test/Validation: Detournay and Cheng 1988 show results obtained by mode 2 in Fig. 2 that should be considered to validate the present code
    p0 = 1e6
    pi = 0.

    ri = 0.1
    r = np.arange(ri,3.*ri,0.01*ri)

    # Note that the normalized stress sig_tt/(eta*p0) is independent of these parameters as mentioned in Detournay and Cheng 1988
    # Therefore, they are all set to 1. except the Poisson ratio is set to 0. Any other values can be considered.
    c = 1.
    alpha = 1.
    nu = 0.
    G = 1.

    K = 2.*G*(1+nu)/3./(1.-2.*nu)
    M11 = K + 4.*G/3.
    eta = alpha*(1-2.*nu)/2./(1-nu)  # eta = alpha*G/M11

    tstar = [0.01, 0.1, 1, 10, 1000]
    for ts in tstar:
        t = ts * ri**2./c
        p,sig_rr,sig_tt = inTime_mode2(t, r, ri, p0, pi, c, alpha, G, M11)
        plt.plot(r/ri,sig_tt/(eta*p0),label=str(ts))

    plt.xlabel('r/Rhole')
    plt.ylabel('sig_tt/(eta*p0)')
    plt.ylim(-0.2,2)
    plt.xlim(1,3)
    plt.legend()
    plt.title('vs Fig.2 in Detournay and Cheng 1988\n')
    plt.show()

def vsDetournay1988_mode3():
    #Test/Validation: Detournay and Cheng 1988 show results obtained by mode 3 in Figs. 5, 6, 9 and 10 that should be considered to validate the present code
    fig = plt.figure(figsize=[12,12])

    theta = 0.
    nu = 0.2
    nuu = 0.4
    B = 0.8 # Skempton
    theta_r = 0. # Detournay and Cheng 1988 do not consider sig_xy, then theta_r = 0

    # Note that the normalized pore pressure p/S0 is independent of S0, G, kappa, ri. Therefore, they are all set to 1. Any other values can be considered.
    S0 = 1.
    ri = 1.
    G = 1.
    kappa = 1.

    r = np.arange(ri,3.*ri,0.005*ri)

    c = 2.*kappa*B**2.*G*(1.-nu)*(1.+nuu)**2./9./(1.-nuu)/(nuu-nu)

    alpha = 3.*(nuu-nu)/B/(1.-2.*nu)/(1.+nuu)
    K = 2.*G*(1+nu)/3./(1.-2.*nu)
    Ku = 2.*G*(1+nuu)/3./(1.-2.*nuu)
    M = (Ku - K)/(alpha**2.)
    M11 = K + 4.*G/3.
    M12 = K - 2.*G/3.

    tstar = [1e-4, 1e-3, 1e-2, 1e-1, 1.]
    for ts in tstar:
        t = ts * ri**2./c

        p, sig_rr, sig_tt, sig_rt = inTime_mode3(t, r, ri, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)

        # Note: Fig.6 in in Detournay and Cheng 1988 does not present the results of mode 3 only
        # The term - S0*np.cos( 2.*(theta - theta_r) ) need to be added to sig_tt to fit with fig.6 in Detournay and Cheng 1988. 
        # It corresponds to the initial shear stress, see Eq.21b in Abousleiman and Cui 1998

        sig_tt += - S0*np.cos( 2.*(theta - theta_r) ) 

        # Note: Fig.9 in in Detournay and Cheng 1988 does not present the results of mode 3 only
        # The term +S0*np.cos( 2.*(theta - theta_r) ) need to be added to sig_rr to fit with fig.9 in Detournay and Cheng 1988. 
        # The term - S0*np.sin(2.*(theta-theta_r)) need to be added to sig_rt to fit with fig.9 in Detournay and Cheng 1988. 
        # It corresponds to the initial shear stress, see Eq.21d in Abousleiman and Cui 1998
        sig_rr += S0*np.cos( 2.*(theta - theta_r) ) 
        sig_rt += - S0*np.sin(2.*(theta-theta_r))

        # Compute deviatoric and mean stress, see Eqs. 57 and 58 in Detournay and Cheng 1988
        S = 1./2.*((sig_tt-sig_rr)**2.+4.*sig_rt**2.)**0.5
        P = -1./2.*(sig_tt + sig_rr)

        plt.subplot(221)
        plt.plot(r/ri,p/S0,label=str(ts))
        plt.subplot(222)
        plt.plot(r/ri,sig_tt/S0,label=str(ts))
        plt.subplot(223)
        plt.plot(P/S0,S/S0,label=str(ts))
        plt.subplot(224)
        plt.plot((P-p)/S0,S/S0,label=str(ts))

    plt.subplot(221)
    plt.xlabel('r/Rhole')
    plt.ylabel('p/S0')
    plt.ylim(0,1.5)
    plt.xlim(1,2.5)
    plt.legend()
    plt.title('vs Fig.5 in Detournay and Cheng 1988\n')

    plt.subplot(222)
    plt.xlabel('r/Rhole')
    plt.ylabel('sig_tt/S0')
    plt.ylim(-4,-1.5)
    plt.xlim(1,1.5)
    plt.title('vs Fig.6 in Detournay and Cheng 1988\n')

    plt.subplot(223)
    plt.xlabel('P/S0')
    plt.ylabel('S/S0')
    plt.ylim(0.5,2.)
    plt.xlim(0,2.)
    plt.title('vs Fig.9 in Detournay and Cheng 1988\n')

    plt.subplot(224)
    plt.xlabel('P_eff Terzaghi/S0')
    plt.ylabel('S/S0')
    plt.ylim(0.5,2.)
    plt.xlim(0,2.)
    plt.title('vs Fig.10 in Detournay and Cheng 1988\n')


    #plt.subplots_adjust(bottom=0, right=1., top=1.)
    plt.show()

def vsDetournay1988_mode123():
    #Test/Validation: Detournay and Cheng 1988 show results obtained for mixed modes 1, 2 and 3 in Fig. 11 that should be considered to validate the present code
    fig = plt.figure(figsize=[15,7])

    theta = 0.
    nu = 0.2
    nuu = 0.4
    B = 0.8 # Skempton
    theta_r = 0. # Detournay and Cheng 1988 do not consider sig_xy, then theta_r = 0

    # The normalized loading are given Detournay and Cheng 1988, p.179 just after Eq.61
    PP0 = 7.
    pw = 4.
    p0 = 4.
    pi = 4. # pi = pw, no cake filter
    S0 = 1.

    # Any value can be taken for following parameters as they are canceled in the calculations
    ri = 1.
    G = 1.
    kappa = 1.



    c = 2.*kappa*B**2.*G*(1.-nu)*(1.+nuu)**2./9./(1.-nuu)/(nuu-nu)
    alpha = 3.*(nuu-nu)/B/(1.-2.*nu)/(1.+nuu)
    K = 2.*G*(1+nu)/3./(1.-2.*nu)
    Ku = 2.*G*(1+nuu)/3./(1.-2.*nuu)
    M = (Ku - K)/(alpha**2.)
    M11 = K + 4.*G/3.
    M12 = K - 2.*G/3.

    # Fixed time
    r = np.arange(ri,3.*ri,0.005*ri)
    tstar = [1e-4, 1e-3, 1e-2, 1e-1, 1.]
    for ts in tstar:
        t = ts * ri**2./c
        p, sig_rr, sig_tt, sig_rt = inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)

        # Compute deviatoric and mean stress, see Eqs. 57 and 58 in Detournay and Cheng 1988
        S = 1./2.*((sig_tt-sig_rr)**2.+4.*sig_rt**2.)**0.5
        P = -1./2.*(sig_tt + sig_rr)

        plt.subplot(111)
        plt.plot((P-p)/S0,S/S0,label='tstar='+str(ts)+', r varies')



    # Fixed radius
    rho = [1., 1.05, 1.1, 1.2]
    list_r = [val*ri for val in rho]
    tstar = 10**np.arange(-4,0,0.1)


    for r in list_r:
        t = tstar * ri**2./c
        p, sig_rr, sig_tt, sig_rt = inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)

        # Compute deviatoric and mean stress, see Eqs. 57 and 58 in Detournay and Cheng 1988
        S = 1./2.*((sig_tt-sig_rr)**2.+4.*sig_rt**2.)**0.5
        P = -1./2.*(sig_tt + sig_rr)

        plt.subplot(111)
        plt.plot((P-p)/S0,S/S0,'--',label= 'r/Rhole = '+str(r)+', t varies')

    plt.subplot(111)
    plt.xlabel('P_eff Terzaghi/S0')
    plt.ylabel('S/S0')
    plt.xlim(3.,5.)
    plt.ylim(1.,6.)
    plt.legend()#loc='top left', bbox_to_anchor=(1,1))
    plt.title('vs Fig.11 in Detournay and Cheng 1988\n')

    #plt.subplots_adjust(bottom=0, right=1., top=1.)
    plt.show()


def vsAbousleiman1998():
    # Test/Validation
    ri = 0.1 #borehole radius
    theta = 90. * np.pi / 180. # converted degree to radial

    # Five poroelastic properties are required for an isotropic problem
    E = 20.6e9 # in-plane Young's modulus
    nu = 0.189 # in-plane Poisson's ratio
    M = 15.8e9 # Biot's modulus
    Ks = 48.2e9 # Bulk modulus of the solid phase
    kappa = 1e-16 # (m2/Pa/s) ratio between the intrinsic permeability and the dynamic viscosity of fluid

    # Anisotropy: three additional parameters are needed for a transversly isotropic problem
    # Out-plane shear modulus is not required because only stresses are calculated
    nE = 1.0 # ratio between the in-plane and out-plane Young's moduli, =1 for the isotropic case
    #nnu varies, =1 for the isotropic case
    nkappa = 1.0 # ratio between the in-plane and out-plane permeability, =1 for the isotropic case

    # Loading: in-situ stress, in-situ pore pressure, borehole pore pressure and mud pressure
    pi = 0.0
    pw = 0.0
    p0 = 9.8e6
    Sx_p = 29e6
    Sy_p = 20e6
    Sz_p = 25e6
    phi_x = 30.*np.pi/180.
    phi_z = 60.*np.pi/180.
   
    fig = plt.figure(figsize=[15,15])

    r = np.arange(ri,1.5*ri,0.005*ri)

    list_nnu = [0.5,1.,2.]
    for nnu in list_nnu:
        
        t = 0.001 * 24 * 3600 # 0.001 day
        p, sig_rr, sig_tt, sig_rt, sig_zz, tau_rz, tau_tz = inTime(t, r, theta, ri, 0, pi, pw, p0, Sx_p, Sy_p, Sz_p, phi_x, phi_z, E, nu, M, Ks, kappa, nE, nnu, nkappa, 1., 1., 1.)
        
        plt.subplot(221)
        plt.plot(r/ri,p/1e6,label='n_nu='+str(nnu)+', t='+str(round(t/24./3600.,3))+' (day)')
        plt.subplot(222)
        plt.plot(r/ri,-(sig_rr+p)/1e6, label='n_nu='+str(nnu)+', t='+str(round(t/24./3600.,3))+' (day)')
        plt.subplot(223)
        plt.plot(r/ri,-(sig_tt+p)/1e6, label='n_nu='+str(nnu)+', t='+str(round(t/24./3600.,3))+' (day)')
        plt.subplot(224)
        plt.plot(r/ri,-(sig_zz+p)/1e6, label='n_nu='+str(nnu)+', t='+str(round(t/24./3600.,3))+' (day)')


        t = 1. * 24 * 3600 # 1 day
        p, sig_rr, sig_tt, sig_rt, sig_zz, tau_rz, tau_tz = inTime(t, r, theta, ri, 0, pi, pw, p0, Sx_p, Sy_p, Sz_p, phi_x, phi_z, E, nu, M, Ks, kappa, nE, nnu, nkappa, 1., 1., 1.)
        plt.subplot(221)
        plt.plot(r/ri,p/1e6, '--',label='n_nu='+str(nnu)+', t='+str(round(t/24./3600.,3))+' (day)')
        plt.subplot(222)
        plt.plot(r/ri,-(sig_rr+p)/1e6, '--',label='n_nu='+str(nnu)+', t='+str(round(t/24./3600.,3))+' (day)')
        plt.subplot(223)
        plt.plot(r/ri,-(sig_tt+p)/1e6, '--',label='n_nu='+str(nnu)+', t='+str(round(t/24./3600.,3))+' (day)')
        plt.subplot(224)
        plt.plot(r/ri,-(sig_zz+p)/1e6, '--', label='n_nu='+str(nnu)+', t='+str(round(t/24./3600.,3))+' (day)')


    plt.subplot(221)
    plt.xlabel('r/Rhole')
    plt.ylabel('Pore pressure (MPa)')
    plt.xlim(1,1.5)
    plt.ylim(-2,12)
    plt.title('vs Fig.3 in Abousleiman and Cui 1998\n')

    plt.subplot(222)
    plt.xlabel('r/Rhole')
    plt.ylabel('-(sig_rr+p) (MPa)')
    plt.xlim(1,1.5)
    plt.ylim(-4,12)
    plt.legend()#loc='top left', bbox_to_anchor=(1,1))
    plt.title('vs Fig.4 in Abousleiman and Cui 1998\n')

    plt.subplot(223)
    plt.xlabel('r/Rhole')
    plt.ylabel('-(sig_tt+p) (MPa)')
    plt.xlim(1,1.5)
    plt.ylim(24,52)
    plt.title('vs Fig.5 in Abousleiman and Cui 1998\n')

    plt.subplot(224)
    plt.xlabel('r/Rhole')
    plt.ylabel('-(sig_zz+p) (MPa)')
    plt.xlim(1,1.5)
    plt.ylim(10,28)
    plt.title('vs Fig.6 in Abousleiman and Cui 1998\n')

    #plt.subplots_adjust(bottom=0, right=1., top=1.)
    plt.show()

def vsAbousleiman1998_continue():
	# Test/Validation
	ri = 0.1 #borehole radius
	theta = 90. * np.pi / 180. # converted degree to radial

	# Five poroelastic properties are required for an isotropic problem
	E = 20.6e9 # in-plane Young's modulus
	nu = 0.189 # in-plane Poisson's ratio
	M = 15.8e9 # Biot's modulus
	Ks = 48.2e9 # Bulk modulus of the solid phase
	kappa = 1e-16 # (m2/Pa/s) ratio between the intrinsic permeability and the dynamic viscosity of fluid

	# Anisotropy: three additional parameters are needed for a transversly isotropic problem
	# Out-plane shear modulus is not required because only stresses are calculated
	#nE varies, =1 for the isotropic case
	nnu = 1.0 # =1 for the isotropic case
	nkappa = 1.0 # ratio between the in-plane and out-plane permeability, =1 for the isotropic case

	# Loading: in-situ stress, in-situ pore pressure, borehole pore pressure and mud pressure
	pi = 0.0
	pw = 0.0
	p0 = 9.8e6

	# For inclined borehole, the in-situ stress must be rotated to the local coordinates of the borehole
	# The solutions of Abousleiman and Cui 1998 are restricted to the case where the borehole is oriented in the direction of the material anisotropy
	S = stressRotation(29e6,20e6,25e6,30.*np.pi/180.,60.*np.pi/180.)
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


	fig = plt.figure(figsize=[7,7])

	r = np.arange(ri,1.5*ri,0.005*ri)

	list_nE = [0.5,1.,2.]
	for nE in list_nE:
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


		t = 0.001 * 24 * 3600 # 0.001 day
		p, sig_rr, sig_tt, sig_rt = inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)
		sig_zz,tau_rz,tau_tz = inTime_outPlane(r, ri,theta,p0,Sx,Sy,Sz,Sxz,Syz,sig_rr,sig_tt,p,nu_p,alpha,alpha_p)

		plt.subplot(221)
		plt.plot(r/ri,p/1e6,label='n_nu='+str(nnu)+', t='+str(round(t/24./3600.,3))+' (day)')
		plt.subplot(222)
		plt.plot(r/ri,-(sig_rr+p)/1e6, label='E/Ep='+str(nE)+', t='+str(round(t/24./3600.,3))+' (day)')
		plt.subplot(223)
		plt.plot(r/ri,-(sig_tt+p)/1e6, label='E/Ep='+str(nE)+', t='+str(round(t/24./3600.,3))+' (day)')
		plt.subplot(224)
		plt.plot(r/ri,-(sig_zz+p)/1e6, label='E/Ep='+str(nE)+', t='+str(round(t/24./3600.,3))+' (day)')


		t = 1. * 24 * 3600 # 1 day
		#p, sig_rr, sig_tt, sig_rt, sig_zz, tau_rz, tau_tz = inTime(t, r, theta, ri, 0, pi, pw, p0, Sx_p, Sy_p, Sz_p, phi_x, phi_z, E, nu, M, Ks, kappa, nE, nnu, nkappa, 1., 1., 1.)
		p, sig_rr, sig_tt, sig_rt = inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)
		sig_zz,tau_rz,tau_tz = inTime_outPlane(r, ri,theta,p0,Sx,Sy,Sz,Sxz,Syz,sig_rr,sig_tt,p,nu_p,alpha,alpha_p)

		plt.subplot(221)
		plt.plot(r/ri,p/1e6, '--',label='n_nu='+str(nnu)+', t='+str(round(t/24./3600.,3))+' (day)')
		plt.subplot(222)
		plt.plot(r/ri,-(sig_rr+p)/1e6, '--',label='E/Ep='+str(nE)+', t='+str(round(t/24./3600.,3))+' (day)')
		plt.subplot(223)
		plt.plot(r/ri,-(sig_tt+p)/1e6, '--',label='E/Ep='+str(nE)+', t='+str(round(t/24./3600.,3))+' (day)')
		plt.subplot(224)
		plt.plot(r/ri,-(sig_zz+p)/1e6, '--', label='E/Ep='+str(nE)+', t='+str(round(t/24./3600.,3))+' (day)')


	plt.subplot(221)
	plt.xlabel('r/Rhole')
	plt.ylabel('Pore pressure (MPa)')
	plt.xlim(1,1.5)
	plt.ylim(-2,12)
	plt.title('vs Fig.3 in Abousleiman and Cui 1998\n')

	plt.subplot(222)
	plt.xlabel('r/Rhole')
	plt.ylabel('-(sig_rr+p) (MPa)')
	plt.xlim(1,1.5)
	plt.ylim(-4,12)
	plt.legend()#loc='top left', bbox_to_anchor=(1,1))
	plt.title('vs Fig.4 in Abousleiman and Cui 1998\n')

	plt.subplot(223)
	plt.xlabel('r/Rhole')
	plt.ylabel('-(sig_tt+p) (MPa)')
	plt.xlim(1,1.5)
	plt.ylim(24,52)
	plt.title('vs Fig.5 in Abousleiman and Cui 1998\n')

	plt.subplot(224)
	plt.xlabel('r/Rhole')
	plt.ylabel('-(sig_zz+p) (MPa)')
	plt.xlim(1,1.5)
	plt.ylim(10,28)
	plt.title('vs Fig.6 in Abousleiman and Cui 1998\n')

	#plt.subplots_adjust(bottom=0, right=1., top=1.)
	plt.show()



def vsAbousleiman1998_continue_bu():
    # Test/Validation
    ri = 0.1 #borehole radius
    theta = 90. * np.pi / 180. # converted degree to radial

    # Five poroelastic properties are required for an isotropic problem
    E = 20.6e9 # in-plane Young's modulus
    nu = 0.189 # in-plane Poisson's ratio
    M = 15.8e9 # Biot's modulus
    Ks = 48.2e9 # Bulk modulus of the solid phase
    kappa = 1e-16 # (m2/Pa/s) ratio between the intrinsic permeability and the dynamic viscosity of fluid

    # Anisotropy: three additional parameters are needed for a transversly isotropic problem
    # Out-plane shear modulus is not required because only stresses are calculated
    #nE varies, =1 for the isotropic case
    nnu = 1.0 # =1 for the isotropic case
    nkappa = 1.0 # ratio between the in-plane and out-plane permeability, =1 for the isotropic case

    # Loading: in-situ stress, in-situ pore pressure, borehole pore pressure and mud pressure
    pi = 0.0
    pw = 0.0
    p0 = 9.8e6

    # For inclined borehole, the in-situ stress must be rotated to the local coordinates of the borehole
    # The solutions of Abousleiman and Cui 1998 are restricted to the case where the borehole is oriented in the direction of the material anisotropy
    S = stressRotation(29e6,20e6,25e6,30.*np.pi/180.,60.*np.pi/180.)
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


    fig = plt.figure(figsize=[7,7])

    r = np.arange(ri,1.5*ri,0.005*ri)

    list_nE = [0.5,1.,2.]
    for nE in list_nE:
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


        t = 0.001 * 24 * 3600 # 0.001 day
        p, sig_rr, sig_tt, sig_rt = inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)
        sig_zz,tau_rz,tau_tz = inTime_outPlane(r, ri,theta,p0,Sx,Sy,Sz,Sxz,Syz,sig_rr,sig_tt,p,nu_p,alpha,alpha_p)

        plt.subplot(111)
        plt.plot(r/ri,-(sig_zz+p)/1e6, label='n_nu='+str(nnu)+', t='+str(round(t/24./3600.,3))+' (day)')


        t = 1. * 24 * 3600 # 1 day
        p, sig_rr, sig_tt, sig_rt = inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)    
        sig_zz,tau_rz,tau_tz = inTime_outPlane(r, ri,theta,p0,Sx,Sy,Sz,Sxz,Syz,sig_rr,sig_tt,p,nu_p,alpha,alpha_p)

        plt.subplot(111)
        plt.plot(r/ri,-(sig_zz+p)/1e6, label='n_nu='+str(nnu)+', t='+str(round(t/24./3600.,3))+' (day)')


    plt.subplot(111)
    plt.xlabel('r/Rhole')
    plt.ylabel('-(sig_zz+p) (MPa)')
    plt.xlim(1,1.5)
    plt.ylim(10,28)
    plt.title('vs Fig.7 in Abousleiman and Cui 1998\n')
    plt.legend()
    #plt.subplots_adjust(bottom=0, right=1., top=1.)
    plt.show()

def vsCui1997():
    # Test/Validation
    ri = 0.1 #borehole radius
    theta = 84.4 * np.pi / 180. # converted degree to radial

    # Five poroelastic properties are required for an isotropic problem
    E = 20.6e9 # in-plane Young's modulus
    nu = 0.189 # in-plane Poisson's ratio
    M = 15.8e9 # Biot's modulus
    Ks = 48.2e9 # Bulk modulus of the solid phase
    kappa = 1e-16 # (m2/Pa/s) ratio between the intrinsic permeability and the dynamic viscosity of fluid

    # Anisotropy: three additional parameters are needed for a transversly isotropic problem
    # Out-plane shear modulus is not required because only stresses are calculated
    nE =1.0 # =1 for the isotropic case
    nnu = 1.0 # =1 for the isotropic case
    nkappa = 1.0 # ratio between the in-plane and out-plane permeability, =1 for the isotropic case

    # Loading: in-situ stress, in-situ pore pressure, borehole pore pressure and mud pressure
    pi = 0.0
    pw = 0.0
    p0 = 10e6

    # For inclined borehole, the in-situ stress must be rotated to the local coordinates of the borehole
    # The solutions of Abousleiman and Cui 1998 are restricted to the case where the borehole is oriented in the direction of the material anisotropy
    S = stressRotation(29e6,20e6,25e6,0.*np.pi/180.,70.*np.pi/180.)
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


    fig = plt.figure(figsize=[7,7])

    r = np.arange(ri,3.*ri,0.005*ri)

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


    t = 1.3*60. # 1.3 min
    p, sig_rr, sig_tt, sig_rt = inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)
    #sig_zz,tau_rz,tau_tz = inTime_outPlane(r, ri,theta,p0,Sx,Sy,Sz,Sxz,Syz,sig_rr,sig_tt,p,nu_p,alpha,alpha_p)

    plt.subplot(111)
    plt.plot(r/ri,p/1e6, label='t='+str(round(t/60.,1))+' (min)')


    t = 21.6*60. # 21.6min
    p, sig_rr, sig_tt, sig_rt = inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)    
    #sig_zz,tau_rz,tau_tz = inTime_outPlane(r, ri,theta,p0,Sx,Sy,Sz,Sxz,Syz,sig_rr,sig_tt,p,nu_p,alpha,alpha_p)

    plt.subplot(111)
    plt.plot(r/ri,p/1e6, label='t='+str(round(t/60.,1))+' (min)')


    plt.subplot(111)
    plt.xlabel('r/Rhole')
    plt.ylabel('p (MPa)')
    plt.xlim(1,3.)
    plt.ylim(0,14.)
    plt.legend()
    plt.title('vs Fig.4 in Cui et al. 1997\n')

    #plt.subplots_adjust(bottom=0, right=1., top=1.)
    plt.show()

def vsCui1997_continue():
    # Test/Validation
    ri = 0.1 #borehole radius
    theta = 5.7 * np.pi / 180. # converted degree to radial

    # Five poroelastic properties are required for an isotropic problem
    E = 20.6e9 # in-plane Young's modulus
    nu = 0.189 # in-plane Poisson's ratio
    M = 15.8e9 # Biot's modulus
    Ks = 48.2e9 # Bulk modulus of the solid phase
    kappa = 1e-16 # (m2/Pa/s) ratio between the intrinsic permeability and the dynamic viscosity of fluid

    # Anisotropy: three additional parameters are needed for a transversly isotropic problem
    # Out-plane shear modulus is not required because only stresses are calculated
    nE =1.0 # =1 for the isotropic case
    nnu = 1.0 # =1 for the isotropic case
    nkappa = 1.0 # ratio between the in-plane and out-plane permeability, =1 for the isotropic case

    # Loading: in-situ stress, in-situ pore pressure, borehole pore pressure and mud pressure
    pi = 0.0
    pw = 0.0
    p0 = 10e6

    # For inclined borehole, the in-situ stress must be rotated to the local coordinates of the borehole
    # The solutions of Abousleiman and Cui 1998 are restricted to the case where the borehole is oriented in the direction of the material anisotropy
    S = stressRotation(29e6,20e6,25e6,0.*np.pi/180.,70.*np.pi/180.)
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


    fig = plt.figure(figsize=[7,7])

    r = np.arange(ri,3.*ri,0.005*ri)

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


    t = 1.3*60. # 1.3 min
    p, sig_rr, sig_tt, sig_rt = inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)
    #sig_zz,tau_rz,tau_tz = inTime_outPlane(r, ri,theta,p0,Sx,Sy,Sz,Sxz,Syz,sig_rr,sig_tt,p,nu_p,alpha,alpha_p)

    plt.subplot(111)
    plt.plot(r/ri,-(sig_rr+p)/1e6, label='t='+str(round(t/60.,1))+' (min)')


    t = 21.6*60. # 21.6min
    p, sig_rr, sig_tt, sig_rt = inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G, M11, M12, kappa)    
    #sig_zz,tau_rz,tau_tz = inTime_outPlane(r, ri,theta,p0,Sx,Sy,Sz,Sxz,Syz,sig_rr,sig_tt,p,nu_p,alpha,alpha_p)

    plt.subplot(111)
    plt.plot(r/ri,-(sig_rr+p)/1e6, label='t='+str(round(t/60.,1))+' (min)')


    plt.subplot(111)
    plt.xlabel('r/Rhole')
    plt.ylabel('-(sig_rr+p) (MPa)')
    plt.xlim(1,3.)
    plt.ylim(-4,14.)
    plt.legend()
    plt.title('vs Fig.5 in Cui et al. 1997\n')

    #plt.subplots_adjust(bottom=0, right=1., top=1.)
    plt.show()


def vsWang1994():
    ri=0.1
    T1 = 100. 
    p1 = T1/(6.)*1e6

    #Westerly Granite, see Tab 1 in Wang and Papamichos 1994

    # Notations by Cheng (2016)
    alpha_e = 0.39e6 #equivalent to parameter cp,PaC
    c =0.22e-4 #m2/s
    kappaT = 0.5*c/0.23 #equivalent to parameter c0, the ratio given in Tab1 of Wang and Papamichos 1994 is for c/c0, not c0/c; Also, a factor 0.5 was added to match the results given by the Figs of Wang and Papamichos.

    # Additional parameters to calculate stresses and displacement
    alpha_d = 2.48e6 #N/m2/K, thermal dilatation coefficient
    Ku = 23.3e9
    nuu = 0.274
    nu = 0.25
    eta = 0.04 # eta = alpha * (1-2*nu)/2/(1-nu)


    G = 3.*Ku*(1.-2.*nuu)/2./(1.+nuu)
    E = 2.*G*(1.+nu)
    alpha = eta / (1-2*nu) * 2.* (1-nu)
    K = 2.*G*(1.+nu)/3./(1.-2.*nu)
    M11 = K + 4./3.*G
    M = (Ku-K)/alpha**2.
    kappa = c/( M*M11/(M11+alpha**2.*M) )
    Ks = K/(1.-alpha)


    fig = plt.figure(figsize=[15,7])
    list_tstar = [0.01,0.1,1.]
    for iTime in range(len(list_tstar)):
        tstar = list_tstar[iTime]
        t = tstar*ri**2./c

        r = np.arange(ri,4.*ri,0.01*ri)
        #p, sig_rr, sig_tt, sig_rt, sig_zz, tau_rz, tau_tz                                = inTime(t, r, theta, ri, Ti, pi, pw, p0, Sx_p, Sy_p, Sz_p, phi_x, phi_z, E, nu, M, Ks, kappa, nE, nnu, nkappa, alpha_d, kappaT, alpha_e)
        p, sig_rr, sig_tt, sig_rt, sig_zz, tau_rz, tau_tz                                 = inTime(t, r, 0.,    ri, T1, p1, 0., 0., 0.,   0.,   0.,   0.,    0.,    E, nu, M, Ks, kappa, 1., 1.,  1., alpha_d, kappaT, alpha_e)
        p_thermal, sig_rr_thermal, sig_tt_thermal, sig_rt_thermal, sig_zz, tau_rz, tau_tz = inTime(t, r, 0.,    ri, T1, 0., 0., 0., 0.,   0.,   0.,   0.,    0.,    E, nu, M, Ks, kappa, 1., 1.,  1., alpha_d, kappaT, alpha_e)


        if (iTime  == 0):
            line_type1 = 'k'
            line_type2 = 'k--'
        if (iTime  == 1):
            line_type1 = 'r'
            line_type2 = 'r--'
        if (iTime  == 2):
            line_type1 = 'b'
            line_type2 = 'b--'

        plt.plot(r/ri,p/p1, line_type1, linewidth=2, label=r'$tstar$='+str(tstar))
        plt.plot(r/ri,p_thermal/p1, line_type2, linewidth=2, label=r'Pure thermal loading: $tstar$='+str(tstar))

        plt.ylabel('Normalized Pore Pressure')
        plt.xlabel('Normalized Radial Distance')
        plt.title('vs Fig.3a in Wang and Papamichos 1994')
        plt.xlim(1, 4)
        plt.ylim(0,1.2)
        plt.legend()

    plt.show()




# Run the tests
'''
testStressRotation()
testMode1()
vsDetournay1988_mode2()
vsDetournay1988_mode3()
vsDetournay1988_mode123()
vsAbousleiman1998()
'''
vsAbousleiman1998_continue()
'''
vsCui1997()
vsCui1997_continue()
vsWang1994()
'''
