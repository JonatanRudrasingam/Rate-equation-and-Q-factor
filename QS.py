
#%%

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from astropy import constants as const

#Ved Q-switching ændre man på sin losses efter et stykke tid

c = const.c.value
alpha = 10**(-15)
S = alpha

laser = 'Nd:YAG'
#starname     = 'KWHya'

re = 2
#re = 4

Qsw = False
#Qsw = True
Title = 'Wp = 900'
if (laser == 'Nd:YAG' and re == 2):
    #Active medium
    Ntot = 10**23         
    sigma = 2.8*10**(-23)
    A = 1/(230*10**(-6))
    
    #Resonator
    losses = 0.8
    L = 1 #Length
    tau = (2*L)/(c*(1-losses))
    
    #Laser parameters
    Ith = 1/(c*sigma*tau)
    Wth = (Ith*A)/(Ntot-Ith)
    Wp = 900*Wth
    
    #Solve_ivp solves the problem:
    y0 = [0, 1]
    
    tf = 1*10**(-5)
    
    #Colour parameters
    colour0 = 'pink'
    colour1 = 'limegreen'

if (laser == 'Nd:YAG' and re == 4):
    global N1B, N2B, N3B
    
    #Properties
    OPL = 0.3277 #m (Optical path length)
    n = 1.82 # (Index of YAG)
    Am = 4.2429*10**(-3)*10**(-3) #m^2 (Mode area)
    Vm = 2.5457*10**(-2)*10**(-6) #m^3 (Mode volume)
    N0 = 3.513*10**(18) # (Population of ground state)
    sigma21 = 6.5*10**(-19)*10**(-3) #m^2 (Stimulated Emission Cross Section Area)
    CWPUMP = 2
    hv808 = 2.4601*10**(-19)
    P = CWPUMP/hv808
    
    #Resonator
    losses = 0.2
    L = 1 #Length
    taucav = (2*L)/(c*(losses))
    
    #Initial population
    N1B = 1.39767*10**(14)
    N2B = 3.556*10**(-6) #0
    N3B = 5.661*10**(-8) #0
    
    #Solve_ivp solves the problem
    y0 = [N1B, N2B, N3B, 0]
    
    tf = 4*10**(-5)
    
    #Time constants
    tau10 = 11*10**(-9)
    tau20 = 395*10**(-6)
    tau21 = 550*10**(-6)
    tau30 = 50*10**(-6)
    tau32 = 450*10**(-12)
    tau2 = (tau20**(-1)+tau21**(-1))**(-1)
    tau3 = (tau30**(-1)+tau32**(-1))**(-1)
    tauRT = (2*OPL)/c


def RE2(t, y):
    global losses_Q, tau, Qsw
    
    I = y[0]
    p = y[1]
    
    if t > 1e-5/2 and Qsw == True:
        losses_Q = 0.9
        tau = (2*L)/(c*(losses_Q))
    if t > 1.5e-5/2 and Qsw == True:
        tau = (2*L)/(c*(losses))
        Qsw = False
    
    dIdt = Wp*(Ntot-I) - sigma*c*p*I-A*I
    dpdt = p*(sigma*c*I-1/tau)+alpha*A*I
    
    return dIdt, dpdt


def RE4(t, y):
    
    N1 = y[0]
    N2 = y[1]
    N3 = y[2]
    I21 = y[3]
    
    dN1dt = I21*sigma21*(N2-N1)+N2/tau21-(N1-N1B)/tau10
    dN2dt = -I21*sigma21*(N2-N1)-(N2-N2B)*(1/tau20+1/tau21)+N3/tau32
    dN3dt = P-(N3-N3B)*(1/tau30+1/tau32)
    dI21dt = I21*sigma21*(N2-N1)*c/(n*Vm)-I21/taucav+S
    
    return dN1dt, dN2dt, dN3dt, dI21dt

#Time settting output:
time_step = 0.1
time_interval = np.linspace(0,tf,10000)
t = time_interval

time_init = 0
time_end = tf

if (re == 2):
    RE = RE2
elif (re == 4):
    RE = RE4

mysol = solve_ivp(RE, (time_init, time_end), y0=y0, t_eval=time_interval)
plt.semilogy(mysol.t,mysol.y[0])


#Output is collected:
if (re == 2):
    fig, ax = plt.subplots(2)
    ax[0].plot(mysol.t,mysol.y[0],label='Inversion', color = colour0)
    ax[0].title.set_text(Title)
    ax[1].plot(mysol.t,mysol.y[1],label='Photon density', color = colour1)
    ax[1].set_xlabel('Time (s)')
    ax[0].set_ylabel('$\\frac{1}{m^3}$',fontsize=14)
    ax[1].set_ylabel('$\\frac{1}{m^3}$',fontsize=14)
    ax[0].legend()
    ax[1].legend()
    plt.show()
    
if (re == 4):
    fig = plt.figure()
    fig.subplots_adjust(top = 0.8)

    ax1 = fig.add_axes([0.15,1,1,1])
    line1 = ax1.plot(mysol.t,mysol.y[0],label='Population 1')

    ax2 = fig.add_axes([0.15,-0.2,1,1])
    line2 = ax2.plot(mysol.t,mysol.y[1],label='Population 2')
    
    ax3 = fig.add_axes([0.15,-1.4,1,1])
    line3 = ax3.plot(mysol.t,mysol.y[2],label='Population 3')

    ax4 = fig.add_axes([0.15,-2.6,1,1])
    line4 = ax4.plot(mysol.t,mysol.y[3],label='Inversion')

    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    plt.show()
