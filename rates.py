#!/usr/bin/python3

#Calculates maximim pericenter distance for gravitational radiation capture
#Cross section, timescale for 2body encounters
#and number of encounters per relaxation time and cluster lifetime (~100 relaxation times)

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

#N: total number of objetcs,
#m: Mass of a single object in solar masses
#rmin: Smallest cluster radius in parsecs
#rmax: The bigger cluster radius in parsecs
#b: [0.1,1] passage distance in maximum pericenter distances

def encounters(N,m,rmin,rmax,b):
    pi=math.pi
    G=6.67408E-11       # [m³/kgs²]
    pipi=2*pi
    Msol=1.989E30       # [kg]
    AU=1.49597e11       # [m]
    c=3E8               # [m/s]
    c2=c**2             # [m²/s²]
    r_sol=6.957E8       # m
    N=1E6               # Numero de bh
    pc= 3.086e+16       # parsec en m
    s_yr=365*24*60*60   # seconds in a year


    def fmt(n):
            decimals='{0:.2f}'.format(n) #Just for format
            return decimals


    lnn=[]          #Density List
    lnenclife=[]    #Number of encounters in lifetime List
    lencratetot=[]  #Encounters rate List
    lencrelaxt=[]   #Encounters in a relaxation time List
    ld=[]           #d List
    lclustersize=[] #Cluster size List
    llifetime=[]    #Lifetime cluster List
    ltrelax=[]      #Relaxation time for cluster List
    lvdisp=[]       #Dispersion velocity List
    lr_perimax=[]   #Maximum pericenter distance List
    lS=[]           #Capture cross section List
    ltmerger=[]     #Merger time at given maximum pericenter List
    la=[]           #Semimajor axis of the newly formed binary List
    lecc=[]         #Eccentricity List of the newly formed binary List
    lh=[]           #Hardness of the newly formed binary List
    lb=[]           #impact_parameter list
    LE_fin=[]       #Final Energy list
    lrp=[]          #Passage distance list

    m=m*Msol    #[kg]  All stars eqAUl mass
    m_2=m_1=m
    M=N*m       #[kg]  Total Mass
    m_red=(m_1*m_2)/(m_1+m_2)  #[kg] Reduced mass
    r_sch=2*G*m/c2     #Schwarzschild Radius

    #Cluster size, rmin and rmax given in terminal
    clustersize=np.arange(rmin,rmax,0.01)

    for r in clustersize:
        R=r*pc       #[m]
        n=N/((4/3)*pi*R**3)   #[particles/m³]

        # Velocity Dispersion for a self-gravitating cluster
        vdisp=((G*N*m)/R)**(1/2)        #[m/s]
        vrel=vdisp
        #vrel=vdisp=20E3                 #[m/s]
        lvdisp.append(vdisp/1000)       #[km/s]


                    ###############################################
                    # Maximum pericenter for GW radiation capture #
                    ###############################################
                            # if r>r_perimax then no capture
        cte1=((85*pi*(2**(1/2)))/12)**(2/7)
        cte2=((G**(7/2))/(c**5))**(2/7)
        mmv=(((m_1*m_2)*((m_1+m_2)**(3/2)))/(vrel**2))**(2/7)

        r_perimax= cte1*cte2*(mmv)   #[m]
        lr_perimax.append(r_perimax/1000) #[km]

        d=r_perimax/r_sch
        ld.append(d) #maximum pericenter in Schwarzschild radius [r_schw]

        #impact_parameter b
        rp=b*r_perimax     #[m] Distance passage

        E_ini=-(vrel**2)*m_red/2 #Initial energy  [J= kgm²s⁻²]

                            ###########################################
                            # Dissipated energy for a parabolic orbit #
                            #                 t rp (DeltaE)           #
                            ###########################################

        de1=-85*pi*(G**(7/2))*((m_1**2)*(m_2**2))*((m_1+m_2)**(1/2))
        de2=12*(2**(1/2))*(c**5)*(rp**(7/2))
        DeltaE=de1/de2              #[J= kgm²s⁻²]

        E_fin=E_ini+DeltaE
        LE_fin.append(E_fin) #Energy of the newly formed binary [J]
        lrp.append(rp/r_sch) #passage distance in shcwarszchild radius

                    ###################################
                    # Semimajor axis and eccentricity #
                    ###################################

#Semimajor axis of the newly formed binary
# Semimajor axis in terms of r_perimax and rp
        sma_0=-(G*(m_1*m_2)/(2*E_ini)) #[m]
        a=sma_0*(-1 + (r_perimax/rp)**(7/2))**-1 #[m] a=(-G*(m_1*m_2))/(2*E_fin)
        la.append(a/AU) #Semimajor axis in AU

        #eccentricity
        ecc=1-rp/a
        lecc.append(ecc)

        #Merger time

        ca=5*(c**5)*(a**4)
        cte3=(73/24)*(ecc*2)
        cte4=(37/96)*(ecc**4)
        fe=((1-ecc**2)**(7/2))/(1+cte3+cte4)
        Gm=256*(G**3)*(m_1*m_2)*(m_1+m_2)
        tmerger=(ca/Gm)*fe    #[s]
        ltmerger.append(tmerger)  #[s] ~10^6 s... a few days

                ##################
                # Cross Sections #
                ##################

        #trelax Cluster
        trelax=(N*R)/(math.log(0.1*N)*vrel) #[s]
        ltrelax.append(trelax)
        #lifetime cluster = 100 trelax
        tlife=100*trelax                          # [s]
        llifetime.append(tlife)                   # [s]

        #Capture cross section for 2 body at r_perimax.
        S = pi*(r_perimax**2)*(1+((G*(m_1+m_2))/((vrel**2)*r_perimax)))    # [m²]
        lS.append(S)    # [m²]

        #timescale for  2 body passage within r_perimax
        timescale_bhencounter= 1/(n*vrel*S)       # [s]
        #encounter rate: 1/timescale_bhencounter
        encounter_rate= n*vrel*S                  # [1/s]
        #Encounter rate for the whole Cluster
        encounter_rateTot= encounter_rate*N       # [1/s]
        lencratetot.append(encounter_rateTot)     # [1/s]


                    ########################################
                    #Encounters during lifetime and t relax#
                    ########################################

        encrelaxt=trelax*encounter_rateTot  #number of encounters in trelax
        lencrelaxt.append(encrelaxt)

        n_encountlife=tlife*encounter_rateTot   #number of encounters in tlife
        #n_encountlife=(100*(n**3)*(R**7)*S)/(math.log(0.1*(n*(r*pc)**3)))
        lnenclife.append(n_encountlife)
        #n_encountlife=(100*(n**3)*(R**7)*S)/(math.log(0.1*(n*(r*pc)**3)))

        lnn.append(n*(pc**3))
        lclustersize.append(r)

    print ('Results for minimum, maximum cluster size\n' )
    print('Passage distance in Schwarszchild radius:', fmt(lrp[0]),',', fmt(lrp[-1]))
    print ('Encounter rate in cluster 1/t[Gyr⁻¹]=',fmt(lencratetot[0]*(s_yr*1E9)), ',',fmt(lencratetot[-1]*(s_yr*1E9)))

    print ('Encounters in t relax',fmt(lencrelaxt[0]), ',',fmt(lencrelaxt[-1]))
    print ('Encounters in lifetime',fmt(lnenclife[0]), ',',fmt(lnenclife[-1]))

    print('tmerger rin [days]', fmt(ltmerger[0]/(60*60*24)),'to', fmt(ltmerger[-1]/(60*60*24)))

    print('vdisp From', fmt(lvdisp[0]),'to', fmt(lvdisp[-1]), ' [km/s]')

    print('rperimax[km]:', fmt(lr_perimax[0]),'to', fmt(lr_perimax[-1]))

    print('a[AU]', fmt(la[0]),'to', fmt(la[-1]))
    print('e=', lecc[0],'to', lecc[-1])

###############
#    PLOTS    #
###############

    fig0=plt.figure()

    ax = fig0.add_subplot(2, 1, 1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.xlabel('Cluster Size [pc]')
    plt.ylabel('# encounters in cluster lifetime')
    plt.title('Number of Encounters in lifetime vs Cluster Size')
    plt.plot(lclustersize,lnenclife)

    plt.subplots_adjust(hspace=0.8)

    ax = fig0.add_subplot(2,1,2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.xlabel('Cluster Size [pc]')
    plt.ylabel('# encounters in a relaxation time')
    plt.title('Number of Encounters in t_relax vs Cluster Size')
    plt.plot(lclustersize,lencratetot)



    fig1=plt.figure()
    ax = fig1.add_subplot(3,1,1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.ylabel('Max Pericenter [km]')
    plt.xlabel('Relative velocity [km/s]')
    plt.title('Maximum pericenter Distance vs Relative Velocity')
    plt.plot(lvdisp,lr_perimax)

    plt.subplots_adjust(hspace=0.8)

    ax = fig1.add_subplot(3,1,2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.ylabel('Merger time [s]')
    plt.xlabel('Relative Velocity [km/s]')
    plt.title('Tmerger vs Relative velocity')
    plt.plot(lvdisp, ltmerger)

    ax = fig1.add_subplot(3,1,3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.ylabel('S[m²]')
    plt.xlabel('Relative Velocity [km/s]')
    plt.title('Cross Section (S) vs Relative Velocity')
    plt.plot(lvdisp,lS)
    fig0.savefig('./fig0.png')



    fig1.savefig('./fig1.png')

    plt.close(fig0)
    plt.close(fig1)

print ('encounters(N,m,rmin,rmax,b)')
print ('N: total number of objetcs \nm: mass of a single object in M_sol \nrmin: Smallest cluster radius in parsecs \nrmax: The bigger cluster radius in parsecs \nb: [0.1,1] passage distance in maximum pericenter distances')
print ('Example:  encounters(1E6,10,0.1,5.1,1)')
