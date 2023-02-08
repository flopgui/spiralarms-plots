import numpy as np
import astropy.units as u
##############################################
##############################################
# our standard cyl2cart are for right-handed. negative phi and Vphi for clockwise rotations
# origin of phi at X<0
# general rotation of coordinate system by an angle phi
def rot(phi,X,Y):
    s=np.sin(phi)
    c=np.cos(phi)  
    X2=X*c-Y*s
    Y2=X*s+Y*c
    return X2,Y2
##############################################
# from cartesian to cylindric
def cart2cyl(X,Y,Z,VX,VY,VZ):
    phi=np.arctan2(-Y,-X)
    R=np.sqrt(X*X+Y*Y)
    VR,Vphi=rot(-(np.pi+phi),VX,VY)
    return R,phi,Z,VR,Vphi,VZ
# from cylindric to cartesian
def cyl2cart(R,phi,Z,VR,Vphi,VZ):
    X=-R*np.cos(phi)
    Y=-R*np.sin(phi)
    VX,VY=rot(np.pi+phi,VR,Vphi)
    return X,Y,Z,VX,VY,VZ
##############################################
##############################################
#time factor to convert units of frecuencies in km/s/kpc to Gyr
# needs to be multiplied to the frecuency in usual units to obtain Gyr
# or devided to times in Gyr to obatin frecuency units inside cos and sin
#ft=2.*np.pi*1000.*3.086*10.**(13.)/(31556926*10.**(9.)) factor 2pi was wrong !!!
#ft=1000.*3.08567758*10.**(13.)/(31556926*10.**(9.)) # OK BUT WE DO IT WITH ASTROPY IN THE END, value slighlty different
#0.9778131051167658
ft=(1*u.km/u.second/u.kpc).to(1./u.Gyr).value
#1.0227121650456952
##############################################
##############################################
def Vf(r,Vc=200.,n=-10.,hr=8.):
    if n<=99:V=Vc*(r/hr)**(1./n)
    if n>99:V=r*0.+Vc
    return V
def omegaf(r,Vc=200.,n=-10.,hr=8.):
    return Vf(r,Vc,n,hr)/r
def kappaf(r,Vc=200.,n=-10.,hr=8.):
    if n<=99:
        k2=2.*(1.+1./n)*(Vc**2./hr**2.)*(r/hr)**(2./n-2.)
        k=np.sqrt(k2)
    if n>99:k=np.sqrt(2)*Vc/r
    return k
def gammaf(r,Vc=200.,n=-10.,hr=8.):
    if n<=99:g=2.*omegaf(r,Vc,n,hr)/kappaf(r,Vc,n,hr)
    if n>99:g=np.sqrt(2)
    return g
##############################################
##############################################
def orbit(Rg,A,phi0,psi,kappa,omega,gamma,t):
    
    #print((Rg,A,psi))
    
    #R=Rg-A*Rg*np.sin(kappa*t+psi)#struck
    #VR=-A*Rg*kappa*np.cos(kappa*t+psi)   #struck 
    #phi=phi0+omega*t+gamma*A*(np.cos(kappa*t+psi)-np.cos(psi))#struck
    R=Rg+A*Rg*np.cos(kappa*t+psi)
    VR=-A*Rg*kappa*np.sin(kappa*t+psi)    
    phi=phi0+omega*t-gamma*A*(np.sin(kappa*t+psi)-np.sin(psi))#defined positive in the sense of rotation 
    Vphi=(Rg**2.)*omega/R#defined positive in the sense of rotation 
    phi2=phi-omega*t    
    
    #phi=phi0+omega*t-gamma*A*(np.sin(kappa*t+psi))#defined positive in the sense of rotation 
    #vphi=omega*rg-gamma*A*kappa*rg*np.cos(kappa*t+psi)#wrong!!

    return R,VR,phi,Vphi,phi2
##############################################
##############################################

def impact2d(R0,phi0,D=20.,Dve=10.):#Struck et al. 2011

    #print('in')
    DVR=Dve*R0*np.cos(2.*phi0)/D
    DVphi=-Dve*R0*np.sin(2.*phi0)/D
    #print('vrvphi',DVR,DVphi)

    return DVR,DVphi

def impact(R0,phi0,rper=1.,Dve=10.,epsilon=1.,m=0,Vc=200,n=-10.,hr=8.):
    dist=(R0**2.+rper**2.-2.*R0*rper*np.cos(phi0))**(0.5)
    sinalpha=rper*np.sin(phi0)/dist
    alpha=np.arcsin(sinalpha)
    
    DV=Dve*(dist/epsilon)**(-m)    
    DVR=DV*np.cos(alpha)
    DVphi=-DV*np.sin(alpha)
    #print(DVR)
    
    return DVR,DVphi
    

##############################################
##############################################    
def IC(R0,DVR,DVphi,Vc=200,n=-10.,hr=8.):

    Rg=R0*(1.+((R0/hr)**(-1./n))*DVphi/Vc)**(n/(n+1.)) #ARREGLAR AIXO 
    #print(Rg)
    if n>99.:Rg=R0*(1+DVphi/Vc)
    #print(Rg)
    #delta=(Rg-R0)/Rg #struck
    delta=-(Rg-R0)/Rg

    kappa=kappaf(Rg,Vc,n,hr)
    omega=omegaf(Rg,Vc,n,hr)
    gamma=gammaf(Rg,Vc,n,hr)
    
    #A
    A=np.sqrt(delta**2.+(DVR/(Rg*kappa))**2.)
    #print(A)
    
    #psi
    #sinpsi=delta/A#struck
    #cospsi=-DVR/(A*Rg*kappa)#struck
    sinpsi=-DVR/(A*Rg*kappa)
    cospsi=delta/A
    #psi=np.arcsin(sinpsi) 
    #print('psi1',np.rad2deg(psi))
    #psi=np.arccos(cospsi) 
    #print('psi2',np.rad2deg(psi))
    psi=np.arctan2(sinpsi,cospsi) 
    #print('psi2',np.rad2deg(psi))
    
    omega0=omegaf(R0,Vc,n,hr)
    #phase=R0*omega0+DVphi-Rg*omega+(Rg-R0)*gamma*kappa
    #print(Vphi00)
    #Vphi00=0.
    #print(Vphi00)
    return Rg,A,psi,kappa,omega,gamma
    
#for flat rotation curve. no need to run impact; needs to be revised!!!!!!!!!
def IC0(R0,phi0,t,rper=2.,Dve=10.,Vc=200.):
    
    dist=(R0**2.+rper**2.-2.*R0*rper*np.cos(phi0))**(0.5)
    sinalpha=rper*np.sin(phi0)/dist
    alpha=np.arcsin(sinalpha)
    
    #psi
    DVc=Dve/Vc
    tanpsi=np.sqrt(2)*np.tan(alpha)/(1-DVc*np.sin(alpha))
    psi=np.arctan(tanpsi) 
    #print(psi)
    
    #A
    At1=np.sin(alpha)/(1.-DVc*np.sin(alpha))**2.
    At2=(np.cos(alpha)/np.sqrt(2))**2.
    A=DV*np.sqrt(At1+At2)
    #print(A)
    
    Rg=R0*(1.-DVc*np.sin(alpha))
    
    kappa=kappaf(Rg,Vc=Vc,n=100)
    omega=omegaf(Rg,Vc=Vc,n=100)
    gamma=gammaf(Rg,Vc=Vc,n=100)
    #print(r0,rg,kappa,A*rg)
    
    return Rg,A,psi,kappa,omega,gamma
     