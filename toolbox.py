from scipy.integrate import odeint
import pickle
from pylab import *
from scipy.fftpack import dst
import matplotlib
import matplotlib.animation as animation
from scipy.optimize import *
import numpy.linalg

def xPos(I,iter1):
    return -I['dx']*(I['nx']-1.)/2. + iter1%I['nx']*I['dx']

def yPos(I,iter1):
    return -I['dy']*(I['ny']-1.)/2. + int(iter1/I['nx'])%I['ny']*I['dy']

def zPos(I,iter1):
    if iter1<I['n']:
        z0=-I['L'] + int(iter1/(I['nx']*I['ny']))*I['dz']
        if I['curved_lay']==0:
            return z0
        else:
            r2=xPos(I,iter1)**2+yPos(I,iter1)**2
            zR=I['zR']
            a=(-9.*r2*z0/4. + z0**3 + 9.*z0*zR**2 + 0.5*sqrt(0.5*(3*r2-2.*z0**2+6.*zR**2)**3+(4.5*r2*z0-2.*(z0**3+9.*z0*zR**2))**2))**(1./3.)
            return z0/3.+(-1.5*r2+z0**2-3.*zR**2)/(3.*a)+a/3.
    else:
        z0= -I['L'] + int(iter1/I['nx']/I['ny'])%I['nz']*I['dz'] 
        if I['curved_lay']==0:
            return -z0
        else:
            r2=xPos(I,iter1)**2+yPos(I,iter1)**2
            zR=I['zR']
            a=(-9.*r2*z0/4. + z0**3 + 9.*z0*zR**2 + 0.5*sqrt(0.5*(3*r2-2.*z0**2+6.*zR**2)**3+(4.5*r2*z0-2.*(z0**3+9.*z0*zR**2))**2))**(1./3.)
        return -(z0/3.+(-1.5*r2+z0**2-3.*zR**2)/(3.*a)+a/3.)

def GaussianBeam(omega0, lambda0,r,z):
    zR=pi*omega0**2/lambda0
    omegaz=omega0*sqrt(1+(z/zR)**2)
    Rz=z*(1+(zR/z)**2)
    psiz=math.atan(z/zR)
    return (omega0/omegaz)*exp(-r**2/omegaz**2)*exp(-1j*psiz)*exp(2*pi*1j*(z+r**2/(2*Rz))/lambda0)

def GreensTensor(lambda0, r,z):
    R=sqrt(r**2+z**2)
    k0R=2.*pi*R/lambda0
    return exp(1j*k0R)*(1.+1j/k0R-1./k0R**2-(1.+3*1j/k0R-3./k0R**2)*r**2/(2*r**2+2*z**2))/(4.*pi*R)
    
def ParseElectricField(I, S, zmin, zmax, nbz, xmin, xmax, nbx, y, c):
    zvec=linspace(zmin,zmax,nbz);
    xvec=linspace(xmin,xmax,nbx);
    
    EF=zeros((nbz,nbx))+1j*zeros((nbz,nbx))
    for iterz in range(nbz):
        for iterx in range(nbx):
            for itern in range(I['n']*I['N']):
                rdif=sqrt((y-yPos(I,itern))**2+(xvec[iterx]-xPos(I,itern))**2)
                zdif=zvec[iterz]-zPos(I,itern)
                EF[iterz,iterx]+=GreensTensor(I['lambda0'],rdif,zdif)*S['c'][itern]/I['lambda0']**2
                
    return abs(EF)           
    

def dphi(phir, t, I, S):
    nh = I['N']*I['n']
    phic = phir[:nh+2] + 1j * phir[nh+2:]    
    
    dphic=1j*zeros(2*I['n']+2)
    dphic[0:I['n']] = -1j * (S['h0']*get_Omega1(t, I, S))*phic[2*I['n']]
    dphic[2*I['n']] = -1j * (conj(S['h0'])*get_Omega1(t, I, S)).dot(phic[0:I['n']])
    dphic[I['n']:2*I['n']] = -1j * (S['h1']*get_Omega2(t, I, S))*phic[2*I['n']+1]
    dphic[2*I['n']+1] = -1j * (conj(S['h1'])*get_Omega2(t, I, S)).dot(phic[I['n']:2*I['n']])
    dphic[2*I['n']] += -1j * (-get_Omega1(t, I, S)**2*S['Stark'])*phic[2*I['n']]
    dphic[2*I['n']+1] += -1j * (-get_Omega2(t, I, S)**2*S['Stark'])*phic[2*I['n']+1]
    dphic[0:nh] += -1j * S['h2'].dot(phic[0:nh])
    dphir = zeros(2*(nh+2))
    dphir[:nh+2] = real(dphic)
    dphir[nh+2:] = imag(dphic)
    return dphir


def get_Omega1(t, I, S):
    gamma =  exp(I['gammaeff_max'] * t) / (2 - exp(I['gammaeff_max'] * t)) * (t < 0) + 1 * (t >= 0)
#    gamma=1
    return sqrt(gamma)


def get_Omega2(t, I, S):
    return get_Omega1(- t, I, S)


def evol(I):
    S = {}
    ## 3 HAMILTONIANS
    c1=1j * np.zeros((I['N']*I['n']))
    c2=1j * np.zeros((I['N']*I['n']))
    for itern in range(I['N']*I['n']):
        if itern<I['n']:
            c1[itern] = I['c'][itern]
        elif itern<2*I['n']:  
            c2[itern] = I['c'][itern]   
    
    nh = I['N']*I['n'] #number of non-excited states
    
    if I['full_inv']==1:
        Hnh = 1j * np.zeros((nh,nh))

        # Detuning and decays of emitting states        
        for iter1 in range(I['N']*I['n']):
                x1=xPos(I,iter1)
                y1=yPos(I,iter1)
                z1=zPos(I,iter1)
                Hnh[iter1,iter1]= - I['Delta'] - 1j*I['gamma']/2.
                for iter2 in range(iter1):
                    x2=xPos(I,iter2)
                    y2=yPos(I,iter2)
                    z2=zPos(I,iter2)
                    Jhop=1.5*I['gamma']*I['lambda0']*GreensTensor(I['lambda0'],sqrt((x2-x1)**2+(y2-y1)**2),z2-z1)
                    J=-real(Jhop)
                    Gamma=2.*imag(Jhop)

                    Hnh[iter1,iter2]+= J-1j*Gamma/2.
                    Hnh[iter2,iter1]+= J-1j*Gamma/2.
        
        S['selfE1']= (linalg.solve(Hnh,c1)).dot(conj(c1))
        S['selfE2']= (linalg.solve(Hnh,c2)).dot(conj(c2))
        S['gammaR']= (linalg.solve(Hnh,c1)).dot(conj(c2))
        S['gammaL']= (linalg.solve(Hnh,c2)).dot(conj(c1))
        
    else:
        Hnh_inv = 1j * np.zeros((nh,nh))
      
        for iter1 in range(I['N']*I['n']):
                x1=xPos(I,iter1)
                y1=yPos(I,iter1)
                z1=zPos(I,iter1)
                Hnh_inv[iter1,iter1]= - 1j*I['gamma']/(2.*I['Delta']**2)
                for iter2 in range(iter1):
                    x2=xPos(I,iter2)
                    y2=yPos(I,iter2)
                    z2=zPos(I,iter2)
                    Jhop=1.5*I['gamma']*I['lambda0']*GreensTensor(I['lambda0'],sqrt((x2-x1)**2+(y2-y1)**2),z2-z1)
                    J=-real(Jhop)
                    Gamma=2.*imag(Jhop)

                    Hnh_inv[iter1,iter2]+= (J-1j*Gamma/2)/I['Delta']**2.
                    Hnh_inv[iter2,iter1]+= (J-1j*Gamma/2)/I['Delta']**2.
        
        S['selfE1']= (Hnh_inv.dot(c1)).dot(conj(c1))
        S['selfE2']= (Hnh_inv.dot(c2)).dot(conj(c2))
        S['gammaR']= (Hnh_inv.dot(c1)).dot(conj(c2))
        S['gammaL']= (Hnh_inv.dot(c2)).dot(conj(c1))
    
    S['c']=c1-c2
    return S

def GaussianOpt(omega0,I):
    I['zR']=pi*omega0**2/I['lambda0']
    I['c']= 1j*np.zeros((I['N']*I['n']))
    for itern in range(I['N']*I['n']):
        I['c'][itern] = GaussianBeam(omega0, I['lambda0'],sqrt(xPos(I,itern)**2+yPos(I,itern)**2),zPos(I,itern))        
    S=evol(I)
    return 1-abs(S['gammaR'])**2/(2*imag(S['selfE1'])*2*imag(S['selfE2']))

def FullOpt(c,I):
    I['c']=c[:I['N']*I['n']]+1j*c[I['N']*I['n']:2*I['N']*I['n']]
    S = evol(I)
    return 1-abs(S['gammaR'])**2/(2*imag(S['selfE1'])*2*imag(S['selfE2']))

def mySolver(I0,gaussian_opt,curved_lay,full_inv):

    I=I0.copy()
    
    I['curved_lay']=curved_lay
    I['full_inv']=full_inv
    
    omega0= sqrt(I['L']*I['lambda0']/pi) # waist
    I['zR']=pi*omega0**2/I['lambda0']
           
    if gaussian_opt:
        Fopt=minimize(GaussianOpt,omega0, args=I, method='SLSQP')
        I['zR']=pi*Fopt.x**2/I['lambda0']
        I['c']= 1j*np.zeros((I['N']*I['n']))
        for itern in range(I['N']*I['n']):
            I['c'][itern] = GaussianBeam(Fopt.x, I['lambda0'],sqrt(xPos(I,itern)**2+yPos(I,itern)**2),zPos(I,itern))
        
        S=evol(I)
        S['omega0']=Fopt.x
        
    else:
        c= np.zeros((2*I['N']*I['n']))
        for itern in range(I['N']*I['n']):
            field = GaussianBeam(omega0, I['lambda0'],sqrt(xPos(I,itern)**2+yPos(I,itern)**2),zPos(I,itern)) 
            c[itern] = real(field)
            c[itern+I['N']*I['n']] = imag(field)
        Fopt=minimize(FullOpt, c,args=I, method='SLSQP')
        I['c']=Fopt.x[:I['N']*I['n']]+1j*Fopt.x[I['N']*I['n']:2*I['N']*I['n']]
        S = evol(I)
        S['omega0']=omega0
    
    #print(I['curved_lay'])
    return S,I
