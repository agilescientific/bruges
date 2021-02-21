'''
===================
rpt.py
===================
Rock physics templates
Uses rock physics models defined in rpm.py
For Voigt-Reuss-Hill averages I define my own function here but I should
really change it and use the bounds already defined in bruges.
@author: Alessandro Amato del Monte
'''
import numpy as np
import matplotlib.pyplot as plt
from bruges.rockphysics import rpm
from bruges.rockphysics import fluidsub

def vrh(volumes,k,mu):
    '''
    Calculates Voigt-Reuss-Hill bounds.
    INPUT
    volumes: array with volumetric fractions
    k: array with bulk modulus
    mu: array with shear modulus
    OUTPUT
    k_u, k_l: upper (Voigt) and lower (Reuss) average for k
    mu_u, mu_l: upper (Voigt) and lower (Reuss) average for mu
    k0, mu0: Hill average of k and mu
    '''
    f=np.array(volumes).T
    k=np.resize(np.array(k),np.shape(f))
    mu=np.resize(np.array(mu),np.shape(f))
    ax=0 if f.ndim==1 else 1
    k_u = np.sum(f*k,axis=ax)
    k_l = 1./np.sum(f/k,axis=ax)
    mu_u = np.sum(f*mu,axis=ax)
    mu_l = 1./np.sum(f/mu,axis=ax)
    k0 = (k_u+k_l)/2.
    mu0 = (mu_u+mu_l)/2.
    return k_u, k_l, mu_u, mu_l, k0, mu0


def rpt(model='soft',vsh=0.0,fluid='gas',phic=0.4,Cn=8,P=10,f=1,cement='quartz'):
    if cement=='quartz':
        Kc, Gc = 37, 45
    elif cement=='calcite':
        Kc, Gc = 76.8, 32
    elif cement=='clay':
        Kc, Gc = 21, 7
    phi=np.linspace(0.1,phic-.1,6)
    sw=np.linspace(0,1,5)
    (Khc, Dhc) = (Kg, Dg) if fluid == 'gas' else (Ko,Do)
    K0,G0 = vrh([vsh, 1-vsh],[Ksh,Kqz],[Gsh,Gqz])[4:]
    D0 = vsh*Dsh+(1-vsh)*Dqz
    if model=='soft':
        Kdry, Gdry = rpm.softsand(K0,G0,phi,phic,Cn,P,f)
    elif model=='stiff':
         Kdry, Gdry = rpm.stiffsand(K0,G0,phi,phic,Cn,P,f)
    elif model=='cem':
         Kdry, Gdry = rpm.contactcement(K0,G0,phi,phic,Cn,Kc,Gc,scheme=2)
    elif model=='crit':
         Kdry, Gdry = rpm.critpor(K0,G0,phi,phic)

    xx=np.empty((phi.size,sw.size))
    yy=np.empty((phi.size,sw.size))

    for i,val in enumerate(sw):
        Kf = vrh([val,1-val],[Kb,Khc],[999,999])[1]
        Df = val*Db+(1-val)*Dhc
        vp,vs,rho,_= fluidsub.vels(Kdry,Gdry,K0,D0,Kf,Df,phi)
        xx[:,i]=vp*rho
        yy[:,i]=vp/vs

    plt.figure(figsize=(10,6))
    plt.plot(xx, yy, '-ok', alpha=0.3)
    plt.plot(xx.T, yy.T, '-ok', alpha=0.3)
    for i,val in enumerate(phi):
        plt.text(xx[i,-1],yy[i,-1]+.01,'$\phi={:.02f}$'.format(val), backgroundcolor='0.9')
    plt.text(xx[-1,0]-100,yy[-1,0],'$S_w={:.02f}$'.format(sw[0]),ha='right', backgroundcolor='0.9')
    plt.text(xx[-1,-1]-100,yy[-1,-1],'$S_w={:.02f}$'.format(sw[-1]),ha='right', backgroundcolor='0.9')
    plt.xlabel('Ip'), plt.ylabel('Vp/Vs')
    plt.xlim(xx.min()-xx.min()*.1,xx.max()+xx.max()*.1)
    plt.ylim(yy.min()-yy.min()*.1,yy.max()+yy.max()*.1)
    plt.title('RPT (N:G={0}, fluid={1})'.format(1-vsh, fluid))
