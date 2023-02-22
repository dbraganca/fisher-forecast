import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.special import legendre
import seaborn as sns
from pathlib import Path
import sys
from classy import Class
import copy
import itertools
import pickle

kk = np.logspace(-5, 0, 200)
#the formula for omega_cdm is such that we change Omegam by dOm
"""def cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,dlnAs=0,dns=0,dh=0,dOm=0,dob=0,dmnu=0):
    return {'omega_b': obback*(1+dob),
          'omega_cdm':Omback*(1+dOm)*hback**2-obback*(1+dob)-mnuback*(1+dmnu)/93.14+2*Omback*hback**2*dh,
          'h': hback*(1+dh),
          'ln10^{10}A_s': lnAsback+dlnAs,
          'n_s': nsback*(1+dns),
          'output': 'mPk,mTk',
          'P_k_max_h/Mpc': 1,
          'z_pk': zpk,
          'N_ncdm':1,
          'N_ur': 2.0328,
          'm_ncdm':mnuback*(1+dmnu),#0.15#0.06
          'T_ncdm': 0.71611
         }"""

def cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,dlnAs=0,dns=0,dh=0,dOm=0,dob=0,dmnu=0,dOk=0):
    return {'omega_b': obback*(1+dob),
          'omega_cdm':Omback*(1+dOm)*hback**2-obback*(1+dob)-mnuback*(1+dmnu)/93.14+2*Omback*hback**2*dh,
          'h': hback*(1+dh),
          'ln10^{10}A_s': lnAsback+dlnAs,
          'n_s': nsback*(1+dns),
          'output': 'mPk,mTk',
          'P_k_max_h/Mpc': 1,
          'z_pk': zpk,
          'N_ncdm':1,
          'N_ur': 2.0328,
          'm_ncdm':mnuback*(1+dmnu),#0.15#0.06
          'T_ncdm': 0.71611,
          'Omega_k':Okback+dOk
         }

def getD(dic):
    M = Class()
    M.set(dic)
    M.compute()
    zpk = dic['z_pk']
    
    # k in h/Mpc
    
    # P(k) in (Mpc/h)**3
    dd= M.scale_independent_growth_factor(dic['z_pk'])
    return  dd


def getpk(dic):
    M = Class()
    M.set(dic)
    M.compute()
    # k in h/Mpc
    zpk = dic['z_pk']
    
    # P(k) in (Mpc/h)**3
    Pk = np.array([M.pk(ki*M.h(), zpk)*M.h()**3 for ki in kk])
    return  Pk

def getf(dic):
    M = Class()
    M.set(dic)
    M.compute()
    zpk = dic['z_pk']
    
    # k in h/Mpc
    
    # P(k) in (Mpc/h)**3
    ff= M.scale_independent_growth_factor_f(dic['z_pk'])
    return  ff

def getT(dic):
    M = Class()
    M.set(dic)
    M.compute()
    zpk = dic['z_pk']

    Dg = M.scale_independent_growth_factor(zpk)/((1+100)*M.scale_independent_growth_factor(100))
    Possionpre =  (3.e5)**2/(1.5*100.**2*M.Omega0_m())

    
    tdic = M.get_transfer(zpk)
    kpsi= tdic['k (h/Mpc)']
    Tpsi = tdic['psi']/tdic['psi'][0]
    return kpsi,Tpsi*kpsi**2*Possionpre*Dg


def getalldics(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,dlnAs=0.1,dns=0.02,dh=0.02,dOm=0.04,dob=0.02,dmnu=0.25,dOk=0.002):
    back = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,0,0,0)
    Asp = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,dlnAs,0,0,0,0,0,0)
    Asl = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,-dlnAs,0,0,0,0,0,0)
    nsp = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,dns,0,0,0,0,0)
    nsl = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,-dns,0,0,0,0,0)
    hp = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,dh,0,0,0,0)
    hl = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,-dh,0,0,0,0)
    Omp = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,dOm,0,0,0)
    Oml = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,-dOm,0,0,0)
    obp = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,dob,0,0)
    obl = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,-dob,0,0)
    nup = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,0,dmnu,0)
    nul = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,0,-dmnu,0)
    Okp = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,0,0,dOk)
    Okl = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,0,0,-dOk)
    allp =  np.array([back,Asp,Asl,nsp,nsl,hp,hl,Omp,Oml,obp,obl,nup,nul,Okp,Okl])
    return allp


def getalldicsforf(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,dlnAs=0.1,dns=0.02,dh=0.02,dOm=0.04,dob=0.02,dmnu=0.25,dOk=0.002):
    back = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,0,0,0)
    hp = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,dh,0,0,0,0)
    hl = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,-dh,0,0,0,0)
    Omp = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,dOm,0,0,0)
    Oml = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,-dOm,0,0,0)
    obp = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,dob,0,0)
    obl = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,-dob,0,0)
    nup = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,0,dmnu,0)
    nul = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,0,-dmnu,0)
    Okp = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,0,0,dOk)
    Okl = cosmicdicn(zpk,lnAsback,nsback,hback,Omback,obback,mnuback,Okback,0,0,0,0,0,0,-dOk)
    allp =  np.array([back,Omp,Oml,obp,obl,hp,hl,nup,nul,Okp,Okl])
    return allp


def getallpks(dics):
    allp=[]
    allp.append(kk)
    for i in range(len(dics)):
        allp.append(getpk(dics[i]))
    return np.array(allp)

def getallTs(dics):
    allT=[]
    for i in range(len(dics)):
        ks,Ta = getT(dics[i])
        allT.append(ks)
        allT.append(Ta)
    return np.array(allT)

def getallfs(dics):
    allf=[]
    for i in range(len(dics)):
        allf.append(getf(dics[i]))
    return np.array(allf)

names = np.array(['fid','lnAs_pl','lnAs_min','ns_pl','ns_min','h_pl','h_min','Om_pl','Om_min','ob_pl','ob_min','nu_pl','nu_min','Ok_pl','Ok_min'])

def savepk(pks,folder):
    os.mkdir(folder)
    for i in range(len(pks)-1):
        np.savetxt(os.path.join(folder, 'pk_'+names[i] + '.dat'), np.transpose(np.array([pks[0],pks[i+1]])))

def savepkandf(pks,fs,folder):
    os.mkdir(folder)
    for i in range(len(pks)-1):
        np.savetxt(os.path.join(folder, 'pk_'+names[i] + '.dat'), np.transpose(np.array([pks[0],pks[i+1]])))
    np.savetxt(os.path.join(folder, 'allfs.dat'),fs)


def savepkandfandT(pks,Ts,fs,folder):
    os.mkdir(folder)
    for i in range(len(pks)-1):
        np.savetxt(os.path.join(folder, 'pk_'+names[i] + '.dat'), np.transpose(np.array([pks[0],pks[i+1]])))
        np.savetxt(os.path.join(folder, 'Tk_'+names[i] + '.dat'), np.transpose(np.array([Ts[2*i],Ts[2*i+1]])))
    np.savetxt(os.path.join(folder, 'allfs.dat'),fs)

def allsaved(dics,dicsf,folder):
    fs = getallfs(dicsf)
    pks = getallpks(dics)
    savepkandf(pks,fs,folder)

def allsavedT(dics,dicsf,folder):
    fs = getallfs(dicsf)
    pks = getallpks(dics)
    Ts = getallTs(dics)
    savepkandfandT(pks,Ts,fs,folder)


