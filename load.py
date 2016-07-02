import sys
import sdf
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

sigma=16
xx=np.arange(-40,41)
gaussf=np.exp(-0.5*xx**2/sigma**2)
gaussf=gaussf/sum(gaussf)
n_files=140
n_files=n_files+1
step=1

time = 0

max_fwd=np.arange(0,n_files+1)
max_bwd=np.arange(0,n_files+1)

for n in range(0, n_files):
    num='{:04}'.format(n)
    fname = "../Data/"+num+".sdf"
    l = sdf.read(fname)
    print(fname)
    
    q_e=1.602e-19
    m_e=9.109e-31
    c=3.0e8
    lambda1=800.0e-9
    lambda2=920.0e-9
    epsilon0=8.8541e-12
    
    
    ey_fwd=l.Electric_Field_Ey.data*q_e/m_e/c/(2*np.pi*c/lambda1)
    bz_fwd=l.Magnetic_Field_Bz.data*q_e/m_e/(2*np.pi*c/lambda1)
    ey_bwd=l.Electric_Field_Ey.data*q_e/m_e/c/(2*np.pi*c/lambda2)
    bz_bwd=l.Magnetic_Field_Bz.data*q_e/m_e/(2*np.pi*c/lambda2)
    
    bz_fwd=0.5*(bz_fwd[0:len(bz_fwd)-1]+bz_fwd[1:len(bz_fwd)])
    ey_fwd=ey_fwd[1:len(ey_fwd)]
    bz_bwd=0.5*(bz_bwd[0:len(bz_bwd)-1]+bz_bwd[1:len(bz_bwd)])
    ey_bwd=ey_bwd[1:len(ey_bwd)]

    x_tuple=l.Grid_Grid.data
    for m in x_tuple:
        x = m
    x1 = 1e6*0.5*(x[1:len(x)-1]+x[2:len(x)]);
    
    fwd=(ey_fwd+bz_fwd)/2
    fwd=fwd/(q_e/m_e/c/(2*np.pi*c/lambda1))
    bwd=(ey_bwd-bz_bwd)/2
    bwd=bwd/(q_e/m_e/c/(2*np.pi*c/lambda2))
    
    
    fwd_intensity=0.5*epsilon0*c/10000*fwd**2
    bwd_intensity=0.5*epsilon0*c/10000*bwd**2
    
    fwd_intensity=np.convolve(abs(fwd_intensity),gaussf,'same')
    bwd_intensity=np.convolve(abs(bwd_intensity),gaussf,'same')
    
    
    max_fwd[step]=max(fwd_intensity)
    max_bwd[step]=max(bwd_intensity)
    step=step+1
    tt = 't = '+str(time)+'s'
    
    
    fig1=plt.figure()
    axes = fig1.add_axes([0.1, 0.1, 0.8, 0.8]) 
    axes.set_xlabel('x (micrometer)',fontsize=16)
    axes.set_ylabel('I (W/cm^2)',fontsize=16)
    axes.plot(x1,fwd_intensity, x1, bwd_intensity)
    axes.set_ylim([0,3e15])
    axes.set_title(tt, family="serif", fontsize=16)
    fig1.savefig(num+"_intensity.pdf")
    fig1.clf()
    plt.close()

    time = time + 50e-15
    
    
t=np.linspace(0,14e-12,step)
fig2=plt.figure()
axes = fig2.add_axes([0.1, 0.1, 0.8, 0.8]) 
axes.set_xlabel('x (micrometer)',fontsize=16)
axes.set_ylabel('I (W/cm^2)',fontsize=16)
axes.plot(t,max_bwd, t, max_fwd);
fig2.savefig("vyvoj.pdf")
