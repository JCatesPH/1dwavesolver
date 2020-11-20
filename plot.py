# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

initarr = pd.read_csv('initialdistribution.csv')
# %%
initarr.plot('x',['Re[E]','Im[E]'])

Emax = np.amax(np.abs(initarr['Re[E]'].to_numpy()))
#Emax = 3e9

# %%
for i in range(0,10000,200):
    tmp = pd.read_csv('data/Ep_' + str(i) + '.csv')
    xarr = tmp['x'].to_numpy()
    Ereal = tmp['Re[E]'].to_numpy()
    Eimag = tmp['Im[E]'].to_numpy()
    Earr = Ereal + 1j * Eimag
    Ek = np.fft.fft(Earr)
    freq = np.fft.fftfreq(Ek.shape[-1])
    
    fig, axs = plt.subplots(1,2,figsize=(12,4))
    axs[0].plot(xarr,Ereal,xarr,Eimag)
    axs[0].set_ylim(-3*Emax,3*Emax)
    axs[0].legend(['Re[E(x)]','Im[E(x)]'])
    ti = "{:.1E} s".format(i*1e-15)
    axs[0].set_title(r'|E| for $t=$' + ti)

    axs[1].plot(freq, np.real(Ek), freq, np.imag(Ek))
    axs[1].set_title(r'$E(k)$ for $t=$' + ti)
    axs[1].legend(['Re[E(k)]','Im[E(k)]'])

    fig.savefig('figs/Ep_' + str(i) + '.png')

   

# %%
