# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

initarr = pd.read_csv('initialdistribution.csv')
# %%
initarr.plot('x',['Re[E]'])

#%%
#initarr.plot('x',['Im[E]'])

Emax = np.amax(np.abs(initarr['Re[E]'].to_numpy()))
#Emax = 3e9

# %%
for i in range(0,10000,200):
    tmp = pd.read_csv('data/Ep_' + str(i) + '.csv')
    xarr = tmp['x'].to_numpy()
    Ereal = tmp['Re[E]'].to_numpy()
    
    fig, axs = plt.subplots(1,1,figsize=(6,4))
    axs.plot(xarr,Ereal)
    axs.set_ylim(-3*Emax,3*Emax)
    ti = "{:.1E} s".format(i*1e-15)
    axs.set_title(r'|E| for $t=$' + ti)

    fig.savefig('figs/Ep_' + str(i) + '.png')

   

# %%
