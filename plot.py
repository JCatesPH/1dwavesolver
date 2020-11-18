# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

initarr = pd.read_csv('initialdistribution.csv')
# %%
initarr.plot('x',['Re[u]','Im[u]'])
# %%
E_20 = pd.read_csv('data/Ep_20.csv')
E_20.plot('x',['Re[E]','Im[E]'])

# %%
E_480 = pd.read_csv('data/Ep_480.csv')
E_480.plot('x',['Re[E]','Im[E]'])
# %%
for i in range(0,1000,100):
    tmp = pd.read_csv('data/Ep_' + str(i) + '.csv')
    ax = tmp.plot('x',['Re[E]','Im[E]'])
    fig = ax.get_figure()
    fig.savefig('figs/Ep_' + str(i) + '.png')