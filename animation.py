#%%
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')

#%%

initarr = pd.read_csv('initialdistribution.csv')
initarr.plot('x',['Re[E]'])
Emax = np.amax(np.abs(initarr['Re[E]'].to_numpy()))
xmin = np.amin(initarr['x'].to_numpy())
xmax = np.amax(initarr['x'].to_numpy())

#%%
fig = plt.figure()
ax = plt.axes(xlim=(xmin, xmax), ylim=(-3*Emax, 3*Emax),
    xlabel="x [m]", ylabel="Re[E(x,t)]")
line, = ax.plot([], [], lw=3) 
text = ax.text(0.5*xmax, -2.6*Emax, '')
#%%
def init():
    line.set_data([], [])
    text.set_text('')
    return line, text
def animate(i):
    tmp = pd.read_csv('data/Ep_' + str(i) + '.csv')
    x = tmp['x'].to_numpy()
    E = tmp['Re[E]'].to_numpy()
    line.set_data(x, E)

    tstr = "t = {:.3E} s".format(i*1e-13)
    text.set_text(tstr)
    return line, text

anim = FuncAnimation(fig, animate, init_func=init,
                            frames=range(0,10000,200), interval=200, 
                            blit=True)

#%%
anim.save('wave.gif', savefig_kwargs={'transparent' : False}, dpi=50)
