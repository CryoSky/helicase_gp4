#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import matplotlib.tri as tri

b=pd.read_csv('pmf-mask-t6.dat', header=None, delim_whitespace=True, comment='#')

b.columns = ['beta', 'alpha', 'pmf','prob']
a=b[b.pmf<9999]

#b.pmf[b.pmf>9999] = np.max(a.pmf) + 1
#levels=np.arange(0.0, np.max(a.pmf) + 0.4, 0.4)
#b.pmf[b.pmf>9999] = 20
b.pmf[b.pmf>20] = 20
levels=np.arange(0.0, 20, 0.5)


plt.figure(figsize=(9,7.5))
ax = plt.gca()
CS = ax.tricontourf(b.beta, b.alpha, b.pmf, levels, cmap=plt.get_cmap('jet'))
cbar = plt.colorbar(CS, ax=ax)
cbar.ax.tick_params(labelsize=18)

plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)

plt.xlabel("Q_beta", fontsize=24)
plt.ylabel("Q_alpha", fontsize=24)



#y1=np.max(a.tc)+50
#if y1>350: y1=350
#
#plt.ylim(0,y1)
##plt.show()
plt.savefig('pmf-p.png')

