# Written by PJ Stanley 2017
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def power_func(ws):
    cin = 4.
    cout = 25.
    rws = 11
    rp = 10

    if ws <= cin:
        p = 0.
    elif ws >= cout:
        p = 0.
    elif ((ws >= cin) and (ws <= rws)):
        p = rp*((ws-cin)/(rws-cin))**3
    elif ((ws >= rws) and (ws <= cout)):
        p = rp

    return p

color = (0,0.4470,0.7410)
ws = np.linspace(0.,30.,1000)
p = np.zeros(1000)
for i in range(1000):
    p[i] = power_func(ws[i])

font = {'family' : 'normal',
        'size'   : 10}

matplotlib.rc('font', **font)

fig, ax = plt.subplots(figsize=(3.5, 2.5))
ax.plot(ws,p,linewidth=2,color=color)
ax.set_ylim(-0.1,10.1)
ax.set_xlim(0, 30)
ax.set_xlabel('V (m/s)',fontsize=10)
ax.set_ylabel('P (MW)',fontsize=10)

# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

plt.tight_layout()
plt.savefig('power_curve.pdf',transparent=True)
plt.show()
