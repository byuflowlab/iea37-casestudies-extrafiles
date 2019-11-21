from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
from math import radians
import scipy.io


font = {'size' : 8}

mpl.rc('font', **font)  # pass in the font dict as kwargs

"""colors"""
sky = '#375E97'
sunset = '#FB6542'
sunflower = '#FFBB00'
grass = '3F681C'

"""load data"""
Denver = scipy.io.loadmat('data/Denver_data.mat')
denver_dir =  Denver['windDirections'][0]
for i in range(len(denver_dir)):
    denver_dir[i] = np.deg2rad(denver_dir[i])
denver_freq =  Denver['windFrequencies'][0]
denver_speed =  Denver['windSpeeds'][0]
denver_dir = np.append(denver_dir,denver_dir[0])
denver_freq = np.append(denver_freq,denver_freq[0])
denver_speed = np.append(denver_speed,denver_speed[0])
denver_dir += np.deg2rad(270.)
denver_dir = denver_dir*-1.

Amalia = scipy.io.loadmat('data/Amalia_data.mat')
amalia_dir =  Amalia['windDirections'][0]
for i in range(len(amalia_dir)):
    amalia_dir[i] = np.deg2rad(amalia_dir[i])
amalia_freq =  Amalia['windFrequencies'][0]
amalia_speed =  Amalia['windSpeeds'][0]
amalia_dir = np.append(amalia_dir,amalia_dir[0])
amalia_freq = np.append(amalia_freq,amalia_freq[0])
amalia_speed = np.append(amalia_speed,amalia_speed[0])
amalia_dir += np.deg2rad(270.)
amalia_dir = amalia_dir*-1.

Redding = scipy.io.loadmat('data/Redding_data.mat')
redding_dir =  Redding['windDirections'][0]
for i in range(len(redding_dir)):
    redding_dir[i] = np.deg2rad(redding_dir[i])
redding_freq =  Redding['windFrequencies'][0]
redding_speed =  Redding['windSpeeds'][0]
redding_dir = np.append(redding_dir,redding_dir[0])
redding_freq = np.append(redding_freq,redding_freq[0])
redding_speed = np.append(redding_speed,redding_speed[0])
redding_dir += np.deg2rad(270.)
redding_dir = redding_dir*-1.

nDirections = len(denver_dir)
width = (2*np.pi) / nDirections
bottom = 0

fig, ((ax1, ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3, 2, figsize=[5.,6.], subplot_kw=dict(projection='polar'))
"""frequencies"""
ax1.plot(denver_dir, denver_freq, '-', color=sunflower, linewidth=3,zorder=0)
ax1.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'])
ax1.set_rgrids([0.0125,0.025], angle=-35.)
ax1.set_yticklabels(['1.25%','2.5%'],horizontalalignment='center')
ax1.set_title('wind rose',fontsize=8)


fig.text(0.03,0.86,'denver',rotation=90)
plt.tight_layout()

#plt.savefig('wind_roses.pdf')
plt.show()
