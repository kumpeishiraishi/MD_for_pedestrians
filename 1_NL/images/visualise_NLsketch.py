import numpy as np
import matplotlib.pyplot as plt

a = np.load('./conf.npz')
conf = a['position']
sigma = a['sigma']
N = sigma.size
Lbox = (N/1.0)**(1/2)

Lfig = 5
fig, ax = plt.subplots(
    figsize=(Lfig, Lfig),
    # dpi=480
)
points_whole_ax = Lfig * 72

points_radius = 2.3e-2*sigma / 1.0 * points_whole_ax
ax.scatter(conf[:,0], conf[:,1], s=points_radius**2, 
           color='darkgray', linewidth=0)
#########################################################################################
index = 100
vec = conf - conf[index,:]
vec -= Lbox * np.floor(vec/Lbox + 0.5)
distance = np.linalg.norm(vec, axis=1)

thres_large = 6
index = distance<thres_large
ax.scatter(conf[index,0], conf[index,1], s=points_radius[index]**2, 
           color='mediumseagreen', linewidth=0)

thres_small = 4
index = distance<thres_small
ax.scatter(conf[index,0], conf[index,1], s=points_radius[index]**2, 
           color='teal', linewidth=0)

#########################################################################################

#########################################################################################
index = 100
points_radius = 2.6e-2*sigma / 1.0 * points_whole_ax
ax.scatter(conf[index,0], conf[index,1], s=points_radius[index]**2, 
           color='red', linewidth=0)
for thres in [thres_large, thres_small]:
    my_ls = 'solid'
    if thres>thres_small:
        my_ls = 'dotted'
    ax.add_patch(
        plt.Circle((conf[index,0], conf[index,1]), radius=thres, fill=False, ls=my_ls, lw=3)
    )
#########################################################################################
ax.set_xlim(-Lbox/2, Lbox/2)
ax.set_ylim(-Lbox/2, Lbox/2)
ax.set_aspect('equal', adjustable='box')

ax.set_xticks([])
ax.set_yticks([])

ax.set_facecolor('lemonchiffon')

plt.tight_layout(pad=0.25)
plt.savefig('NL.png')
plt.close()
