import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plotxlim = (0.5,2)
plotylim = (0.5,2)

#### plotting density - use KDE

## get all pieces together, then colors and naming consistent in plots
G5_Out = np.array(['GLUT5_out', 'Out Open', 'green', "Greens"])
G5_OutOcc = np.array(['GLUT5_out_occ', 'Out Occ.', 'grey', "Greys"])
G5_Occ = np.array(['GLUT5_occ', 'Occluded', 'red', "Reds"])
G5_InOcc = np.array(['GLUT5_in_occ', 'In Occ.', 'orange', "Oranges"])
G5_In = np.array(['GLUT5_in', 'In Open', 'blue', "Blues"])
all_sims = np.array([G5_Out, G5_OutOcc, G5_Occ, G5_InOcc, G5_In])

fig, ax = plt.subplots()

def plot_gate_dist(name, col, color_map):
    ec = np.loadtxt('../gate_dists/extracellular/%s.EC.xvg' %name)[:,1]
    ic = np.loadtxt('../gate_dists/intracellular/%s.IC.xvg' %name)[:,1]
    
    sns.kdeplot(ic,ec, shade=True, shade_lowest=False, alpha=0.5,legend=True, cbar=False, ax=ax, cmap=color_map)
   # b= ax.scatter(ic,ec, color=col, label=name, alpha = 1)   #needs to be after so that the spots come on top

start_EC = []
start_IC = []

for sim in all_sims:
    plot_gate_dist(sim[0], sim[2], sim[3])
    start_EC.append(np.loadtxt('../gate_dists/extracellular/%s.EC.starting_str.xvg' % sim[0])[1])
    start_IC.append(np.loadtxt('../gate_dists/intracellular/%s.IC.starting_str.xvg' % sim[0])[1])
    
plt.ylabel('Extracellular gate distance (nm)')
plt.xlabel('Intracellular gate distance (nm)')


colors = ["green", "grey", "red", "orange", "blue"]
texts = ["Out Open", "Out occ.", "Occ.", "In Occ.", "In Open"]

plt.scatter(start_IC, start_EC, color = colors)

patches = [ plt.plot([],[], marker="o", ms=10, ls="", mec=None, color=colors[i], alpha = 0.5, 
            label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]   #seaborn legends are weird... 
plt.legend(handles=patches)



plt.savefig('gate_distances_KDE.png', dpi =500)
plt.clf()




### plotting as a function of time - use contour
fig, axs = plt.subplots(5)
for n, sim in enumerate(all_sims):
    name = sim[0]
    ec = np.loadtxt('../gate_dists/extracellular/%s.EC.xvg' %name)[:,1]
    ic = np.loadtxt('../gate_dists/intracellular/%s.IC.xvg' %name)[:,1]
    times = np.loadtxt('../gate_dists/intracellular/%s.IC.xvg' %name)[:,0] / 1000
    t = axs[n].scatter(ic, ec, c=times, cmap = sim[3], alpha = 0.5)
    axs[n].set_xlim(plotxlim)
    axs[n].set_ylim(plotylim)
    
    axs[n].set_title(sim[0])
    fig.colorbar(t, ax = axs[n])
    
axs[2].set_ylabel('Extracellular gate (nm)')
#plt.ylabel('Extracellular gate (nm)')
plt.xlabel("intracellular gate (nm)")
fig.set_figheight(15)
plt.tight_layout()
plt.savefig('gate_distance_time.png', dpi = 1000)









