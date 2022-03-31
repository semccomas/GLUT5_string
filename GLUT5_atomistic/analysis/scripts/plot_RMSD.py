import numpy as np
import matplotlib.pyplot as plt

G_in = np.loadtxt('../input_f/GLUT5_in.rmsd.xvg')[:,1]
G_in_occ = np.loadtxt('../input_f/GLUT5_in_occ.rmsd.xvg')[:,1]
G_occ = np.loadtxt('../input_f/GLUT5_occ.rmsd.xvg')[:,1]
G_out_occ = np.loadtxt('../input_f/GLUT5_out_occ.rmsd.xvg')[:,1]
G_out = np.loadtxt('../input_f/GLUT5_out.rmsd.xvg')[:,1]


a = 0.8
plt.plot(G_out, color = '#5912A3', alpha = a, label = 'Outward open')
plt.plot(G_out_occ, color = '#0059E6', alpha = a, label = 'Outward occluded')
plt.plot(G_occ, color = '#0096F5', alpha = a, label = 'Fully occluded')
plt.plot(G_in_occ, color = '#00A396', alpha = a, label = 'Inward occluded')
plt.plot(G_in, color = '#007359', alpha = a, label = 'Inward open')

plt.xlabel('Time (ns)')
plt.ylabel('RMSD (nm)')
plt.xlim(0)
plt.ylim(0)
plt.legend()
plt.savefig('../images_figures/RMSD_all.png', dpi = 1000)

