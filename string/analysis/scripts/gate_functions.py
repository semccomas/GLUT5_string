
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################

## get gate distances as an array
## gate_EC and gate_IC shape should be : gate_EC = [(TM1 start,TM1 end), (TM7 start, TM7 end)],
## and  gate_IC = [(TM4 start, TM4 end), (TM10 start, TM10 end)]

#quick refs: GLUT5: gate_EC = [(30,37), (289,295)], gate_IC = [(136,145), (386,394)]
# GLUT1: gate_EC = [(29,37), (288,295)], gate_IC = [(137,146), (385,394)]
# PfHT: gate_EC = [(43,51), (311,318)], gate_IC = [(145,154), (409,418)]


def make_gate_arr(md_uni, gate_EC, gate_IC):
    from MDAnalysis.analysis import distances
    import numpy as np

    gate_EC_dists = []
    gate_IC_dists = []
    
    for timestep in md_uni.trajectory:
        tm1 = md_uni.select_atoms('resid %i-%i' %(gate_EC[0][0], gate_EC[0][1])).center_of_mass()
        tm7 = md_uni.select_atoms('resid %i-%i' %(gate_EC[1][0], gate_EC[1][1])).center_of_mass()
        tm4 = md_uni.select_atoms('resid %i-%i' %(gate_IC[0][0], gate_IC[0][1])).center_of_mass()
        tm10 = md_uni.select_atoms('resid %i-%i' %(gate_IC[1][0], gate_IC[1][1])).center_of_mass()    


        gate_EC_dists.append(float(distances.distance_array(tm1, tm7)))
        gate_IC_dists.append(float(distances.distance_array(tm4, tm10)))
    print("returning EC gate, IC gate dists")
    gate_EC_dists = np.array(gate_EC_dists)
    gate_IC_dists = np.array(gate_IC_dists)
    return gate_EC_dists, gate_IC_dists








#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################


# The idea of this function is that you have several conditions that you would like to compare
## Just plot distance over time and have both graphs in one plot

## an example:
# make_gate_plot(EC_list = [PfHT_EC, apo_EC, GLUT1_EC],\
#               IC_list = [PfHT_IC, apo_IC, GLUT1_IC],\
#               label_list = ['PfHT', 'PfHT_apo', 'GLUT1'],\
#               color_list = ['red', 'green', 'blue'],\
#               figname='../image_figs/gate_dist.GLUT1_vs_PfHT')

def gate_dist_over_time(EC_list, IC_list, label_list, color_list, figsize = (13,6), ylim = (8,17), figname = None, trending_line = None):
    import matplotlib.pyplot as plt 
    import numpy as np

    markersize = 8
    fig, (ax1, ax2) = plt.subplots(nrows = 2, sharex = True, figsize = figsize)

    ## ec gate
    for n, EC in enumerate(EC_list):
        ax1.plot(EC, label = label_list[n], color = color_list[n])
        if trending_line:
            ax1.plot(np.arange(len(EC)), np.ones(len(EC))*np.mean(EC), color = 'black', alpha = 0.5, linewidth = 1)
    ax1.set_title("Extracellular gate")
    ax1.set_xlim(0)
    ax1.set_ylim(ylim)
    ax1.set_ylabel('Distance (A)')


    ## ic gate
    for n, IC in enumerate(IC_list):
        ax2.plot(IC, label = label_list[n], color = color_list[n])
        if trending_line:
            ax2.plot(np.arange(len(IC)), np.ones(len(IC))*np.mean(IC), color = 'black', alpha = 0.5, linewidth = 1)
    ax2.set_title("Intracellular gate")
    ax2.set_xlim(0) 
    ax2.set_ylim(ylim)
    ax2.set_xlabel('Time (ns)')
    ax2.set_ylabel('Distance (A)')

    
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])

    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])

    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25),
              fancybox=True, shadow=True, ncol=2)

    if figname:
         plt.savefig(f'{figname}.png', dpi = 400)
    else:
         plt.show()




#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################


## This function will load a trajectory, measure gate dist, and plot IC vs EC gate in scatter plot for one condition
## This is specific to GLUT5 at the moment!! (b/c gate dists)


## an example:
#indir = '../../state_by_state_running/targeted_MD'
#uni_to_gate_scatter(trajdir = indir, trajname = 'efflux_apo_gate_CV', ext = '10.string.pdb', color = 'blue')


def uni_to_gate_scatter(trajdir, trajname, ext, color, a = 0.4, scatter = True, label = None, uni_top = None, highlight_frames = None):
    import MDAnalysis as mda
    import matplotlib.pyplot as plt
    import numpy as np
    
    
    ## LOADING DATA
    if ext.endswith("gro") or ext.endswith("pdb"):
        traj = mda.Universe(f'{trajdir}/{trajname}.{ext}')
    if ext.endswith('xtc'):
        traj = mda.Universe(uni_top, f'{trajdir}/{trajname}.{ext}')

    ## do gate calcs
    EC, IC = make_gate_arr(traj, gate_EC = [(30,37), (289,295)], gate_IC = [(136,145), (386,394)])
    
    ## plot, don't show immediately so you can stack values
    plot_gate_scatter(EC, IC, label = label, color_list = color, show = False, scatter = scatter)






#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################

## This will plot the EC and IC gates as the usual scatter plot. You can either specify a list of colors, or a colormap
## which will automatically color each item in the list

# this input should be a numpy array, either 1d or 2d (for several conditions)
## eg:
#d = np.array((gate_EC_dists, gate_IC_dists))
#s = np.array((gate_IC_dists, gate_EC_dists))
#plot_gate_scatter(s, d, colormap = 'viridis', label = ['r','d']) 




def plot_gate_scatter(EC_arr, IC_arr, scatter = True, label = None, colormap = None, color_list = None, title = None, show = None, lims = None):
    import matplotlib.pyplot as plt
    import numpy as np
 
    if colormap:
        cmap_vals = plt.cm.get_cmap(colormap)
        color_vals = cmap_vals(np.linspace(0,1, np.shape(IC_arr)[0]))
        #color_vals = color_vals[::-1]
    else:
        color_vals = color_list

    ## need to have separate treatment for multiple gate dists. 
    ### Seems a bit unneccessary but for coloring this is useful
    if IC_arr.ndim > 1:
        for n in range(np.shape(IC_arr)[0]):
            plt.plot(IC_arr[n]/10.0, EC_arr[n]/10.0, marker = 'o', color = color_vals[n], label = label[n]) 

 
   ## now plot for only one pair of dists
    else:
            if colormap and scatter:
                plt.scatter(IC_arr/10.0, EC_arr/10.0, marker = 'o', c = color_vals, cmap = colormap, label = label)
            elif colormap and not scatter:
                print("can't plot plt.plot with colormap, set scatter = True!")
            elif not colormap and scatter:
                plt.scatter(IC_arr/10.0, EC_arr/10.0, marker = 'o', color = color_vals, label = label) 
            else:
               plt.plot(IC_arr/10.0, EC_arr/10.0, marker = 'o', color = color_vals, label = label)


    if not lims:
        plt.xlim(0.9, 1.82)
        plt.ylim(0.74, 1.7)

    plt.xlabel("Intracellular gate (nm)")
    plt.ylabel("Extracellular gate (nm)")
    plt.title(title)
    plt.legend()
    if show:
        plt.show()



#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################


## this will just plot classic RMSD's from MDA in a way that you can compare them

def compare_RMSD(traj_list, label_list = None, color_list = None):
    import MDAnalysis.analysis.rms
    import matplotlib.pyplot as plt
    
    rmsd_list = []
    for traj in traj_list:
        RMSD = MDAnalysis.analysis.rms.RMSD(traj, center = True)
        RMSD.run(0)
        rmsd = RMSD.rmsd.T
        rmsd_list.append(rmsd)

    fig = plt.figure(figsize = (15,5))

    for n, rmsd in enumerate(rmsd_list):
        time = rmsd[1]
        if label_list and color_list:
	        plt.plot(rmsd[2], label = label_list[n], color = color_list[n])
        elif label_list and not color_list:
                plt.plot(rmsd[2], label = label_list[n])
        elif not label_list and color_list:
                plt.plot(rmsd[2], color = color_list[n])
        else:
                plt.plot(rmsd[2])

    plt.legend()
    plt.xlim(0)
    plt.ylim(0)
    
    plt.xlabel('Iteration #')
    plt.ylabel('RMSD (A)')








