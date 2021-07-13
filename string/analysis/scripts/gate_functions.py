
## get gate distances as an array
## gate_EC and gate_IC shape should be : gate_EC = [(TM1 start,TM1 end), (TM7 start, TM7 end)],
## and  gate_IC = [(TM4 start, TM4 end), (TM10 start, TM10 end)]

#quick refs: GLUT5:
# GLUT1: gate_EC = [(29,37), (288,295)], gate_IC = [(137,146), (385,394)]
# PfHT: gate_EC = [(43,51), (311,318)], gate_IC = [(145,154), (409,418)]


def make_gate_arr(md_uni, gate_EC, gate_IC):
    from MDAnalysis.analysis import distances

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
    return gate_EC_dists, gate_IC_dists










# The idea of this function is that you have several conditions that you would like to compare

## an example:
# make_gate_plot(EC_list = [PfHT_EC, apo_EC, GLUT1_EC],\
#               IC_list = [PfHT_IC, apo_IC, GLUT1_IC],\
#               label_list = ['PfHT', 'PfHT_apo', 'GLUT1'],\
#               color_list = ['red', 'green', 'blue'],\
#               figname='../image_figs/gate_dist.GLUT1_vs_PfHT')

def gate_dist_over_time(EC_list, IC_list, label_list, color_list, figsize = (13,6), ylim = (8,17), figname = None):
    markersize = 8
    fig, (ax1, ax2) = plt.subplots(nrows = 2, sharex = True, figsize = figsize)

    ## ec gate
    for n, EC in enumerate(EC_list):
        ax1.plot(EC, label = label_list[n], color = color_list[n])
    ax1.set_title("Extracellular gate")
    ax1.set_xlim(0)
    ax1.set_ylim(ylim)
    ax1.set_ylabel('Distance (A)')


    ## ic gate
    for n, IC in enumerate(IC_list):
        ax2.plot(IC, label = label_list[n], color = color_dir[label_list[n]])
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






## This function will plot one condition, as the classic KDE plot for GLUT5 string project format, eg. IC gate on x axis
## and EC gate on y axis. 

## an example:
#indir = '../../state_by_state_running/targeted_MD'
#plot_comparison(trajdir = f'{indir}/influx_apo_all_heavy', trajname = 'OutOpen-InOpen', label = 'influx_all',\
#               color = 'green', skip25 = True)

def plot_comparison(trajdir, trajname, color, a = 0.4, label = None, skip25 = None, highlight_frames = None):
    if skip25:
        traj = mda.Universe(f'{trajdir}/{trajname}/{sim_ref_dict[trajname][0]}', f'{trajdir}/{trajname}/{trajname}.skip25.xtc')
    else:
        traj = mda.Universe(f'{trajdir}/{trajname}/{sim_ref_dict[trajname][0]}', f'{trajdir}/{trajname}/{trajname}.xtc')
        
    gate_EC_dists = []
    gate_IC_dists = []
    for timestep in traj.trajectory:
        tm1,tm7,tm4,tm10 = get_tm_COM(traj)
        gate_EC_dists.append(float(distances.distance_array(tm1, tm7)) / 10)
        gate_IC_dists.append(float(distances.distance_array(tm4, tm10))/ 10) #keep in nm
    gate_EC_dists = np.array(gate_EC_dists)
    gate_IC_dists = np.array(gate_IC_dists)

    if not highlight_frames:
        plt.scatter(gate_IC_dists, gate_EC_dists, label=label, alpha = a, color = color)
    
    if highlight_frames:
        plt.scatter(gate_IC_dists[highlight_frames], gate_EC_dists[highlight_frames],\
                    label=label, alpha = 1, color = color, s=100, edgecolor = 'black')

    plt.xlim(0.9,1.82)
    plt.ylim(0.74, 1.7)
    plt.ylabel("Extracellular gate distance (nm)")
    plt.xlabel("Intracellular gate distance (nm)")
    plt.legend()
















