from blues.analysis import msm, tools, cluster, population, plots
import os, glob, traceback, pickle
from argparse import ArgumentParser
import pyemma.coordinates as coor
from pyemma import config

import numpy as np


def main(fpath, molid, args):
    prefix = 't4-tol'
    #prefix = os.path.join(fpath, molid, molid)
    #trajfiles = ['test.nc']
    #topfiles = ['test.pdb']
    trajfiles = glob.glob(prefix + '*-centered.nc')
    topfiles = glob.glob(prefix + '*-centered.pdb')
    #jsonfile = os.path.join(prefix+'-acc_its.json')

    #Pre-process trajectory files, check for unbound ligands
    trajfiles = tools.check_bound_ligand(trajfiles, topfiles[0])

    #Select features to analyze in trajectory
    feat = coor.featurizer(topfiles[0])
    lig_atoms = feat.select("resname LIG and not type H")
    feat.add_selection(lig_atoms)
    inp = coor.source(trajfiles, feat)

    #Define some input parameters
    dt = 8  #Frames in MD taken every 8ps == 0.008 ns
    lag_list = np.arange(1, 40, 5)  #Give a range of lag times to try

    #Initailize object to assign trajectories to clusters
    data = msm.ConstructMSM(inp)

    #Select the apprioriate lagtime from the implied timescale plots
    #data.plotImpliedTimescales(data.Y, dt, lag_list, outfname=prefix)

    #Selecting a lagtime of 150ps, every BLUES iteration is 176ps
    lagtime = 150

    #Discretize the trajectories
    #Estimate our Markov model from the discrete trajectories
    data.getMSM(data.Y, dt, lagtime, fixed_seed=1)

    #Analyze our assigned clusters by Silhouette method
    fbm = cluster.FindBindingModes(data)

    #Get the optimimal number of clusters by silhouette score
    n_clusters = fbm.getNumBindingModes(range_n_clusters=range(2, 10), dim=2, outfname=prefix)

    #Assign trajectory frames to PCCA centers
    traj_cluster_labels = fbm.assignTrajFramesToPCCACenters(n_clusters)

    #Draw samples from the PCCA metastable distributions
    pcca_outfiles = fbm.savePCCASamples(n_clusters, outfname=prefix)

    #Every iteration stores 5 frames, every frame stores 0.0352ns
    bmp = population.BindingModePopulation(molid, n_clusters, time_per_frame=1, frames_per_iter=5)

    traj_populations = bmp.calcPopulation(traj_cluster_labels)
    bmp.barplotPopulation(n_clusters, traj_populations, bmp.orig_populations, outfname=prefix)

    for trj, lab in zip(trajfiles, traj_cluster_labels):
        bmp.plotTrajPopulation(trj, lab, outfname=prefix)


parser = ArgumentParser()
parser.add_argument('-f', '--file_path', dest='fpath', type=str, help='parent directory of BLUES simluations')
parser.add_argument('-o', '--output', dest='outfname', type=str, default="blues", help='Filename for output DCD')
parser.add_argument('-m', '--molid', dest='molid', type=str, help='molecule ID')
parser.add_argument('--show_progress_bars', action='store_true', default=False, dest='show_progress_bars')
args = parser.parse_args()

try:
    if args.show_progress_bars:
        config.show_progress_bars = True
    else:
        config.show_progress_bars = False

    main(args.fpath, args.molid, args)
except Exception as e:
    print('\nERROR!!!', args.molid)
    with open('%s.err' % args.molid, 'w') as errfile:
        errfile.write('### %s \n' % args.molid)
        errmsg = traceback.format_exc()
        errfile.write(errmsg)
        print('\n' + str(errmsg))
