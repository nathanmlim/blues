from blues.analysis import msm, tools, cluster, population
import os, glob, traceback, pickle
from argparse import ArgumentParser
import pyemma.coordinates as coor
import numpy as np
import fnmatch
from blues import utils


def main(fpath, molid, args):
    #prefix = 't4-tol/t4-tol'
    prefix = 'test'
    #trajfiles = utils.get_data_filename('blues.analysis', 'tests/data/run03-centered.dcd')
    trajfiles = ['test.nc']
    pdbfile = 'test.pdb'
    #pdbfile = utils.get_data_filename('blues.analysis', 'tests/data/run03-centered.pdb')
    #Select features to analyze in trajectory
    feat = coor.featurizer(pdbfile)
    lig_atoms = feat.select("resname LIG and not type H")
    feat.add_selection(lig_atoms)
    #atom_index = np.array([
    #    2634, 2635, 2636, 2637, 2638, 2639, 2640, 1605, 1622, 1638, 1658, 1675, 1692, 1700, 1714, 1728, 1735, 1751,
    #    1761, 1768, 1788
    #])
    #paired_index = feat.pairs(atom_index)
    #feat.add_distances(paired_index)
    inp = coor.source(trajfiles, feat)

    #Initailize object to assign trajectories to clusters
    data = msm.ConstructMSM(inp)

    #Select the apprioriate lagtime from the implied timescale plots
    lag_list = np.arange(1, 16, 5)  #Give a range of lag times to try
    #data.plotImpliedTimescales(data.Y, 30, lag_list)

    #Frames in MD taken every 30ps
    #Discretize the trajectories
    #Estimate our Markov model from the discrete trajectories
    data.getMSM(data.Y, dt=30, lagtime=200)
    #with open('/home/nathanlim/SCRATCH/blues/blues/analysis/tests/data/t4-tol-msm.pkl', 'wb') as outf:
    #    pickle.dump(data,outf)

    #Analyze our assigned clusters by Silhouette method
    #with open('/home/nathanlim/SCRATCH/blues/blues/analysis/tests/data/t4-tol-msm.pkl', 'rb') as inf:
    #    data = pickle.load(inf)
    fbm = cluster.FindBindingModes(data)

    #Get the optimimal number of clusters by silhouette score
    n_clusters = fbm.getNumBindingModes(range_n_clusters=range(2, 5), outfname=prefix)
    #with open('/home/nathanlim/SCRATCH/blues/blues/analysis/tests/data/silhouette_pcca.pkl', 'wb') as outf1:
    #    pickle.dump(fbm.silhouette_pcca,outf1)
    #Draw samples from the PCCA metastable distributions
    pcca_outfiles = fbm.savePCCASamples(n_clusters, outfname=prefix)
    leaders, leader_labels = fbm.selectLeaders(pcca_outfiles, n_clusters, n_leaders_per_cluster=3, outfname=prefix)


parser = ArgumentParser()
parser.add_argument('-f', '--file_path', dest='fpath', type=str, help='parent directory of BLUES simluations')
parser.add_argument('-o', '--output', dest='outfname', type=str, default="blues", help='Filename for output DCD')
parser.add_argument('-m', '--molid', dest='molid', type=str, help='molecule ID')
parser.add_argument('--show_progress_bars', action='store_true', default=False, dest='show_progress_bars')
args = parser.parse_args()

try:
    main(args.fpath, args.molid, args)
except Exception as e:
    print('\nERROR!!!', args.molid)
    with open('%s.err' % args.molid, 'w') as errfile:
        errfile.write('### %s \n' % args.molid)
        errmsg = traceback.format_exc()
        errfile.write(errmsg)
        print('\n' + str(errmsg))
