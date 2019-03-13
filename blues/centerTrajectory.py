import mdtraj as md


def centerTrajectory(traj, outfname):
    #Recenter/impose periodicity to the system
    anchor = traj.top.guess_anchor_molecules()[0]
    imgd = traj.image_molecules(anchor_molecules=[anchor])
    #Let's align by the protein atoms in reference to the first frame
    prt_idx = imgd.top.select('protein')
    superposed = imgd.superpose(reference=imgd, frame=0, atom_indices=prt_idx)
    #Save out the aligned frames
    superposed.save_netcdf(outfname + '-centered.nc')
    superposed[0].save_pdb(outfname + '-centered.pdb')


traj = md.load('test.nc', top='test.pdb')
centerTrajectory(traj, 't4-tol')
