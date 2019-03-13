The way the system is set up is that the ethylene is restrained in the center of two particles. Since the ethylene hydrogens are charged, there is some orientational bias that causes populations of roughly 75:25 (at 200K).


eth.prmtop/eth.inpcrd  - ethylene with the charges on the hydrogen modifed
sqB.pdb - Just nitrogen atoms to act as particle centers (each with slightly different charges). Used in creating the system
tsystem.py - script used to create the system. If you just want to use the output of this, you can use the serialized system
test_sytem.xml - Serialized XML of the system
test_system.pdb - pdb of the stest system
output.dcd - example output of the test_system
aplot.py - plotting script to get the populations
