Scripts and input files for:

Yong Wang, Elena Papaleo, Kresten Lindorff-Larsen

Mapping transiently formed and sparsely populated conformations on a complex energy landscape

https://doi.org/10.7554/eLife.17505

-----
Path_Moil_workflow.sh is the main script to generate an initial path connecting the major and the minor state of T4L L99A based on the SDP (steepest descent path) module in Moil11 software. This script depends on a few other scripts, such as:
Moil2Charmm.pl used for converting the format of outputed path from Moil to MD force field (here CHARMM),
pth2pdb.sh for converting Moil .pth to pdb structures,
pdb2multimodel.sh for separating the path into single pdb structures per frame.

plumed-PathCV-WTMetaD.dat shows the plumed settings used for path CV driven well-temperated metadynamics with adaptive width gaussians.
