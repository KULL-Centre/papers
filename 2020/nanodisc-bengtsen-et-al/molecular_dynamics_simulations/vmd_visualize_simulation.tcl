
#===============================================================
# Tone Bengtsen
# Creation Date: 03-19-2020
# Purpose: script that
#        1) loads simulation
#        2) aligns all frames 
#        3) setup settings for displaying them different 
# run by: vmd -e vmd_visualize_simulation.tcl OR  vmd -e vmd_visualize_simulation.tcl -args topology.pdb all_simulations.xtc
#===============================================================



#=============================================
# load molecule files 
#=============================================


if { [llength $argv] == 2} {
    set topol_file [lindex $argv 0]
    set traj_file [lindex $argv 1]
    mol load pdb $topol_file
    mol addfile $traj_file waitfor all


} else {
    mol load pdb topology.pdb
    mol addfile all_simulations.xtc waitfor all 
}



#=============================================
# remove rotation and translation by alignment 
#=============================================

# calculate transformation matrix
set frame0 [atomselect 0 "name CA and not resid 234 to 243 and not resid 55 to 65" frame 0]

set nr_frames [molinfo 0 get numframes]
for {set i 0} {$i < $nr_frames} {incr i} {
    set traj [atomselect 0 "name CA and not resid 234 to 243 and not resid 55 to 65" frame $i]
    set traj_align [measure fit $traj $frame0]
    set whole_traj [atomselect 0 all frame $i]
    $whole_traj move $traj_align
}

#=============================================
# set selections 
#=============================================

# nanodisc helices  
atomselect macro Helix1 {resid 44 to 65}
atomselect macro Helix2 {resid 66 to 87}
atomselect macro Helix3 {resid 88 to 99}
atomselect macro Helix4 {resid 100 to 120}
atomselect macro Helix6 {resid 143 to 164}
atomselect macro Helix7 {resid 165 to 187}
atomselect macro Helix8 {resid 188 to 208}
atomselect macro Helix9 {resid 209 to 219}
atomselect macro Helix10 {resid 220 to 243}
atomselect macro lip_tails {lipid and noh and  (not name C13 C14 C15 O14 O13 O22 O32)}
atomselect macro lip_heads {lipid and  name N}


#=============================================
# change representation structures
#=============================================

# background color 
color Display Background white 

#change repr for simulation of protein only
mol delrep 0 0
mol representation NewCartoon
mol color ColorID 0 
mol selection {protein}
mol addrep 0

#change repr for simulation of  lipids
mol representation DynamicBonds  2.300000 0.300000 12.000000
mol color ColorID 15 
mol selection {lip_tails}
mol addrep 0

#change repr for simulation of  lipids nitrogens
mol representation VDW
mol color ColorID 0
mol selection {lip_heads}
mol addrep 0


#=============================================
# set selections 
#=============================================

mol color ColorID 3
mol selection {Helix1}
mol representation NewCartoon
mol addrep 0
#mol showrep 0 3 off


mol color ColorID 4
mol selection {Helix2}
mol representation NewCartoon
mol addrep 0
#mol showrep 0 1 off

mol color ColorID 5
mol selection {Helix3}
mol representation NewCartoon
mol addrep 0
#mol showrep 0 1 off


mol color ColorID 6
mol selection {Helix4}/
mol representation NewCartoon
mol addrep 0
#mol showrep 0 1 off

mol color ColorID 7
mol selection {Helix6}
mol representation NewCartoon
mol addrep 0
#mol showrep 0 1 off

mol color ColorID 8
mol selection {Helix7}
mol representation NewCartoon
mol addrep 0
#mol showrep 0 1 off

mol color ColorID 9
mol selection {Helix8}
mol representation NewCartoon
mol addrep 0
#mol showrep 0 1 off

mol color ColorID 10
mol selection {Helix9}
mol representation NewCartoon
mol addrep 0
#mol showrep 0 1 off

mol color ColorID 11
mol selection {Helix10}
mol representation NewCartoon
mol addrep 0
#mol showrep 0 1 off




