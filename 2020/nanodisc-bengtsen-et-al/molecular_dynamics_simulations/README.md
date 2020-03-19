# Molecular dynamics simulation files

## Visualization
See simulation visualised by using the visualization program VMD (better than pymol for large simulations). Run:
``` 
vmd -e vmd_visualize_simulation.tcl
```

## Methods: 
    - ff: CHARMM36-m 
    - programme: GROMACS 5.0.7
    - ensemble: NVT
    - temp: 303.15 K
    - pressure: 1 bar
    - 2 simulations, concatenated for analysis
    - see mdout.mdp for more info


## Files: 
    - topology.pdb: 
            pdb structure of lipids and MSP proteins, excluding solvent. 
    - all_simulations.xtc
        concatenation of simulation_1 and simulation_2 into one long and processed to have continous timesteps (-settime in gmx trjcat). 
    - vmd_visualize_simulation.tcl
        file to help visualise the simulation in the visualization program VMD (better than pymol for large trajectories). Run by: vmd -e vmd_visualize_simulation.tcl
    - gromacs_run_files/simulation_1.xtc
        First simulation excluding solvent. Outputted every 1 ns, first 10ns discarded. Processed to avoid possible PBC jumps. 
    - gromacs_run_files/simulation_2.xtc
        Second simulation excluding solvent. Outputted every 1 ns,  first 10ns discarded. Processed to avoid possible PBC jumps.·
    - all_simulations.xtc
        concatenation of simulation_1 and simulation_2 into one long and processed to have continous timesteps (-settime in gmx trjcat). 
    - gromacs_run_files/input_run.tpr
        Gromacs production run input file. 
    - gromacs_run_files/index.ndx 
        Gromacs index file to handle only protein and lipids. 
    - gromacs_run_files/mdout.mdp
        All simulation settings 
    - gromacs_run_files/topologies/ 
        force-field and gromacs itp files
    - clustering_for_visualisation/
        QT clustering of concatenated simulations 
        
