==============================================================================
This is a README file for jackknife.py

Written By: Amanda Buyan

Please cite pubs.acs.org/doi/abs/10.1021/ct500003g
==============================================================================

The jackknife script is for use on GROMACS trajectories where there are two TM
helices that dimerize.  The general idea is to figure out when two properties
(the crossing angle between the two helices, as well as the residues
contacting each other) converge.  

The principle comes from the jackknife method (Quenouille, 1949, extended by
Tukey in 1959).  The idea is to compute the estimate of the statistic of
interest (in this case, crossing angle or contacting residues) by leaving one
or more of the samples out of the calculation, to see how the statistic would
change.  This was further extended by Wilke, 2012, who extended the idea by
iterating over all possible subset sizes to determine the minimum sample size 
needed to detect patches in fMRI scans of patients.  I have adapted this for 
TM helices.  This script works for both homo- and heterodimers.

The analysis consists of the following:

1. Boxplots (min, first quartile, median, third quartile, max) of the
distribution of the crossing angle values for both right-handed and
left-handed helices 

2. Frequency plots of the contacts.  This is done by taking all of the frames
of all of the simulations in one subset, averaging the distance between
residues for the contact matrix.  Then, the closest N contacts are chosen
(these are the contacts with the lowest average distances), given a 1 if it
met this criteria, and a 0 if it did not meet this criteria.  This is done for
all subsets in the simulation.  The boolean matrices are then added together,
and the frequencies of each residue-residue contact that happens are apparent. 

==============================================================================
Installation Prerequisites
==============================================================================

Before running the jackknife script, there are a few things that need to be
installed:

1. MDAnalysis (includes NumPy and SciPy)
2. Subsets
3. ctypes
4. tables
5. gnuplot

To install each of these individually, use the following instructions:

1. MDAnalysis (includes NumPy and SciPy)

    Website for installation (depending on your version of Ubuntu):

        https://code.google.com/p/mdanalysis/wiki/InstallRecipes

    This website should contain all of the instructures for installation 

2. subsets

    You have to download this from this website:

        https://github.com/EricBurnett/EnumeratedSubsets/tree/master/python

    Keep this script in the same directory as jackknife.py

3. ctypes

    For Linux:

        sudo apt-get install python-ctypes

4. tables

    For Linux:

        sudo apt-get install libhdf5-serial-dev
        sudo apt-get install python-tables

5. gnuplot

    For Linux:

        sudo apt-get install gnuplot

==============================================================================
Preprocessing of files
==============================================================================

The files should be preprocess so they contain only the two TM helices, 
and that the trajectory is centred on one of them.  To do this, run the 
following commands in GROMACS for each directory (change residue numbers 
for your particular system):

echo -e 'r 1-27\nq\n' | make_ndx -f md.tpr -o TM1.ndx

echo -e '16\n1\n' | trjconv -f md.xtc -s md.tpr -n TM1.ndx -pbc mol -center -o
md_wrap.xtc

echo -e '1\n' | editconf -f md.gro -n TM1.ndx -o protein.gro

To run them for multiple directories, run the following command:

for i in `seq 1 50`; do cd $i;echo -e 'r 1-27\nq\n' | make_ndx -f md.tpr -o
TM1.ndx;echo -e '16\n1\n' | trjconv -f md.xtc -s md.tpr -n TM1.ndx -pbc mol
-center -o md_wrap.xtc;echo -e '1\n' | editconf -f md.gro -n TM1.ndx -o
protein.gro;cd ../;done

==============================================================================
Example Usage
==============================================================================

To get the help screen for jackknife.py, type

python jackknife.py -h

and all of the options will be displayed, as well as the defaults.

***All input files MUST have the same name***

Example of usage:

I have two transmembrane helices, with 33 residues per helix (and 78
particles in each helix, making the selection for the first helix 1-77 and the
second helix 78-154), run in MARTINI2.2.  The number of directories I have is
50, timestep is 200ps, total time of the simulation is 1 microsecond, and the
bounds of the crossing angle is -60 and 60.  I want the crossing angle
binsize to be 1 degree, number of subsets to be 1000, and the number of closely 
contacting residues to be 10. Lastly, the name of my gro files are protein.gro 
and my xtc files are md_wrap.xtc (like above).  The command you would run is:

python jackknife.py -g protein.gro -x md_wrap.xtc -res 1 77 78 154 -numdirs 50
-dstart 1 -dend 49 -ts 200 -totaltime 1000000 -bounds -60 60 -binsize 1
-numsubsets 1000 -topres 10

-dstart is the starting size of the subsets (1)
-dend is the ending size of the subsets (49 - N-1, where N is the total number
of simulations you have)

==============================================================================
Notes
==============================================================================

-All options should be explained in the help screen
-Other suggestions are welcome, contributions from these suggestions are
encouraged
