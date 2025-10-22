**CM script Manual**

**By Sven Weerdenburg, 29 April 2022**

**Install instructions:**

I have used Python 3.7, I think you will need at least Python 3.6 in
order to run the code. (The code uses f-Strings which were introduced in
Python 3.6) Since we need some modules for the code to run, I installed
Python 3.7 locally. I am not sure if it is really needed, but that’s how
I got it to work at least.

**Installing Python & modules locally on Snellius**

I do not fully remember how I installed Python locally, but I know I did
it with Anaconda.

Load modules:

*module load 2021*

*module load Anaconda3/2021.05*

conda install -c anaconda python=3.7
--prefix=’/home/\$USER/.local/python’

Python3.7 should then be installed in your .local directory.

Once Python is installed, we need to install 3 modules:

- **Abipy:** <https://abinit.github.io/abipy/installation.html>

- **Joblib** <https://joblib.readthedocs.io/en/latest/installing.html>

- **Tqdm:** <https://tqdm.github.io/>

**Install Abipy (version: 0.9.1)**

Abipy is absolutely needed for the code to run, because this module can
extract information from the \*GSR.nc file.

To install issue the following commands:

*conda config --add channels conda-forge*

conda install abipy --prefix=”/home/\$USER/.local/python”

That should be all. The installer might ask you to update & install some
other modules too. Just answer yes when prompted.

**Install Joblib (version: 1.0.1)**

This module makes parallel computations possible. It should be possible
to install it with pip (see joblib install guide), but I got it working
by downloading the tarball from: <https://pypi.org/project/joblib/>

Once the tarball is unpacked in Snellius, enter the joblib directory
that was made and issue:

*python setup.py install --prefix*=”/home/\$USER/.local/python”

Again, if prompted to update other modules, just answer yes.

**  
**

**Install tqdm (version: 4.59.0)**

The tqdm module adds a progress bar to the slurm file while the
calculations are running. I highly recommend using tqdm, but if there
are problems installing tqdm you can use the code without by deleting
line 21 and 384, and by removing the word tqdm from the cmscript.py
file.

You can install tqdm with:

*conda install -c conda-forge tqdm
--prefix=*”/home/\$USER/.local/python”

**  
**

**Using the cmscript:**

Now that all modules are installed, the cmscript can be run. This part
of the manual will explain the input that is required and what will be
output.

**Abinit DFT calculation**

The cmscript.py script does carrier multiplication calculations from the
output files of Abinit. It is therefore important that the proper input
variables are set for the Abinit DFT calculations. This sub-section will
explain what needs to be defined in the Abinit input file, so that the
output files can be analyzed with the cmscript.

First, the cmscript requires a full Brillouin Zone calculation. This
means that the k-point grid that will be generated is spanned over the
entire Brillouin Zone and not just the irreducible part. To do this make
sure you set the input variable:

kptopt 3

cmscript.py requires the energies in eV. The script should be able to
automatically convert energies from Hartree to eV, but this has not been
tested. So just to be sure set the following input variable:

enunit 1

so that all energies will be given in eV.

Another thing that is important to note, is that the results of the
cmscript depend heavily on the number of k-points that are included in
the DFT calculations. A 4x4x4 grid will have converged much less than a
16x16x16 grid. It is therefore good practice to do DFT calculations for
multiple k-point grids, and then analyze the results for all the grids
with cmscript so that the results can be compared.

**Input**

In order to run cmscript.py you will first need to define the input.
This is done in the cmscript.py file itself, in lines 30 - 62. This
subsection will describe what all input variables mean.

**run_cm_calculations:** If this is set to True, the cmscript will
calculate the number of carrier multiplication transitions. If set to
False, the script will not do the calculations but will load the files
and variables into memory. Only set this to False for testing and
modifying the script.

**abinit_file:** Here you need to define the absolute path to a
\*\_GSR.nc file. This is an output file of Abinit and contains all the
information that the cmscript needs.

**TrueBG:** Here you define the experimental bandgap. This will be used
to calculate the required offset for the scissor operator. Note that
this is given in eV.

**Emax:** The maximum photon energy for which you want calculate the
N<sub>CM</sub>.

This will be used to reduce the number of bands that are included in the
calculations. If Emax = 4, then the code will not include bands whose
minimal energy is less than 4 eV. Based on TrueBG, the minimum energy
will be calculated with:  
Emin = Emax - TrueBG  
Any valence band whose maximum energy is below Emin, will not be
included in the calculations.

$$|E_{1i} - E_{1f} + E_{2i} - E_{2f}| < etol$$

**logfilename:** Name of the log file. This file will contain
information about the input and progress.

**save_kfile** True if you want to save the integer numbers representing
the k-points of all transitions set this to True. This will require much
more computation time and memory storage. It is recommended to leave
this at False.

**save_CMfile** Saves the k-point and energy of the initial and final
state of the electron. Also by default set to False for same reasons as
save_kfile.

If you set save_kfile or save_CMfile to True, you need to define a
filename for **kfilename** and **CMfilename**. Note: The band index
saved in CMfilename is not the real band index. This is because the
script uses a reduced form of the eigenvalue matrix (red_energies
variable). In order to convert the reduced band index to the real band
index, use the red2e function. For example red2e(bandindex), the output
is the real bandindex.

**remove_degenerate**: This variable is set to True by default. It
removes degenerate energy bands from calculations, which saves
computation time. At the end the N<sub>CM</sub> is corrected with a
factor to account for degenerate bands.

**Running the cmscript**

After defining the input, you can run the script. In this package an
example is given in the directory \example_cmscript\\. I will explain
how to run cmscript.py at the hand of this example.

Running csmcript requires three files:

- cmscript.py

- runjob.sh

- The Abinit \*\_GSR.nc output file.

All the required files are included in the example directory. To run the
example, you will need to change the **abinit_file** variable to the
absolute path of the \*\_GSR.nc file that is also included in the
example directory. Furhtermore, you will need to change the path to
Python 3.7. This means changing the line:

/home/\$USER/.local/bin/python-modules/bin/python3.7 cmscript.py

Optionally, you can also change the number of processors in line

\#SBATCH -c 16

For 500 k-points and more I recommend using 128 processors.

Once all is done, issue sbatch runjob.sh to run the cmscript.

Once the job is submitted with sbatch runjob.sh a log.txt file should be
generated soon. For the CM transitions calculations you can track
progress in the slurm file. The tqdm module will print a progress bar.

The duration of the calculations depend heavily on the number of
k-points, Emax and etol. A 16x16x8 k-point grid for MoTe<sub>2</sub>
took me \~5 days and 12x12x6 took \~ 8 hours on 128 processors.

**Output files**

**Ncm.csv** This contains the number of CM transitions for each kpoint
and reduced band index. The row represents a kpoint and column the
reduced band index. Again, to convert reduced band index to real band
index use the red2e(bandindex) function. To obtain the energy of a
specific row & column, simply use red_energies\[row, column\].

reduced_energies.csv file that contains the energies that are used. Each
row represents a kpoint. The column represents reduced band index

log.txt The log file that contains information about the calculations.

(karray_4x4x4.csv) see kfilename in this manual

(CM_initial_states_4x4x4_0.85eV.csv) see CMfilename in this manual

**Analyzing the Ncm.csv file**

After the N<sub>CM</sub> calculations, the Ncm.csv file can be analyzed.
In /example_analysis/ a Jupyter Notebook file is included, which was
used to make several figures.

To run the example analysis, open up Jupyter Notebook and open the
analysis.ipynb file in the /example_analysis/ directory on your
computer. This example uses the files in the
/example_analysis/12x12x6_0.90eV_30Ha/ directory. These files were
generated by running the cmscript on the output of an Abinit DFT
calculation on a Full Brillouin zone 12x12x6 grid, with a scissor
operator to match an experimental indirect gap of 0.90 eV, and a 30
Hartree cut-off energy. Running the three cells should result in 4
figures.

The first is the density of transitions. This is similar to the density
of states, only gives the y-axis the density of CM transitions for
specific photon energies. Better quantitative results are shown in the
second figure, which is the average number of carrier multiplication
transitions, normalized with respect to the number of k-points, plotted
as function of photon energy. The third figure is the quantum yield as
function of photon energy. For this plot the experimental data of
MoTe<sub>2</sub> was required, which is included in the
/example_analysis/Literature Data/ directory. The fourth figure is the
R<sub>CM</sub> / R<sub>cool</sub> rate ratio, which was obtained from
fitting the calculated data to the experimental data.

Note that the third and fourth figure can only be produced if
experimental data is provided to the code. If not, then only the first
and second figure can be generated.

To understand what the analysis script does, have a look at the
comments.
