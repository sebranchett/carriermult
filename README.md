# Charge Carrier Multiplication

[Sven Weerdenburg](https://orcid.org/0009-0001-0623-5551) created the initial contents of this repository as part of his Master's thesis entitled "A Theoretical Study on Charge Carrier Multiplication in MoTe2".

Please see Sven's excellent [documentation](./cmscript_manual.md).

If you use this code in your research, please cite the following publication:
[New Theoretical Model to Describe Carrier Multiplication in Semiconductors: Explanation of Disparate Efficiency in MoTe<sub>2</sub> versus PbS and PbSe.](https://doi.org/10.1021/acs.jpcc.4c00383)

Sven's work has been extended here to read input from Quantum ESPRESSO (as an alternative to ABINIT) and an example of using the code on [DelftBlue](https://www.tudelft.nl/dhpc/system) is provided.
If you use DelftBlue, please cite as described [here](https://doc.dhpc.tudelft.nl/delftblue/Citing-DHPC/).

## Usage
A complete analysis consists of three steps:
1. Perform a DFT calculation to obtain the band structure and wavefunctions
2. Calculate the number of carrier multiplication transitions using the `cmscript.py` script
3. Analyze the results using the `analysis.ipynb` Jupyter notebook

Depending on the size of your system, steps 1 and 2 can take a significant amount of time and are best performed on a computing cluster.
The files `Ncm.csv` and `reduced_energies.csv` can then be transferred to your local machine for analysis in step 3.

## Installation
The Python code has a number of dependencies, which can be installed with [Conda](https://docs.conda.io/en/latest/).
```bash
conda env create -f environment.yml
```
This will create a conda environment named `carriermult`.

You will need to activate the conda environment each time you start a new terminal session, and in your batch script.
```bash
conda activate carriermult
```

Depending on your setup, you may need to do the following to use the `carriermult` environment in Jupyter notebooks:
```bash
python -m ipykernel install --name carriermult --display-name "carriermult"
```

## Warning
This code is provided as-is, without any guarantees or warranty.
So far, it has been used on the following systems:
- MoTe2 (as described in the publication and using Quantum ESPRESSO on DelftBlue)
- PbS (as described in the publication)
- PbSe (as described in the publication)
