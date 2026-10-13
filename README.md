# Charge Carrier Multiplication

[Sven Weerdenburg](https://orcid.org/0009-0001-0623-5551) created the initial contents of this repository as part of his Master's thesis entitled "A Theoretical Study on Charge Carrier Multiplication in MoTe2".

Please see Sven's excellent [documentation](./cmscript_manual.md).

[Aniol Anglada Trias](https://es.linkedin.com/in/aniol-anglada-trias-04b9983b2) extended the functionality to allow for calculations on crystal supercells for the study of defects, taking into account intra-bandgap defect states. Aniol did this as part of his B.Sc. thesis project.

If you use this code in your research, please cite the following publication:
[New Theoretical Model to Describe Carrier Multiplication in Semiconductors: Explanation of Disparate Efficiency in MoTe<sub>2</sub> versus PbS and PbSe.](https://doi.org/10.1021/acs.jpcc.4c00383)

Sven's work has been extended here to read DFT input from [Quantum ESPRESSO](https://www.quantum-espresso.org/), via [Yambo](https://www.yambo-code.eu/), as an alternative to ABINIT. You can also read in Yambo GW data. An example of using the code on [DelftBlue](https://www.tudelft.nl/dhpc/system) is provided.
If you use DelftBlue, please cite as described [here](https://doc.dhpc.tudelft.nl/delftblue/Citing-DHPC/).

## Usage
A complete analysis consists of three steps:
1. Perform a DFT calculation to obtain the band structure in the full Brillouin zone. Optionally, perform a GW calculation
2. Calculate the number of carrier multiplication transitions using the `cmscript.py` script
3. Analyze the results starting from the `analysis.ipynb` Jupyter notebook

In step 2, you will need to provide the input parameters in the file `cm_input.yaml`.

In step 3, you will need to pay close attention to paths of input files and graph titles and labels. You will also need to provide experimental data for comparison.

Depending on the size of your system, steps 1 and 2 can take a significant amount of time and are best performed on a computing cluster.
The files `Ncm.csv` and `reduced_energies.csv` can then be transferred to your local machine for analysis in step 3.

## Installation
The Python code has a number of dependencies, which can be installed with [uv](https://docs.astral.sh/uv/). Make sure you have uv on both your laptop and on DelftBlue, see the [uv installation instructions](https://docs.astral.sh/uv/getting-started/installation/) and the [DelftBlue Documentation](https://doc.dhpc.tudelft.nl/delftblue/Python/#uv).

Clone this repository and move into the cloned directory. Then create a virtual environment with the Python dependencies defined in `pyproject.toml`.

To create the uv virtual environment on your laptop:
```bash
uv sync
```

To create the uv virtual environment on DelftBlue type:
```bash
module load 2026 cpu python
uv venv --python "$(which python)"
uv sync
```

This sets up the virtual environment. You will need to **activate the virtual environment every time** you start a new terminal session, and in your batch script.

To activate the virtual environment on a Windows Git Bash system, use the following command:

```bash
source .venv/Scripts/activate
```
On macOS or Linux, use the following command:

```bash
source .venv/bin/activate
```

In your SLURM script on DelftBlue:
```
source /path/to/project/.venv/bin/activate
```

Depending on your setup, you may need to do the following to use the virtual environment in Jupyter notebooks on your laptop:
```bash
uv run ipython kernel install --user --name=carriermult --display-name "carriermult"
```

## Warning
This code is provided as-is, without any guarantees or warranty.
So far, it has been used on the following systems:
- MoTe2 (as described in the publication and using Quantum ESPRESSO on DelftBlue)
- PbS (as described in the publication)
- PbSe (as described in the publication)
- WSe2

> [!WARNING]
> Since the introduction of the DelftBlue 2026 stack, the Abinit example is no longer working.
