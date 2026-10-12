# install uv once per computer (DelftBlue/laptop) in the home directory
curl -LsSf https://astral.sh/uv/install.sh | sh
# on DelftBlue, get the right version of python
module load 2026 cpu python
# on DelftBlue, tell uv which version of python to use
uv venv --python "$(which python)"
# per directory, tell uv to set up an environment with the packages in pyproject.toml
uv sync  # or: uv pip install -r pyproject.toml
# uv run python my_script.py  # to run commands in the environment without activating it
# source .venv/bin/activate  # in interactive work
# source /path/to/test_project/.venv/bin/activate  # in slurm job
