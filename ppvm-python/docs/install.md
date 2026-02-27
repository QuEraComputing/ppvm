# Installation

## Requirements

- Python 3.10 or later

## Install from source

Clone the repository and install with [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/QuEraComputing/ppvm
cd ppvm
uv sync --project ppvm-python
```

This compiles the Rust core and installs the package into a local virtual
environment. The first build takes a minute; subsequent builds are cached.

## Using the package

Activate the environment or prefix commands with `uv run`:

```bash
uv run --project ppvm-python python -c "from ppvm import PauliSum; print(PauliSum(initial_terms=['ZZ']))"
```
