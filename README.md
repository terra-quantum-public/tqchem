<div align="center">
  <h1>tqchem</h1>
  <a href="https://github.com/terra-quantum-public/tqchem">Repository</a>
  &middot; 
  <a href="https://prefix.dev/channels/terraquantumag/packages/tqchem">tqchem on Prefix.dev</a>
</div>

## Table of Contents

1. [Introduction to tqchem](#introduction-to-tqchem)
2. [Getting Started](#getting-started)
3. [Installation](#installation)
4. [Running an Exemplary Experiment](#running-an-exemplary-experiment)
5. [Documentation](#documentation)

# Introduction to tqchem

The `tqchem` library provides sophisticated tools for the manipulation and optimization of molecular systems using internal BAT coordinates. With this tool, managing conformer optimization and molecular structure data conversion is made efficient and straightforward.

Utilize the following functionalities after installation:

- **TTconf**: High accuracy conformer optimization.
- **Generate**: Quickly convert SMILES strings to .xyz files and rapidly generate conformers using RDKit.

No more need to manually manage molecular setups—let `tqchem` handle it. From CLI or any preferred environment:

- Optimize molecular structures.
- Convert molecular data formats.
- Manage and visualize molecular data efficiently.

# Getting Started

## System Requirements

To run `tqchem`, ensure your system meets the following requirements:

- Modern operating system (macOS 12.0+, Windows 10+, or Linux)
- Python 3.9 or above
- Conda for environment management

## Quick Start

To quickly set up your environment and start using `tqchem`, follow these steps:

```bash
conda create -n my_env_name python=3.9  # create the environment
conda activate my_env_name             # activate your environment
conda install tqchem -c conda-forge -c https://repo.prefix.dev/terraquantumag  # install tqchem
tqchem ttconf -h                       # for ttconf help and options
```

# Installation

To install `tqchem` using `conda`, the following command will set up all necessary dependencies:

```bash
conda install tqchem -c conda-forge -c https://repo.prefix.dev/terraquantumag
```

# Running an Exemplary Experiment

After setting up `tqchem`, running computations is straightforward. For example, to run a conformer optimization with a SMILES string:

```sh
tqchem ttconf COCCNC
```

# Documentation

For comprehensive documentation and further resources, visit our [Documentation](https://tqchem-docs.terraquantum.io).

