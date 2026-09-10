# Installation Guide

## **Table of Contents**

1. [Prerequisites for Installation](#prerequisites-for-installation)
2. [Installation Steps](#installation-steps)
3. [Important Notes](#important-notes)
4. [Support](#support)

## **Prerequisites for Installation**

1. **Operating System**:
   Linux (x86-64) or macOS. **Windows is not supported**, because CREST, a required
   dependency, has no Windows build. Testing covers Linux x86-64 and macOS on Apple
   silicon; Intel macOS should work but is untested.

2. **Python Version**:
   Ensure you have Python **3.11 or 3.12** installed.

3. **Conda Environment Manager**:
   You must have one of the following environment managers installed:
    - [Conda](https://docs.conda.io/projects/conda/en/stable/)
    - [Anaconda](https://www.anaconda.com/download)
    - [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
    - [Mamba](https://github.com/mamba-org/mamba)

   If you do not have any of the above installed, download and set up your preferred environment manager before proceeding.

4. **License Requirement**:
    - If you have an interest in using TQChem, go to the following link to request access:
      [https://terraquantum.swiss/tqchem/enabling-technology/](https://terraquantum.swiss/tqchem/enabling-technology/).
    - Usage is limited to **non-commercial academic purposes**, as specified in the **TQChem Academic Use License Agreement**.

5. **License Key**:
   Once obtained, set the key as an environment variable:
   ```bash
   export TQCHEM_LICENSE_KEY='<YOUR_LICENSE_KEY>'
   ```

## **Installation Steps**

Create a dedicated environment and install TQChem into it:

```bash
conda create -n tqchem python=3.12
conda activate tqchem
conda install tqchem -c conda-forge -c https://repo.prefix.dev/terraquantumag --override-channels
```

Use a fresh environment rather than an existing one. TQChem pins several compiled
chemistry packages, and resolving them alongside an established environment often
downgrades them silently.

`--override-channels` keeps the `defaults` channel out of the solve. Mixing `defaults`
with `conda-forge` is the most common way this installation goes wrong; alternatively
set `channel_priority: strict` in your `.condarc`.

## **Important Notes**

- **Licensing on first import:**
  The `TQCHEM_LICENSE_KEY` environment variable is checked when the package is
  imported. There is no interactive prompt: the key is validated and, the first
  time you use TQChem on a given machine, that machine is registered
  automatically.

  On the first import from a new machine you will see:

   ```
   Validating license for machine id: <your-machine-id>
   Activated machine: <your-machine-id>
   License is valid, have fun using tqchem
   ```

  and on every import after that:

   ```
   Validating license for machine id: <your-machine-id>
   License is valid, have fun using tqchem
   ```

- **If the key is missing**, the import fails with
  `Please provide a valid License key with the env var TQCHEM_LICENSE_KEY`.

- **If you reach your machine activation limit**, the import fails with
  `Machine limit exceeded, please contact support@terraquantum.swiss`.

- Validation happens online, so the machine running TQChem needs network access.

- The academic license strictly limits usage to **non-commercial academic purposes**.

## **Support**

- For technical issues or license-related inquiries, contact:
  **[support@terraquantum.swiss](mailto:support@terraquantum.swiss)**

Enjoy working with TQChem!
