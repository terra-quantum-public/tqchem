# Installation Guide

## **Table of Contents**

1. [Prerequisites for Installation](#prerequisites-for-installation)
2. [Installation Steps](#installation-steps)
3. [Important Notes](#important-notes)
4. [Support](#support)

## **Prerequisites for Installation**

1. **Python Version**:
   Ensure you have Python **>=3.10** installed.

2. **Conda Environment Manager**:
   You must have one of the following environment managers installed:
    - [Conda](https://docs.conda.io/en/latest/miniconda.html)
    - [Anaconda](https://www.anaconda.com/products/distribution#download-section)
    - [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
    - [Mamba](https://github.com/mamba-org/mamba)

   If you do not have any of the above installed, download and set up your preferred environment manager before proceeding.

3. **License Requirement**:
    - If you have an interest in using TQChem, go to the following link to request access:
      [https://terraquantum.swiss/tqchem-request-access](https://terraquantum.swiss/tqchem-request-access).
    - Usage is limited to **non-commercial academic purposes**, as specified in the **TQChem Academic Use License Agreement**.

4. Once obtained, set the license key as an environment variable:
   ```bash
   export TQCHEM_LICENSE_KEY='<YOUR_LICENSE_KEY>'
   ```

## **Installation Steps**

Run the following command to install TQChem using **conda**:

```bash
conda install tqchem -c conda-forge -c https://repo.prefix.dev/terraquantumag
```

Ensure you have Conda or one of its variants (e.g., Anaconda, Miniconda, or Mamba) installed and set up properly in your system to execute this command.

## **Important Notes**

- **Acceptance of EULA:**
  Upon the first import of the TQChem library in your project, if the EULA (End User License Agreement) has not yet been accepted, you will be presented with the following prompt:

   ```
   Validating license for machine id: <your-machine-id>
   Contract acceptance is not registered, please proceed with contract acceptance
   Do you accept the license at /path/to/tqchem/LICENSE.txt? Y(es)/n(o):
   ```

    - To proceed, press **Enter** or type **"Y"** (Yes). This is required to continue using the library.
    - Upon success, the activation of your machine will occur automatically.

  When the machine activation is successful, you will see the following message:

   ```
   Validating license for machine id: <your-machine-id>
   Activated machine: <your-machine-id>
   License is valid, have fun using tqchem
   ```

  If you do not accept the license, you will be unable to use the library.
- **Environment Variable:**
  The `TQCHEM_LICENSE_KEY` environment variable is automatically verified during the package import. Make sure it is properly set according to the steps provided above.

- The academic license strictly limits usage to **non-commercial academic purposes**.

## **Support**

- For technical issues or license-related inquiries, contact:
  **[support@terraquantum.swiss](mailto:support@terraquantum.swiss)**
- **Machine Activation Limit:**
  If you reach the limit of machine activations allowed by your license, please contact our support team via the email provided above to resolve the issue.

Enjoy working with TQChem!
