Welcome to tqchem's documentation!
==================================
``tqchem`` is a python package developed by Terra Quantum AG for molecular modelling and simulation.
More specifically, it provides algorithms to find conformers of molecules, perform docking simulations,
predict spectra, and provides interfaces to various other packages.

The non-commercial version of ``tqchem`` comes with limited functionality.


Installation
============

You can install tqchem from conda via

.. code-block:: shell

    conda install tqchem -c conda-forge -c https://repo.prefix.dev/terraquantumag


You need a license key to perform simulations with tqchem.
You can request non-commercial access via `this request form <https://terraquantum.swiss/tqchem-request-access>`_.

Usage
=====

``tqchem`` can be used as an SDK from python via its application programming interface (API) or via its command line interface (CLI).
Usage examples for the API can be found in the `Examples <examples.html>`_ section and references in the `API references <autoapi/index.html>`_.
For examples of the CLI see below.

Command Line Interface (CLI)
============================

After installation, ``tqchem`` can be executed as a standalone application.
To get an overview of all the functionalities, simply type

.. code-block:: shell

   tqchem --help

To directly get started, you can generate a molecule file with instantly generated conformers using

.. code-block:: shell

   tqchem generate --smiles 'CCOC' --n_conformers 10

where ``--n_conformer`` determines the number of conformers.
The result will be a file called ``COOC.xyz`` with heuristic conformers.

To run ``ttconf`` on these starting conformers, run

.. code-block:: shell

   tqchem ttconf "CCOC.xyz" --method gfn2-xtb --n_sweeps 8 --rank 3

See the CLI help message (``tqchem ttconf --help``) for a list of all command options.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   examples
   autoapi/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
