# Example notebooks

These notebooks show how to use `tqchem` from Python. A good entry point is
`conformer_search.ipynb`, which walks through a conformer search with the TTConf
algorithm; `conformer_coordinate_screening.ipynb` covers the interactive coordinate
viewer.

## What you already have

`numpy`, `plotly`, `tblite` and `py3Dmol` are installed with `tqchem`, so nothing in the
notebooks needs them added separately.

## Running the notebooks

You need a notebook interface, which `tqchem` does not install:

```bash
conda install jupyterlab   # into the environment tqchem is installed in
```

```bash
pixi add jupyterlab
```

Then start it in this directory:

```bash
jupyter lab
```

The notebooks are also stored with their outputs, so they can be read without running
anything.

## Extra tools for `create_vmd_movie.ipynb`

That one notebook renders a trajectory into an animated GIF and calls two external
programs to do it. The others need nothing beyond the above.

**ImageMagick** provides the `magick` command:

```bash
conda install imagemagick
```

```bash
pixi add imagemagick
```

**VMD** has to be installed separately. conda-forge carries it only for Linux x86-64, at
version 1.9.3:

```bash
conda install vmd     # Linux only
```

On macOS, download it from
[the VMD site](https://www.ks.uiuc.edu/Research/vmd/) and edit the path at the top of the
notebook's last cell to match your installation — it is currently set to a macOS Apple
silicon install and will not match yours.

Note that the notebook calls both programs through Jupyter's `!` shell syntax, which does
not report a failure. If either is missing, the notebook runs to the end and simply writes
no movie. If you get no `movie.gif`, check that `magick` and VMD are both callable before
looking anywhere else.
