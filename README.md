# q2 Python package

The `q2` code allows you to use the 2019 MOOG version (in its SILENT mode) to calculate elemental abundances of stars and/or determine their atmospheric parameters using the standard techniques of iron line excitation/ionization equilibrium. It also allows you to calculate other fundamental stellar parameters such as mass and age using isochrones. A tutorial is available <a href="https://github.com/astroChasqui/q2_tutorial">here</a>.

Installation
------------
Note that the `q2` package requires Python 3.7 or later.
The new version of `q2` can only be installed via pip, simply try:

```bash
pip install qoyllur-quipu
```

If you have installed the old version of `q2`, you must delete it from your HOME directory and also remove its PATH from bashrc (.bash_profile for Mac OS). Once you installed q2, .


Future updates will be performed via:

```bash
pip install qoyllur-quipu --upgrade
```

Quickstart
----------
Find spectroscopic parameters of a sample of stars using the Sun as the reference star (strict line-by-line differential analysis):

```python
import q2
data = q2.Data('stars.csv', 'lines.csv')
sp = q2.specpars.SolvePars(grid='marcs')
q2.specpars.solve_all(data, sp, 'solution.csv', 'Sun')
```

Measure elemental abundances of a sample of stars with respect to the solar abundances (line-by-line):

```python
species_ids = ['CI', 'OI', 'BaII']
q2.abundances.get_all(data, species_ids, 'abundances.csv', 'Sun')
```

MOOG
----

q2 needs a recent version of <a href="http://www.as.utexas.edu/~chris/moog.html">MOOGSILENT</a> (>2003) installed in your system. It has been fully tested only with the FEB 2017 version.

Data
----

You wont be able to use all the package features without the associated data, particularly the model atmospheres and isochrone grids. Use <a href="http://www.astrochasqui.com/share/q2Data.zip">this link</a> to download the data and make sure it is extracted into the q2/Data folder.

Preferred citation
------------------

<a href="https://doi.org/10.1051/0004-6361/201424244">Ramirez et al. 2014, A&A, 572, A48</a>

Contact info
------------

Ivan Ramirez (iramirez@tacomacc.edu)
