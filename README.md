q2 Python package
=================

q2 allows you to use MOOG (in its SILENT mode) to calculate elemental abundances of stars and/or determine their atmospheric parameters using the standard techniques of iron line excitation/ionization equlibrium. It also allows you to calculate other fundamental stellar parameters such as mass and age using isochrones. A tutorial is available <a href="https://github.com/astroChasqui/q2_tutorial">here</a>.

Examples
--------

Find spectroscopic parameters of a sample of stars using the Sun as the reference star (strict line-by-line differential analysis):

```python
import q2
data = q2.Data('stars.csv', 'lines.csv')
sp = q2.specpars.SolvePars(grid='marcs')
q2.specpars.solve_all(data, sp, 'solution.csv', 'Sun')
```

Measure elemental abundances of a sample of stars with respect to the solar abudances (line-by-line):

```python
species_ids = ['CI', 'OI', 'BaII']
q2.abundances.all(data, species_ids, 'abundances.csv', 'Sun')
```

MOOG
----

q2 needs a recent version of <a href="http://www.as.utexas.edu/~chris/moog.html">MOOGSILENT</a> (>2003) installed in your system.

Data
----

You wont be able to use all the package features without the associated data, particularly the model atmospheres and isochrone grids. Use <a href="http://www.astrochasqui.com/projects/astro/share/q2Data.tar.gz">this link</a> to download the data and make sure it is extracted into the q2/Data folder.

Contact info
------------

Ivan Ramirez (ivan@astro.as.utexas.edu)
