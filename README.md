q2 Python package
=================

q2 allows you to use MOOG (in its SILENT mode) to calculate elemental abundances of stars and/or determine their atmospheric parameters using the standard techniques of iron line excitation/ionization equlibrium.

Examples
--------

```python
import q2
data = q2.Data('stars.csv', 'lines.csv')
sp = q2.specpars.SolvePars(grid='marcs')
q2.specpars.solve_all(data, sp, 'solution.csv', 'Sun')
```

MOOG
----

Folder Moog contains files that you should replace in your MOOG installation in order to use q2. Copy them into your computer and re-compile MOOG.

Data
----

You wont be able to use all the package features without the associated data, particularly the model atmospheres and isochrone grids. These files are large and I haven't yet figured out the best way of sharing them. In the meantime send me an email and ask me for a zip or tar file with the Data folder.

Contact info
------------

Ivan Ramirez (ivan@astro.as.utexas.edu)
