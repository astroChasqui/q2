# q2 Python package

<a href="http://www.as.utexas.edu/~chris/moog.html">
    <img src="https://img.shields.io/badge/powered%20by-MOOG-orange"/>
</a>
<a href="https://matplotlib.org/3.1.1/index.html">
    <img src="https://img.shields.io/badge/powered%20by-matplotlib-orange"/>
</a>
<a href="https://docs.python.org/3/library/os.html">
    <img src="https://img.shields.io/badge/powered%20by-os-orange"/>
</a>
<a href="https://numpy.org/">
    <img src="https://img.shields.io/badge/powered%20by-numpy-orange"/>
</a>
<a href="https://scipy.org/">
    <img src="https://img.shields.io/badge/powered%20by-scipy-orange"/>
</a>
<a href="https://bokeh.org/">
    <img src="https://img.shields.io/badge/powered%20by-bokeh-orange"/>
       
The `q2` code allows you to use the 2019 <a href="http://www.as.utexas.edu/~chris/moog.html">MOOG</a> version (in its SILENT mode) to calculate elemental abundances of stars and/or determine their atmospheric parameters using the standard techniques of iron line excitation/ionization equilibrium. It also allows you to calculate other fundamental stellar parameters such as mass and age using isochrones. A tutorial is available <a href="https://github.com/astroChasqui/q2_tutorial">here</a>.

Installation
------------
The new version of `q2` can only be installed via pip. Simply try:

```bash
pip install qoyllur-quipu
```

If you have installed the old version of `q2`, you must delete it from your HOME directory and also remove its PATH from bashrc (.bash_profile for Mac OS). Once you have installed `q2`, open a Jupyter Notebook Environment (or iPython) and import `q2`.

```python
import q2
```

By importing `q2`, the latest version of <a href="http://www.as.utexas.edu/~chris/moog.html">MOOG</a> (2019) will begin to install. The only thing you need to do is declare the kind of machine you are using, i.g., 'pcl' for linux (Linux Mint, Ubuntu, etc) and 'mac' for Mac Os. That's all folks. This process is done only once. Note that the `q2` package requires Python 3.7 or later. 

Future releases will be upgraded via:

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

Author
------
- Ivan Ramirez (iramirez@tacomacc.edu)

Maintainers
-----------
- <a href="https://github.com/ramstojh">Jhon Yana Galarza</a> 
- <a href="https://github.com/kaykeigh">Kayleigh Meneghini</a>

Preferred citation
------------------

If you make use of this code, please cite <a href="https://doi.org/10.1051/0004-6361/201424244">Ramirez et al. 2014, A&A, 572, A48</a>. The BibTeX entry for the paper is:
```bibtex
@ARTICLE{Ramirez2014,
       author = {{Ram{\'\i}rez}, I. and {Mel{\'e}ndez}, J. and {Bean}, J. and {Asplund}, M. and {Bedell}, M. and {Monroe}, T. and {Casagrande}, L. and {Schirbel}, L. and {Dreizler}, S. and {Teske}, J. and {Tucci Maia}, M. and {Alves-Brito}, A. and {Baumann}, P.},
        title = "{The Solar Twin Planet Search. I. Fundamental parameters of the stellar sample}",
      journal = {\aap},
     keywords = {stars: abundances, stars: fundamental parameters, planetary systems, Astrophysics - Solar and Stellar Astrophysics},
         year = 2014,
        month = dec,
       volume = {572},
          eid = {A48},
        pages = {A48},
          doi = {10.1051/0004-6361/201424244},
archivePrefix = {arXiv},
       eprint = {1408.4130},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2014A&A...572A..48R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
