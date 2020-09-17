[![Travis CI w/ Logo](https://img.shields.io/travis/grburgess/nazgul/master.svg?logo=travis)](https://travis-ci.org/grburgess/nazgul)  [![codecov](https://codecov.io/gh/grburgess/nazgul/branch/master/graph/badge.svg)](https://codecov.io/gh/grburgess/nazgul)
[![PyPi Downloads](http://pepy.tech/badge/nazgul)](http://pepy.tech/project/nazgul)
[![PyPI version fury.io](https://badge.fury.io/py/nazgul.svg)](https://pypi.python.org/pypi/nazgul/)
[![Documentation Status](https://readthedocs.org/projects/nazgul/badge/?version=latest)](https://nazgul.readthedocs.io/?badge=latest)
![GitHub contributors](https://img.shields.io/github/contributors/grburgess/nazgul)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)


# nazgul

![alt text](https://raw.githubusercontent.com/grburgess/nazgul/master/logo.png)

## What is this?
Nazgul is a framework for performing GRB localization via fitting non-parametric models to their data time-series and computing the the time delay between them. It is currentrly built upon the magic of Stan and implements a parallel version of non-stationary Random Fourier Features. The idea is get away from heuristic methods such as cross-correlation which do not have a self-consistent statitical model. 

The idea is that satellites throughout the Sol system observe gamma-ray bursts at different times due to the finite speed of light. This creates a time delay in their observed light curves which can be used to triangulate the gamma-ray burst position on the sky. These triangulation create annuli or rings on the sky which Nazgul searchs for so that it, in the darkness, it can bind them to a location on the sky. (We are nerds, get over it).

![alt text](https://raw.githubusercontent.com/grburgess/nazgul/master/idea.png)


The heriarchical model is shown below and details can be found in [link here](). If you find the method and/or code useful in your research we ask that you please cite the paper. 

![alt text](https://raw.githubusercontent.com/grburgess/nazgul/master/model.png)



The sister program to simulate time-delayed light curves is [pyIPN](https://github.com/grburgess/pyipn) and can be used to generate time-delayed light curves for algorithm testing. 

## Installation
```bash
pip install nazgul
```

# Note!

Nazgul is still under heavy development! Expect problems!




