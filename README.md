# paircorrelation2d

paircorrelation2d is a Python module to compute the 2D pair correlation function ([radial distribution function](https://en.wikipedia.org/wiki/Radial_distribution_function)) g(r) for a set of points, corrected to take account of the boundary effects.

## Installation

Simply put the script in your working directory (or in any directory you added to your python path).

## Usage

### For the impatient

```python
import matplotlib.pyplot as plt
from paircorrelation2d import pcf2d

#Assuming array_positions is a Nx2 numpy array containing the 2D coordinates of N points
#Assuming bins_distances is a Mx1 numpy array containing bin edges defining the values of r for which g(r) is going to be computed

[g_of_r,r]=pcf2d(array_positions, bins_distances)

plt.plot(r,g_of_r)
```
### For the patient

See [example.ipynb](https://cameleon.univ-lyon1.fr/aberut/paircorrelation2d/-/blob/master/example.ipynb) for a detailed presentation.

## Requirements

- [numpy](https://numpy.org/) (tested with ver 1.19.1)
- [matplotlib](https://matplotlib.org/index.html) (tested with ver 3.3.1)
- [shapely](https://shapely.readthedocs.io/en/latest/manual.html) (tested with ver 1.7.0)

## License

[MIT](https://choosealicense.com/licenses/mit/)
