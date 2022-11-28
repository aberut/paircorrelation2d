# paircorrelation2d

paircorrelation2d is a Python module to compute the 2D pair correlation function ([radial distribution function](https://en.wikipedia.org/wiki/Radial_distribution_function)) g(r) for a set of points, corrected to take account of the boundary effects.

A more detailed explanation on the different type of corrections to take account of boundary effects can be found in Natsuda Klongvessa (Ong)'s PhD Thesis, section 3.3.1 (pages 51-52), available online: [Study of Dense Assemblies of Active Colloids : collective Behavior and Rheological Properties](https://tel.archives-ouvertes.fr/tel-03282035/document).

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

See [example.ipynb](./example.ipynb) for a detailed presentation.

## Requirements

- [numpy](https://numpy.org/) (>=1.19)
- [matplotlib](https://matplotlib.org/index.html) (>=3.3)
- [shapely](https://shapely.readthedocs.io/en/latest/manual.html) (>=1.8)

## License

[MIT](https://choosealicense.com/licenses/mit/)
