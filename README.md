# GLOBAGRIM Global <sub> Atmospheric Grid Point Model <\sub>

This model is based on the GLOBAGIM of the CLISAP (more information can be found via: https://www.clisap.de/research/b:-climate-manifestations-and-impacts/crg-dynamical-systems/models/globagim/)

## Installation
Create a conda environment with
```

conda env create -f environment.yml
```
This installs GLOBAGRIM and all dependencies.
Then, activate the environment with:
```
conda activate globagrim_env
```
---
If you already have ` numpy, xarray, netcdf4, matplotlib, cartopy, jupyter` installed, you can also add the GLOBAGRIM package to your current environment by executing:
```
pip install git+git://github.com/MarcoKulueke/globagrim_project.git#egg=globagrim
```
---

## Example Usage

Import the model from globagrim:

```
from globagrim import model
```

You can run the model with the following parameters:
- `NJ` Number of west-east circles of longitude
- `NK` Number of north-south circles of latitude
- `NL` Number of vertical levels 
- `TF` Total integration time in hours
- `DT` Integration time step in seconds
- `IEXP` Experiment Number (`1` low pressure system over the pacific, `2`stream over montain in North America)

All parameters are optional. If you do not pass any parameters the default parameters will be used.

Below you can see an example call:
```
model.run(NJ=40, NK=80, TF=1)
```

You can visualize the results with: `model.plot(<variable name>, <model time step>)`
```
model.plot('PSG', 40)
```
Again, all parameters are optional. You can specify the visualization with the following parameters:
- `var_name` String variable name
- `time_step` Integer model time step
- `center` Touple `(<lon>,<lat>)` for center of visualization
