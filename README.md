# Global Atmospheric Grid Point Model (GLOBAGRIM)

This model is based on the GLOBAGIM of the CLISAP (more information can be found via: https://www.clisap.de/research/b:-climate-manifestations-and-impacts/crg-dynamical-systems/models/globagim/)

## Installation
```
pip install git+git://github.com/MarcoKulueke/globagrim_project.git#egg=globagrim
```

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
- `IEXP` Experiment Number

All parameters are optional. If you do not pass any parameters the default parameters will be used.

Below you can see an example call:
```
model.run(NJ=40, NK=80, TF=1)
```

You can visualize the results with `model.plot(<variable name>, <model time step>)`
```
model.plot('PSG', 40)
```
