# Global Atmospheric Grid Point Model (GLOBAGRIM)

This model is based on the GLOBAGIM of the CLISAP (more information can be found via: https://www.clisap.de/research/b:-climate-manifestations-and-impacts/crg-dynamical-systems/models/globagim/)

## Installation
```
pip install git+git://github.com/MarcoKulueke/globagrim_project.git#egg=globagrim
```

## Usage
```
from globagrim import run
run.run(NJ=20, NK=10, TF=1) # 20 grid points in W-E, 10 in N-S, 1 h integration time
```
