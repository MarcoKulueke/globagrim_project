from setuptools import setup, find_packages

setup(
    name="globagrim",
    version="0.1",
    packages=find_packages(exclude=["tests*"]),
    license="MIT",
    description="Global Atmospheric Grid Point Model",
    long_description=open("README.md").read(),
    install_requires=["numpy", "xarray", "netcdf4", "matplotlib", "cython", "geos", "pyshp", "pyproj", "six", "cartopy", "matplotlib"],
    url="https://github.com/MarcoKulueke/globagrim_project",
    author="Marco Kulüke",
    author_email="marco.kulueke@gmail.com",
)
