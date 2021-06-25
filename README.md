# PyThermo
## Installation
### Installation with pip
The package may be installed via pip with the following Shell command:
```console
pip install pythermo
```

### Manual installation
First, download this GitHub repository as a .zip file, and save the contents in a directory. Open a shell environment and type out the following
```console
cd <installation directory>
```
where <installation directory> is the directory which the .zip file has extracted to. Then run the following command to install the package
```console
python setup.py install
```

##Project Structure, Versioning, and Documentation

Project structure follows python-guide.org recommendations.

Docstring format follows Google style recoomendations.

Versioning for publishing to PyPI follows the "major.minor.patch" format based on https://semver.org/ recommendations.

* major version - when you make incompatible API changes,
* minor version - when you add functionality in a backwards-compatible manner, and
* patch version - when you make backwards-compatible bug fixes.

The Markdown cheat sheet is a useful reference for keeping documentation up to date.


## General Information
After installation, the package is imported in python by the following line of code.
```python
from pythermo import pythermo as pt
```
