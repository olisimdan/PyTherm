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

## Testing
To run unit tests, open a shell environment and run the following command in the top-level directory.
```console
python -m unittest discover -v
```

## Project Structure, Versioning, and Documentation

Project structure follows python-guide.org recommendations.

Docstring format follows Sphinx markdown recommendations. This allows for syntax documentation via Sphinx.

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

## Learning Material
Learning material is available in the repository folders "Learning Material" and "docs". Within "Learning Material" you will find Jupyter Notebooks and sample python codes. In the "docs" folder you will find user guide to help set up the package, and a syntax page where the syntax of all functions are described. Even after installing with pip it is recommended to download the repository as a .zip file directly from GitHub in order to have access to these folders.