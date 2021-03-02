import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pythermo",
    version="1.0.0",
    author="Xiaodong Liang & Daniel Qvistgaard",
    author_email="s164067@student.dtu.dk",
    description="Thermodynamic Modelling Tools for Phase Equilibria",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/olisimdan/PyTherm",
    packages=['pythermo'],
    setup_requires=['numpy'],
    install_requires=['numpy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
   # include_package_data=True,
    package_data = {
        'pythermo' : ['dll_files/xThermoWrapper.F90'],
        'pythermo' : ['dll_files/*.dll'],
    },
    data_files = None,
)