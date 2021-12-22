import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyeosim",  # Replace with your own username
    version="0.1",
    author="Joseph T. Fennell",
    author_email="info@joefennell.org",
    description="Python Earth Observation Simulator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/joe-fennell/pyeosim/",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        'xarray',
        'pyyaml',
        'scipy',
        'numpy',
        'dask',
        'py6s',
        'pandas',
        'xlrd',
        'netcdf4',
        'matplotlib'
    ],
    scripts=[],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
