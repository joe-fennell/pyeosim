import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyeosim",  # Replace with your own username
    version="0.0.2",
    author="Joseph T. Fennell",
    author_email="info@joefennell.org",
    description="Python Earth Observation Simulator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/joe-fennell/pyeosim/",
    packages=setuptools.find_packages(),
    package_dir={'georis': 'georis/'},
    package_data={'georis': ['*.tif', '*.yaml', '*.shp', '*.shx', '*.cpg',
                             '*.dbf', '*.prj', '*.yaml']},
    include_package_data=True,
    install_requires=[
          'xarray',
          'pyyaml',
          'scipy',
          'numpy',
          'dask',
          'Click'
      ],
    scripts=[],
    # entry_points='''
    #     [console_scripts]
    #     georis=georis.cli:georis
    # ''',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
