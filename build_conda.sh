git clone git@github.com:joe-fennell/pyeosim.git
cd pyeosim
conda-build . --no-anaconda-upload
conda install pyeosim --use-local
