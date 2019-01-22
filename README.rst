MMODES
#########################
MMODES (Metabolic Models based Ordinary Differential Equations Simulation) provides an straight-forward framework to run dynamic simulations of metabolic models consortia using cobrapy and scipy. The general schema of the script was taken from `DAPHNE <https://github.com/libretro/daphne/tree/master/daphne>`__, though it was simplified and rebuilt for a general case. It also provides I/O utilities to generate COMETS-like files and to manage lenguages of models extracellular metabolites.

Installation
-------------
On a Linux shell::
    cd path_to_mmodes/mmodes; sudo python3 setup.py install

Example script simulation
-------------
An example with one strain can be found on "example.py"
I'd recommend to work on a virtualenv since the package uses the latest python3 versions of the modules required.
More information about virtualenvs `here <https://www.configserverfirewall.com/ubuntu-linux/create-python-virtualenv-ubuntu/>`__.
