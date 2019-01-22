MMODES
#########################
MMODES (Metabolic Models based Ordinary Differential Equations Simulation) provides an straight-forward framework to run dynamic simulations of metabolic models consortia using cobrapy and scipy. The scipy ODE implementation of the script was inspired by `DAPHNE <https://github.com/libretro/daphne/tree/master/daphne>`__, although this code was built for a general case. It also provides I/O utilities to generate and parse COMETS-like files and to manage lenguages of models extracellular metabolites.

Installation
-------------
On a Linux shell:

.. code:: bash

    cd path_to_mmodes/mmodes
    sudo python3 setup.py install

Example script simulation
-------------
An example with one strain can be found on "example.py"
I'd recommend to work on a virtualenv since the package uses the latest python3 versions of the modules required.
More information about virtualenvs `here <https://www.configserverfirewall.com/ubuntu-linux/create-python-virtualenv-ubuntu/>`__.
