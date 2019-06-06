
.. image:: build_the_docks/logo_f.svg
   :width: 400px
   :alt: MMDOES logotype
   :align: center

MMODES
======
MMODES (Metabolic Models based Ordinary Differential Equations Simulation) provides an straight-forward framework to run dynamic simulations of metabolic models consortia using cobrapy and scipy.
It implements a version of dFBA built on top of cobrapy and scipy. The scipy ODE implementation of the script was inspired by `DAPHNE <https://github.com/libretro/daphne/tree/master/daphne>`__,
although this code was built for a general case. It also provides I/O utilities to generate and parse COMETS-like files and to manage "languages" of GEM models.


Documentation
~~~~~~~~~~~~~
Documentation can be found on https://mmodes.readthedocs.io/en/latest/.

Installation
~~~~~~~~~~~~
I'd recommend to work on a virtualenv since the package uses the latest python3 versions of the modules required.
More information about virtualenvs `here <https://www.configserverfirewall.com/ubuntu-linux/create-python-virtualenv-ubuntu/>`_.
Clone this repository and install it via setup.py. On a bash shell:

.. code:: bash
    git clone https://github.com/carrascomj/mmodes.git # or ssh
    cd path_to_mmodes/mmodes
    sudo python3 setup.py install

Docker
~~~~~~
A `Docker <https://www.docker.com/get-started>`_ container is under development.

Example script simulation
~~~~~~~~~~~~~~~~~~~~~~~~~
An example with one strain can be found on "example.py".
