
.. image:: build_the_docks/logo_f.svg
   :width: 400px
   :alt: MMDOES logotype
   :align: center

What is MMODES?
===============

MMODES (Metabolic Models based Ordinary Differential Equations Simulation) provides
an straight-forward framework to run dynamic simulations of metabolic models communities
using `COBRApy <https://opencobra.github.io/cobrapy/>`_ and `Scipy <https://www.scipy.org/>`_, `numpy <https://www.numpy.org/>`_.

It implements a version of dFBA. The scipy ODE
implementation of the script was inspired by `DAPHNE <https://github.com/libretro/daphne/tree/master/daphne>`__,
although this code was built for a general case. It also provides I/O utilities
to generate and parse COMETS-like files and to manage "languages" of GEM models.


Documentation
~~~~~~~~~~~~~
User documentation can be found on https://mmodes.readthedocs.io/.

Installation
~~~~~~~~~~~~
Working on a **virtualenv** is highly recommended since the package uses the latest python3 versions of the modules required.
More information about Python virtualenv's `here <https://virtualenv.pypa.io/en/stable/>`_.

MMODES is on PyPi, therefore, given that **pip** in installed:

.. code:: bash

    pip upgrade pip --user
    pip3 install mmodes --user


Another option is **building from source**. Clone this repository and install it via *setup.py*:

.. code:: bash

    git clone https://github.com/carrascomj/mmodes.git # or ssh
    cd path_to_mmodes/mmodes
    sudo python3 setup.py install

On both cases, **uninstalling** can be accomplished via pip:

.. code:: bash

    pip3 uninstall mmodes --user # if user install
    sudo pip3 uninstall mmodes # if superuser install (from source)

Docker
~~~~~~
A `Docker <https://www.docker.com/get-started>`_ container is under development.

Example script simulation
~~~~~~~~~~~~~~~~~~~~~~~~~
An example with one strain can be found on "example.py".
