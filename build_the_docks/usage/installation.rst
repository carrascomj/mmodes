Installing MMODES
=================

A Python version >= 3.6 is required in order to install MMODES (install `Python 3 <https://www.python.org/downloads/>`_).
To avoid messing up with versions, working on a `virtualenv <https://virtualenv.pypa.io/en/stable/>`_ might be a good idea.

PIP Install
~~~~~~~~~~~
It's recommended to install MMODES via pip. This way, the latest stable version
of the package is guaranteed. Naturally, in order to install MMODES through pip,
pip tool is required (if it isn't installed, check the `pip's installation instructions <https://pip.pypa.io/en/stable/installing/>`_).
Usually, the *--user* flag is required. On a bash shell:

.. code:: bash

    pip3 install mmodes --user

Build from SOURCE
~~~~~~~~~~~~~~~~~
Otherwise, you can build from the :ref:`GitHub repository`_. `COBRApy <https://opencobra.github.io/cobrapy/>`_
version used is 15.2. `Scipy <https://www.scipy.org/>`_, `numpy <https://www.numpy.org/>`_,
`matplotlib <https://matplotlib.org/>`_ and `dill <https://pypi.org/project/dill/>`_
are also required. A cobra version >= 14 should also work, albeit not being guaranteed.

.. code:: bash

    git clone https://github.com/carrascomj/mmodes.git # or ssh
    cd path_to_mmodes/mmodes
    sudo python3 setup.py install

DOCKER Install
~~~~~~~~~~~~~~
A docker image is currently under development.

Uninstall
~~~~~~~~~
Uninstalling can be accomplished via *pip*:

.. code:: bash

    pip3 uninstall mmodes --user # if user install
    sudo pip3 uninstall mmodes # if superuser install (from source)
