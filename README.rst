mmodes
Metabolic Models based Ordinary Differential Equations Simulation
-------------
The package provides an straight-forward framework to run dynamic simulations
of metabolic models consortia using coprapy and scipy.

Installation
-------------
On a Linux shell:
$ cd path_to_mmodes/mmodes
$ sudo python3 setup.py install

Example script simulation
-------------
On python:
import mmodes
# 1) instantiate Consortium
cons = mmodes.Consortium(stcut = -100, mets_to_plot = ["glc_D[e]", "pi[e]"])
# 2) add models
cons.add_model("ModelsInput/BGN4_eu.xml", 0.001)
cons.add_model("ModelsInput/Fae_ok.xml", 0.001)
# 3) instantiate dMetabolites
# coming soon...
# 4) instantiate media
cons.media = cons.set_media(media)
# 5) Run it
cons.run(verbose=True)
