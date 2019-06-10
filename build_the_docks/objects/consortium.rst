.. _consortiumob:

The Consortium Object
=====================

The Consortium Class is the central axis of the MMODES package.
It carries all the required parameters to run a dynamic simulation of a microbial
community.

Most common parameters
~~~~~~~~~~~~~~~~~~~~~~

.. py:class:: Consortium(max_growth = 10, v = 1, stcut = 1e-8, title = "draft_cons", mets_to_plot = [], work_based_on = "id", manifest = "", comets_output = False)

 :param if max_growth: is the maximum biomass that a GEM model is allowed to reach (*default = 10*)
 :param float v: is the volume of the modeled space. Amounts will be transformed to  concentrations using this parameters (*default = 1*). Units are arbitrary, but L are used by convention.
 :param float stcut: is the limit of biomass flux (growth increment) where the simulation is considered to have reached a stable state and stops. Turn to a negative number to keep the simulation running (*default = 1e-8*).
 :param str title: of the generated plot of the simulation
 :param list mets_to_plot: are the metabolites to be later plotted.
 :param str work_based_on: *(="id" | "name")* is a REALLY important parameter. It indicates whether extracellular metabolite names or ids should be used to communicate models and understand the medium. One should use the attribute (id or name) that is consistent among all the GEM models (just consistency on the extracellular metabolites is required) (*default = "id"*)
 :param str manifest: if a non-empty string is provided, it will output a fluxes TSV file to this path (default="").
 :param bool comets_output: whether output (including fluxes) should be written in COMETS-like format.

 COMETS-like output can be used to visualize reaction fluxes in `VisAnt <http://visant.bu.edu/>`_.

Setting models
~~~~~~~~~~~~~~
dModel objects are containers of COBRA model objects, with some
features to compute the multi-strain simulation. The method to add a GEM model to
the :py:class:`Consortium` is the following:

.. py:method:: Consortium.add_model(mod_path, volume_0, solver = "glpk", method = "pfba", dMets = {}, limit = False)

 :param string mod_path: path to the model in SBML, MAT or JSON format.
 :param float volume_0: initial concentration (usually, g/L) of biomass.
 :param dictionary dMets: of metabolite.id : dMetabolites.
 :param maximum limit: biomass values that is allowed for this particular dModel (*default = False*, no limitation).

More information about the *limit* parameter and the dModel in general can be found on :doc:`dModel`.

Setting media
~~~~~~~~~~~~~
Consortium.media is a simple dictionary of metabolite ids/names : concentrations.
It can be passed as a simple dictionary. Also,
a handy method is provided to read from a JSON file object (with the structure of a dictionary).

.. py:method:: Consortium.set_media(media, concentration = False)

   Adds **media** as the Consortium medium object. If concentration is True, values of the dictionary
   will be converted to concentrations (using the volume parameter).

.. py:method:: Consortium.media_from_json(jfile, concentration = False)

   Uploads medium from a **jfile** path. This path corresponds to a JSON file which contains a dictionary.

.. note::
    Concentration units are arbitrary, although the convention dictates *mmol/L* for metabolites and *g/L* for biomass.
    Take into account consistency among units when instantiating the medium. Also, time is assumed to be in hours.

Running the Community simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. py:method:: Consortium.run(maxT=10, integrator='vode', stepChoiceLevel=(0., 0.5, 1000.), verbose = False, outf = "plot.tsv", outp = "plot.png", plot = True, actualize_every = float('-inf'))

   Starts the community simulation, solving the system of ODE's.

   :param float maxT: in time simulation units (hours), simulation will stop when it reaches this parameter (*default = '10'*).
   :param str integrator: *('vode' 'dopri5' 'fea'  'rk4')* type of ODE integrator (*default = 'vode'*).
   :param str stepChoiceLevel: (0, max time step, max number of time-steps) for *vode* and (time-step, 0, max number of time-steps) for the rest of integrators (*default = 0., 0.5, 100*).
   :param bool verbose: a verbose simulation will show a progress bar and the reason of exiting the simulation (*default = False*).
   :param str outf: path where the output will be generated (*default= plot.tsv*).
   :param str outp: path where the plot will be generated (*default= "plot.png"*).
   :param bool plot: whether to generate the plot (*default= True*).

Assigning a *maxT* parameter doesn't guarantee to reach that time, since simulation
could be terminated when it reaches the maximum number of steps (in *stepChoiceLevel*),
when it's stabilized (*stcut* in the :py:class:`Consortium`) or when some simulation
exceeds that maximum growth (*max_growth* in the :py:class:`Consortium`, different
from *limit* in :doc:`dModel`).

| Once the simulation is finished, the output could be later generated:

.. code:: python3

    from mmodes.vis import plot_comm
    plot_comm(cons) # cons is a Consortium object which has already run

On this point, other metabolites could've been plotted changing *mets_to_plot*
attribute of :py:class:`Consortium`.

.. warning::
    Please, take into account that the results will be *appended* to **outf** and
    the plot will be generated from this path. Thus, make sure a file with this
    name doesn't exist before a simulation is started.

Adding perturbations
~~~~~~~~~~~~~~~~~~~~

Metabolites and perturbations are added with the following method:

.. py:method:: Consortium.add_mets(pert, concentration = False)

 :param dict pert: same format as media. Additionally, keys corresponding to model ID's can be used to add biomass to a model.
 :param bool concentration: whether amounts in pert should be transformed to concentration units.

Once the method has been called, the same Consortium can run again, simulating a perturbation.

| Putting it all together:

.. code:: python3

    from mmodes import Consortium
    from mmodes.vis import plot_comm

    cons = Consortium()
    cons.add_models(mod_path = "path_to_some_model_file.mat", volume_0 = 0.001)
    cons.add_models(mod_path = "path_to_some_other_model_file.xml", volume_0 = 0.0012)
    cons.media = cons.media_from_json(jfile = 'some_dict_file.json')
    cons.run(plot = False)
    cons.add_mets({'glc[e] : 0.02'})
    cons.run(plot = False)
    plot_comm(cons)
