.. _consortiumob:

The Consortium Object
=====================

The Consortium Class is the central axis of the MMODES package.
It carries all the necessary parameters to run a dynamic simulation of a microbial
community.

Most common parameters
~~~~~~~~~~~~~~~~~~~~~~

.. py:class:: Consortium(max_growth = 10, v = 1, stcut = 1e-4, title = "draft_cons", mets_to_plot = [], work_based_on = "id", manifest = "", comets_output = False)

 * **v** is the volume of the modeled space. Amounts will be transformed to  concentrations using this parameters (*default = 1*). Units are arbitrary, but L are used by convention.
 * **stcut** is the limit of biomass flux (growth increment) where the simulation is considered to have reached a stable state and stops. Turn to a negative number to keep the simulation running (*default = 1e-4*).
 * **title** of the generated plot of the simulation
 * **mets_to_plot** are the metabolites to be later plotted.
 * **work_based_on** *(="id" | "name")* is a REALLY important parameter. It indicates whether extracellular metabolite names or ids should be used to communicate models and understand the medium. One should use the attribute (id or name) that is consistent among all the GEM models (just consistency on the extracellular metabolites is required) (*default = "id"*)
 * **manifest**: if a non-empty string is provided, it will output a fluxes TSV file to this path (default="").
 * **comets_output**: bool, whether output (including fluxes) should be written in COMETS-like format.

 COMETS-like output can be used to visualize reaction fluxes in `VisAnt <http://visant.bu.edu/>`_.

Setting models
~~~~~~~~~~~~~~
dModel objects are containers of COBRA model objects, with some
features to compute the multi-strain simulation. The method to add a GEM model to
the :py:class:`Consortium` is the following:

.. py:method:: Consortium.add_model(mod_path, volume_0, solver = "glpk", method = "pfba", dMets = {}, limit = False)

* **mod_path**, string: path to the model in SBML, MAT or JSON format.
* **volume_0**, float: initial concentration (usually, g/L) of biomass.
* **dMets**: dictionary of metabolite.id : dMetabolites.
* **limit**: maximum biomass values that is allowed for this particular dModel (*default = False*, no limitation).

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

**NOTE**: Concentration units are arbitrary, although the convention dictates *mmol/L* for metabolites and *g/L* for biomass.
Take into account consistency among units when instantiating the medium. Also, time is assumed to be in hours.

Running the Community simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. py:method:: Consortium.run(maxT=10, integrator='vode', stepChoiceLevel=(0., 0.5, 1000.), verbose = False, outf = "plot.tsv", outp = "plot.png", plot = True)

   Starts the community simulation, solving the system of ODE's.

* **maxT**, in time simulation units (hours), simulation will stop when it reaches this parameter (*default = '10'*).
* **integrator** *('vode' 'dopri5' 'fea'  'rk4')* type of ODE integrator (*default = 'vode'*).
* **stepChoiceLevel** (0, max time step, max number of time-steps) for *vode* and (time-step, 0, max number of time-steps) for the rest of integrators (*default = 0., 0.5, 100*).
* **verbose**: bool, a verbose simulation will show a progress bar and the reason of exiting the simulation.
* **outf** string, path where the output will be generated.
* **outp** string, path where the plot will be generated.
* **plot** bool, whether to generate the plot.

Assigning a *maxT* parameter doesn't guarantee to reach that time, since simulation
could be stop when it reaches the maximum number of steps (in *stepChoiceLevel*)
or when it's stabilized (*stcut* in the :py:class:`Consortium`).

| Once the simulation is finished, the output could be later generated:

.. code:: python3

    from mmodes.vis import plot_comm
    plot_comm(cons) # cons is a Consortium object which has already run

Adding perturbations
~~~~~~~~~~~~~~~~~~~~

Metabolites and perturbations are added with the following method:

.. py:method:: Consortium.add_mets(pert, concentration = False)

* **pert** is a dictionary with the format of media. Additionally, keys corresponding to model ID's can be used to add biomass to a model.
* **concentration** whether amounts in pert should be transformed to concentration units.

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
