.. _consortiumob:

The Consortium Object
=====================

The Consortium Class is the central axis of the MMODES package.
It carries all the necessary parameters to run a dynamic simulation.

Most common parameters
~~~~~~~~~~~~~~~~~~~~~~

.. py:class:: Consortium(max_growth = 10, death_rate = 0, v = 1, stcut = 1e-4, title = "draft_cons", mets_to_plot = [], work_based_on = "id", manifest = "", comets_output = False)

 * **v** is the volume of the modeled space. It's important to work with concentrations (*default = 1*).
 * **stcut** is the limit of biomass flux (growth increment) where the simulation is considered to have reached a stable state and stops. Turn to a negative number to keep the simulation running (*default = 1e-4*).
 * **title** is just to write it in a possible automatically generated plot of the simulation
 * **mets_to_plot** are the metabolites to be later plotted.
 * **work_based_on** *(="id" | "name")* is a REALLY important parameter. It indicates whether extracelullar metabolite names or ids should be used to communicate models and understand the medium. One should use the attribute (id or name) that is consistent among all the GEM models (just consistency on the extracellular metabolites is required) (*default = "id"*)
 * **manifest**: if a non-empty string is provided, it will output a fluxes tsv to this path (default="").

Setting models
~~~~~~~~~~~~~~
dModel objects are containers of COBRA model objects, with some
features to compute the multi-strain simulation. More information
can be found on _dModel_obj. The method to add GEM models to the :py:class:`Consortium`
is the following:

.. py:method:: Consortium.add_model(mod_path, volume_0, solver = "glpk", method = "pfba", dMets = {}, limit = False)

* **mod_path**, string: path to the model in SBML, MAT or JSON format.
* **volume_0**, float: initial concentration (usually, g/L) of biomass.
* **dMets**: dictionary of metabolite.id : dMetabolites
* **limit**: maximum biomass values that is allowed for this particular dModel.

More information can be found on the :doc:`dModel`.

Setting media
~~~~~~~~~~~~~
Consortium.media is a simple dictionary of metabolite ids/names : concentrations.
It can be passed as a simple dictionary Also,
a handy method is provided to read from a JSON file object (with the form of a dictionary).

.. py:method:: Consortium.set_media(media, concentration = False)

   Adds **media** as the Consortium medium object. If concentration is True, values of the dictionary
   will be converted to concentrations (using the volume parameter).

.. py:method:: Consortium.media_from_json(jfile, concentration = False)

   Uploads medium from a **jfile** path. This path corresponds to a JSON file which contains a dictionary.
