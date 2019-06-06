.. dModel_obj:

The dModel Object
=================

The dModel object expands the COBRApy Model object with some features. It acts as
a container, not as a subclass of it.

Adding models to the Consortium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Usually, models are added with Consortium methods. Please refer to :doc:`consortium` to see how.
It's important to note that when a dModel is added, it will be optimized to check if
it's operative. The model.id cobra attribute will be used as dModel.id.

Furthermore, the growth increment of the dModel will be taken from de biomass function. It doesn't
need to be the objective funtion that will be optimized during the simulation, MMODES
will take the first reaction that starts with "biomass". This allows some functionalities
but can be troublesome when more than one biomass reaction is present in the model.

Once added to the Consortium, dModels are accesible from the attribute models.
This attribute is a dictionary of dModel.id : dModel object.

Limiting growth of the model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
At the moment of adding the model, **limit** is a parameter that can be tuned:

* if *False* (default), no limitations will be applied;
* if *True*, model won't grow, but the solution of fluxes will be used to update the medium and added to the output files;
* if a *numeric* value is passed, the model will be allowed to grow until this amount is reached, then it will behave as if *limit* = True.

Death rate
~~~~~~~~~~
death_rate is an attribute that can be added once the dModel is instantiated
(or added to the Consortium). This parameter is **incompatible** with *limit*.
