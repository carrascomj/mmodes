#!/usr/bin/env python3

# All the classes, modules required and exceptions of the consortium package.
# The package provides an straight-forward framework to run dynamic simulations
# of metabolic models consortia using coprapy and scipy.

# Author: Jorge Carrasco Muriel
# e-mail: jorge.cmuriel@alumnos.upm.es
# Date of first version: 30/12/2018

import os
import re
import cobra
import random
import warnings
import numpy as np
from scipy.integrate import ode
from copy import deepcopy as dcp
from mmodes.io import Manifest
from mmodes.io import load_model
from mmodes.vis import plot_comm
from . import _fsolvers
from mmodes.io.models_in import ProBar


class NoBiomassException(Exception):
    '''
    NoBiomassException is raised when a dynamicModel.reactions attribute does not
    contain 'biomass' as key
    Biomass is needed to follow growth!
    '''
    pass

class NotIntegratorError(Exception):
    '''
    When integrator is False... check Consortium.run() parameters
    '''
    pass

class InfeasibleSolution(Exception):
    '''
    If method (pFBA or FBA) reports Infeasible when dModel is instantiated
    '''
    pass

class Volume():
    '''
    Volume object that will be created inside a dModel object to contain biomass info.
    '''
    def __init__(self, volume_0, model):
        self.bm = self.get_bm_reaction(model) # biomass reaction cobra object
        self.q = volume_0 # biomass concentration
        return

    def get_bm_reaction(self, model):
        '''
        Look for the biomass reaction
        INPUT -> COBRA model object
        OUTPUT -> biomass rection cobra object
        '''
        # There are some strange biomass exchange reactions, so it doesn't properly
        # work sometimes.
        for reac in model.reactions:
            if reac.id.lower().find("biomass") != -1 and reac not in model.exchanges:
                bm = reac.id
                break
        else:
            raise NoBiomassException(f"No biomass reaction found in {model.id}")
            # f-strings require >= python 3.6!
        return bm

    def dupdate(self):
        '''
        Update value of biomass concentration
        '''
        self.q += self.bm.value*self.q
        return self.bm.value*self.q

class dModel():
    '''
    Expanded COBRA model object that Consortium object will manage to run the
    dynamic simulations.
    '''
    def __init__(self, mod_path, volume_0, solver = "glpk", method = "pfba", dMets = {}, work_based_on = "id", limit = False):
        self.path = mod_path if type(mod_path) is str else "Unknown"
        self.model = load_model(mod_path) # cobra model
        self.volume = Volume(volume_0, self.model) # volume object
        self.identifier = work_based_on.lower()
        self.exchanges = self.get_exchange_reactions() # dict met:exchange reactio
        self.dMets = dMets
        self.model.solver = solver
        if solver == "glpk":
            self.model.solver.timeout = 3 # avoid "glpk" hangs
        self.stuck = False
        self._death_rate = 0
        self.limit = limit
        self.method = method
        self.free_media()
        self.check_feasible()
        return

    def __str__(self):
        return "MMODES model object of {0}, path in {1} and actual biomass of {2}".format(self.model.id, self.path, self.volume.q)

    @property
    def death_rate(self):
        '''
        Private property just to apply the limit feature. If true -> dBM = 0
        '''
        return self._death_rate

    @death_rate.setter
    def death_rate(self, value):
        '''
        If value is a COBRA model and self.limit is true, the biomass reaction
        flux will be set to make dBM 0. Else, it's a normal setter.
        If limit is numeric, it won't affect death_rate until biomass of the
        model (g) reaches this limit value.
        '''
        try:
            self._death_rate = float(value)
        except TypeError:
            try:
                if self.limit:
                    if type(self.limit) == bool: # never growing
                        self._death_rate = value.reactions.get_by_id(self.volume.bm).flux
                    elif (type(self.limit == float) or type(self.limit == float)):
                        q = self.volume.q
                        if q >= limit: # it's reached the maximum biomass
                            vBM = value.reactions.get_by_id(self.volume.bm).flux
                            self._death_rate = vBM - 1 + limit/q # so dBM = q - limit
            except AttributeError:
                return

    def update(self):
        '''
        Updates value of biomass concentration
        '''
        return self.volume.dupdate()

    def add_biom(self, added):
        '''
        Add biomass to Volume
        '''
        self.volume.q += added
        return

    def get_exchange_reactions(self):
        '''
        model.exchanges cobra method isn't properly working. This function
        provides a fix as a dictionary.
        OUPUT -> dictionary, metabolite: exchange reactions
        '''
        dic_out = {}
        for reac in self.model.exchanges:
                for met in reac.metabolites:
                    if met.id.find("biomass") == -1 and met.compartment in ["e", "e0", "ExtraCellular", "extracellular"]: #strange metabolite in some models
                        if self.identifier == "name":
                            dic_out[met.name] = reac.id
                        else:
                            dic_out[met.id] = reac.id
        return dic_out

    def free_media(self):
        '''
        Eliminates boundaries associated with cobra.model.medium object. The
        simulation fix bounds, not medium.
        '''
        for reac in self.exchanges.values():
            # I don't know if all exchanges reactions are included in medium.
            self.model.medium[reac] = 1000
        return

    def opt(self, mod = ""):
        '''
        Solves FBA or pFBA
        INPUT: mod, cobra model, so it's able to work as context
        '''
        # print("Biomass of "+mod.id+" "+str(mod.reactions.Biomass.flux)) ############ DEBUG ############
        if mod == "":
            mod = self.model
        if self.method == "pfba":
            try:
                cobra.flux_analysis.pfba(mod)
            except: # TODO: specify the exception
                mod.optimize()
                self.stuck = True
        elif self.method == "fba":
            mod.optimize()
        return

    def check_feasible(self):
        with self.model as mod:
            if self.method == "pfba":
                try:
                    cobra.flux_analysis.pfba(mod)
                except: # TODO: specify the exception
                    raise InfeasibleSolution(f"pFBA is infeasiable for model {mod.id}. You may want to change the 'method' param to 'fba'.")
            elif self.method == "fba":
                try:
                    mod.optimize()
                except:
                    raise InfeasibleSolution(f"FBA is infeasible for model {mod.id}. You may want to check your model.")

    def add_dMet(self, met):
        '''
        Adds a dMetabolite to self.dMets
        INPUT -> met, dMetabolite object
        '''
        if self.identifier == "name":
            self.dMets[met.name] = met
        else:
            self.dMets[met.id] = met
        return

class dMetabolite():
    def __init__(self, id, concentration_fed = 0, Vmax = 20, Km = 0.01, each_fed = 0):
        self.id = id # as in models
        self.concentration_fed = concentration_fed # concentration to refresh
        self.Vmax = Vmax
        self.Km = Km
        self.each_time = each_fed

    def refresh(self, t):
        '''
        Computes times when the metabolite is added
        INPUT -> t, iteration time of the solver
        '''
        # TODO: since iteration won't always work with fixed step time, I guess
        # it's necessary to change each_time to an interval ????
        ref = True
        if self.each_time == 0: # avoid ZeroDivisionError
            ref = False
        elif t % self.each_time != 0:
            ref = False
        return ref

class Consortium():
    '''
    Main object were simulations are computed. It keeps dModels in a dictionary
    of id: <dModel>, media in a dictionary of met.id : <COBRA metabolite>, mets
    with specified experimental Michaellis Menten in a dictionary of met.id :
    <dMetabolites> and other parameters of the simulation.
        -> dinamicpFBA is where dModels are updated and optimized to compute biomass.
        -> run() is the main method, which controls the simulation.
        *-> plot_comm() is a function in "vis" dependency required to plot the output
    '''
    def __init__(self, media = {}, models = {}, max_growth = 10, v = 1,
    mets_def = [], stcut = 1e-8, title = "draft_cons", mets_to_plot = [],
    work_based_on = "id", manifest = "", comets_output = False):
        self.identifier = "name" if work_based_on.lower() in ["name", "names"] else "id"
        self.models = models
        self.v = v # volume is in liters
        self.max_growth = max_growth
        self.__media = {}
        self.media = self.set_media(media, concentration = True) # dict met.id : concentration
        self.mets_def = self.set_dMetabolites(mets_def) # dict of dMetabolite objects
        self.T = [0.]
        self.stopDFBA = (False, "not running yet!")
        self.stcut = stcut
        self.title = title
        self.mets_to_plot = mets_to_plot
        self.orgs_to_plot = ""
        self.manifest = manifest
        self.comets = comets_output
        return

    def __str__(self):
        echo = f"Consortium MMODES object with {self.v} volume and {len(self.models)} models: "
        for mod in self.models.values():
            echo += "\n"+mod.__str__()
        return echo

    def set_dMetabolites(self, mets):
        '''
        Initializes self.mets_def as dictionary of met.id: dMetabolite obj
        INPUT -> mets, list of dMetabolite objects
        OUTPUT -> dict
        '''
        if mets == []:
            return {}
        dict_mets = {}
        for met in mets:
            dict_mets[met.id] = met
        return dict_mets

    def cobrunion(self):
        '''
        Union of the models as a list
        OUPUT -> list of union of extracellular metabolites
        '''
        models = self.models.values()
        ex_mets = set()
        for model in models:
            for met in model.model.metabolites:
                if not met.compartment:
                    print(f"{met.id} with name {met.name} haven't got an specified compartment attribute.")
                elif met.compartment in ['e', 'e0', 'ExtraCellular', 'extracellular']:
                    ex_mets.add(met.__getattribute__(self.identifier))
        return ex_mets

    def add_model(self, mod_path, volume_0, solver = "glpk", method = "pfba", dMets = {}, limit = False):
        '''
        Adds a model to the consortium
        '''
        if type(limit) != bool:
            limit *= v # concentration to mass
        mod = dModel(mod_path, volume_0 * self.v, solver, method, dMets, work_based_on = self.identifier, limit = limit)
        if mod.model.id not in self.models:
            self.models[mod.model.id] = mod
        else:
            self.models[mod.model.id+"_"+str(len(self.models))] = mod
        # If metabolite already in another model, refresh value will be replaced
        self.mets_def.update(mod.dMets)
        self.media = self.set_media(self.media)
        return

    def add_mets(self, pert, concentration = False):
        '''
        Adds a perturbation to media
        INPUT -> pert, dict met.id : concentration/quantity
        '''
        v = self.v if concentration else 1
        for met in pert:
            if met in self.media:
                self.media[met] += float(pert[met])*v
            elif met in self.models:
                # pert[met] could be negative too
                self.models[met].add_biom(v*float(pert[met]))
            else:
                # TODO: use a custom warning class MediaWarning
                warnings.warn(f"\nMetabolite {met} wasn't added to media.")
        return

    def set_media(self, media, concentration = False):
        '''
        Set self.media, as extracellular metabolites in all models updated
        with initial concentrations
        INPUT -> media: dictionary, inital concentrations in media,
            concentration: bool, True in 1st call, to multiply per volume
        OUPUT -> media0, dict of all extracellular metabolites initialized
        '''
        # It's important to keep previous settings of media even if some metabolites
        # aren't added to the medium. If a new model is later appended, some of the
        # "incorrect" metabolites (not in current GEMs ex space) could be effective.
        self.__media.update(media)
        media0 = {}
        if self.models:
            v = self.v if concentration else 1
            ex_mets = self.cobrunion()
            for met in ex_mets:
                if met in self.__media:
                    media0[met] = float(self.__media[met])*v
                else:
                    media0[met] = 0
        return media0

    def medium_from_json(self, jpath, concentration = False):
        ''' Read medium from JSON file'''
        import json
        with open(jpath) as jdata:
            self.media = self.set_media(json.load(jdata), concentration)

    def get_manifest(self, manifest):
        if manifest:
            if not self.comets:
                return Manifest(models = self.models, media = self.media, fflux = manifest,
                fmedia = self.output, comets = False)
            else:
                return Manifest(models = self.models, media = self.media, fman = manifest, comets = True)
        else:
            return ""

    def mm(self, c, v, k):
        '''
        Michaelis Menten computation of lower bound. It assumes uptake is always
        computed as lower bound in exchange reactions.
        INPUTS -> floats: c (amount), v (concentration/time), k (concentration).
        OUTPUT -> float (amount/time)
        '''
        return (-c*v)/(c/self.v+k)

    def dinamicpFBA(self, t, log_texts = False):
        '''
        f in ODEsolver, computes dBM and dC
        INPUT -> t: timestep, just to compute refresment of media
               log_texts: ignore this parameter.
        OUPUT -> numpy array ([dBM/dt]+[dCj/dt]])
        '''
        dBMdt = {k: 0 for k in self.models}
        ran_mods_id = [k for k in dBMdt]
        random.shuffle(ran_mods_id) # to randomly iterate over models
        copy_media = dcp(self.media) # to transform
        dMedia = dcp(self.media) # to set dCdt
        for key in dMedia:
            dMedia[key]=0
        # 1) Add refresh concentrations
        if self.mets_def:
            for met in self.mets_def:
                if self.mets_def[met].refresh(t):
                    copy_media[met] += self.mets_def[met].concentration_fed
                    dMedia[met] = self.mets_def[met].concentration_fed
        for mID in ran_mods_id:
            # use model as a context, to restore possible changes in bounds later.
            # I guess it breaks backwards compatibility with cobra versions but it's lighter.
            org = self.models[mID]
            if not org.volume.q:
                # Biomass of organism = 0, don't waste resources here
                dBMdt[mID] = 0
                if self.manifest:
                    self.manifest.write_fluxes(org.model, self.T[-1])
                continue
            with org.model as mod:
                # 2) Updates lower bounds
                for met in copy_media:
                    if met not in org.exchanges:
                        continue
                    if copy_media[met] == 0:
                        mod.reactions.get_by_id(org.exchanges[met]).lower_bound = 0
                    elif met not in org.dMets:
                        # leave them as defined in cobra model
                        if copy_media[met] < abs(mod.reactions.get_by_id(org.exchanges[met]).lower_bound*org.volume.q):
                            mod.reactions.get_by_id(org.exchanges[met]).lower_bound = -copy_media[met]/org.volume.q
                    else:
                        # or dMetabolites has specific Vmax, Km
                        mod.reactions.get_by_id(org.exchanges[met]).lower_bound = self.mm(copy_media[met],org.dMets[met].Vmax, org.dMets[met].Km)
                # 3) Run pFBA and update volume (it will be properly updated later)
                try:
                    org.opt(mod)
                    if self.manifest:
                        self.manifest.write_fluxes(mod, self.T[-1])
                    if mod.reactions.get_by_id(org.volume.bm).flux > 0:
                        org.death_rate = mod
                        dBMdt[mID] = org.volume.q*mod.reactions.get_by_id(org.volume.bm).flux - org.death_rate*org.volume.q
                    else:
                        # wrong behaviour
                        dBMdt[mID] = 0
                    # 4) Updates media with fluxes (just for the next organisms bounds, media
                    # will be updated later in the ODE execution, self.update_true_ode)
                    for met, reac in org.exchanges.items():
                        copy_media[met] += mod.reactions.get_by_id(reac).flux*org.volume.q
                        dMedia[met] += mod.reactions.get_by_id(reac).flux*org.volume.q
                except AttributeError: # this workarounds infeasible solutions with gurobi
                    print(f"\nOrganism {mod.id} couldn't be optimized in time {self.T[-1]}. Setting growth to 0.\n")
                    if self.manifest:
                        #self.manifest.write_fluxes(load_model(org.path), self.T[-1]) # quite dirty...
                        self.manifest.write_fluxes(org_model, self.T[-1]) # TODO: see if works
                    dBMdt[mID] = 0

        # 5) returns ODE system solution[t]
        dQdt = [dBMdt[mod] for mod in sorted(dBMdt)] + [dMedia[met] for met in sorted(dMedia)]
        self.dQdt = dQdt # it will be used by other methods (is_stable, e.g.)
        return dQdt

    def update_true_ode(self, sol):
        '''
        Once an iteration of ODE solver is done, it updates media with true
        values from sol.
        INPUT -> sol, list of dQdt solutions in a time step of ODE solution
        '''
        s = list(sol)
        i = 0
        for mod in sorted(self.models):
            if s[0] > 0:
                self.models[mod].volume.q = s.pop(0)
            else:
                self.models[mod].volume.q = 0
                del(s[0])
            if self.models[mod].volume.q > self.max_growth:
                self.stopDFBA = (True, "\nMax growth reached at "+str(self.T[-1])+" hours")
            elif self.models[mod].stuck:
                self.stopDFBA = (True, "\npFBA infeasible at "+str(self.T[-1])+" hours")
        for met in sorted(self.media):
            if s[0] > 0:
                self.media[met] = s.pop(0)
            else:
                self.media[met] = 0
                del(s[0])
        return

    def get_conditions(self):
        ''' Get initial conditions to start the integrator'''
        t = self.T[-1]
        q = [self.models[mod].volume.q for mod in sorted(self.models)] + [self.media[met] for met in sorted(self.media)]
        return q,t

    def is_stable(self):
        '''
        Check if the slope is less than near-0 value
        '''
        biomassesdt = self.dQdt[0:len(self.models)] # it would be more correct to take solver output
        for b in biomassesdt:
            if b>self.stcut:
                break
        else:
            self.stopDFBA = (True, "\nStationary state has been reached.\n")
        return

    def run(self, maxT=10, integrator='vode', stepChoiceLevel=(0., 0.5, 1000.), verbose = False, outf = "plot.tsv", outp = "plot.png", plot = True, actualize_every = float("-inf")):
        '''
        Solves systems of ODEs while updating values.
        INPUTS -> maxT = maxtime condition
                integrator = type of solver of scipy.integrator.ode/odeint class
                stepChoiceLevel=(0, max T step, max N steps) or (0,endValue, nSteps)
                                depending on the 'integrator' parameter
                verbose = bool, to output stuff
                outf = path to tsv to write
                outp = path to png to write
                plot = bool, if plot should be generated (call to function)
                actualize_every = float, interval time when outputs are not written
        1) Def f as one of object methods
        2) while no convergence:
            3) advances one step of ODE solver
            4) updates values in object of biomasses and media concentrations
        '''
        def dqdt(t, y, mod, log_texts = False):
            '''
            Return the dQdt system of ODEs defined from the flux solutions
            '''
            return mod.dinamicpFBA(t, log_texts) # "mod" will be "self"

        self.stopDFBA = (False, "\nRunning...")
        self.output = outf
        self.outplot = outp
        iotimes = 0
        curr_act = actualize_every
        # 1. Set parameters of solver
        integratorSet = False
        nMaxSteps = stepChoiceLevel[2]
        if integrator.upper() == "FEA":
            solver = _fsolvers.FEA(f = dqdt, dt = stepChoiceLevel[0], mod = self)
            q0, t0 = self.get_conditions()
            solver.set_initial_value(q0, t0)
            step = 0
        elif integrator.upper() in ["RK", "RUNGEKUTTA", "RK4"]:
            solver = _fsolvers.RungeKutta4(f = dqdt, dt = stepChoiceLevel[0], mod = self)
            q0, t0 = self.get_conditions()
            solver.set_initial_value(q0, t0)
            step = 0
        else:
            # as in DAPHNE, https://github.com/QTB-HHU/daphne_ecoli-diauxie (Succurro et al., 2018)
            solver = ode(dqdt).set_integrator(integrator)
            if integrator in ['dopri5', 'lsoda']:
                nMaxSteps -= 1
                # In this case: stepChoiceLevel=(0,endValue, nSteps)
                grid_t = np.linspace(stepChoiceLevel[0], stepChoiceLevel[1], stepChoiceLevel[2])
                grid_dt = grid_t[1] - grid_t[0]
                solver.set_integrator(integrator, nsteps=1, max_step=grid_dt)
                integratorSet = True
            else:
                # maybe vode should be the only option.
                solver.set_integrator(integrator, min_step=stepChoiceLevel[0], max_step=stepChoiceLevel[1], method = 'bdf', order = 5)
                integratorSet = True
            if integratorSet:
                # set the parameters of the differential function dqdt: model and verbosity
                solver.set_f_params(self, verbose)
                # suppress Fortran-printed warning
                solver._integrator.iwork[2] = -1
                warnings.filterwarnings("ignore", category=UserWarning)
                q0, t0 = self.get_conditions()
                solver.set_initial_value(q0, t0)
                step = 0
                substeps = 0
                eventIdx = 1
            else:
                raise NotIntegratorError("ODE Parameters weren't properly supplied.")
        if verbose:
            bar = ProBar(maxT)
        manifest = ""
        if self.manifest and isinstance(self.manifest, str): # self.manifest is a string to be intitialized
            manifest = self.get_manifest(self.manifest)
        elif self.manifest: # Consortium has already been run with manifest
            manifest = dcp(self.manifest)
        # 2. ODE solver loop.
        self.manifest = False
        while not self.stopDFBA[0] and self.T[-1] < maxT and step < nMaxSteps:
            # 2.1. Advances in solver
            step += 1
            write = self.T[-1] >= iotimes*curr_act or iotimes == 0
            if write:
                self.manifest = manifest
                iotimes += 1
            else:
                self.manifest = False
            solver.integrate(maxT, step=True)
            # 2.2. Update object consortium parameters
            self.update_true_ode(solver.y)
            # write media before appending time
            if write: self.write_plot_tsv()
            self.T.append(solver.t)
            # 2.3. Check if stationary phase, else continue
            self.is_stable()
            if verbose:
                bar.progress(self.T[-1])
        self.manifest = manifest
        # 3. Plot and final message, to check why it was finished
        if plot:
            plot_comm(self)
        if verbose:
            if self.stopDFBA[0]:
                print(self.stopDFBA[1])
            else:
                print("\nMax time or steps reached.\n")
        return

    def write_plot_tsv(self):
        '''
        Records the output of the solver as a TSV
        INPUTS -> out_file: path to output file
        '''
        if not os.path.isfile(self.output):
            # write header
            with open(self.output, "w") as f:
                line1 = ""
                i = 1
                for mod in sorted(self.models):
                    line1 += mod+"\t" if mod != "" else "biomass"+str(i)+"\t"
                    i += 1
                self.orgs_to_plot = line1.split(sep = "\t")[:-1]
                for met in sorted(self.media):
                    line1 += met+"\t"
                f.write("time"+"\t"+line1[:-1]+"\n")
        with open(self.output, "a") as f:
            line = ""
            for mod in sorted(self.models):
                line += str(self.models[mod].volume.q)+"\t"
            for met in sorted(self.media):
                line += str(self.media[met])+"\t"
            f.write(str(self.T[-1])+"\t"+line[:-1]+"\n")
        if self.manifest:
            self.manifest.write_media()
            self.manifest.write_biomass()
