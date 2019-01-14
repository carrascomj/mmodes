#!/usr/bin/env python3

# All the classes, modules required and exceptions of the consortium package.
# The package provides an straight-forward framework to run dynamic simulations
# of metabolic models consortia using coprapy and scipy.

# Author: Jorge Carrasco Muriel
# e-mail: jorge.cmuriel@alumnos.upm.es
# Date of 1st version: 30/12/2018
# Date of last update: 14/01/20189

import sys
import os
import re
import cobra
import random
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import ode, odeint # "integrate" module is not importable!
from copy import deepcopy as dcp
from mmodes.io import Manifest


class NoBiomassException(Exception):
    """
    NoBiomassException is raised when a dynamicModel.reactions attribute does not
    contain 'biomass' as key
    Biomass is needed to follow growth!
    """
    pass
class NotIntegratorError(Exception):
    """ When integrator is False... check Consortium.run parameters"""
    pass

class InfeasibleSolution(Exception):
    """ If method (pFBA or FBA) reports Infeasible when dModel is instantiated"""
    pass

class Volume():
    def __init__(self, volume_0, model):
        self.bm = self.get_bm_reaction(model) # biomass reaction cobra object
        self.q = volume_0 # biomass concentration
        return

    def get_bm_reaction(self, model):
        """ Look for the biomass reaction
        INPUT: model
        OUTPUT: biomass rection cobra object """
        # There are some strange models that have a biomas exchange reaction,
        # so it doesn't properly works sometimess.
        for reac in model.reactions:
            if reac.id.lower().find("biomass") != -1:
                bm = reac.id
                break
        else:
            raise NoBiomassException("No biomass reaction found in "+str(model.id))
        return bm

    def dupdate(self):
        """ Update value of biomass concentration
        INPUT: volume """
        self.q += self.bm.value*self.q
        return self.bm.value*self.q

class dModel():
    def __init__(self, mod_path, volume_0, solver = "glpk", method = "pfba", dMets = {}):
        self.path = mod_path
        self.model = self.load_model() # cobra model
        self.volume = Volume(volume_0, self.model) # volume object
        self.exchanges = self.get_exchange_reactions() # dict met:exchange reactio
        self.dMets = dMets
        self.model.solver = solver
        self.stuck = False
        self.method = method
        self.free_media()
        self.check_feasible()
        return

    def load_model(self):
        """ Loads a matlab or sbml model
        INPUT: path
        OUTPUT: cobra model object """
        mod = self.path
        try:
            model_obj = cobra.io.read_sbml_model(mod)
        except:
            model_obj = cobra.io.load_matlab_model(mod)
        return model_obj

    def update(self):
        """ Updates value of biomass concentration
        INPUT: volume """
        return self.volume.dupdate()

    def get_exchange_reactions(self):
        """ model.exchanges cobra method isn't properly working. This function
        provides a fix as a dictionary.
        OUPUT -> dictionary, metabolite: exchange reactions"""
        dic_out = {}
        for reac in self.model.reactions:
            if reac.id.find("EX_") != -1:
                for met in reac.metabolites:
                    if met.id.find("biomass") == -1: #strange metabolite in some models
                        dic_out[met.id] = reac.id
        return dic_out

    def free_media(self):
        """ Eliminates boundaries associated with cobra.model.medium object. The
        simulation fix bounds, not medium."""
        for reac in self.exchanges.values():
            # I don't know if all exchanges reactions are included in medium.
            self.model.medium[reac] = 1000
        return

    def opt(self, mod = ""):
        """ Solves FBA or pFBA
        INPUT: mod, cobra model, so it's able to work as context """
        if mod == "":
            mod = self.model
        if self.method == "pfba":
            try:
                cobra.flux_analysis.pfba(mod)
            except:
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
                except:
                    raise InfeasibleSolution("pFBA is infeasiable for model "+mod.id+". You may want to change the method param to 'fba'.")
            elif self.method == "pfba":
                try:
                    mod.optimize()
                except:
                    raise InfeasibleSolution("FBA is infeasiable for model "+mod.id+". You may want to check your model.")
    def add_dMet(self, met):
        """ Adds a dMetabolite to self.dMets
        INPUT -> met, dMetabolite object"""
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
        """ Computes times when the metabolite is added
        INPUT -> t, iteration time of the solver"""
        # TODO: since iteration won't always work with fixed step time, I guess
        # it's necessary to change each_time to an interval ????
        ref = True
        if self.each_time == 0: # avoid ZeroDivisionError
            ref = False
        elif t % self.each_time != 0:
            ref = False
        return ref

class Consortium():
    def __init__(self, media = {}, models = {}, max_growth = 10, death_rate = 0, v = 1, timeStep = 0.1, mets_def = [], defaultKm = 0.01, defaultVmax = 20, stcut = 1e-4, title = "draft_cons", mets_to_plot = [], manifest = ""):
        self.models = models
        self.media = self.set_media(media, concentration = True) # dict met.id : concentration
        self.max_growth = max_growth
        self.death_rate = death_rate
        self.mets_def = self.set_dMetabolites(mets_def) # dict of dMetabolite objects
        self.Km = defaultKm
        self.Vmax = defaultVmax
        self.T = [0.]
        self.v = v # volume is in liters
        self.stopDFBA = (False, "not running yet!")
        self.stcut = stcut
        self.title = title
        self.mets_to_plot = mets_to_plot
        self.orgs_to_plot = ""
        self.manifest = manifest
        return

    def set_dMetabolites(self, mets):
        """ Initializes self.mets_def as dictionary of met.id: dMetabolite obj
        INPUT -> mets, list of dMetabolite objects
        OUTPUT -> dict """
        if mets == []:
            return {}
        dict_mets = {}
        for met in mets:
            dict_mets[met.id] = met
        return dict_mets

    def cobrunion(self):
        """union of the models as a list
        OUPUT -> list of union of extracellular metabolites"""
        models = self.models.values()
        ex_mets = set()
        for model in models:
            for met in model.model.metabolites:
                if met.id.find("[e]") != -1:
                    ex_mets.add(met.id)
        return ex_mets

    def add_model(self, mod_path, volume_0, solver = "glpk", method = "pfba", dMets = {}):
        """ Adds a model to the consortium """
        mod = dModel(mod_path, volume_0 * self.v, solver, method, dMets)
        self.models[mod.model.id] = mod
        # If metabolite already in another model, refresh value will be replaced
        self.mets_def.update(mod.dMets)
        self.media = self.set_media(self.media)
        return

    def add_mets(self, pert, concentration = False):
        """ Adds a perturbation to media
        INPUT -> pert, dict met.id : concentration"""
        # TODO: maybe it'd be nice to do it as in set_media throwing a warning
        # if metabolite isn't in media (and won't be added to it).
        if concentration:
            pert = {k: float(pert[k])*self.v+self.media[k] for k in pert}
        else:
            pert = {k: float(pert[k])+self.media[k] for k in pert}
        self.media.update(pert)
        return

    def set_media(self, media, concentration = False):
        """ Set self.media, as extracellular metabolites in all models updated
        with initial concentrations
        INPUT -> media, dict of inital concentrations in media,
            concentration -> bool, True in 1st call, to multiply per volume
        OUPUT -> media0, dict of all extracellular metabolites initialized """
        # TODO: work on units
        media0 = {}
        if self.models: # empty sequences are false...
            ex_mets = self.cobrunion()
            for met in ex_mets:
                if met in media:
                    media0[met] = float(media[met])
                    if concentration:
                        media0[met] = media0[met]*self.v
                else:
                    media0[met] = 0
        return media0

    def get_manifest(self, manifest):
        if manifest:
            return Manifest(models = self.models, media = self.media, fman = manifest)
        else:
            return ""

    def mm(self, c, v=None, k=None):
        """ Michaellis menten computation of lower bound """
        if v == None:
            v = self.Vmax
        if k == None:
            k = self.Km
        return (-c*v)/(k*c)

    def dinamicpFBA(self, t, log_texts = False):
        """ f in ODEsolver, computes dBM and dC
        INPUT -> t: timestep, just to compute refresment of media
               dC
        OUPUT -> numpy array ([dBM/dt]+[dCj/dt]])"""
        # TODO: should write a COMETS-like manifest with fluxes, etc.
        # The main reason to preserve format is to run in VisAnt, therefore it
        # isn't urgent at all.
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
            # use model as a context, to restore possible changes in bounds later
            # I guess it breaks backwards compatibility with cobra versions but it's lighter.
            org = self.models[mID]
            with org.model as mod:
                # 2) Updates lower bounds
                for met in copy_media:
                    if met not in org.exchanges:
                        continue
                    if copy_media[met] == 0:
                        mod.reactions.get_by_id(org.exchanges[met]).lower_bound = 0
                    elif met not in org.dMets:
                        # set as default
                        # or let it be as defined in cobra model
                        if copy_media[met] < mod.reactions.get_by_id(org.exchanges[met]).lower_bound:
                            mod.reactions.get_by_id(org.exchanges[met]).lower_bound = copy_media[met]
                    else:
                        # dMetabolites has specific Vmax, Km
                        mod.reactions.get_by_id(org.exchanges[met]).lower_bound = self.mm(copy_media[met],org.dMets[met].Vmax, org.dMets[met].Km)
                # 3) Run pFBA and update volume (it will be updated properly later)
                org.opt(mod)
                if self.manifest:
                    self.manifest.write_fluxes(mod, self.T[-1])
                dBMdt[mID] = org.volume.q*mod.reactions.get_by_id(org.volume.bm).flux - self.death_rate*org.volume.q
                # 4) Updates media with fluxes (just for the next organisms bounds, media
                # will be updated later in the ODE execution, self.update_true_ode)
                for met, reac in org.exchanges.items():
                    copy_media[met] += mod.reactions.get_by_id(reac).flux*org.volume.q
                    dMedia[met] += mod.reactions.get_by_id(reac).flux*org.volume.q
        # 5) returns ODE system solution[t]
        dQdt = [dBMdt[mod] for mod in sorted(dBMdt)] + [dMedia[met] for met in sorted(dMedia)]
        self.dQdt = dQdt # it will be used by other methods (is_stable, e.g.)
        return dQdt

    def update_true_ode(self, sol):
        """ Once an iteration of ODE solver is done, it updates media with true
        values from u_new
        INPUT -> sol, list of dQdt solutions in a time step of ODE solution"""
        s = list(sol)
        i = 0
        for mod in sorted(self.models):
            if s[0] > 0:
                self.models[mod].volume.q = s.pop(0)
            else:
                self.models[mod].volume.q = 0
                del(s[0])
            if self.models[mod].volume.q > self.max_growth:
                self.stopDFBA = (True, "\nMax growth reached at "+str(self.T[-1])+"hours")
            elif self.models[mod].stuck:
                self.stopDFBA = (True, "\npFBA infeasible at "+str(self.T[-1])+"hours")
        for met in sorted(self.media):
            if s[0] > 0:
                self.media[met] = s.pop(0)
            else:
                self.media[met] = 0
                del(s[0])
        if self.manifest:
            self.manifest.write_media()
            self.manifest.write_biomass()
        return

    def get_conditions(self):
        t = self.T[-1]
        q = [self.models[mod].volume.q for mod in sorted(self.models)] + [self.media[met] for met in sorted(self.media)]
        return q,t

    def is_stable(self):
        """ Check if the slope is less than near-0 value """
        biomassesdt = self.dQdt[0:len(self.models)] # it would be more correct to take solver output
        for b in biomassesdt:
            if b>self.stcut:
                break
        else:
            self.stopDFBA = (True, "\nStationary state has been reached.\n")
        return

    def run(self, maxT=10, integrator='dopri5', stepChoiceLevel=(0., 10., 50.), verbose = False, outf = "plot.tsv", outp = "plot.png", plot = True):
        """Solves systems of ODEs while updating values.
        INPUTS -> maxT = maxtime condition
                integrator = type of solver of scipy.integrator.ode/odeint class
                stepChoiceLevel=(0, max T step, max N steps) or (0,endValue, nSteps)
                                depending on the 'integrator' parameter
                verbose = bool, to output stuff
        1) Def f as one of object methods
        2) while no convergence:
            3) advances one step of ODE solver
            4) updates values in object of biomasses and media concentrations"""
        def dqdt(t, y, mod, log_texts = False):
            '''
            Return the dQdt system of ODEs defined from the flux solutions
            '''
            return mod.dinamicpFBA(t, log_texts) # "mod" will be "self"

        self.stopDFBA = (False, "\nRunning...")
        self.output = outf
        self.outplot = outp
        # 1. Set parameters of solver
        integratorSet = False
        nMaxSteps = stepChoiceLevel[2]
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
            solver.set_integrator(integrator, min_step=stepChoiceLevel[0], max_step=stepChoiceLevel[1], method = 'bdf', order = 5)
            integratorSet = True
        if integratorSet:
            #set the parameters of the differential function dqdt: model and verbosity
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
        if self.manifest and isinstance(self.manifest, str):
            self.manifest = self.get_manifest(self.manifest)

        # 2. ODE solver loop.
        while not self.stopDFBA[0] and self.T[-1]<maxT and step < nMaxSteps:
            # 2.1. Advances in solver
            step += 1
            solver.integrate(maxT, step=True)
            # 2.2. Update object consortium parameters
            self.T.append(solver.t)
            self.update_true_ode(solver.y)
            # 2.3. Check if stationary phase, else continue
            self.is_stable()
            self.write_plot_tsv()
            if verbose:
                bar.progress(self.T[-1])
        # 3. Final message. Check why it was finished
        if plot:
            self.plot_comm()
        if verbose:
            if self.stopDFBA[0]:
                print(self.stopDFBA[1])
            else:
                print("\nMax time or steps reached.\n")
            # We're pretty polite.
            #print("Thank you for using my script!")
        return

    def write_plot_tsv(self):
        """ Records the output of the solver as a TSV
        INPUTS -> out_file: path to output file """
        if not os.path.isfile(self.output):
            # write header
            with open(self.output, "w") as f:
                line1 = ""
                i = 1
                for mod in sorted(self.models):
                    line1 += mod+"\t" if mod != "" else "biomass"+i+"\t"
                    i += 1
                self.orgs_to_plot = line1.split(sep = "\t")[:-1]
                for met in sorted(self.media):
                    line1 += met+"\t"
                f.write("time"+"\t"+line1+"\n")
        with open(self.output, "a") as f:
            line = ""
            for mod in sorted(self.models):
                line += str(self.models[mod].volume.q)+"\t"
            for met in sorted(self.media):
                line += str(self.media[met])+"\t"
            f.write(str(self.T[-1])+"\t"+line+"\n")

    def plot_comm(self):
        '''Plots the concentration of given microorganisms and metabolites
        INPUTS -> path: string, path to the file to be used as input
                output: string, path where the plot will be generated
                organisms: list of strings to label microorganisms
                metabolites: list of strings to label metabolites
        OUTPUTS -> returns nothing, generates a plot in "outfile'''
        title = self.title
        path = self.output
        output = self.outplot
        # some colors to plot...
        colors_hex_2 = "#ff8c00", "#a020f0", "#00ff00", "#ffff00", "#7a8b8b", "#cd5555", "#1e90ff", "#787878", "#ff7f50", "#000000"
        # or nominal data color scheme found in http://geog.uoregon.edu/datagraphics/color/Cat_12.txt
        colors_hex = "#ff7f00", "#32ff00", "#19b2ff", "#654cff", "#e51932", "#000000", "#ffff32", "#ff99bf", "#ccbfff", "#a5edff", "#b2ff8c", "#ffff99", "#ffbf7f"
        to_plot=pd.read_csv(path, sep='\t', header = 0)
        for col in to_plot:
            if col not in self.mets_to_plot and col not in self.orgs_to_plot and col != "time":
                del to_plot[col]
        to_plot = to_plot / self.v # return concentrations
        fig, ax1 = plt.subplots()
        t = to_plot.time
        i = 0
        for col in to_plot:
            if 0 < i < len(self.orgs_to_plot) + 1:
                s1 = to_plot[col]
                ax1.plot(t, s1, colors_hex[i], alpha = 0.8)
            elif i >= len(self.orgs_to_plot) + 1:
                if i == len(self.orgs_to_plot) + 1:
                    ax1.set_xlabel("time(h)")
                    ax1.set_ylabel('organisms (g/L)')
                    ax1.tick_params('y')
                    ax2 = ax1.twinx()
                s2 = to_plot[col]
                ax2.plot(t, s2, colors_hex[i])
                # Maybe I should put the axis on black...
                ax2.set_ylabel('metabolites (mmol/L)')
                ax2.tick_params('y')
            i += 1
        plt.title(title)
        ax1.legend(loc=6)
        ax2.legend(loc=2)
        fig.tight_layout()
        plt.savefig(output)
        return

class ProBar():
    """ Just a simple object progress bar to output on screen """
    def __init__(self, n, blength = 20):
        self.n = n
        self.blength = 20
        self.progress(0)

    def progress(self, pro):
        pro = pro/self.n
        block = int(round(self.blength*pro))
        bar = "\rRunning... [" + "#"*block + "-"*(self.blength-block) + "] " + '{:.2f}'.format(pro*100) + "%"
        if pro == 1:
            bar += " DONE!\n"
        # I don't really know if this properly works on every IDE.
        sys.stdout.write(bar)
        sys.stdout.flush()