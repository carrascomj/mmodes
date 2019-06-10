#!/usr/bin/python3

# Experiment class, subclass of Consortium.
# Inputs are mainly supplied as files, perturbations and medium are
# read from a JSON file.

# @author: Jorge Carrasco Muriel
# Created: 09/06/2019

import json, os
# import datatable as dt
# from datatable import f
import numpy as np
from .core import Consortium
from mmodes.vis import plot_comm
import warnings
import random

# GLOBAL VARIABLES
working_set = 0 # to then check datatable installed packages
dt = 0 # datatable package
f = 0 # FrameProxy object of datatable package
pd = 0 # pandas package
re = 0 # re built-in module
# datatable module is a recent h2o module coming from the R package. It improves
# perfomance, but, since it's not that common as pandas, It's preferable to avoid
# putting it as a requirement when it's not used by the main modules

class NoModelsException(Exception):
    '''
    NoModelsException is raised when no models could be added to the Experiment
    object or models_dir wasn't supplied
    '''
    pass

class NoMedium(Exception):
    '''
    NoMedium is raised when no medium could be added to the Experiment
    object or models_dir wasn't supplied
    '''
    pass

class Experiment(Consortium):
    '''
    Experiment subclass of Consortium. In order to ease the configuration of experiments
    from files, the Experiment class implements four keypoints:
        1. Models are loaded from a directory path, supplied as input
        2. Medium is loaded from a path, as a JSON:
          - JSON is a list of dictionaries.
          - Dictionaries contains metabolites as keys and amounts/concentrations
              as values.
          - First dictionary is medium. The rest are perturbations.
          - Names of perturbations are input
       3. Filtering of the output.
       4. Save output.
    '''
    def __init__(self, medium_path = "", models_dir = "", rand_biomasses = [0.0001,0.0005], perturbations = [], max_growth = 10, v = 1,
    mets_def = [], stcut = 1e-8, title = "draft_cons", mets_to_plot = [],
    work_based_on = "id", manifest = "", comets_output = False, lp = "fba", solver = "glpk"):
        self.identifier = "name" if work_based_on.lower() in ["name", "names"] else "id"
        self.v = v # volume is in liters
        self.models = {}
        self.media = {}
        self._Consortium__media = {}
        self.mets_def = self.set_dMetabolites(mets_def) # dict of dMetabolite objects
        self.manifest = manifest
        self.comets = comets_output
        self.max_growth = max_growth
        self.stopDFBA = (False, "not running yet!")
        self.mets_to_plot = mets_to_plot
        self.stcut = stcut
        self.title = title if title else 'Some_Experiment'
        self.orgs_to_plot = ""
        self.T = [0.]
        self.path_to_models = models_dir
        self.read_models(models_dir, rand_biomasses, lp, solver)
        self.pers_state = self.read_exp_media(medium_path, perturbations) # list of dict of MEDIA: dict + PERTURBATION: name + time = ""
        self.tsv_filter = None
        return

    def __str__(self):
        echo = f"Consortium MMODES object with {self.v} volume and {len(self.models)} models in {self.path_to_models}: "
        for mod in self.models.values():
            echo += "\n"+mod.__str__()
        pers = ', '.join([per["PERTURBATION"] for per in self.pers_state])
        echo += "\nPERTURBATIONS: " + pers + "\n\n"
        return echo

    def read_models(self, models_dir, rand_biomasses, lp, solver):
        if models_dir == "":
            raise NoModelsException('models_dir parameter is required!')
        self.models_dir = models_dir if models_dir.endswith('/') else models_dir + '/'
        b_rand = lambda: random.uniform(*rand_biomasses)
        for modf in os.listdir(self.models_dir):
            if (modf.endswith("xml") or modf.endswith("mat") or modf.endswith("json")) and os.path.isfile(self.models_dir+modf):
                try:
                    self.add_model(self.models_dir+modf, b_rand(), solver = solver, method = lp)
                except:
                    warnings.warn(f'{modf} in {self.models_dir} is not a model and will be ignored')
        if not self.models:
            raise NoModelsException(f'No model could be added to the community from directory {models_dir}')
        return

    def read_exp_media(self, medium_path, perturbations):
        if not medium_path:
            raise NoMedium('medium_path parameter is required!')
        try:
                with open(medium_path) as jdata:
                    gen_media = json.load(jdata)
        except json.decoder.JSONDecodeError:
            raise NoMedium(f'{medium_path} is not correctly formatted as JSON')
        nm = len(gen_media)
        np = len(perturbations)
        if nm != np and nm < np:
            warnings.warn(f"\nNumber of perturbations names > perturbations! Please, check your input")
        elif nm > np:
            warnings.warn(f"\nNumber of perturbations names < perturbations! Filling with generic names...")
            for i in range(nm-np):
                perturbations.append(f'generic_perturbation{i}')
        out_media = [{} for _ in gen_media]
        for g in range(len(gen_media)):
            # restructure the input to a dictionary of tree keys
            out_media[g] = {}
            out_media[g]["MEDIA"] = gen_media[g] # original input
            out_media[g]["PERTURBATION"] = perturbations[g] # name, identifying perturbation
            out_media[g]['time'] = 0 # to save time at when perturbation was added
        return out_media

    def _tsv_filter_dt(self, inplace = False, equif = True):
        '''
        Function that filters medium and fluxes TSV files based on perturbation times.
                inplace: bool, whether overwrite input paths (default False);
                v: float, volume magnitude to obtain medium concentrations;
                equif: bool, whether write an additional fluxes filtered file,
                        with 100 equidistant points (default True)
        OUTPUT -> it returns None, writes 2(3) TSV files
        '''
        def equidistant(df, n):
            sample = np.linspace(df.nrows-1,1,n).astype('int')
            sample.sort()
            return df[sample, :]

        dfs = []
        if not self.output:
            warnings.warn("In filtering: Medium parameter wasn't supplied, it won't be generated.")
        else:
            dfs.append([dt.fread(self.output), self.output, 0])
            if self.v != 0:
                for i in range(1,dfs[0][0].ncols-1):
                    dfs[0][0][:,i] = dfs[0][0][:,f[i]/self.v]
        if not self.manifest:
            warnings.warn("In filtering : Fluxes parameter wasn't supplied, it won't be generated.")
        else:
            dfs.append([dt.fread(self.manifest), self.manifest, 1])

        for log, path, n in dfs:
            log = log[f.Perturbations != "FALSE", :]
            if n != 0 and equif:
                log_equif = equidistant(log,100) # take 100 equidistant rows
                log_equif.to_csv("equi_" + path)
                del(log_equif)
                # TODO: I don't know how to implement a condroll with datatable
                # We aren't currentyly using it, anyway
            if inplace:
                log.to_csv(path)
            else:
                log.to_csv("filtered_" + path)
        return

    def _tsv_filter_pd(self, inplace = False, equif = True):
        '''
        Function that filters medium and fluxes TSV files based on perturbation times.
                inplace: bool, whether overwrite input paths (default False);
                v: float, volume magnitude to obtain medium concentrations;
                equif: bool, whether write an additional fluxes filtered file,
                        with 100 equidistant points (default True)
        OUTPUT -> it returns None, writes 2(3) TSV files
        '''
        def equidistant(df,n):
            sample = np.linspace(df.shape[0]-1,1,n).astype('int')
            sample.sort()
            return df.iloc[sample].reset_index(drop = True)

        def condroll(df, cond, n = 1, forward = True, extend_sel_fluxes = 1):
            '''
            Function that extends rows selected by condition n positions.
            INPUTS -> df, pandas DataFrame object with condition applied to rows;
                      cond, bool selection of dataframe rows;
                      n, number of rows to extend;
                      forward (bool), extend rows selection forward, backwards or 'both'.
            OUTPUT -> pandas DataFrame
            '''
            n += 1
            inx = df[cond].index.repeat(n)
            ninx = np.empty(0)
            if not forward:
                inx = inx - (n-1)
            nrow = df.shape[0]
            for i in range(0,len(inx),n):
                for j in range(n):
                    if 0 <= inx[i]+j < nrow: # row index isn't out of bounds
                        ninx = np.append(ninx, inx[i] + j)
            if forward in ['both', 'b']:
                # natural join, sort by time and reset index
                return pd.merge(df.iloc[ninx,:],
                                condroll(df, cond, n = n-1, forward = False),
                                how = "outer").sort_values(by = "time")#.reset_index(drop = True)
            else:
                return df.iloc[ninx,:]

        dfs = []
        if not self.output:
            warnings.warn("In filtering: Medium parameter wasn't supplied, it won't be generated.")
        else:
            dfs.append([pd.read_csv(self.output), self.output, 0])
            if self.v != 0:
                for i in range(1,dfs[0][0].shape[1]):
                    dfs[0][0].iloc[:,1:-1] = dfs[0][0].iloc[:,1:]/v
        if not self.manifest:
            warnings.warn("In filtering: Fluxes parameter wasn't supplied, it won't be generated.")
        else:
            dfs.append([pd.read_csv(self.manifest), self.manifest, 1])

        for log, path, n in dfs:
            if n == 0:
                log = log[log['Perturbations'] != "FALSE"].dropna(axis = 1, how = 'all').reset_index(drop = True)
            else:
                if equif:
                    log_equif = equidistant(log,100) # take 100 equidistant rows
                    log_equif.to_csv("equi_" + path, sep = "\t", index = False)
                    del(log_equif)
                log = condroll(log, log['Perturbations'] != "FALSE",
                    n = extend_sel_fluxes, forward = 'both').dropna(axis = 1, how = 'all').reset_index(drop = True)
            if inplace:
                log.to_csv(path, sep = "\t", index = False)
            else:
                log.to_csv("filtered_" + path, sep = "\t", index = False)
        return

    def log(self):
        '''
        Writes information of consortium object to file
        '''
        global re
        if not re:
            import re
        logf = '../simulations.txt'
        p = re.compile(r'#+ SIMULATION (\d+) #+')
        if os.path.isfile(logf): # parse last simulation number
            with open(logf) as l:
                for line in l.readlines():
                    num_sim = p.search(line)
                    if num_sim:
                        head = " SIMULATION "+str(int(num_sim.group(1))+1)+" "
        else:
            head = " SIMULATION 1 "
        lines = '{:{fill}{align}{width}}'.format(head,
            fill = '#',
            align = '^',
            width = 30) + "\n"
        lines += self.__str__()
        with open(logf, "a") as l:
            l.write(lines)
        return

    def run_experiment(self, integrator = 'vode', verbose = True, stepChoiceLevel = (),
        outp = "Some_Experiment.png", intervl = 10, filter = False, equif = True,
        inplace_filter = False, plot = True, actualize_every = float('-inf')):
        '''
        Main function that runs sequentally dFBA's, one per perturbation.
        '''
        it = 1
        print(intervl)
        for i in range(len(self.pers_state)):
            if it == 1:
                # 1) instantiate media
                self.media = self.set_media(self.pers_state[i]["MEDIA"], True)
            else:
                for k in self.pers_state[i]["MEDIA"]: # not needed at all...
                    if self.pers_state[i]["MEDIA"][k] == 0:
                        self.media[k] = 0
                self.add_mets(self.pers_state[i]["MEDIA"], True)
            # 2) run it
            self.pers_state[i]["time"] = self.T[-1] # to then filter output by perturbations
            if verbose:
                print(f'\n\033[1;32;40m Iteration {it}: Loading {self.pers_state[i]["PERTURBATION"]}!\033[0m')
            self.run(verbose=verbose, plot=False, maxT = intervl+self.T[-1], integrator = "FEA",
            stepChoiceLevel=stepChoiceLevel, outp = outp, actualize_every = actualize_every)
            it += 1
            if i < len(self.pers_state) - 1:
                self.rewrite_with_per(self.pers_state[i+1]["PERTURBATION"])
            else:
                self.rewrite_with_per("END")
        # 3) processing
        if plot: plot_comm(self)
        if filter:
            self.init_filter()
            self.tsv_filter(equif = equif, inplace = inplace_filter)
        self.log()
        return

    def init_filter(self):
        '''
        Check if datatable is installed. Else, uses pandas to manage the filtering.
        '''
        global working_set, pd, dt, f
        if not dt and not pd:
            if not working_set:
                from pkg_resources import working_set
                if 'datatable' in [d.project_name for d in working_set]:
                    import datatable as dt
                    from datatable import f
                    self.tsv_filter = self._tsv_filter_dt
                else:
                    import pandas as pd
                    self.tsv_filter = self._tsv_filter_pd
        else:
            if dt:
                self.tsv_filter = self._tsv_filter_dt
            else:
                self.tsv_filter = self._tsv_filter_pd
        return

    def write_plot_tsv(self):
        '''
        Overwrites Consortium function to write a Perturbations columns.
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
                f.write("time\t"+line1[:-1]+"\tPerturbations"+"\n")
        with open(self.output, "a") as f:
            line = ""
            for mod in sorted(self.models):
                line += str(self.models[mod].volume.q)+"\t"
            for met in sorted(self.media):
                line += str(self.media[met])+"\t"
            if self.T[-1] == 0:
                f.write(str(self.T[-1])+"\t"+line[:-1]+f'\t{ self.pers_state[0]["PERTURBATION"] }\n')
            else:
                f.write(str(self.T[-1])+"\t"+line[:-1]+f"\tFALSE\n")
        if self.manifest:
            self.manifest.write_media()
            self.manifest.write_biomass()
        return

    def rewrite_with_per(self, per):
        '''
        Replace last FALSE in Perturbations' column with per
        '''
        with open(self.output) as old:
            with open('tmp.mmodes', 'w') as tmp:
                last_line = False
                for line in old.readlines():
                    if last_line:
                        tmp.write(last_line)
                    last_line = line
                tmp.write(last_line.replace('FALSE', per))
        os.rename('tmp.mmodes', self.output)
        return
