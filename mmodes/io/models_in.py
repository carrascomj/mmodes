#!/usr/bin/env python3

# I/O utilities to manage and translate models to one 'lenguage'.
#   0. Utility functions:
#       0.1. find_models
#       0.2. load_model
#   change_model_prefix
#   1. Tranlation worflow extracellular ids to BiGG
#   2. Translation workflow extracellular names/ids from pool (model based)

# Author: Jorge Carrasco Muriel
# e-mail: jorge.cmuriel@alumnos.upm.es
# Date of 1st version: 30/12/2018

import cobra
import os, re, sys
from warnings import filterwarnings as wfilt
# global variables that refers to not yet loaded packages
requests = 0
np = 0
cobrunion = 0
csv = 0
pickle = 0

############# UTILY FUNCTIONS USED BY OTHER MODULES IN THE PACKAGE #############
def find_models(dir_models, just_path = False):
    '''
    Find models in "dir_models" path
    INPUTS -> dir_models: path to models directory (string)
            just_path: bool, True if just outputs a string list with the paths.
    OUPUT -> list of model objects
    '''
    path_models = []
    for file in os.listdir(dir_models):
        # Find all sbml or matlab models and exclude comets models
        if (file.find("xml") != -1 or file.find("mat") != -1) and file.find("xml.cmt") == -1:
            path_models.append(dir_models+file)
    if just_path:
        return path_models
    models = []
    strain_number = 1
    for mod in path_models:
        models.append(load_model(mod))
        # model.id attribute will be used to name files
        if models[-1].id == "":
            models[-1].id = "strain_" + str(strain_number)
        strain_number += 1
    return models

def load_model(model_path):
    '''
    Function that loads models
    '''
    if hasattr(model_path, "metabolites"):
        # A cobra model was passed as argument
        model = model_path
    else:
        wfilt("ignore", category=UserWarning) # although sbml warnings aren't warnings...
        try:
            model = cobra.io.load_json_model(model_path)
        except:
            try:
                model = cobra.io.read_sbml_model(model_path)
            except:
                model = cobra.io.load_matlab_model(model_path)
    return model

class ProBar():
    '''
    Just a simple progress bar object to output on screen
    '''
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
################################################################################


def change_model_prefix(model_path, save_model=None, suffix = r'__91__(\w)__93__'):
    '''
    Function that translate model extracellular metabolites to other 'lenguages'
    '''
    mod = load_model(model_path)
    new_model = cobra.Model()
    new_model.id = mod.id
    p = re.compile(r'(.+)'+suffix)
    pb = re.compile(r'biomass', re.I)
    pdash = re.compile(r'DASH_')
    for i in mod.reactions:
        for reactant in i.metabolites:
            reactant.id=re.sub(pdash, r'',reactant.id)
            reactant.id=re.sub(p,r'\1[\2]',reactant.id)
        for react in i.reaction:
            react=re.sub(p,r'\1[\2]',react)
        new_model.add_reaction(i)
    for reac in new_model.reactions:
        if pb.search(reac.id):
            new_model.objective = [new_model.reactions.get_by_id(reac.id)]
            break
    if save_model != None:
        cobra.io.write_sbml_model(new_model,save_model)
    return new_model

############################# TRANSLATION WORKFLOW #############################
# 1. BiGG equivalences (name -> id)
def download_file(url):
    '''
    Function that downloads a file from url
    '''
    local_filename = url.split('/')[-1]
    if not os.path.exists(local_filename):
        global requests
        if not requests:
            import requests
        # NOTE the stream=True parameter below
        with requests.get(url, stream=True) as r:
            with open(local_filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk: # filter out keep-alive new chunks
                        f.write(chunk)
                        # f.flush()
    return local_filename


def lev_distance(seq1,seq2):
    '''
    Computes Levenshtein Distance
    INPUTS -> seq1 and seq2: strings to be compared.
    OUPUT -> int, Levenshtein distance (score).
    '''
    len1 = len(seq1) + 1
    len2 = len(seq2) + 1
    mat = np.zeros((len1, len2))
    for i in range(len1):
        mat[i,0] = i
    for i in range(len2):
        mat[0,i] = i
    for x in range(1, len1):
        for y in range(1, len2):
            if seq1[x-1] == seq2[y-1]:
                mat [x,y] = min(
                    mat[x-1, y] + 1,
                    mat[x-1, y-1],
                    mat[x, y-1] + 1
                )
            else:
                mat [x,y] = min(
                    mat[x-1,y] + 1,
                    mat[x-1,y-1] + 1,
                    mat[x,y-1] + 1
                )
    return (mat[len1 - 1, len2 - 1])

def write_eqs(mod_mets, outp):
    '''
    Interactive function, tries to match every metabolite in mod_mets with
    BiGG database based on name atribute. It accounts for perfect match or
    best similarity match.
    INPUTS -> mod_mets: list of COBRA metabolite objects
            output: string, path to output
    OUTPUT -> dict, met.name: id BiGG
    '''
    global np
    global csv
    if not np:
        import numpy as np
    if not csv:
        import csv
    print("Comparing with BiGG...\n")
    # Download all metabolites in BiGG to dict
    p_num = re.compile(r'^[\d]+$')
    p_name = re.compile(r'^([\da-zA-Z_\(\)\,\-\: ]+)_[CHONSPW(?:Fe)(?:Co)\d]+$')
    localf = download_file("http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt")
    metabolites = {}
    with open(localf, "r") as dbf:
        for line in dbf.readlines():
            linearr = line.split("\t")
            metabolites[linearr[2]] = linearr[1] # name -> id
    memod = {}
    # Load memo file if output has already been generated
    if os.path.isfile(outp):
        with open(outp) as f:
            reader = csv.DictReader(f, fieldnames = ["id", "name"], delimiter = "\t")
            memod = {row["name"]: row["id"] for row in reader}
    bar = ProBar(len(mod_mets))
    # Match them with model metabolites.
    dreacs = {}
    leave = ""
    with open(outp, "a") as f:
        for met in mod_mets:
            if met.name in memod:
                dreacs[met.name] = memod[met.name]
                bar.progress(mod_mets.index(met))
                continue
            elif met.name in metabolites:
                matched = metabolites[met.name]
            else:
                # Separate name from formula
                is_match = re.search(p_name, met.name)
                if is_match:
                    st_name = is_match.group(1).lower()
                elif met.name.endswith("_"):
                    st_name = met.name[:-1].lower()
                else:
                    st_name = met.name.lower()
                best = []
                best_distance = 100000
                for seq in metabolites:
                    candidate = (lev_distance(st_name, seq.lower()),seq)
                    if best_distance > candidate[0]:
                        best = [candidate[1]]
                        best_distance = candidate[0]
                    elif best_distance == candidate[0]:
                        best.append(candidate[1])
                if len(best) == 1 and best_distance < 10:
                    matched = metabolites[best[0]]
                else:
                    # Let the user decide
                    while True:
                        print()
                        print("Several matches were found for", met.name, "with score", str(best_distance))
                        for i in range(len(best)):
                            print(i, ". ", metabolites[best[i]], " -> ", best[i], sep = "")
                        choose = input("Which hit would you like to choose? [0-"+str(len(best)-1)+"]: ")
                        if p_num.match(choose):
                            matched = metabolites[best[int(choose)]]
                            break
                        else:
                            ided = input("Non numerical value provided, accept user provided id?[Y/N]: ")
                            if ided.lower() not in ["no", "n", "nein", "nao", "ez", "non"]:
                                matched = choose
                                break
                            else:
                                leave = input("Leave unmatched?[Y/N]: ")
                                if leave.lower() not in ["no", "n", "nein", "nao", "ez", "non"]:
                                    matched = met.id
                    print()
            line = matched+"\t"+met.name
            if leave.lower() in ["no", "n", "nein", "nao", "ez", "non"]:
                # to DEBUG
                line += "\tUNMATCHED"
            f.write(line + "\n")
            dreacs[met.name] = matched
            bar.progress(mod_mets.index(met))
    print("Searching BiGG finished.\n")
    return dreacs

# 2. Main call to translation with Levenshtein method.
def translate_from_BiGG(mod, outp = False, memo = "metabolites_modelxBiGG.tsv"):
    '''
    Function that translates model extracellular metabolites to BiGG nomenclature
    based on name of metabolites.
    INPUTS -> mod: string, path to sbml model
              output: string, path to output or False to don't write output
    OUPUT -> model translated
    '''
    model = load_model(mod)
    ex_mets = []
    for met in model.metabolites:
        if met.compartment == 'e' or met.compartment == 'e0':
            ex_mets.append(met)
    dcomps = write_eqs(ex_mets, memo)
    for met in ex_mets:
        model.metabolites.get_by_id(met.id).id = dcomps[met.name]+"_e"
    if outp:
        cobra.io.write_sbml_model(model, outp)
    return model
################################################################################

######################### MODEL POOL BASED TRANSLATION #########################
def process_all_formulas(allmets, based_on = "id"):
    '''
    Function that prepares metabolites to be processed in terms of formula.
    INPUTS -> metabolites: iterable with names/ids/object of metabolites to be
                           processed;
              based_on: "id" or "string".
    OUTPUT -> dictionary with metabolites name/id mapped to formula pieces.
    '''
    map_met = {}
    # I'd expect every item in input to be of the same class
    metabolites = list(allmets)
    if type(metabolites[1]) == str:
        p_name = re.compile(r'(.+)_([CHONSPKW(?:Fe)(?:Cd)(?:Zn)(?:Na)(?:Hg)(?:Co)(?:Ca)(?:Mg)(?:Mg)(?:Fru)(?:Glc)\d]*)$')
        for met in metabolites:
            is_match = re.search(p_name, met)
            try:
                st_name = is_match.group(1).lower()
                st_formula = is_match.group(2)
                map_met[st_name] = (break_formula(st_formula), met)
            except:
                print(f"Something's wrong with metabolite { met }")
                sys.exit(1)
    else:
        piece = "name" if based_on.lower() in ["name", "names"] else "id"
        for met in metabolites:
            map_met[met.__getattribute__(piece).lower()] = (break_formula(met.formula), met.__getattribute__(piece))
    return map_met

def break_formula(f, debug = False):
    '''
    Break a chemical formula into pieces to be easily compared later.
    '''
    if f == None:
        return {}
    fspl = {}
    if debug:
        print(f)
    for i in f:
        if i.isdigit():
            if curr:
                if fspl[curr] == "1":
                    fspl[curr] = i
                else:
                    fspl[curr] += i
            else:
                fspl[i] = "1"
        elif not i.istitle():
            if curr:
                fspl[curr+i] = fspl[curr]
                del(fspl[curr])
                curr = curr+i
            else:
                fspl[i] = "1"
        else:
            curr = i
            fspl[i] = "1"
    fspl = {k : float(fspl[k]) for k in fspl}
    return fspl

def score_formula1(f1,f2):
    '''
    Score the SIMILARITY of 2 formulas.
    Inputs are list of characters; output is float or false.
    '''
    if f1 == [] or f2 == []:
        return False
    elif abs(len(f1) - len(f2)) >= 2: # if theres that much difference, don't even compare
        score = -1000
    else:
        score = 0
        for piece1 in f1:
            if piece1 not in f2:
                score -= 5
            elif f1[piece1] == f2[piece1]:
                score += 3
            else:
                if piece1 == "C":
                    score -= abs(f1[piece1] - f2[piece1]) # penalize a lot differences in carbons
                else:
                    score -= abs(f1[piece1] - f2[piece1])/3
        for piece2 in f2:
            if piece2 not in f1:
                score -= 5
    return score

def score_formula(f1,f2):
    '''
    Binary strict version
    '''
    score = -1000
    if len(f1) != len(f2) or not f1 or not f2:
        return -1000
    for piece1 in f1:
        if piece1 not in f2:
            break
        elif f1[piece1] != f2[piece1]:
            break
    else:
        score = 1000
    return score

def get_metabolites_pool(models_scheme, based_on = "id"):
    '''
    Generate a set of extracellular metabolites names/ids
    '''
    piece = "name" if based_on.lower() in ["name", "names"] else "id"
    if hasattr(models_scheme, "cobrunion"):
        # if a Consortium object is provided, just call the apposite method
        return models_scheme.cobrunion
    elif type(models_scheme) not in [tuple, list]:
        # if a COBRA model object is provided, return a set ex metabolties
        return {met.__getattribute__(piece) for met in models_scheme}
    elif models_scheme == 1:
        return {met.__getattribute__(piece) for met in models_scheme[0]}
    else:
        return cobrunion(models_scheme, based_on = piece)

# TODO: make this more compact, combining both "write_eqs" functions
def write_eqs_models(allmets, db, outp):
    '''
    Interactive function, tries to match every metabolite in allmets with set of
    metabolites as "lenguage". It accounts for perfect match on formula and filter
    by Levenshtein distance among them.
    INPUTS -> allmets: list of strings (names or ids of COBRA metabolites) to
                        be translated;
            db: list of strings, names or ids of COBRA metabolites in the correct
                        lenguage;
            outp: string, path to file where output will be written and whose content
                        will be used as memoization speed-up, if possible;
            output: string, path to output.
    OUTPUT -> dict, met.name (model to be translated): met.name/id (scheme models)
    '''
    global np
    global csv
    if not np:
        import numpy as np
    if not csv:
        import csv

    dcomps = {} # output returned
    memod = {} # info from memofile
    # Load memo file if output has already been generated
    if os.path.isfile(outp):
        with open(outp) as f:
            reader = csv.DictReader(f, fieldnames = ["id", "name"], delimiter = "\t")
            memod = {row["name"]: row["id"] for row in reader}
    else:
        with open(outp, "w") as f:
            f.write("\t".join(["Metabolite", "Match", "Automatic"])+"\n")

    p_num = re.compile(r'^[\d]+$')
    print("Preparing input...")
    mod_mets = process_all_formulas(allmets)
    metabolites = process_all_formulas(db)
    print("Translating...")
    bar = ProBar(len(mod_mets))
    leave = ""
    f = open(outp, "a")
    i = 0
    for met in mod_mets:
        auto = True # this will be appended to file as a column
        if met in memod:
            dcomps[mod_mets[met][1]] = memod[met]
            bar.progress(i)
            i+=1
            continue
        elif met in metabolites:
            matched = met
        else:
            # first, compare with formula
            best_for = []
            best_score = -100
            for seq in metabolites:
                candidate = (score_formula(metabolites[seq][0], mod_mets[met][0]), metabolites[seq][1])
                if candidate[0] >= best_score:
                    best_for = [candidate[1]]
                    best_score = candidate[0]
                elif candidate[0] == best_score:
                    best_for.append(candidate[1])
            best_distance = 100000
            best = []
            # second, if more than one match was retrieved, filter by Levenshtein distance
            if len(best_for) > 1:
                for seq in best_for:
                    candidate = (lev_distance(seq.lower(), met.lower()), seq)
                    if best_distance > candidate[0]:
                        best = [candidate[1]]
                        best_distance = candidate[0]
                    elif best_distance == candidate[0]:
                        best.append(candidate[1])
            else:
                best = best_for
            if len(best) == 1 and best_score > len(mod_mets[met][0])*2:
                matched = best[0]
            elif best_score == -100:
                matched = "UNMATCHED"
            else:
                # Let the user decide
                auto = False
                while True:
                    print()
                    print(f"Several matches were found for {met} with score {int(best_distance)}")
                    for i in range(len(best)):
                        print(i, ". ", best[i], sep = "")
                    choose = input(f"Which hit would you like to choose? [0-{len(best)-1}]: ")
                    if p_num.match(choose):
                        matched = best[int(choose)]
                        break
                    else:
                        ided = input("Non numerical value provided, accept user provided id?[Y/N]: ")
                        if ided.lower() not in ["no", "n", "nein", "nao", "ez", "non"]:
                            matched = choose
                            break
                        else:
                            leave = input("Leave unmatched?[Y/N]: ")
                            if leave.lower() in ["no", "n", "nein", "nao", "ez", "non"]:
                                matched = met
                            else:
                                matched = "UNMATCHED"
                            break
                print()
        line = "\t".join([mod_mets[met][1], matched, str(auto)])
        f.write(line + "\n")
        dcomps[mod_mets[met][1]] = matched
        bar.progress(i)
        i+=1
    f.close()
    return dcomps

def translate_from_models(model, models_scheme, based_on = "id", scheme_string = False, memo = "", save_model = False):
    '''
    Translate model's extracellular metabolites (names/ids) to the lenguage of
    other models supplied.
    INPUTS -> model, COBRA model object to be translated;
            models_scheme: list of [cobra models or strings] or Consortium object;
            based_on (id or name): attribute to be translated;
            scheme_string: bool, if models_scheme is a list of strings;
            memo: string, path to memoization file;
            save_model : False or string, path to model to be written as file.
    OUTPUT -> COBRA model object translated
    '''
    global picle
    if not pickle:
        import dill as pickle
    based_on = "name" if based_on.lower() in ["names", "name"] else "id"
    print("Loading model metabolites...")
    # take ExtraCellular metabolites from model to translate
    if not os.path.isfile("mets_to_translate.p"):
        model = load_model(model)
        mets_to_translate = [met for met in model.metabolites if met.compartment in ["e", "e0"]]
        with open("mets_to_translate.p", 'wb') as f:
            pickle.dump(mets_to_translate, f)
    else: # memo if you can!
        with open("mets_to_translate.p", 'rb') as f:
            mets_to_translate = pickle.load(f)
    # take extracelullar union of ExtraCellular metabolites to use as lenguage
    print("Working on the pool...")
    if scheme_string:
        metabolites = set(models_scheme)
    else:
        metabolites = get_metabolites_pool(models_scheme, based_on)
    # call interactive translating function
    dcomps = write_eqs_models(metabolites, mets_to_translate, memo) # WARNING: there's a tweak here JUST TO GENERATE THE OUTPUT MEMO
    # translate the model with the dictionary generated
    for met in model.metabolites:
        if met.name in dcomps:
            if dcomps[met.name] != "UNMATCHED":
                setattr(model.metabolites.get_by_id(met.id), based_on, dcomps[met.name])
    if save_model: # save to file
        cobra.io.load_matlab_model(model, save_model)
    return model
