#!/usr/bin/env python3

# I/O utilities to manage and translate models to one 'lenguage'.

# Author: Jorge Carrasco Muriel
# e-mail: jorge.cmuriel@alumnos.upm.es
# Date of 1st version: 30/12/2018

import cobra
import os, re
from warnings import filterwarnings as wfilt
requests = 0
np = 0

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
        wfilt("ignore", category=UserWarning)
        try:
            model = cobra.io.read_sbml_model(model_path)
        except:
            try:
                model = cobra.io.load_matlab_model(model_path)
            except:
                model = cobra.io.load_json_model(model_path)
    return model
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

############################ TRANLASTAION WORKFLOW #############################
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
def healthy_translate(mod, outp = False, memo = "metabolites_modelxBiGG.tsv"):
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
