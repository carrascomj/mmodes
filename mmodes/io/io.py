#!/usr/bin/env python3

# I/O utilities to generate COMETS-like files and translate models to one 'lenguage'.
# TODO: it needs to be separated into different files.

# Author: Jorge Carrasco Muriel
# e-mail: jorge.cmuriel@alumnos.upm.es
# Date of 1st version: 30/12/2018

import cobra
from decimal import Decimal
import os
import re
from warnings import filterwarnings as wfilt
np = 0
requests = 0

class ImplementationError(Exception):
    '''
    If Manifest was instantiated without models, raise (Implementation error,
    it should be removed after it's properly integrated with mmodes
    '''
    pass

class Manifest():
    def __init__(self, models = [], media = {}, fman = "COMETS_manifest.txt", fflux = "flux_log_template.txt", fbiom = "total_biomass_log_template.txt", fmedia = "media_log_template.txt"):
        self.fflux = fflux
        self.fbiom = fbiom
        self.fmedia = fmedia
        self.fman = fman
        self.models = self.get_models(models)
        self.media = self.get_media(media)
        self.T = 0
        self.curr_t = 0
        self.first = 0
        self.write_manifest()

    def write_manifest(self):
        with open(self.fman, "w") as f:
            # should it write a layout?
            f.write("LayoutFileName: no_lay.txt\n")
            for mod in sorted(self.models):
                f.write("ModelFileName: "+self.models[mod][0].path+"\n")
            f.write("FluxFileName: "+self.fflux+"\n")
            f.write("MediaFileName: "+self.fmedia+"\n")
            f.write("TotalBiomassFileName: "+self.fbiom+"\n")
        return

    def get_models(self, models):
        '''
        Getter of models attribute
        INPUT -> models: dict
        OUPUT -> models_out: dictionary of model.id : [cobra_model, 0, num_model]
        '''
        if not models:
            raise ImplementationError("You haven't passed models to Manifest!")
        elif isinstance(models, dict):
            models_out = {}
            num_mod = 1
            for mod in sorted(models):
                models_out[mod] = [models[mod], 0, str(num_mod)]
                num_mod += 1
        return models_out

    def get_media(self, media):
        '''
        Getter of media
        INPUT -> media: dictionary of met.id: value/concentration
        OUPUT -> media: dict
        '''
        if isinstance(media, dict):
            return media
        else:
            raise ImplementationError("Media has to be a dictionary")

    def write_biomass(self):
        '''
        Write biomass log file
        '''
        with open(self.fbiom, "a") as f:
            f.write(str(self.T)+"\t"+"\t".join([str(self.models[k][0].volume.q) for k in sorted(self.models)])+"\n")
        return

    def write_media(self):
        '''
        Write media, log file
        '''
        if not os.path.isfile(self.fmedia):
            st_list = ["'"+k+"'" for k in sorted(self.media)]
            with open(self.fmedia, "w") as f:
                f.write("media_names = { "+", ".join(st_list)+"};\n")
        with open(self.fmedia, "a") as f:
            line = ""
            i = 1
            for met in sorted(self.media):
                line += "media_"+str(self.T)+"{"+str(i)+"} = sparse(zeros(1,1));\n"
                if float(self.media[met]) != 0:
                    line += "media_"+str(self.T)+"{"+str(i)+"}(1, 1) = "+str('%.2E' % Decimal(self.media[met]))+";\n"
                i+=1
            f.write(line)
        return

    def write_fluxes(self, model, t):
        '''
        Writes fluxes log file
        INPUT -> t, float representing time
                model, cobra model object
        '''
        # lines have got the format "fluxes{t}{1}{1}{num strain} = [1E0 1.3E21 ...]"
        if self.curr_t != t:
            # since COMETS is fixed-step and solver calls dinamicpFBA several times,
            # management of time value is a bit tricky
            self.T += 1
            self.curr_t = t
        if self.T != self.models[model.id][1] or self.first < len(self.models):
            # 1st call is the only one being taken
            self.first += 1
            self.models[model.id][1] = self.T
            with open(self.fflux, "a") as f:
                if model.solver.status: # model has been optimized at least 1 time
                    f.write("fluxes{"+str(self.T)+"}{1}{1}{"+self.models[model.id][2]+"} = ["+" ".join(['%.2E' % Decimal(reac.flux) for reac in model.reactions])+" ];\n")
                else:
                    f.write("fluxes{"+str(self.T)+"}{1}{1}{"+self.models[model.id][2]+"} = ["+" ".join(['%.2E' % Decimal(0) for reac in model.reactions])+" ];\n")

        return

########################## READ FROM COMETS FUNCTIONS ##########################

def cobrunion(models=[], AGORA = True):
    '''
    union of the models as a list
    INPUT -> list of model_objects
    OUPUT -> list of union of extracellular metabolites
    '''
    ex_mets = set()
    for model in models:
        for met in model.metabolites:
            if met.compartment == 'e' or met.compartment == 'e0':
                ex_mets.add(met.id)
    return ex_mets

def json_parser(path):
    '''
    Generator that parse a JSON file which contains an array of hashes.
    INPUT -> string, path to JSON file
    OUTPUT -> dict, yield a media
    '''
    with open(path) as json_file:
        json_data = json.load(json_file)
    for i in range(len(json_data)):
        yield json_data[i]

def to_comets(model, dir_models):
    '''
    Write a COMETS model file based on a model object
    INPUT -> model object,
             string of directory path where models can be found
    OUPUT -> creates a "xml.cmt" file,
             return a path of model
    '''
    if cobra.__version__ == "0.13.4":
        old_v = False
    else:
        old_v = True
    # Open output file:
    sbmlInputFile = dir_models + model.description + ".xml.cmt"
    with open(sbmlInputFile, mode='w') as f:
        # Print the S matrix
        f.write("SMATRIX  "+str(len(model.metabolites))+"  "+str(len(model.reactions))+"\n")
        for x in range(len(model.metabolites)):
            for y in range(len(model.reactions)):
                if (model.metabolites[x] in model.reactions[y].metabolites):
                    coeff=model.reactions[y].get_coefficient(model.metabolites[x])
                    f.write("    "+str(x+1)+"   "+str(y+1)+"   "+str(coeff)+"\n")
        f.write("//\n")

        # Print the bounds
        f.write("BOUNDS  -1000  1000\n");
        for y in range(len(model.reactions)):
            lb=model.reactions[y].lower_bound
            up=model.reactions[y].upper_bound
            f.write("    "+str(y+1)+"   "+str(lb)+"   "+str(up)+"\n")
        f.write("//\n")

        # Print the objective reaction
        f.write('OBJECTIVE\n')
        if not old_v:
            obj=str(model.objective.expression.as_leading_term()).split(sep='+')
            try:
                obj[1][5:]
                obj_rr=obj[1][5:]
            except IndexError:
                obj=str(model.objective.expression.as_leading_term()).split(sep='-')
                obj_rr=obj[0][4:-1]
            for y in range(len(model.reactions)):
                if (model.reactions[y].id == obj_rr):
                    indexObj=y+1
        else:
            for y in range(len(model.reactions)):
                if (model.reactions[y] in model.objective):
                    indexObj=y+1
        f.write("    "+str(indexObj)+"\n")
        f.write("//\n")

        # Print metabolite names
        f.write("METABOLITE_NAMES\n")
        for x in range(len(model.metabolites)):
            f.write("    "+model.metabolites[x].id+"\n")
        f.write("//\n")

        # Print reaction names
        f.write("REACTION_NAMES\n")
        for y in range(len(model.reactions)):
            f.write("    "+model.reactions[y].id+"\n")
        f.write("//\n")

        # Print exchange reactions
        f.write("EXCHANGE_REACTIONS\n")
        for y in range(len(model.reactions)):
            if (model.reactions[y].id.find('EX_')==0):
                f.write(" "+str(y+1))
        f.write("\n//\n")
    return sbmlInputFile

def write_layout(files_model, ex_mets, media, biomasses, outfile = "Consortium_layout.txt", maxCycles = 1000, timeStep = 0.1,  spaceWidth = 0.05, deathRate = 0, numRunThreads = 2, maxSpaceBiomass = 10, biom_name = "total_biomass_log_template.txt", media_name = "media_log_template.txt", flux_name = "flux_log_template.txt"):
    '''
    Writes the comets layout
    INPUTS -> files_model: list of strings containing model files,
              ex_mets: list of extracellular metabolites,
              media: dictionary of metabolite:concentration
              dictionary of files_model:initial biomass
              ... some of the COMETS parameters, that are described in http://www.bu.edu/segrelab/comets-parameters/
    OUPUTS -> creates a COMETS layout in "Consortium_layout.txt
    '''
    # It's important to notice that name of model files are assigned in terms of
    # model.description. If models are saved without it, the assignment is ASCII
    # ordered.
    with open(outfile, "w") as f:
        # Print lines 1-4
        lines4 = "model_file\t" + "\t".join(files_model)
        lines4 += "\n\tmodel_world\n\t\tgrid_size\t1\t1\n\t\tworld_media\n"
        f.write(lines4)

        # Print extracellular metabolites with concentrations
        for met in ex_mets:
            if met not in media:
                f.write("\t\t\t"+met+"\t"+"0.0\n")
            else:
                f.write("\t\t\t"+met+"\t"+str(media[met])+"\n")

        # Print end of world_media
        f.write("\t\t//\n\t\tmedia\n\t\t//\n\t//\n\tinitial_pop\n")

        # Print initial biomasses
        biomass_line = "\t\t0\t0"
        for file in files_model:
            biomass_line += "\t" + biomasses[file]
        f.write(biomass_line + "\n\t//\n")

        # Print parameters
        params = "parameters\n    maxCycles = "+str(maxCycles)+"\n"
        params += "    timeStep = "+str(timeStep)+"\n"
        params += "    spaceWidth = "+str(spaceWidth)+"\n"
        params += "    deathRate = "+str(deathRate)+"\n"
        params += "    maxSpaceBiomass = "+str(maxSpaceBiomass)+"\n"
        params += "    numRunThreads = "+str(numRunThreads)+"\n"
        params += "    pixelScale = 100\n"
        params += "    allowCellOverlap = true\n"
        params += "    saveslideshow = false\n"
        params += "    writetotalbiomasslog = true\n"
        params += "    totalbiomasslogname = "+biom_name+"\n"
        params += "    writemedialog = true\n"
        params += "    medialogname = "+media_name+"\n"
        params += "    writefluxlog = true\n"
        params += "    fluxlogname = "+flux_name+"\n"
        params += "    uselognametimestamp = false\n"
        params += "    exchangestyle = Monod Style\n"
        params += "    defaultKm = 0.01\n"
        params += "    defaultVmax = 20\n//\n"
        f.write(params)
    return

def write_cmt_script(outfile):
    with open ("comets_script_template", "w") as comets:
        comets.write("load_layout  "+str(outfile))
    # grant execution permissions (maybe should be granted when its created)
    # mode = ((os.stat("comets_script_template_tmp").st_mode) | 0o555) & 0o7777
    # os.chmod("comets_script_template_tmp", mode)

def make_lay(dir_models = "", file_media = "", outfile = "Consortium_layout.txt",  maxCycles = 1000, timeStep = 0.1,  spaceWidth = 0.05, deathRate = 0, numRunThreads = 2, maxSpaceBiomass = 10):
    '''
    To use as a module
    INPUTS -> dir_models: directory where models can be found
              file_media: path to media file, default as dir_models/media.json
              ... some of the COMETS parameters, that are described in http://www.bu.edu/segrelab/comets-parameters/
    OUPUTS -> creates a COMETS layout in "Consortium_layout.txt",
              returns the name of the outfile
    '''
    if dir_models[-1] != "/":
        dir_models += "/"
    if file_media == "":
        file_media = dir_models+"media.json"
    # Instantiate models in a list
    models = find_models(dir_models)
    # Union of extracellular metabolites to include in layout
    ex_mets = cobrunion(models)
    files_model = []
    for model in models:
        files_model.append(to_comets(model, dir_models))
    # Generate comets layout file
    media_generator = json_parser(file_media)
    media = next(media_generator) # dictionary of metabolite:concentration
    biomasses = next(media_generator) # dictionary of files_model:initial biomass
    write_layout(files_model, ex_mets, media, biomasses, outfile, maxCycles, timeStep, spaceWidth, deathRate, numRunThreads, maxSpaceBiomass)
    write_cmt_script(outfile)
    return outfile

def plot_parser(suffix='', mets = [], endCycle=int(100)):
    '''
    Parser that translates the output of comets to a tsv easy-to-plot file
    INPUTS -> suffix: of the media_log file
            mets: list of metabolites to be later plotted
            endCycle: not necessary at all...
    OUPUTS -> returns the names of the file generated and the file to be later ploted
            generates a TSV
    '''
    mets_title = ""
    num_mets=[]
    with open("media_log_"+suffix+".txt") as media:
        lines = [line for line in media.readlines()]
    # 1st line as list
    media_names = re.sub(r'(?:media_names)|(?:\[e\])|=|{|}|\'| |\n|;', r'',lines[0]).split(sep=",")
    for met in mets:
        mets_title = mets_title + "_" + met
        # extract the number of metabolites
        num_mets.append(media_names.index(met)+1)
    outFile = "biomass_vs"+mets_title+"_"+suffix+".txt"
    plotFile = "biomass_vs"+mets_title+"_plot.pdf" # plotFile path will be used just once; no need to differ among iterations
    # Compile the regex patterns just once
    patterns = [re.compile(r'media\_\d+\{'+str(num)+r'\}\(1, 1\) = ([\d\.E\-]+)') for num in num_mets]
    patterns_0 = [re.compile(r'media\_\d+\{'+str(num)+r'\} = s') for num in num_mets]
    concs_to_print = [ [] for i in range(len(patterns))]
    for i in range(1,len(lines)):
        for p in range(len(patterns)):
            match = patterns[p].search(lines[i])
            match_0 = patterns_0[p].search(lines[i])
            if match_0:
                concs_to_print[p].append("0")
            if match:
                concs_to_print[p].pop() # delete last 0
                concs_to_print[p].append(match.group(1))
                #concs_to_print[p].append(match.groups(1)[0])
    # Now, concs_to_print is a list of lists containing [metabolites][time]
    with open(outFile, "w") as FHO:
        with open("total_biomass_log_"+suffix+".txt") as FHI:
            i = 0
            for line_in in FHI.readlines():
                if i > endCycle:
                    break
                line_add = ""
                for conc in concs_to_print:
                    line_add = line_add + "\t" + conc[i]
                # "[-1]" so it doesn't include \n
                FHO.write(line_in[:-1]+line_add+"\n")
                i += 1
    # Returns the paths to be managed by other functions
    return outFile, plotFile

def get_media_composition(path_to_media = "media_log_template.txt", suff = "[e]"):
    '''
    Parser that translates the information generated by COMETS of the last cycle
    media compostion of the last iterarion to a python dictionary.
    INPUTS -> path_to_media: string, path to media_log_file,
            suff: string that indicates compartment on each metabolites of the models,
                  such as '[e]' or '_e' or '__91__e__93__'
    OUPUTS -> met_dict: dictionary, metabolites as keys and concentrations as values
    '''
    with open(path_to_media) as mediaf:
        line1 = False
        lines = ""
        met_dict = {}
        for line in mediaf.readlines():
            if not line1:
                line1 = re.sub(r'(?:media_names)|(?:\[e\])|=|{|}|\'| |\n|;', r'',line).split(sep=",")
            else:
                lines += line
                lastline = line
    # This line could be substituted by another argument of the function. For now, it's just taking the last cycle computed
    # that should have been imposed in the layout
    end_cycle = re.search(r'media_(\d+)',lastline).group(1)
    for i in range(1,len(line1)+1):
        conc = re.search(r'media\_'+str(end_cycle)+r'\{'+str(i)+r'\}(?:\(1, 1\)) = (.+);', lines)
        #conc = 0,re.search(r'media\_'+str(end_cycle)+r'\{'+str(i)+r'\}(?:\(1, 1\)) = (.+);', lines).group(1)  # python <=3.5
        if not conc:
            conc = 0
        else:
            conc = conc.group(1)
        met_dict[line1[i-1]+suff] = conc
    return met_dict

def parse_flux(ff, t = 10, st = 2):
    '''
    Select line from a COMETS flux_log file
    INPUTS -> ff: str, path to flux_log COMETS file
            t: int, time to extract
            st: int, number of strain
    OUTPUT -> string, fluxes separated by ' '
    '''
    p = re.compile(r'^fluxes\{'+str(t)+r'\}\{\d+\}\{\d+\}\{'+str(st)+r'\} = \[(.+) ?\]')
    with open(ff) as f:
        for line in f.readlines():
            match = p.search(line)
            if match:
                return match.group(1)
                break
        else:
            print("No matched string")
            return ""

def write_flux_line(cha, mod, outp, sep = "\t"):
    '''
    Write tsv with reactions and fluxes from 'cha'
    INPUTS -> cha: str, from comets flux log
            mod: str, path to model (or cobra model object)
            outp: str, path to output
    '''
    # Prepare input
    fluxes = cha.replace(r' ', sep)
    if not os.path.isfile(outp):
        # if 1st call, write header
        model = load_model(mod)
        line1 = "\t".join([reac.id for reac in model.reactions])
        with open(outp, "w") as f:
            f.write(line1+"\n")
    with open(outp, "a") as f:
        f.write(fluxes[:-1]+"\n")
######################## END READ FROM COMETS FUNCTIONS ########################

############################# MANAGEMENT OF MODELS #############################

def find_models(dir_models, just_path = False):
    '''
    Find models in "dir_models" path
    INPUTS -> dir_models: path to models directory (string)
            just_path: bool, True if just outputs a string with the path.
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
        # models_description attribute will be used to name files, maybe it's
        # better to just apply "strain_"...
        if models[-1].description == "":
            models[-1].description = "strain_" + str(strain_number)
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
    # Charge memo file if output has already been generated
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

########################### END MANAGEMENT OF MODELS ###########################
