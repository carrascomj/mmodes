#!/usr/bin/env python3

# I/O utilities to generate COMETS-like files and translate models to one 'lenguage'.

# Author: Jorge Carrasco Muriel
# e-mail: jorge.cmuriel@alumnos.upm.es
# Date of 1st version: 30/12/2018
# Date of last update: 14/01/20189

import cobra
from decimal import Decimal
import os

class ImplementationError(Exception):
    """ If Manifest was instantiated without models, raise (Implementation error,
    it should be removed after it's properly integrated with mmodes """
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
        self.write_manifest()

    def write_manifest(self):
        with open(self.fman, "w") as f:
            # should it write a layout?
            f.write("LayoutFileName: no_lay.txt\n")
            for mod in sorted(self.models):
                f.write("ModelFileName: "+self.models[mod][0].path+"\n")
            f.write("FluxFileName: "+self.fmedia+"\n")
            f.write("MediaFileName: "+self.fflux+"\n")
            f.write("TotalBiomassFileName: "+self.fbiom+"\n")
        return

    def get_models(self, models):
        """ Getter of models attribute
        INPUT -> models: dict
        OUPUT -> models: dictionary of model.id : [cobra_model, 0] """
        if not models:
            raise ImplementationError("You haven't passed models to Manifest!")
        elif isinstance(models, dict):
            models_out = {}
            num_mod = 1
            for mod in sorted(models):
                models_out[mod] = [models[mod],0, str(num_mod)]
                num_mod += 1
        return models_out

    def get_media(self, media):
        """ Getter of media
        INPUT -> media, dictionary of met.id: value/concentration
        OUPUT -> media, dict"""
        if isinstance(media, dict):
            return media
        else:
            raise ImplementationError("Media has to be a dictionary")

    def write_biomass(self):
        """ Write biomass log file """
        with open(self.fbiom, "a") as f:
            f.write(str(self.T)+"\t"+"\t".join([str(self.models[k][0].volume.q) for k in sorted(self.models)])+"\n")
        return

    def write_media(self):
        """ Write media, log file """
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
                    # TODO: check this
                    line += "media_"+str(self.T)+"{"+str(i)+"}(1, 1) = "+str('%.2E' % Decimal(self.media[met]))+";\n"
                i+=1
            f.write(line)
        return

    def write_fluxes(self, model, t):
        """ Writes fluxes log file
        INPUT -> t, float representing time
                model, cobra model object"""
        # lines have got the format "fluxes{t}{1}{1}{num strain} = [1E0 1.3E21 ...]"
        if self.curr_t != t:
            # since COMETS is fixed-step and solver calls dinamicpFBA several times,
            # managing of time value is a bit tricky
            self.T += 1
            self.curr_t = t
        if self.T == self.models[model.id][1]:
            # We're just taking the 1st call of each time value (just because
            # it's easier, maybe it should be implemented as mean value)
            return
        else:
            self.models[model.id][1] = self.T
            with open(self.fflux, "a") as f:
                f.write("fluxes{"+str(self.T)+"}{1}{1}{"+self.models[model.id][2]+"} = ["+" ".join(['%.2E' % Decimal(reac.flux) for reac in model.reactions])+" ];\n")
        return

# Some util functions to manage data

def find_models(dir_models):
    """Find models in "dir_models" path
    INPUT -> path to models directory (string)
    OUPUT -> list of model objects"""
    path_models = []
    for file in os.listdir(dir_models):
        # Find all sbml or matlab models and exclude comets models
        if (file.find("xml") != -1 or file.find("mat") != -1) and file.find("xml.cmt") == -1:
            path_models.append(dir_models+file)
    models = []
    strain_number = 1
    for mod in path_models:
        try:
            models.append(cobra.io.load_matlab_model(mod))
        except:
            models.append(cobra.io.read_sbml_model(mod))
        # models_description attribute will be used to name files, maybe it's
        # better to just apply "strain_"...
        if models[-1].description == "":
            models[-1].description = "strain_" + str(strain_number)
        strain_number += 1
    return models

def cobrunion(models=[], AGORA = True):
    """union of the models as a list
    INPUT -> list of model_objects
    OUPUT -> list of union of extracellular metabolites"""
    ex_mets = set()
    for model in models:
        for met in model.metabolites:
            if met.id.find("[e]") != -1:
                ex_mets.add(met.id)
    return ex_mets

def json_parser(path):
    """ Generator that parse a JSON file which contains an array of hashes.
    INPUT -> string, path to JSON file
    OUTPUT -> dict, yield a media """
    with open(path) as json_file:
        json_data = json.load(json_file)
    for i in range(len(json_data)):
        yield json_data[i]

def to_comets(model, dir_models):
    """Write a COMETS model file based on a model object
    INPUT -> model object,
             string of directory path where models can be found
    OUPUT -> creates a "xml.cmt" file,
             return a path of model"""
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
    """Writes the comets layout
    INPUTS -> files_model: list of strings containing model files,
              ex_mets: list of extracellular metabolites,
              media: dictionary of metabolite:concentration
              dictionary of files_model:initial biomass
              ... some of the COMETS parameters, that are described in http://www.bu.edu/segrelab/comets-parameters/
    OUPUTS -> creates a COMETS layout in "Consortium_layout.txt"""
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
    """To use as a module
    INPUTS -> dir_models: directory where models can be found
              file_media: path to media file, default as dir_models/media.json
              ... some of the COMETS parameters, that are described in http://www.bu.edu/segrelab/comets-parameters/
    OUPUTS -> creates a COMETS layout in "Consortium_layout.txt",
              returns the name of the outfile"""
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

def translate_model(model_path, save_model=None, suffix = r'__91__(\w)__93__'):
    """ Function that translate model extracellular metabolites to other models"""
    mod = cobra.io.read_sbml_model(model_path)
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
