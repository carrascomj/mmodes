#!/usr/bin/python3

# Example script for mmodes package.
# Author: Jorge Carrasco Muriel
# e-mail: jorge.cmuriel@alumnos.upm.es

import json
import mmodes

def main():
    with open("ModelsInput/media.json") as jdata:
        media = json.load(jdata)[0]
    # 1) instantiate Consortium
    # if "manifest" contains a non-empty string, it will generate COMETS-like ouput
    cons = mmodes.Consortium(stcut = 1e-7, mets_to_plot = ["glc_D[e]", "ac[e]"], v = 1, manifest = "COMETS_manifest.txt")
    # 2) instantiate dMetabolites
    # for instance, https://www.ncbi.nlm.nih.gov/pubmed/18791026?dopt=Abstract
    glc = mmodes.dMetabolite(id = "glc_D[e]", Km = 14.8, Vmax = 0.13)
    # 3) add model
    # model from AGORA. A single strain to agilize example.
    cons.add_model("ModelsInput/BGN4_eu.xml", 0.001, solver = "glpk", method = "pfba", dMets = {glc.id: glc}) # Bifidobacterium_adolescentis_ATCC_15703.xml (euavg)
    # 4) instantiate media
    cons.media = cons.set_media(media, True)
    # 5) run it, plotting the output
    cons.run(maxT = 10, outp = "plot_example.png", outf = "tab_example.tsv", verbose=True)
    # 6) print stuff on screen
    for mod in cons.models:
        print(mod, cons.models[mod].volume.q, sep = " -> ")
    print("Glucose", str(cons.media["glc_D[e]"]), sep = " -> ")
    print("Acetate", str(cons.media["ac[e]"]), sep = " -> ")
    print()

if __name__ == "__main__":
    main()
