#!/usr/bin/python3

# Functions to visualize and analize the output generated.

# Author: Jorge Carrasco Muriel
# e-mail: jorge.cmuriel@alumnos.upm.es
# Date of first version: 14/03/2019

import os
import numpy as np
from copy import deepcopy as dc
pd = 0 # pandas
plt = 0 # matplolib.pyplot


def plot_comm(cons, color_set = "tableau20"):
    '''
    Plots the concentration of given microorganisms and metabolites
    INPUT -> cons: MMODES consortium object already simulated
    OUTPUT -> returns nothing, generates a plot in 'cons.outplot'
    '''
    # Call always after cons.runn()!
    if not hasattr(cons, "outplot"):
        print("Consortium object must run before analyzing the output of the run!")
        return
    global plt
    if not plt:
        import matplotlib.pyplot as plt
    global pd
    if not pd:
        import pandas as pd


    title = cons.title
    path = cons.output
    output = cons.outplot
    pallettes = {
        # TODO: maybe, this should be implemented with palletable...
        # some colors to plot...
        "colors_hex_2" : ["#ff8c00", "#a020f0", "#00ff00", "#ffff00", "#7a8b8b", "#cd5555", "#1e90ff", "#787878", "#ff7f50", "#000000"],
        # or nominal data color scheme found in http://geog.uoregon.edu/datagraphics/color/Cat_12.txt
        "colors_hex" : ["#ff7f00", "#32ff00", "#19b2ff", "#654cff", "#e51932", "#000000", "#ffff32", "#ff99bf", "#ccbfff", "#a5edff", "#b2ff8c", "#ffff99", "#ffbf7f"],
        # or tableau20 http://tableaufriction.blogspot.com/2012/11/finally-you-can-use-tableau-data-colors.html
        "tableau20" : ["#1F77B4","#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#17BECF", "#BCBD22", "#AEC7E8", "#98DF8A", "#C5B0D5", "#F7B6D2", "#DBDB8D"]
    }
    colors = pallettes[color_set]
    to_plot=pd.read_csv(path, sep='\t', header = 0)
    # check if metabolites where correctly specified
    mets_ok = False
    for i in cons.mets_to_plot:
        if i in to_plot:
            mets_ok = True
        else:
            print(f"\nMetabolite {i} won't be plotted. Check your spelling of this metabolite.")
    if not mets_ok:
        print("Metabolites weren't properly supplied in 'cons.mets_to_plot'. Plot won't be generated!")
        return False

    for col in to_plot: # leave just selected metabolites and biomasses
        if col not in cons.mets_to_plot and col not in cons.orgs_to_plot and col != "time":
            del to_plot[col]
    to_plot.loc[:, to_plot.columns != 'time'] = to_plot.loc[:, to_plot.columns != 'time'] / cons.v # return concentrations
    fig, ax1 = plt.subplots()
    t = to_plot.time
    i = 0
    for col in to_plot:
        if 0 < i < len(cons.orgs_to_plot) + 1: # plotting biomasses
            s1 = to_plot[col]
            ax1.plot(t, s1, colors[i], alpha = 0.8)
        elif i >= len(cons.orgs_to_plot) + 1: # plotting metabolites
            if i == len(cons.orgs_to_plot) + 1: # set new y axis
                ax1.set_xlabel("time(h)", fontstyle = 'italic')
                ax1.set_ylabel('organisms (g/L)', fontstyle = 'italic')
                ax1.tick_params('y')
                ax2 = ax1.twinx()
            s2 = to_plot[col]
            ax2.plot(t, s2, colors[i], linestyle='--')
            ax2.set_ylabel('metabolites (mmol/L)', fontstyle = 'italic')
            ax2.tick_params('y')
        i += 1
    plt.title(title, loc = 'left')
    N = 6
    ymin, ymax = ax1.get_ylim()
    ax1.set_yticks(np.round(np.linspace(0, ymax, N), 5))
    ymin, ymax = ax2.get_ylim()
    ax2.set_yticks(np.round(np.linspace(0, ymax, N), 2))
    ax1.grid(True, 'major', ls='-', color = 'grey', alpha=.3)
    ax1.legend(loc=6)
    ax2.legend(loc=2)
    fig.tight_layout()
    plt.savefig(output, dpi = 300)
    return

def walkplot(cons):
    '''
    Interactive function to walk through the metabolites in the media. It prints
    a growth plot of the consortium for each four metabolites.
    INPUT -> cons: MMODES consortium object already simulated
    '''
    # Call always after cons.runn()!
    if not hasattr(cons, "outplot"):
        print("Consortium object must run before analyzing the output of the run!")
        return
    global plt
    if not plt:
        import matplotlib.pyplot as plt

    # keep original inputs
    original_metplots = dc(cons.mets_to_plot)
    original_outplot = dc(cons.outplot)
    cons.outplot = "tmp.png"
    mets = [k for k in cons.media]
    # loop over metabolites, plot them and let the user close the window
    for m in range(0, len(mets), 4):
        cons.mets_to_plot = [mets[m], mets[m+1], mets[m+2], mets[m+3]]
        cons.plot_comm()
        plt.show()
        print("Image number", m/4)

    # clean temporary files, return to orginal parameters
    plt.close("all")
    os.remove(cons.outplot)
    cons.mets_to_plot = original_metplots
    cons.outplot = original_outplot
    return

def find_active_mets(cons, sign = "both"):
    '''
    Function that prints to screen metabolites that have changed during the simulation
    sorted by proportional value of change (or absolute value, if initial condition is 0).
    INPUTS -> cons: MMODES consortium object already simulated;
              sign: one of ["both", "+", "-"], it prints both, positive or negative increments.
    OUTPUT -> metabolites ordered by absoulute value of change.
    '''
    # TODO: order by fold change
    global pd
    if not pd:
        import pandas as pd
    def keysort(elem):
        '''Key function to sort the output'''
        if elem[1] == 0:
            return abs(float(elem[3][10:-10]))
        else:
            return abs(float(elem[3][10:-10]))/elem[1]

    # Call always after cons.runn()!
    if not hasattr(cons, "output"):
        print("Consortium object must run before analyzing the output of the run!")
        return

    sim_tsv = pd.read_csv(cons.output, sep = "\t")
    nrow = sim_tsv.shape[0]-1
    table_print = []
    for col in sim_tsv.columns:
        if col == 'time':
            continue
        val0 = sim_tsv[col][0]
        valn = sim_tsv[col][nrow]
        increment = val0-valn
        if increment > 0 and sign != "-":
            increment = "\033[1;32;40m" + str(val0-valn) + "\033[0m"
        elif increment < 0 and sign != "+":
            increment = "\033[1;31;40m" + str(val0-valn) + "\033[0m"
        else:
            continue
        table_print.append([col, val0, valn, increment])
    table_print = sorted(table_print, key = keysort, reverse = True)
    table_print = [("Metabolite", "Initial Value", "Final Value", "Increment")] + table_print
    for tup in table_print:
        for el in tup:
            print(el, end = "\t")
        print("\n")
    return [el[0] for el in table_print]