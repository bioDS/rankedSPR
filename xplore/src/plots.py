__author__ = 'Lena Collienne'


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def plot_hist(d, bins, filehandle = ''):
    # Shows and saves histogram (to specified file)
    df = pd.DataFrame(data = d)
    sns.histplot(data=df, bins = bins, stat = 'density', legend = False)
    if filehandle != '':
        plt.savefig(filehandle)
    plt.tight_layout()
    plt.show()

def plot_dots(d, ylimits = None, filehandle = '', lgnd = False, line = False):
    # Shows and saves values (to specified file)
    if line == False:
        p = sns.scatterplot(data=d, s = 50, legend = lgnd)
    else:
        p = sns.lineplot(data = d, legend = lgnd)
    if ylimits != None:
        p.set(ylim=(ylimits[0],ylimits[1]))
    # p = (p.set_axis_labels("Iteration","Distance").set(ylim=(0,1)))
    if filehandle != '':
        plt.savefig(filehandle)
    plt.tight_layout()
    plt.show()
