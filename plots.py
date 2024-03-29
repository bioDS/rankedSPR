
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_hist(d, xlabel = 'distance', ylabel = 'frequency', filehandle = '', density = True, bwidth = 1):
    # Shows and saves histogram (to specified file)
    plt.clf()
    df = pd.DataFrame(data = d)
    upper_bound = df.max()[0]
    lower_bound = df.min()[0]
    sns.set_theme(font_scale=1.2)
    if density == True:
        sns.histplot(df, palette = ['#b02538'], edgecolor = 'black', alpha=1, binwidth=bwidth, binrange = [ lower_bound - bwidth/2, upper_bound + bwidth/2], stat = 'density', legend = False)
    else:
        sns.histplot(df, palette = ['#b02538'], edgecolor = 'black', alpha=1, binwidth=bwidth, binrange = [ lower_bound - bwidth/2, upper_bound + bwidth/2], stat = 'count', legend = False)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    if filehandle != '':
        plt.savefig(filehandle)
    plt.show()


def plot_dots(d, ylimits = None, filehandle = '', lgnd = False, line = False):
    # Shows and saves values (to specified file)
    if line == False:
        p = sns.scatterplot(data=d, s = 30, legend = lgnd)
    else:
        p = sns.lineplot(data = d, legend = lgnd)
    p.set_xlabel('iteration')
    p.set_ylabel('distance')
    if ylimits != None:
        p.set(ylim=(ylimits[0],ylimits[1]))
    # p = (p.set_axis_labels("Iteration","Distance").set(ylim=(0,1)))
    if filehandle != '':
        plt.savefig(filehandle)
    plt.tight_layout()
    plt.show()
