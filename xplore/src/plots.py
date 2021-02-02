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

def plot_dots(d, filehandle = ''):
    # Shows and saves values (to specified file)
    plt.plot(d, linestyle = 'None', marker = 'o', markersize = 6)
    if filehandle != '':
        plt.savefig(filehandle)
    plt.tight_layout()
    plt.show()
