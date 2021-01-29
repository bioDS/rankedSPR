__author__ = 'Lena Collienne'


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def plot_hist(d, filehandle):
    # Shows and saves histogram (to specified file)
    df = pd.DataFrame(data = d)
    sns.histplot(data=df)
    plt.savefig(filehandle)
    plt.show()

def plot_dots(d, filehandle):
    # Shows and saves values (to specified file)
    plt.plot(distances_rf)
    plt.savefig(filehandle)
    plt.show()
