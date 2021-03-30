__author__ = 'Lena Collienne'


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def plot_hist(d, bins, filehandle = '', density = True):
    # Shows and saves histogram (to specified file)
    df = pd.DataFrame(data = d)
    if density == True:
        p = sns.histplot(data=df, bins = bins, stat = 'density', legend = True)
    else:
        p = sns.histplot(data=df, bins = bins, stat = 'count', legend = True)
    p.set_xlabel('distance')
    p.set_ylabel('frequency')
    if filehandle != '':
        plt.savefig(filehandle)
    plt.tight_layout()
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

def mean_comparison():
    mean_list_caterpillar = np.loadtxt('../simulations/distance_distribution/coalescent/mean_to_caterpillar_distance_repeat_n_16_N_10000_500_iterations.np')
    mean_list_focal = np.loadtxt('../simulations/distance_distribution/coalescent/mean_focal_distance_repeat_n_16_N_10000_500_iterations.np')
    mean_list = np.loadtxt('../simulations/distance_distribution/coalescent/mean_distance_repeat_n_16_N_10000_500_iterations.np')

    d = pd.DataFrame(data = list(zip(mean_list, mean_list_caterpillar, mean_list_focal)), columns = ["mean", "mean to caterpillar", "mean from arbitrary focal"])
    sns.scatterplot(data=d, s = 30, legend = True)
    plt.savefig('../simulations/distance_distribution/coalescent/mean_distance_repeat_n_16_N_10000_500_iterations_all+cat+focal.eps')
    plt.show()

    print('min mean list: ', min(mean_list),'\n max mean list: ', max(mean_list))
    print('max mean caterpillar: ', max(mean_list_caterpillar))
    print('min mean focal: ', min(mean_list_focal))