import numpy as np
import pandas as pd
import argparse, math, sys
import scipy
import scipy.stats
import statsmodels as sm
import statsmodels.stats.multitest
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec

import matplotlib.lines as mlines
import matplotlib.patches as mpatches

# main function
def main():
    
    result = pd.read_table('gcor.txt', delim_whitespace=True)
    ntrait = result.shape[0]
    xvals = np.array(list(range(ntrait)))

    fig, ax0 = plt.subplots(figsize=(15, 4.0))
    ax0.bar(xvals, result['GCOR'], edgecolor='k')

    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.set_xlim(-0.5, ntrait+0.5)
    ax0.set_xticks(xvals)
    ax0.set_xticklabels(result['TRAIT2'].str.replace('_', ' '),
        fontsize=12, rotation=45, ha='right')
   
    ax0.set_ylabel('cross-trait $r_g$', fontsize=16)
    ax0.set_title(r'cross-trait genetic correlation ($r_g$) with COVID-19', fontsize=16)

    plt.savefig('gcor.pdf', bbox_inches='tight')

if(__name__ == '__main__'):
    main()
