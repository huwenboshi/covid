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
    
    result = pd.read_table('ratio.txt', delim_whitespace=True)
    result = result.sort_values(by=['RATIO'], ascending=False)
    result = result.reset_index(drop=True)

    ntissue = result.shape[0]
    xvals = np.array(list(range(ntissue)))

    fig, ax0 = plt.subplots(figsize=(15, 2.0))
    ax0.bar(xvals, result['RATIO'], linewidth=0.75,
        error_kw={'elinewidth':0.75}, yerr=1.96*result['RATIO_SE'],
        edgecolor='k')
    ax0.axhline(y=1.0, color='k', linestyle='--', linewidth=0.75)

    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.set_xlim(-0.75, ntissue)
    ax0.set_ylim(0.95, 1.075)
    ax0.set_xticks(xvals)
    ax0.set_xticklabels(result['TISSUE'].str.replace('_', ' '),
        fontsize=10, rotation=45, ha='right')
   
    ax0.set_ylabel(r'$\chi^2$ enrichment', fontsize=16)
    ax0.set_title(r'enrichment of $\chi^2$ statistic across 53 '
                   'tissue-specifically expressed gene annotations',
                   fontsize=16)

    plt.savefig('ratio.pdf', bbox_inches='tight')

if(__name__ == '__main__'):
    main()
