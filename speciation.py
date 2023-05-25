"""
This code reads a spreadsheet of speciation over a fixed pH range
from PHREEQC and produces a figure that visualizes the HREE/LREE ratio
combinations in all REEs.
"""

import pandas as pd
import matplotlib.pyplot as plt
pd.set_option('display.max_rows', None)

df = pd.read_csv('PHREEQC/ree_speciation.csv')

REEs = ['La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd',
        'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
filtered_REEs = ['La', 'Pr', 'Nd', 'Sm', 'Gd', 'Tb',
                 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
LREEs = ['La', 'Pr', 'Nd', 'Sm']
HREEs = ['Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']

total_count = len(LREEs) * len(HREEs)
count = 1
plt.rcParams["figure.figsize"] = (10, 20)

for hh in HREEs:
    for ll in LREEs:
        plt.subplot(len(HREEs), len(LREEs), count)
        plt.plot(df['pH'], df[f'{hh}CO3+']/df[f'{ll}CO3+'], '-',
                 c='#0d0887', linewidth='3')
        plt.plot(df['pH'], df[f'{hh}(CO3)2-']/df[f'{ll}(CO3)2-'], '-',
                 c='#9c179e', linewidth='3')
        plt.plot(df['pH'], df[f'{hh}3+']/df[f'{ll}3+'], '-',
                 c='#ed7953', linewidth='3')
        plt.xlabel('pH', fontsize=12)
        plt.ylabel(f'{hh}/{ll}', fontsize=12)
        count += 1
"""
plt.legend(['dicarbonate','monocarbonate','free metal ion'], fontsize=12,
            bbox_to_anchor=(0.5, -0.23),fancybox=True, shadow=True, ncol=3)
"""
plt.tight_layout()
plt.savefig("speciation_test.png")
plt.show()
