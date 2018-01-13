### NGS Pool Target Profile Checker ###

import pandas as pd
from Bio.SeqUtils import GC
import numpy as np


from twistdb.db.utils import get_handle
from twistdb.models import *
 
db = get_handle(engine_url="postgresql://twister:PASSWORD@db-qa/twistdb251")


translation_table = str.maketrans('ATGC', 'TACG')
def reverse_complement(sequence):
    """Quickly return reverse complement without using BioPython."""
    return str.translate(sequence, translation_table)[::-1]


tep = db.query(OrderItem).get("OI_5a3afa56975aea0007cbc5c9")

# Alternative to do this in a loop over all NGS TE pools in the DB
# for tep in db.query(NGSPool).all():
    ## do this for each pool


rc_pool = [reverse_complement(x) for x in tep.rc_boosted_pool]


""" Takes a comma seprated string of oligos, splits in to a single column
    Then uses pandas DataFrame to create a data structure out of the oligos with column headings
"""
# print (oligo_str)
# oligo_str = [x.strip() for x in oligos.split(',')]

# seq2 = pd.DataFrame(oligo_str, columns=['oligo'])
seq2 = pd.DataFrame(tep.fwd_boosted_pool + rc_pool, columns=['oligo'])
#print 'seq2', seq2

seq3 = pd.DataFrame(seq2.oligo.value_counts())
#print 'seq3 counts', seq3

seq3.reset_index(level=0, inplace=True)
seq3.columns = ['oligo', 'frequency']
seq3['gc'] = seq3.oligo.apply(GC)
#print 'seq3[gc]', seq3['gc']

""" The next lines of code makes 5% GC bins, then creates a new column in your dataframe that labels each row with its bin
    then calculates the average oligo count frequency by that gc_bin column grouping
"""

### Define Frequency Bins or Target Profiles ###
# This version makes too many bins versus the number of bins that we use operationally
# freq_bins = np.arange(0, 100, 5)

# This creates a  bin structure akin to our default bin structure for boosting
# where anything less than 30% is a single bin and anything more than 75% is a single bin
freq_bins = [30, 35, 40, 45, 50, 55, 60, 65, 70, 75]

# This sets bins as strings - for some reason it worked but I have no idea why
#freq_bins = ['30', '35', '40', '45', '50', '55', '60', '65', '70', '75']

# Divide up each oligo into individual GC bins - each bin is given
# a number value which is awkward labeling but the categorization
# appears to be correct
seq3['gc_bin'] = np.digitize(seq3.gc, bins=freq_bins) 

# We can validate the categorization by looking manually at the number of
# oligos that are in a given bin:
len(seq3.loc[(seq3.gc < 35) & (seq3.gc >= 30),:])
# This gives 6 - which is equal to the number of oligos in the second (because no oligos have <30% GC) bin.
seq3['gc_bin'].value_counts()[1]

len(seq3.loc[(seq3.gc < 35) & (seq3.gc >= 30),:]) == seq3['gc_bin'].value_counts()[1]

# Then we can summarize by grouping by GC bin and looking at summary statistics
# for oligo frequency (i.e. print count) in each bin; these should parallel the
# values in the boosting profile.
# summary_df = seq3.groupby('gc_bin').frequency.mean()

# Target profile for 120-mer boosting
target_profile = pd.Series([1, 1, 1.14, 1.22, 1.29, 1.46, 1.50, 1.55, 1.81, 1.89, 2.44], index=[0,1,2,3,4,5,6,7,8,9,10])

summary_df = pd.concat([seq3.groupby('gc_bin').frequency.mean(),
                        target_profile,
                        seq3['gc_bin'].value_counts()], axis=1)
summary_df.index = ["<30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-100"]
summary_df.columns=['observed_freq', 'target_frequency', 'oligo_count']

summary_df['freq_delta'] = summary_df.observed_freq - summary_df.target_frequency


""" Testing above code script """
#print '\nSeq3 Frequency Chart\n', summary_df
print summary_df


""" Simple nucleotide calculations """
num_of_Ts = oligos.count('T')

num_of_GCs = oligos.count('GC')

num_of_ATs = oligos.count('AT')

num_of_Gs = oligos.count('G')

num_of_As = oligos.count('A')

num_of_Cs = oligos.count('C')


print ('\n')
print "Nucleotide Break Down:", '\n' "A's:", num_of_As, '\n' "T's:", num_of_Ts, '\n' "G's:", num_of_Gs, '\n' "C's:", num_of_Cs, '\n' "GC's:", num_of_GCs, '\n' "AT's:", num_of_ATs
