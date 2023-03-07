from modeller import *
from modeller.automodel import *

env = environ()
env.io.atom_files_directory = ['.']
env.io.hetatm = True

a = automodel(env, alnfile='rat_cow_in.ali', knowns=('cow_in', 'rat_out_N_C_terms_modeling'), sequence='4YBQ:rat_out_fasta', assess_methods=(assess.DOPE))

a.starting_model = 1
a.ending_model = 200
a.make()
