from modeller import *
from modeller.automodel import *

env = environ()
env.io.atom_files_directory = ['.']
env.io.hetatm = True

a = automodel(env, alnfile='rGLUT5_in_occ.improvedPymol.ICH4_ICH5.ali', knowns=('4ja3_in_occ', 'rat_out_ICH4_C_term_modeling'), sequence='4YBQ:rat_out_fasta', assess_methods=(assess.DOPE))


a.starting_model = 1
a.ending_model = 100
a.make()
