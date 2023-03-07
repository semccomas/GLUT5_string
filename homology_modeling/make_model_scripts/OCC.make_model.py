from modeller import *
from modeller.automodel import *

env = environ()
env.io.atom_files_directory = ['.']
env.io.hetatm = True

a = automodel(env, alnfile='rGLUT5_occ.improvedPymol.ICH5.ali', knowns=('PfHT', 'rat_out_C_term_modeling'), sequence='4YBQ:rat_out_fasta', assess_methods=(assess.DOPE))

a.starting_model = 1
a.ending_model = 10
a.make()
