

"""
@author: Nono Saha Cyrille Merleau, email: nonosaha@mis.mpg.de/csaha@aims.edu.gh
"""
from RNA import hamming_distance, fold_compound




class Landscape (object) :

    def __init__ (self, target) :
        self.target = target

    def fitness (self, structure) :
        if '{' in structure :
            str_ = structure.replace('{','(')
            str_ = str_.replace('}',')')
            return 1./(1+hamming_distance(self.target, str_))

        return 1./(1+hamming_distance(self.target, structure))

    def ens_defect(self, sequence) :

        fc = fold_compound(sequence)
        fc.pf()
        fc.bpp()
        ed = fc.ensemble_defect(self.target)
        return ed
