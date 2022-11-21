from PairwiseAlignment import PairwiseAlignment
from MyAlign import MyAlign
from MySeq import MySeq
from SubstMatrix import SubstMatrix

class MultipleAlignment():
    def __init__(self, seqs, alignseq):
        self.seqs = seqs
        self.alignpars = alignseq
    
    def add_seq_alignment (self, alignment, seq):
        res = []
        for i in range(len(alignment.listseqs)+1):
            res.append("")
        cons = MySeq(alignment.consensus(),alignment.al_type)
        self.alignpars.needleman_Wunsch(cons, seq) # AQUÍ
        align2 = self .alignpars.recover_align()
        orig = 0
        for i in range( len (align2)):
            if align2[0,i] == "-":
                for k in range( len (alignment.listseqs)):
                    res[k] += "-"
            else:
                for k in range( len (alignment.listseqs)):
                    res[k] += alignment[k,orig]
                orig+=1
        res[ len (alignment.listseqs)] = align2.listseqs[1]
        return MyAlign(res, alignment.al_type)
    
    def align_consensus(self):
        self .alignpars.needleman_Wunsch(self.seqs[0], self.seqs[1],{},-1)
        res = self .alignpars.recover_align()
        for i in range (2, len(self.seqs)):
            res = self.add_seq_alignment(res, self .seqs[i])
        return res

if __name__ == "__main__":
    s1 = MySeq("ATAGC")
    s2 = MySeq("AACC")
    s3 = MySeq("ATGAC")
    sm = SubstMatrix()
    sm.create_submat(1, -1,"ACGT")
    aseq = PairwiseAlignment(sm, -1)
    ma = MultipleAlignment([s1,s2,s3], aseq)
    al = ma.align_consensus()
    print (al)