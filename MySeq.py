class MySeq:
    """Biological sequence class"""
    def __init__( self , seq, seq_type = "DNA"):
        self .seq = seq
        self .seq_type = seq_type

    def print_sequence( self ):
        print ("Sequence: " + self .seq)

    def get_seq_biotype ( self ):
        return self .seq_type

    def show_info_seq ( self ):
        print ("Sequence: " + self .seq + " biotype: " + self .seq_type)

    def count_occurrences( self , seq_search):
        return self .seq.count(seq_search)

    def __len__(self):
        return len ( self .seq)

    def __str__(self):
        return self .seq_type + ":" + self .seq

    def __getitem__(self, n):
        return self .seq[n]

    def __getslice__(self, i, j):
        return self .seq[i:j]

    def alphabet (self):
        if ( self .seq_type=="DNA"): return "ACGT"
        elif ( self .seq_type=="RNA"): return "ACGU"
        elif ( self .seq_type=="PROTEIN"): return "ACDEFGHIKLMNPQRSTVWY"
        else : return None

    def validate (self):
        alp = self .alphabet()
        res = True
        i=0
        while i < len ( self .seq) and res:
            if self .seq[i] not in alp: res = False
            else : i += 1
        return res

    def transcription (self):
        if ( self .seq_type == "DNA"): return MySeq( self .seq.replace("T","U"), "RNA")
        else: return None
    
    def reverse_comp (self):
        if ( self .seq_type != "DNA"): return None
        comp = ""
        for c in self .seq:
            if (c == "A"): comp = "T" + comp
            elif (c == "T"): comp = "A" + comp
            elif (c == "G"): comp = "C" + comp
            elif (c== "C"): comp = "G" + comp
        return MySeq(comp, "DNA")

    def translate_codon (cod):
        """Translates a codon into an aminoacid using an internal
        dictionary with the standard genetic code."""
        tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "TGT":"C", "TGC":"C",
        "GAT":"D", "GAC":"D",
        "GAA":"E", "GAG":"E",
        "TTT":"F", "TTC":"F",
        "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
        "CAT":"H", "CAC":"H",
        "ATA":"I", "ATT":"I", "ATC":"I",
        "AAA":"K", "AAG":"K",
        "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
        "ATG":"M", "AAT":"N", "AAC":"N",
        "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "CAA":"Q", "CAG":"Q",
        "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
        "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
        "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
        "TGG":"W",
        "TAT":"Y", "TAC":"Y",
        "TAA":"_", "TAG":"_", "TGA":"_"}
        if cod in tc: return tc[cod]
        else : return None
    
    def translate (self, iniPos = 0):
        if (self.seq_type != "DNA"): return None
        seq_aa = ""
        for pos in range(iniPos, len(self.seq) -2, 3):
            cod = self .seq[pos:pos+3]
            seq_aa += MySeq.translate_codon(cod)
        return MySeq(seq_aa, "PROTEIN")
    
#s1 = MySeq("MKKVSJEMSSVPYW", "PROTEIN")
#print (s1)
#print ( len (s1))
#print (s1[4])
#print (s1[2:5])
s1 = MySeq("ATGTGATAAGAATAGAATGCTGAATAAATAGAATGACAT")
s2 = MySeq("MKVVLSVQERSVVSLL", "PROTEIN")
print (s1.validate(), s2.validate())
print (s1)
s3 = s1.transcription()
s3.show_info_seq()
s4 = s1.reverse_comp().translate()
s4.show_info_seq()