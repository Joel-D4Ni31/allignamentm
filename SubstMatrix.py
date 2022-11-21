class SubstMatrix:
    def __init__(self):
        self .alphabet = ""
        self .sm = {}

    def __getitem__(self, ij):
        i, j = ij
        return self .score_pair(i, j)

    def score_pair(self, c1, c2):
        if c1 not in self.alphabet or c2 not in self.alphabet:
            return None
        return self .sm[c1+c2]

    def read_submat_file( self , filename, sep):
        sm = {}
        f = open(filename, "r")
        line = f.readline()
        tokens = line.split()
        ns = len(tokens)
        alphabet = []
        for i in range(0, ns):
            alphabet.append(tokens[i][0])
        for i in range(0, ns):
            line = f.readline();
            tokens = line.split();
            for j in range(0, len(tokens) - 1):
                k = alphabet[i] + alphabet[j]
                sm[k] = int(tokens[j + 1])
        return sm

    def create_submat( self , match, mismatch, alphabet):
        sm = {}
        for c1 in alphabet:
            for c2 in alphabet:
                if (c1 == c2):
                    sm[c1 + c2] = match
                else:
                    sm[c1 + c2] = mismatch
        return sm