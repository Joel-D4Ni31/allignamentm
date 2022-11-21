from MyAlign import MyAlign
from MySeq import MySeq
from SubstMatrix import SubstMatrix

class PairwiseAlignment:
    def __init__( self , sm, g):
        self .g = g
        self .sm = sm
        self .S = None
        self .T = None
        self .seq1 = None
        self .seq2 = None

    def score_pos (self, c1, c2, sm, g):
        if c1 == "-" or c2 == "-": return g
        else: return sm[c1 + c2]

    def score_alin (self, seq1, seq2, sm, g):
        res = 0
        for i in range (len(seq1)):
            res += self.score_pos(seq1[i], seq2[i], sm, g)
        return res

    def max3t (v1, v2, v3):
        if v1 > v2:
            if v1 > v3: return 1
            else: return 3
        else:
            if v2 > v3: return 2
            else: return 3
    
    def needleman_Wunsch (self, seq1, seq2, sm, g):
        S = [[0]]
        T = [[0]]
        ## initialize gaps’ row
        for j in range (1, len (seq2)+1):
            S[0].append(g * j)
            T[0].append(3)
            ## initialize gaps’ column
        for i in range (1, len (seq1)+1):
            S.append([g * i])
            T.append([2])
            ## apply the recurrence relation to fill the remaining of the matrix
        for i in range (0, len (seq1)):
            for j in range( len (seq2)):
                s1 = S[i][j] + self.score_pos(seq1[i], seq2[j], sm, g)
                s2 = S[i][j+1] + g
                s3 = S[i+1][j] + g
                S[i+1].append(max(s1, s2, s3))
                T[i+1].append(self.max3t(s1, s2, s3))
        return (S, T)
    
    def recover_align (self, T, seq1, seq2):
        res = ["", ""]
        i = len (seq1)
        j = len (seq2)
        while i>0 or j>0:
            if T[i][j]==1:
                res[0] = seq1[i-1] + res [0]
                res[1] = seq2[j-1] + res [1]
                i -= 1
                j -= 1
            elif T[i][j] == 3:
                res[0] = "-" + res[0]
                res[1] = seq2[j-1] + res [1]
                j -= 1
            else:
                res[0] = seq1[i-1] + res [0]
                res[1] = "-" + res[1]
                i -= 1
        return res
    
    def smith_Waterman (self, seq1, seq2, sm, g):
        S = [[0]]
        T = [[0]]
        maxscore = 0
        for j in range (1, len (seq2)+1):
            S[0].append(0)
            T[0].append(0)
        for i in range (1, len (seq1)+1):
            S.append([0])
            T.append([0])
        for i in range (0, len (seq1)):
            for j in range( len (seq2)):
                s1 = S[i][j] + self.score_pos (seq1[i], seq2[j], sm, g)
                s2 = S[i][j+1] + g
                s3 = S[i+1][j] + g
                b = max(s1, s2, s3)
                if b <= 0:
                    S[i+1].append(0)
                    T[i+1].append(0)
                else:
                    S[i+1].append(b)
                    T[i+1].append(self.max3t(s1, s2, s3))
                    if b > maxscore:
                        maxscore = b
        return (S, T, maxscore)
    
    def max_mat(mat):
        maxval = mat[0][0]
        maxrow = 0
        maxcol = 0
        for i in range (0, len (mat)):
            for j in range (0, len (mat[i])):
                if mat[i][j] > maxval:
                    maxval = mat[i][j]
                    maxrow = i
                    maxcol = j
        return (maxrow,maxcol)
    
    def recover_align_local (self, S, T, seq1, seq2):
        res = ["", ""]
        i, j = self.max_mat(S)
        while T[i][j]>0:
            if T[i][j]==1:
                res[0] = seq1[i-1] + res [0]
                res[1] = seq2[j-1] + res [1]
                i -= 1
                j -= 1
            elif T[i][j] == 3:
                res[0] = "-" + res[0]
                res[1] = seq2[j-1] + res [1]
                j -= 1
            elif T[i][j] == 2:
                res[0] = seq1[i-1] + res [0]
                res[1] = "-" + res[1]
                i -= 1
        return res