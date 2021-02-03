import sys

"""
Class representing BNX file fragment length distribution.
"""
class BNX:
    def __init__(self, filepath):
        self.path = filepath        
        diffs, dists = self.collectData()
        self.dlen = len(diffs)
        self.p0 = self.distP2(diffs, self.dlen)        
        self.mol_dist = BNX.MolsToDist(dists)

        diffs = None
        dists = None        

    """
    Method collects all fragment lengths in input bnx file.
    """
    def collectData(self):
        f = open(self.path, "r")
        
        lines = f.readlines()        
        diffs = []
        dists = []
        cnt = 0     
        skipped = 0 
        for line in lines:
            #handle line
            if line.startswith("#"):
                skipped += 1
                continue
            
            if cnt % 4 == 0:
                v = line.strip().split("\t")
                m_len = float(v[2])                
                dists.append(m_len)   

            if cnt % 4 == 1:
                values = line.strip().split("\t")
                positions = values[1:]

                for i in range(len(positions)-2):
                    diff = int(float(positions[i+1]) - float(positions[i]))
                    diffs.append(diff)

            cnt += 1

        f.close()

        return (diffs,dists)

    @staticmethod
    def ParseSitesNo(fpath):
        f = open(fpath, "r")
        
        lines = f.readlines()
        sites_no = []
        
        cnt = 0
        for line in lines:
            #handle line
            if line.startswith("#"):
                continue
            
            if cnt % 4 == 0:
                v = line.strip().split("\t")                
                sites = int(v[5])                
                sites_no.append(sites)
            
            cnt += 1

        f.close()

        return sites_no

    """
    Computes probability distribution of fragment lengths.
    """
    def distP(self, diffs, dlen):
        imax = 0

        for diff in diffs:
            if diff > imax:
                imax = diff

        f = (25000)*[0]

        for diff in diffs:
            f[diff] += 1

        p = (imax+1)*[0]

        for i in range(imax+1):
            p[i] = f[i] / dlen

        return p

    def distP2(self, diffs, dlen):
        p = 25000*[0]

        for diff in diffs:
            if diff < 25000:
                p[diff] += 1

        for i in range(25000):
            p[i] /= dlen
        
        return p

    @staticmethod
    def CollectMolLengths(fpath):
        f = open(fpath, "r")
        lines = f.readlines()
        f.close()

        ls = []
        cnt = 0
        for line in lines:
            if line.startswith("#"):
                continue

            l = line.strip()
            if len(l) == 0:
                continue

            if cnt % 4 == 0:
                v = line.split("\t")
                m_len = float(v[2])
                ls.append(m_len)

            cnt += 1

        return ls

    @staticmethod
    def CollectMolLengthsSites(fpath):
        f = open(fpath, "r")
        lines = f.readlines()
        f.close()

        ls = []
        cnt = 0
        for line in lines:
            if line.startswith("#"):
                continue

            l = line.strip()
            if len(l) == 0:
                continue

            if cnt % 4 == 0:
                v = line.split("\t")
                m_len = float(v[2])
                sites_no = int(v[5])      
                ls.append((m_len, sites_no))

            cnt += 1

        lim = 0
        for l,s in ls:
            if l>lim:
                lim = int(l)
        
        f = [0]*(lim+1)
        for l,s in ls:
            f[int(l)] += 1

        p = [0]*lim
        for i in range(lim):
            p[i] += f[i] / len(ls)

        fo = open("mol-lens.csv","w")
        for i in range(lim):
            if i % 375 != 0:
                continue
            fo.write("%d\t%.6f\n" % (i,p[i]))
        fo.close()

        lim = 0
        tot1 = 0
        for l,s in ls:
            if s == 1:                
                tot1 += 1
                if l>lim:
                    lim = int(l)
        
        f = [0]*(lim+1)
        for l,s in ls:
            if s == 1:
                f[int(l)] += 1

        p = [0]*lim
        for i in range(lim):
            p[i] += f[i] / tot1

        fo = open("mol-lens-sites-1.csv","w")
        for i in range(lim):
            if i % 375 != 0:
                continue
            fo.write("%d\t%.6f\n" % (i,p[i]))
        fo.close()        


    @staticmethod
    def MolsToDist(mol_lens):
        h = {}

        for l in mol_lens:
            if l not in h:
                h[l] = 1
            else:
                h[l] += 1

        t = []
        for key in h:
            t.append((key, h[key]))

        t.sort()

        return t

    @staticmethod
    def ExpectedMolsPerGenome(L, t):
        denom = 0
        tot = 0
        for l,f in t:
            tot += f

        for l, f in t:
            p = f/tot
            denom += p*l

        return L/denom

    @staticmethod
    def ToBNX(bnx_path, molecules):
        f = open(bnx_path, "w")
        f.write("#Artificially generated BNX file from Olgen OM - Miss Predictor\n")

        for mol in molecules:
            """
            assume four lines
            1 - header -> only third value is non-empty and expresses molecule len
            2 - pos -> array of positions -> last entry mol len
            3,4 - intensities -> empty
            """
            #line 1
            f.write("0\t\t%d\t\t\t%d\t\n" % (mol.length, len(mol.pos)))
            #line 2
            f.write("1")
            for p in mol.pos:
                f.write("\t%.2f" % p)
            #end line 2 empty lines 3,4
            f.write("\n\n\n")
        
        f.close()        

"""
BNXData - simplified space-delimited format to store optical mapping data
- each row contains entries related to one optical map
- columns:
  1 - molecule id
  2 - molecule length
  3 - runID
  4..n - [label positions]
"""
class BNXData(BNX):
    def collectD(self):
        f = open(self.path, "r")
        
        lines = f.readlines()
        diffs = []
        
        for line in lines:
            values = line.strip().split(" ")
            positions = values[4:]

            for i in range(len(positions)-2):
                diff = int(float(positions[i+1]) - float(positions[i]))

                diffs.append(diff)

        f.close()

        return diffs

    def store(self):
        f = open("dist.csv", "w")

        f.write("Fragment Length\tBNX\n")
        
        for i in range(25000):
            f.write("%d\t%.8f\n" % ((i+1), self.p0[i]))
            
        f.close()

    @staticmethod
    def CollectMolLengths(fpath):
        f = open(fpath, "r")
        lines = f.readlines()
        f.close()

        ls = []
        
        for line in lines:
            l = line.strip()
            if len(l) == 0:
                continue
            
            v = line.split(" ")
            
            m_len = float(v[1])
            
            ls.append(m_len)

        return ls

#b = BNX(sys.argv[1])