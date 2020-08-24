"""
Class responsible for representation of reference genome fragment lengths distribution.
- class assumes human genome with 24 chromosomes.
"""
class OMGenome:
    def __init__(self, g_path):
        self.genome_path = g_path
        self.chrlen = 24*[0]
        self.pos = self.parsePositions()
        self.d = self.convertToDist(self.pos)
        self.dlen = len(self.d)
        self.max_d = self.maxDistance(self.d)
        self.f0 = self.distF(self.d, self.max_d)
        self.p0 = self.distP(self.f0, self.max_d, self.dlen)

    """
    Parse positions from reference cmap file.
    """
    def parsePositions(self):
        f = open(self.genome_path, "r")
        lines = f.readlines()
        f.close()

        pos = [[]]
        last_chr = None
        chr_id = 0
        
        for line in lines:
            if line.startswith("#"):
                continue

            line_values = line.strip().split("\t")
            self.chrlen[int(line_values[0])-1] = int(float(line_values[1]))                

            #end of last chromosome marker.
            if int(line_values[3]) > int(line_values[2]):
                if line_values[0] == "24":
                    break
                pos.append([])
                chr_id += 1
                continue

            pos[chr_id].append(int(float(line_values[5])))

        return pos

    """
    Returns genome size
    """
    def genomeSize(self):
        gs = 0
        for csize in self.chrlen:
            gs += csize
        return gs

    """
    Returns number of labels in the genome.
    """
    def labelsCount(self):
        lc = 0
        for chrpos in self.pos:
            lc += len(chrpos)
        return lc

    """
    Converts positions of restriction sites into distribution of fragment lengths.
    """
    def convertToDist(self, positions):
        d = []
        for chr_pos in positions:
            pos_l = chr_pos[0]

            for i in range(1,len(chr_pos)):
                d.append(chr_pos[i] - pos_l)
                pos_l = chr_pos[i]

        return d

    """
    Finds the largest distance between restriction sites.
    """
    def maxDistance(self, ds):
        maximum = 0
        for d in ds:
            if d > maximum:
                maximum = d
        return maximum

    """
    Counts frequency of fragment lengths.
    """
    def distF(self, dist, maxd):
        f0 = [0]*(maxd+1)
        for d in dist:
            f0[d] += 1
        return f0

    """
    Counts probability distribution of fragment lengths.
    """
    def distP(self, f0, maxd, dlen):
        p0 = [0]*(maxd+1)
        for i in range(maxd+1):
            p0[i] = f0[i]/dlen
        return p0

