"""
Class representing BNX file fragment length distribution.
"""
class BNX:
    def __init__(self, filepath):
        self.path = filepath
        self.diffs = self.collectD()
        self.dlen = len(self.diffs)
        self.p0 = self.distP(self.diffs, self.dlen)

    """
    Method collects all fragment lengths in input bnx file.
    """
    def collectD(self):
        f = open(self.path, "r")

        line = f.readline()
        lines = f.readlines()
        diffs = []
        cnt = 0
        for line in lines:
            #handle line
            if line.startswith("#"):
                line = f.readline()
                continue

            if cnt % 4 == 1:
                values = line.strip().split("\t")
                positions = values[1:]

                for i in range(len(positions)-2):
                    diff = int(float(positions[i+1]) - float(positions[i]))

                    diffs.append(diff)

            cnt += 1

        f.close()

        return diffs

    """
    Computes probability distribution of fragment lengths.
    """
    def distP(self, diffs, dlen):
        imax = 0

        for diff in diffs:
            if diff > imax:
                imax = diff

        f = (imax+1)*[0]

        for diff in self.diffs:
            f[diff] += 1

        p = (imax+1)*[0]

        for i in range(imax+1):
            p[i] = f[i] / dlen

        return p
