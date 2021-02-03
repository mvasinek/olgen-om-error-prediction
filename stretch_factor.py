from genome import *

class StretchFactor(OMGenome):
    def __init__(self, gpath, chrlen, pos, sf):
        if gpath != None:
            super().__init__(gpath)

        self.genome_path = gpath
        
        self.chrlen  = chrlen
        self.pos = pos
        
        self.sf = sf

    def processOne(self, queue):
        p1 = []
        for pos_chr in self.pos:    
            last_val = 0        
            p1_c = []
            for i in range(len(pos_chr)):
                new_val = int(pos_chr[i]*self.sf)

                if last_val + 6 < new_val:
                    p1_c.append(new_val)
                
                #pos_chr[i] = new_val

                #if last_val + 6 >= new_val:
                #    print("sf: ", new_val, self.sf, pos_chr[i], pos_chr[i-1])

                last_val = new_val      

            p1.append(p1_c)     

        queue.put(p1)
