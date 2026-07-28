import bisect, gzip
from collections import defaultdict

tetra_gt_list = ["0/0/0/0", "0/0/0/1", "0/0/1/1", "0/1/1/1", "1/1/1/1"]

hexa_gt_list = ["0/0/0/0/0/0", "0/0/0/0/0/1", "0/0/0/0/1/1", "0/0/0/1/1/1", "0/0/1/1/1/1", "0/1/1/1/1/1", "1/1/1/1/1/1"]
         
ref_encoding_dict = {"A": 0, "C": 1, "G": 2, "T": 3}  
   
class BedIndex:
    def __init__(self, bed_path):
        raw = defaultdict(list)
        f = gzip.open(bed_path, "rt") if bed_path.endswith(".gz") else open(bed_path)
        with f:
            for line in f:
                if not line.strip() or line.startswith(("track","browser","#")): 
                    continue
                chrom, s, e, *_ = line.split("\t")
                s, e = int(s) + 1, int(e) + 1  
                if e > s:
                    raw[chrom].append((s, e))
        self.idx = {}
        for c, iv in raw.items():
            iv.sort()
            merged = []
            for s, e in iv:
                if not merged or s > merged[-1][1]:
                    merged.append([s, e])
                else:
                    merged[-1][1] = max(merged[-1][1], e)
            self.idx[c] = ( [s for s,_ in merged], [e for _,e in merged] )

    def contains(self, chrom, pos0):
    
        if chrom not in self.idx: 
            return False
        starts, ends = self.idx[chrom]
        i = bisect.bisect_right(starts, pos0) - 1
        return i >= 0 and pos0 < ends[i]