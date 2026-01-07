# complexity.py

class Complexity:
    def __init__(self):
        self.depth = 0
        self.cmult = 0
        self.pmult = 0
        self.add = 0
    
    def insert_value(self, depth, cmult, pmult, add):
        self.depth = depth
        self.cmult = cmult
        self.pmult = pmult
        self.add = add
        
    def print_params(self):
        print(f"{'DCPA:':<8}{self.depth}|{self.cmult}|{self.pmult}|{self.add}")
        # print(f"Depth:\t{self.depth}")
        # print(f"CMult:\t{self.cmult}")
        # print(f"PMult:\t{self.pmult}")
        # print(f"Add:\t{self.add}")
        
def attach(c1: Complexity, c2: Complexity, attach_type: str) -> Complexity:
    res = Complexity()
    if attach_type == 'x': # 곱셈
        res.depth = max(c1.depth, c2.depth) + 1
        res.cmult = c1.cmult + c2.cmult + 1
        res.pmult = c1.pmult + c2.pmult
        res.add = c1.add + c2.add
        
    elif attach_type == '+': # 덧셈
        res.depth = max(c1.depth, c2.depth)
        res.cmult = c1.cmult + c2.cmult
        res.pmult = c1.pmult + c2.pmult
        res.add = c1.add + c2.add + 1
        
    return res
    
    
def compare(c1: Complexity, c2: Complexity) -> int:
    if c1.depth != c2.depth:
        return 1 if c1.depth < c2.depth else 2
    if c1.cmult != c2.cmult:
        return 1 if c1.cmult < c2.cmult else 2
    if c1.pmult != c2.pmult:
        return 1 if c1.pmult < c2.pmult else 2
    if c1.add != c2.add:
        return 1 if c1.add < c2.add else 2
    return 0