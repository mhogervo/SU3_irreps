class MultiTraceOperator:
    '''
    Class that describes an abstract multi-trace operator. For example

        myOperator = MultiTraceOperator([(0,1),(0,0,0)])

    encodes the operator Tr[XY]Tr[X^3].
    '''
    def __init__(self,listOfSingleTraces):
        self.descr = tuple(sorted(listOfSingleTraces))
        
        letters = countLetters(self.descr)
        lettersAppearing = set(letters.keys())
        if not {0,1}.issuperset(lettersAppearing):
            raise ValueError(f"The multi-trace operator {listOfSingleTraces} contains operators other than [0,1].")
        self.weight = sum(letters.values())
        self.charge = letters.get(0,0) - letters.get(1,0)

def countLetters(listOfSingleTraces):
    '''
    Plug in a multi-trace operator ( sigma_1, sigma_2, ... ) where sigma_i = (a,b,...)
    is a single-trace operator,
    and count how often the different letters appear.
    '''
    dd, dl = dict(), list(listOfSingleTraces)
    while dl:
        singleTrace = list(dl.pop())
        while singleTrace:
            i = singleTrace.pop(0)
            n_i = dd.get(i,0)
            dd[i] = n_i + 1
    return dd
