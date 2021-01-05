import ROOT

class RawD+3+ata:
    """
    read raw data so as to  waveform analysis
    """

    def __init__(self,files3):
        """
        files: Raw data to be processed
        """
        self.event_tree = ROOT.TChain("HitTree")
        for file in files:
            self.event_tree.Add(file)

    def GetEntry(self,n):
        if n >= self.event_tree.GetEntries() or n<0:
            return False
        self.event_tree.GetEntry(n)
        return True

    def GetEntries(self):
        return self.event_tree.GetEntries()

    def GetAttr(self,attr_name):

        value = getattr(self.event_tree,attr_name,None)
        return value  
 
    def Test(self):
        print("successfully ")
        return True
    