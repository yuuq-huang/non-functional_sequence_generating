from Bio.SeqUtils import MeltingTemp
from time         import time
from past.builtins import xrange
import re
import seqfold
from subprocess import Popen, PIPE,STDOUT
from swscore import swscore

class Checker(object):

    # Setup synthesis assessment variables
    comp_table = str.maketrans('ATGC', 'TACG')
# =============================================================================
#     runs_tuple = (('CCCCCCCCC',            'GGGGGGGGG'),
#                   ('AAAAAAAAAAAAA',        'TTTTTTTTTTTTT'),
#                   ('TCTTCTTCTTCTTCTTCT',   'TCGTCGTCGTCGTCGTCG',   'TCCTCCTCCTCCTCCTCC',   'GGAGGAGGAGGAGGAGGA',   'GGTGGTGGTGGTGGTGGT',
#                    'GGCGGCGGCGGCGGCGGC',   'AATAATAATAATAATAAT',   'AAGAAGAAGAAGAAGAAG',   'AACAACAACAACAACAAC',   'CCACCACCACCACCACCA',
#                    'CCTCCTCCTCCTCCTCCT',   'CCGCCGCCGCCGCCGCCG',   'ATAATAATAATAATAATA',   'ATTATTATTATTATTATT',   'ATGATGATGATGATGATG',
#                    'ATCATCATCATCATCATC',   'AGAAGAAGAAGAAGAAGA',   'AGTAGTAGTAGTAGTAGT',   'AGGAGGAGGAGGAGGAGG',   'AGCAGCAGCAGCAGCAGC',
#                    'ACAACAACAACAACAACA',   'ACTACTACTACTACTACT',   'ACGACGACGACGACGACG',   'ACCACCACCACCACCACC',   'TAATAATAATAATAATAA',
#                    'TATTATTATTATTATTAT',   'TAGTAGTAGTAGTAGTAG',   'TACTACTACTACTACTAC',   'TTATTATTATTATTATTA',   'TTGTTGTTGTTGTTGTTG',
#                    'TTCTTCTTCTTCTTCTTC',   'TGATGATGATGATGATGA',   'TGTTGTTGTTGTTGTTGT',   'TGGTGGTGGTGGTGGTGG',   'TGCTGCTGCTGCTGCTGC',
#                    'TCATCATCATCATCATCA',   'GAAGAAGAAGAAGAAGAA',   'GATGATGATGATGATGAT',   'GAGGAGGAGGAGGAGGAG',   'GACGACGACGACGACGAC',
#                    'GTAGTAGTAGTAGTAGTA',   'GTTGTTGTTGTTGTTGTT',   'GTGGTGGTGGTGGTGGTG',   'GTCGTCGTCGTCGTCGTC',   'GCAGCAGCAGCAGCAGCA',
#                    'GCTGCTGCTGCTGCTGCT',   'GCGGCGGCGGCGGCGGCG',   'GCCGCCGCCGCCGCCGCC',   'CAACAACAACAACAACAA',   'CATCATCATCATCATCAT',
#                    'CAGCAGCAGCAGCAGCAG',   'CACCACCACCACCACCAC',   'CTACTACTACTACTACTA',   'CTTCTTCTTCTTCTTCTT',   'CTGCTGCTGCTGCTGCTG',
#                    'CTCCTCCTCCTCCTCCTC',   'CGACGACGACGACGACGA',   'CGTCGTCGTCGTCGTCGT',   'CGGCGGCGGCGGCGGCGG',   'CGCCGCCGCCGCCGCCGC'),
#                   ('ATATATATATATATATATAT', 'AGAGAGAGAGAGAGAGAGAG', 'ACACACACACACACACACAC', 'TATATATATATATATATATA', 'TGTGTGTGTGTGTGTGTGTG', 'TCTCTCTCTCTCTCTCTCTC',
#                    'GAGAGAGAGAGAGAGAGAGA', 'GTGTGTGTGTGTGTGTGTGT', 'GCGCGCGCGCGCGCGCGCGC', 'CACACACACACACACACACA', 'CTCTCTCTCTCTCTCTCTCT', 'CGCGCGCGCGCGCGCGCGCG'))
# =============================================================================

    def _is_GC_content_pass(self, seq, gc_low, gc_high,start,end):
        """
        >>> synth_obj._is_GC_content_pass(seq='AAAAAAAAAAAAAAAAAAGG', gc_low=0.40, gc_high=0.50,start=15,end=20)
        True
        >>> synth_obj._is_GC_content_pass(seq='AAAAAAAAAAAAAAAGGGGG', gc_low=0.10, gc_high=0.25,start=0,end=20)
        True
        >>> synth_obj._is_GC_content_pass(seq='AAAAAAAAAAAAAAAAGGGG', gc_low=0.10, gc_high=0.25,start=0,end=20)
        True
        """

        gc_count=0.0
        i = start
        ran=end-start
        while i < end:
            if seq[i] in ['G','C']:
                gc_count+=1
            i += 1
        if not gc_low <= (gc_count/ran) <= gc_high:
            return False
        return True


    def _is_melting_temp_pass(self, seq, tm_low, tm_high,start,end):
        """
        >>> synth_obj._is_melting_temp_pass(seq='AAAAAAAAAATTTTTTTTTT', tm_low=36.0, tm_high=73.0,start=0,end=20)
        False
        >>> synth_obj._is_melting_temp_pass(seq='GCGCGCGCGCGCATATATAT', tm_low=36.0, tm_high=73.0,start=0,end=20)
        True
        """
        i=start
        while i<end: 
            sub_seq=seq[start:end]
            #print(MeltingTemp.Tm_NN(sub_seq, dnac1=250.0, dnac2=0.0, saltcorr=7))
            if not tm_low <= MeltingTemp.Tm_NN(sub_seq, dnac1=250.0, dnac2=0.0, saltcorr=7) <= tm_high:
                
                return False
            i+=1
        return True
    
    def _is_forbidden_seq_exist(self,seq,forbid_seq):
        for s in forbid_seq:
            x=re.search(s,seq)
            if x!=None:
                return True
        return False


    def evaluate(self, seq,cons):
        """
        >>> synth_obj.evaluate(seq='AGTCGCAGACAAAAAAAAAAATTTTTTTTTTTCGATGCTAGTCG', constraint={'position':0,'range':20,'tm':[0,0],'GC':[0.3,0.5]})
        False
        >>> synth_obj._is_melting_temp_pass(seq='GCGCGCGCGCGCATATATAT', tm_low=36.0, tm_high=73.0,start=0,end=20)
        True
        """
        gc_low=cons['GC'][0]
        gc_high=cons['GC'][1]
        tm_low=cons['tm'][0]
        tm_high=cons['tm'][1]
        beginning=cons['position']
        ending=beginning+cons['range']
        if gc_low==0 and gc_high ==0:
            gc_low=0
            gc_high=1
        if tm_low==0 and tm_high ==0:
            tm_low=0
            tm_high=len(seq)*4
            
        if self._is_GC_content_pass(seq,gc_low,gc_high,beginning,ending):
            if self._is_melting_temp_pass(seq,tm_low,tm_high,beginning,ending):
                return True
        return False
    
    def inter_fold(self,seq,fcons):
        t=fcons['temp']
        MFE=fcons['inter_MFE']
        FE=seqfold.dg(seq,temp=t)
        #print(FE)
        if FE<MFE:
            return True
        else:
            return False
# =============================================================================
#        RNAfold 
        #inputs='GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC'
#         process = Popen('D:\maker\RNAfold.exe', stdin=PIPE, stdout=PIPE, stderr=STDOUT)
#         command_input=inputs.encode()
#         output = process.communicate(input=command_input)[0]
# =============================================================================
    
    def intra_fold(self,seq,seq_list,fcons):
        MFE=fcons['intra_MFE']
        for s in seq_list:
            inputs=seq+'&'+s
            process = Popen('D:\maker\RNAcofold.exe', stdin=PIPE, stdout=PIPE, stderr=STDOUT)
            command_input=inputs.encode()
            output=process.communicate(input=command_input)[0]
            output=output.decode('utf-8')
            result=output.splitlines()[1]
            FE=float(re.search(' \((.*)\)$',result).group(1))
            #print('MFE',FE)
            if FE<MFE:
                return True
        return False
            
    def SW_score(self,seq,seq_list,max_score):
        #print('begin sw score')
        for s in seq_list:
            sw_score=swscore(seq,s) # you can also choose gap penalties, etc...
            #print('sw score',sw_score)
            if sw_score>max_score:
                return True
        return False
            
        
        
    def check(self,seq,seq_list,constraints):
        if constraints.get('forbid_seq')!=None:
            if self._is_forbidden_seq_exist(seq,constraints['forbid_seq']):
                #print('forbid')
                return False
        if constraints.get('fcons')!=None:
            if self.inter_fold(seq,constraints['fcons']):
                return False
            if seq_list:
                if self.intra_fold(seq,seq_list,constraints['fcons']):
                    return False
        if constraints.get('max_SWscore')!=None:
            #print(constraints['max_SWscore'])
            if self.SW_score(seq,seq_list,constraints['max_SWscore']):
                return False
            
        
        for cons in constraints['cons']:
            if not self.evaluate(seq,cons):
                return False
        return True



if __name__ == '__main__':
    import doctest
    doctest.testmod(extraglobs={'synth_obj': Checker()})
    

    #seq = 'GGGCAGATGTCACTCAGTAGGATTTGATCAGTACGTGCATGATGCTGCGGGCAGATGTCACTCAGTAGGATTTGATCAGTACGTGCATGATGCTGC'*100#*1000
    # seq = 'AGTCGCAGACAAAAAAAAAAATTTTTTTTTTTCGATGCTAGTCG'

    # seq = 'GCGC GACG ATGA CCGG CAGA CCGG TCAT CGTC ATAT'
    # seq = 'GCGCGACGATGACCGGCAGACCGGTCATCGTCATAT'

    # synth_obj = synthesis()
    # t0 = time()
    # print synth_obj.evaluate(seq)
    # print 'Time elapsed: {}'.format(time()-t0)


    # t0 = time()
    # # print synth_obj._rule_GC_content(seq)
    # # print synth_obj._rule_melting_temp(seq)
    # # print synth_obj._rule_slippage(seq)
    # # print synth_obj._rule_runs(seq)
    # # print synth_obj._rule_palindrome(seq)
    # print synth_obj._rule_hairpin(seq)
    # print 'Time elapsed: {}'.format(time()-t0)