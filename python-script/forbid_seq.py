# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 23:41:00 2021

@author: RayVe
"""
def get_forbid_seq():
    ForbSeq={'ATG','TTATNCACA','TGTGANNNNNNTCACANT','CTAG','CTAA','CAAA','CAAG',
             'NGCTNAGCN','GGGNNNNNCCC','GGGG','GAATTC','ACTAGT','TCTAGA','CTGCAG',
             'GGTCTC','ACCTGC','CGTCTC','GAAGAC','AAAAA','TTTTT','GGGGG','CCCCC',
             'ATATAT','ACACAC','AGAGAG','TATATA','TCTCTC','TGTGTG','CACACA','CTCTCT',
             'CGCGCG','GAGAGA','GTGTGT','GCGCGC','SSRCGCC','TATAAA','GGCCAATCT',
             'GGGCGG','RNNAUGG','KHGCGTG','GGWWNDDWWGGDBW','TRATYGNATTA','KNWCWBNBTGTWCY'
             'SCCTGRGGS','GSCYGRGGR','YKNRTATTSK','YACGTG','TGAYGYMW','TGAYGTHM','TGAYKY'
             'KSSTGACGTGG','YGCTGAGTCM','THHHDDWWTGASTMW','RVRGGAAVYRRVR','WTTCNWRGRW',
             'KCACGTGM','GVWCACGTGACHS','TTKAWG','WTTTAY','RTTRWGCAAY','RTTGCDYAAY',
             'TKDBGMAAK','TKSGMAAT','MTKKBGYAATBK','YCCCWRGGGA','YWRAGGTCA',
             'MRAGGTCA','YBRCGTCAY','SRVTGACGTSA','GATTA','CCRSHAGRKGGCRS',
             'CCRSHAGGKGGCRS','RVNDATYRRT','TKRYGTAA','RRVKATTGCA','RMTAATTR',
             'GGCGSSRRR','RGCGCGMAWC','SGCGCGRAA','SGCGSGAAR','GCGCCAAA','RGCGGGAR',
             'RHGKGGGYRK','GVGKGGGCGK','RNVAGGAAR','SMGGAAGY','RVSCGGAAG','MMGGAARY',
             'MGGAARY','SCAYGTG','ACCGGAAGT','SAAGGTCA','TYAAGGTCA','AKGTCAVVVWRNYY',
             'RGGDCWVGSTGMCCY','SMGGAWRY','GGARR','CMGGAAG','MGGAAGY','WGAYAAGATAANM',
             'TTTATKRY','MMGGAAAT','CMGGAAGY','RTGACTCA','TGAGTCAY','TGTTTACWYW','TGTTTACWBW',
             'TGTTTGCYYW','TGTTTACWTAVS','TGTGGATT','WBTGATTGGY','TRTTTATYT','TGTTTATKKTT','YTGTTTW',
             'RYBTGTTTWY','TGTTTABK','TGTTTMY','RTTGTTTATKT','SCGGAAG','CCGGAAGYS',
             'KNNNNNNDVWGATAAS','AGATAAS','AGATARVR','WGATAR','AGATAAGR','RGDACWBTVTGTHCY',
             'KYMSTGATT','RSHSWGATT','TGGGDGGTC','STGGGTGGTC','TGGGTGGTY','KGGNMBCAGCTGCGKCY',
             'KGKTGCCM','ACGTGM','SNNVDCGGACGTW','RTKRYGTAAY','RGTTAWTVWTYDHY','RTTAATRDTTDAY',
             'RKDBCAAAGTYCW','KDNCAAAGKBCR','WYATTGATTW','GAANNWWSBRGMR','SRWBVWKSKVGR',
             'SCAGSTG',	'YNRTTTATG','TWTTATKGS','TGATGGATG','WTGATTRATK','TTGATTRATK',
             'TGGGAR','KGYMAGGGGGCA','RRVNGAAASYGAAASY','RAAARYGAAAS','RRAAASBGAAA',
             'RVGGAASWGRRW','RAAABYRAAW','GAAABHGAAACT','CAGGTGY',
             'RTGAGTCAY','TGAGTCA','KMKCGCGAGRBBBVS','WGGGTGNGG','GGGYGKGGY',
             'CAGGGGGTG','CTTTGWW','AWTTAATTA','WWNTGCTGASWHRKCWHWW','RDNWGCTGASTCAGCA',
             'WWNTGCTGASWYRKCM','YGCTGAS','RHCACGTG','GGGMGGRGSVRSRSSVSSSSSS','SGGCCGGMK'
             'SCCGGAG','CTATTTWKAR','YTATTTWWRR','KCTATTTTTAK','TGACAGBTKWM',
             'TGACAGBTKT','CAYGTG','GKRMCGKGTGCADM','CMGTT','MCRCGTG','CACGTGK',
             'GCAGCTG','CASCTGYY','YWTTSWNWTGYWRW','RRCAGRTGDBM','RTGAY','RYGGAAWDWY',
             'GGAAAW','GGAAAAHY','GGAAARHDW','RTGASTCRGCW','TGGSHNGNHKCCW','YYTGGCWBNNDBYCM'
             'TTAYGTAA','GGRAAKYCCC','RRCCAATSRS','KBAAGTG','TCAAGGA','WWTAAGTAWW',
             'RYWAAGTGSR','TAATTG','WWAAVTAGGTCA','RRGGTCANKRWCCY','RRSTYARNNRAGKTCA',
             'RRGSTCAVVVRAGKWYR','RBYVARAGGTCA','RRSBSARAGGKMR','AAKTCAAAARKMM',
             'RAAGGTCA','YCAAGGBY','RAGKTCAAGKTCA','GCGCVTGCGC','MHVSGATTAR','KGTAACKGT','RRCWTGYYWGRRCWWGYY',
             'RCWDGYHKGRRCWYGY','RCRWGYYKSDRCWWGY','RNBCRVYVRDGCRKRRM','WWGATKGAT',
             'KVMHKTGATTGWWK','YSATTGGY','STMATTAR','YBTGTGGTY','WMWTGVAWATW',
             'ATGYAAAW','RTKSTWATGCWW','YWTTVWHATGCADW','RGGBSAVAGGKCR',
             'WBYRGGBVARAGKKBW','RSWTAGRGWMC','RMAGWGAAAGT','RGAACANHBTGY','TAATTR',
             'RRGARANNBHMCASSTGYCM','RGKTSANVVDVRSKKSR','GBDBRVMNVRAGGTCA','RGGGMTTTSM',
             'KGGRNWTTCC','YCAGSACCDHGGACAGYDSY','GTTGCBWBGSRR','SBGTTGCCRKGGVRMCVS',
             'GTTRCCATGGHR','SVBBDGTYKCYRKG','WAABTAGGTCA','AAWVHAGGTCAG','YTGTGGTYW',
             'TGTGGTK','TGTGGTYW','WBDDRRDSARAGKYHA','WSAGGTCA','GTCAAAGGTCA','WGKCTSDCWSYY',
             'GTCTG','WGKCTVDCWSYY','CAGGTG','YYWTTSTBMTKSWDW','CSCTTTGTTYT','ATTGTT',
             'YYATTGTTY','KGGGMGGRR','SVVVVRRRGGCGGRRSBNVVS','KGGGMGGRS','RRAWGRGGAASTGRDR',
             'AAASMGGAAS','RTSRSSWGW','WGGVSWGR','YBBCCWKHYWYRGBY','WTTGTWW','TTCYGRGAA',
             'TTCYTRGAA','YSVYTTCCHRGAARHYHSHK','RSTTTCNBTTYY','TTCCCRGAA','TTCHHNKRAA',
             'SAAGGHYR','SBSTGRGAA','YKBNNNNNNBWGATAA','KATAWAW','RSWWCAAAGS','RMATTCCW',
             'GGGRHTTTCC','CTTTSAWSY','YCAGCTGYVG','MCAGRTGBY','RRWCAYGTG','GTCACGTG','KKHVRNHBVAGGWSM',
             'RVSYVMBVKSAGGTCA','AAVATGGMBRM','RGTCACGTG','RGGKSABNRRGKDSR','GACGTGTMHHW',
             'CAGGTR','AGGCCBVGS','KGGGWGSK','KGGDGBK','GGGKRVY','GCRCMYAGTAGGBVYBY',
             'KMDDKGMAKKMTGGGWRDKS','SGGGVSRMAGNNVKTTGKBKSC','SGAAAAA','KNWCWBNVTGTHCY',
             'GCCYGRGGC','SCCYSRGGC','RTGAYRTMW','KHHHDDWWTGASTMW','AGGTCA',
             'YBRCRTCAY','SBVVSNVDRGGMGGG','RVSCGGAAGBRS','MMGGAARBS','RGGTCABNBYRHYM',
             'RGGKCAVRSYGHCYY','KBNNNNNNVWGATAAR','WGATAAS','AGATAASR','RGDACWBTNTGTHC',
             'RGDNCAAAGKYCW','KDNCAAAGKBCW','TTTTATKRG','RGRAASWGRRM','RRTGASTCA',
             'RBTCTCGCGAGABBBVS','WWHTGCTGASTCAKS','TGCTGASYCARCMHHW','TGCTGASKHDRSR',
             'RVBCACGTG','SGGGMGGGGSVRRRSVSSSSS','CTATTTWKRR','YWTTSWNWTKYWRWK',
             'YTGGCWBNNDBYC','RRCCAATSRR','RRCCAATSRSVDBBVNVNS','RRCWTGYYWGRRCAWGYY',
             'CWDGYHWSVRCWYGY','RDKYRVYVRRGCRKRRM','GABDGRCRKBNVNNNNVVSSM','WTYTGCATR',
             'YWTTVWHATGCWVWY','WBYRGGDVARAGKKBR','RSWTAGRGWMS','RMAGKGAAAGT',
             'YCAGSACCDHGGACAGYKSC','SVBYDGTYRCYRRG','KKDBBDBTGTGGYYK',
             'WVBRRRDSARAGKYBR','SKCTSDCWSYY','YWTTBWHATGYWRAW','RKGGGYGGRR',
             'RRRGGCGGRRS','WRWGRGGAAGTGRDR','CCWYWYHYRG','TTCYYRGAA','YHVSTTYCNBTTYCY',
             'TTCCTRGAA','KBNNNNNNBWGATAA','RMATTCCWK','GGGRHWTTCC','AAVATGGCSSC',
             'CACGTGAC','RBCACGTGAY','YRRTGCATKVTGGGW'}
    
    iupac_table = {
           'A': 'A',
           'C': 'C',
           'U':'T',
           'G': 'G',
           'T': 'T',
           'R': '[AG]',
           'Y': '[CT]',
           'S': '[GC]',
           'W': '[AT]',
           'K': '[GT]',
           'M': '[AC]',
           'B': '[CGT]',
           'D': '[AGT]',
           'H': '[ACT]',
           'V': '[ACG]',
           'N': '[ATGC]'
       }
    trans_seq=[]
    for seq in ForbSeq:
        iupac_seq=''
        for nt in seq:
            iupac_seq=iupac_seq+iupac_table[nt]
        trans_seq.append(iupac_seq)
    
    return trans_seq
        