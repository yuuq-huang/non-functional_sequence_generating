# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 15:07:14 2021

@author: RayVe
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 12:49:36 2021

@author: RayVe
"""
from past.builtins import xrange
from Bio.SeqUtils import MeltingTemp
from time        import time, sleep
from glob        import glob
import os, psutil
import forbid_seq

from checker import Checker

import random
import uuid
import math
import textwrap
import trie

import utils
import adder
import makerchecks
import berno




#import edlib

class kmerSetArray(object):

    def __init__(self, size):
        self.unique_kmers = dict()
        self.indexd_kmers = [' '] * size

    def __setitem__(self, index, kmer):
        if not self.indexd_kmers[index] == ' ':
            if self.indexd_kmers[index] in self.unique_kmers:
                self.unique_kmers.pop(self.indexd_kmers[index])
        self.unique_kmers[kmer]  = index
        self.indexd_kmers[index] = kmer

    def __contains__(self, kmer):
        if kmer in self.unique_kmers:
            return True
        else:
            return False

    def index(self, kmer):
        if kmer in self.unique_kmers:
            return self.unique_kmers[kmer]
        return None

    def __len__(self):
        return len(self.unique_kmers)

    def __repr__(self):
        return 'kmerSetArray[{}]'.format(', '.join(map(str, self.indexd_kmers)))

    def __del__(self):
        self.unique_kmers = None
        self.indexd_kmers = None

class SeqMaker(object):

    def __init__(self, seed=None):
        self.iupac_table = {
            'A': {'A'},
            'C': {'C'},
            'G': {'G'},
            'T': {'T'},
            'R': {'A', 'G'},
            'Y': {'C', 'T'},
            'S': {'G', 'C'},
            'W': {'A', 'T'},
            'K': {'G', 'T'},
            'M': {'A', 'C'},
            'B': {'C', 'G', 'T'},
            'D': {'A', 'G', 'T'},
            'H': {'A', 'C', 'T'},
            'V': {'A', 'C', 'G'},
            'N': {'A', 'T', 'G', 'C'}
        }
        self.compl_dict = {'A': set(['T']), 'G': set(['C', 'T']), 'C': set(['G']), 'T': set(['A', 'G'])}
        #self.synthesis  = Synthesis()
        self.checker = Checker()
        self.seed       = seed
        self.kmer_db_id = str(uuid.uuid4())
        self.kmer_db    = None
        self.mem        = None
        self.seed = seed if not seed is None else random.random()
        random.seed(self.seed)

    def _get_meta_seq(self, seq):
        meta_seq = []

        for i, iupac_nt in enumerate(seq):

            try:
                meta_seq.append(set(self.iupac_table[iupac_nt]))
            except:
                raise ValueError(' [X] Invalid IUPAC code at index {} in sequence constraint'.format(i))

        return tuple(meta_seq)

    def _reset_candidate_kmer_set(self, candidate, kmer_set, i):
        candidate[i] = '-'
        kmer_set[i]  = ' '

    def _clear_kmers(self, candidate, kmer_set, i, k):
        j = i
        while j >= k:
            kmer_set[j]  = ' '
            j -= 1
        i = j + 1
        return i

    def _freeze_point(self, candidate, tried_set, kmer_set, i, k):
        j = i
        while j > k:
            self._reset_candidate_kmer_set(candidate, kmer_set, j)
            tried_set[j] = set()
            # tried_set[j] = set(candidate[j])
            j -= 1
        else:
            self._reset_candidate_kmer_set(candidate, kmer_set, j)
        i = j
        return i

    # FUNCTION TO BRANCH BACK TO OPTIMAL LOCATION UPON BOUND
    # Note: Compute these values for all position at the beginning 
    # for the given sequence and structure constraints.
    def _get_rollback_index(self, meta_seq, i, homology):
        roll_back_index = None
        # See if any position in last homology places has potential for change
        j = i-1
        while j > i - homology and j > -1:
            # Case ( or . without fixed based in sequence
            if not (len(meta_seq[j]) == 1):
                roll_back_index = j
                break
            j -= 1
        return roll_back_index

    def _roll_back(self, meta_seq, i, homology, candidate, tried_set, kmer_set, rbi=None):
        # Default Case: Explicit roll back
        if not rbi is None:
            roll_back_index = rbi
            self._freeze_point(candidate, tried_set, kmer_set, i, k=roll_back_index)
            return rbi
        
        return i

    def _get_local_roll_back_index(self, candidate, i, local_model_fn, verbose):
        candidate_str = ''.join(candidate[:i+1])
        # Try to evaluate the local_model_fn on candidate_str
        try:
            outcome = local_model_fn(candidate_str)
            if outcome in [True, False]:
                state, index = outcome, i
            else:
                state, index = outcome
            index = int(index) if not index is None else i
        # Handle exception / report exception
        except Exception as e:
            if int(verbose) > 1:
                print(' [x] Local Model fn. failed to evaluate partial path: {}'.format(candidate_str))
                print(' [x] Exception message: {}'.format(e))
                print(' [+] Continuing path finding ... ')
            return None # No intelligence, no need to roll back
        if state:
            return None # No conflict reported, no need to roll back
        elif not state: # Conflict found
            return min(index, i)
        else:
            if int(verbose) > 1:
                print(' [x] Local Model fn. state return type is not boolean: {}'.format(outcome))
                print(' [+] Continuing path finding ... ')
            return None # No valid status retported

    def _get_non_coding_candidate(self,
        meta_seq,
        homology,
        local_model_fn,
        tr,
        verbose,
        jump=False,
        start_seq=None,
        allow_internal_repeat=False):

        # print jump

        # Setup the data structures
        candidate = ['-'] * len(meta_seq) if not start_seq else list(start_seq)
        tried_set = [set() for _ in xrange(len(meta_seq))]
        kmer_set  = kmerSetArray(size=len(meta_seq))

        # Setup indexing
        i               = 0
        roll_back_count = 0

        # Main backtracking code
        while -1 < i < len(meta_seq):
            # Jumping out of iteration
            if jump and roll_back_count == homology:
                print('None_1')
                candidate = None
                break
            
            # Try to build a candidate
            if candidate[i] == '-':

                # Forward phase - A nucleotide may be chosen
                if tr.check_if_end(''.join(candidate[:i])):
                    tried_set[i]=set(tr.get_node(''.join(candidate[:i-1])))
                else:
                    tried_set[i]=tried_set[i]
                if len(tried_set[i])<=len(meta_seq[i]):
                    #Reset roll back count
                    roll_back_count=0
                    if len(tried_set[i])==len(meta_seq[i]):
                        if len(tried_set[i]) == 4 :
                            self._reset_candidate_kmer_set(candidate, kmer_set, i)
                            i-=1
                        else:
                            candidate[i] = ''.join(list(meta_seq[i]))
                    else:
                        candidate[i] = random.choice(list(meta_seq[i]-tried_set[i]))
                        tried_set[i].add(candidate[i])
                         # Backward phase - Nucleotide choices exhausted, so traceback

# =============================================================================
#                 else:
#                     #Update roll_back_count
#                      print(tried_set[i])
#                      print(i)
#                      roll_back_count += 1
#                      print(roll_back_count)
#                      
#                      # Reset and roll back
#                      tr.insert(''.join(candidate[:i-1]))
#                      tried_set[i] = set()
#                      kmer_set[i]  = ' '
#                      if i > 0:
#                          i = self._roll_back(meta_seq, i, homology, candidate, tried_set, kmer_set)
#                          i = self._roll_back(meta_seq, i-1, homology, candidate, tried_set, kmer_set) + 1
#                          i -= 1
#                          continue
# =============================================================================
# =============================================================================
                    

           

            # See if the local model function is violated
# =============================================================================
#             if local_model_fn:
#                 rbi = self._get_local_roll_back_index(candidate, i, local_model_fn, verbose)
#                 if not rbi is None: # Model function violated, and a traceback location was determined
#                     roll_back_count = 0 if rbi < i else roll_back_count
#                     tr.insert(''.join(candidate[:i-1]))
#                     i = self._roll_back(meta_seq, i, homology, candidate, tried_set, kmer_set, rbi)
#                     continue
# =============================================================================

            # Handle internal and shared repeats
            if i >= homology-1:
                # Get the hmer/rmer/mmer
                hmer = ''.join(candidate[i-homology+1:i+1])
                rmer = utils.get_revcomp(hmer)
                mmer = min(hmer, rmer)

                # Are either of these mers seen previously?
                if mmer in self.kmer_db or rmer in self.kmer_db or (not allow_internal_repeat and (hmer in kmer_set or rmer in kmer_set)):
                    if tr.search(''.join(candidate[:i])):
                        #print('repeat')
                        break
                    else:
                       #print('roll back')
                       #print(''.join(candidate[:i]))
                       tr.insert(''.join(candidate[:i]))
                       #print(tr.check_if_end(''.join(candidate[:i])))
                       self._reset_candidate_kmer_set(candidate, kmer_set, i)
                       i-=1
                       
                    continue
                
                else:
                    kmer_set[i] = hmer
                    #print(kmer_set)

            # Roll forward
            i += 1
            
        # Prepare to return candidate
        if candidate is None or '-' in candidate:
            return None
        else:
            return ''.join(candidate)
            tr.insert(''.join(candidate))
        

    def _get_opt_pass_count(self, synth_opt, global_model_fn):
        opt_criteria_count = 1 # Since candidate must at least be non-repetitive
        if synth_opt:
            opt_criteria_count += 1
        if global_model_fn:
            opt_criteria_count += 1
        return opt_criteria_count

    # Diagnostic function -- Shouldn't trigger on experimental changes
    def _is_non_coding_construction_verified(self, meta_seq, candidate):
        i = 0
        while i < len(candidate):
            if not candidate[i] in meta_seq[i]:
                return False
                #if not candidate[j] in self.compl_dict[candidate[i]]:
                    #'{}th Nucleotide {} not compementary to {}th Nucleotide {}'.format(i, candidate[i], j, candidate[j])
                    #return False
            i += 1
        return True




    def _is_checker_verified(self,candidate,seq_list,constraints):
         return self.checker.check(candidate,seq_list,constraints)
    

    def _get_verified_non_coding_candidate(self, 
        homology,
        constraints,
        meta_seq,
        seq_list,
        synth_opt,
        local_model_fn,
        global_model_fn,
        jump_count,
        fail_count,
        tri,
        verbose,
        abortion,
        allow_internal_repeat=False):

        # Setup counts and variables
        current_jump_count = 0
        current_fail_count = 0
        struct_fail_count  = 0
        synth_fail_count   = 0
        model_fail_count   = 0
        seed_seq           = None
        verified_candidate = None
        opt_pass_count     = self._get_opt_pass_count(synth_opt, global_model_fn)
        
        while not verified_candidate:

            # print current_jump_count, current_fail_count

            # Try to get a non-repetitive candidate
            start_seq = self._get_non_coding_candidate(meta_seq,
                homology,
                local_model_fn,
                tri,
                verbose,
                jump=current_jump_count < jump_count,
                start_seq=seed_seq,
                allow_internal_repeat=allow_internal_repeat)
            candidate = start_seq
            
# =============================================================================
#             c.append(candidate)
#             if len(c) != len(set(c)):
#                 print('Fail')
#             else:
#                 print('Pass')
# =============================================================================
            

            # print candidate
            opt_count = 0

            # If valid candidate then process
            if candidate:
                tri.insert(candidate)
                current_jump_count = 0
                opt_count += 1

                # Diagnostic block -- Shouldn't trigger on experimental changes
                if not self._is_non_coding_construction_verified(meta_seq, candidate):
                    current_fail_count += 1
                    raise Exception('Part construction incorrect due to code/logic changes! Report issue.')
                else:
                    pass
                
                # Synthesis optimization            
                if synth_opt:
                    if self._is_checker_verified(candidate,seq_list,constraints):
                        opt_count += 1
                    else:
                        synth_fail_count += 1
                        
                        #print(synth_fail_count)
                
                # Global model optimization
                if global_model_fn:
                    try:
                        outcome = global_model_fn(candidate)
                    except Exception as e:
                        if int(verbose) > 1:
                            print(' [x] Global Model fn. failed to evaluate full path: {}'.format(candidate_str))
                            print(' [x] Exception message: {}'.format(e))
                            print(' [+] Found path accepted ... ')
                            
                    if outcome:
                        opt_count += 1
                    elif not outcome:
                        model_fail_count += 1
                    else:
                        if int(verbose) > 1:
                            print(' [x] Global Model fn. state return type is not boolean: {}'.format(outcome))
                            print(' [+] Found path accepted ... ')
                        model_fail_count += 1

                # Did everything get optimized?
                if opt_count == opt_pass_count:
                    verified_candidate = candidate
                else:
                    current_fail_count += 1

                # Failure count exceeded, terminate
                if current_fail_count == fail_count:
                    break

            # No candidate produced
            else:
                # No jumps made yet no non-repetitive candidate found
                if current_jump_count >= jump_count:
                    break
                # Increase current_jump_count 
                else:
                    current_jump_count += 1
                # Abortion limit reached?
                if abortion and current_jump_count >= jump_count:
                    break

        if int(verbose) > 1:
            print('\n  [seq fails] {}, [struct fails] {}, [synth fails] {}, [global fails] {}, [opt fails] {}'.format(current_jump_count,struct_fail_count,synth_fail_count,model_fail_count,current_fail_count))
        seq_list.append(verified_candidate)
        return verified_candidate, current_jump_count+1, current_fail_count+1

    def _get_non_coding_nrps(self,
        i,
        homology,
        constraints,
        seq,
        seq_list,
        target,
        synth_opt,
        local_model_fn,
        global_model_fn,
        jump_count,
        fail_count,
        verbose,
        abortion,
        allow_internal_repeat=False):

        # Setup structures
        seq         = seq.upper().replace('U', 'T')
        meta_seq    = self._get_meta_seq(seq=seq.upper().replace('U', 'T'))
        
        pre_seq=trie.Trie()
        
        seq_count  = 0
        iter_count = 0
        time_sum   = 0.0
        begin_time = time()
        break_flag = False

        # Setup Bernoulli Success model
        total_jump_trials    = jump_count
        total_jump_successes = 1
        curr_jump_prob       = berno.get_prob(trials=total_jump_trials, success=total_jump_successes)
        curr_jump_trial      = jump_count
        total_fail_trials    = fail_count
        total_fail_successes = 1
        curr_fail_prob       = berno.get_prob(trials=total_fail_trials, success=total_fail_successes)
        curr_fail_trial      = fail_count

        # Stream parts until completion
        while True:
            
            t0 = time()

            # print 'Estimated Jump Count = {}'.format(curr_jump_trial)
            # print 'Estimated Fail Count = {}'.format(curr_fail_trial)

            candidate, curr_jump_trial, curr_fail_trial = self._get_verified_non_coding_candidate(homology,
                constraints,
                meta_seq,
                seq_list,
                synth_opt,
                local_model_fn,
                global_model_fn,
                curr_jump_trial,
                curr_fail_trial,
                pre_seq,
                verbose,
                abortion,
                allow_internal_repeat)

            # print 'Incurred  Jump Count = {}'.format(curr_jump_trial-1)
            # print 'Incurred  Fail Count = {}'.format(curr_fail_trial-1)

            if candidate is None:
                break_flag = True
                print('no candidate verified')

            # Got a canidate -- will try again
            if not break_flag:
                final_candidate = [candidate]
                update_status   = adder.kmerdb_update(self.kmer_db, final_candidate, homology, self.mem)

                if not update_status: # Memory full
                    if verbose:
                        print('[ERROR] Memory Full ... Breaking Loop')
                        yield update_status

                seq_count  += len(final_candidate)
                time_sum   += time()-t0
                iter_count += 1

                if not self.mem:
                    self.kmer_db.commit()

                if verbose:
                    print('[constr] {}, [parts] {}, [{}-mers] {}, [iter {}] {:.2f}s, [avg] {:.2f}s, [total time] {:.2f}h'.format(i,
                        seq_count,
                        homology,
                        len(self.kmer_db),
                        iter_count,
                        time()-t0, time_sum / iter_count,
                        (time() - begin_time) / 3600.0))
                
                yield final_candidate[0]

                # No more parts required
                if seq_count == target:
                    break
            # No more candidates to build
            else:
                yield candidate
                break

            # Update failure limits based on Bernoulli Success model
            total_jump_trials    += curr_jump_trial
            total_jump_successes += 1
            curr_jump_prob       =  berno.get_prob(trials=total_jump_trials, success=total_jump_successes)
            curr_jump_trial      =  berno.get_trials(prob=curr_jump_prob)
            total_fail_trials    += curr_fail_trial
            total_fail_successes += 1
            curr_fail_prob       =  berno.get_prob(trials=total_fail_trials, success=total_fail_successes)
            curr_fail_trial      =  berno.get_trials(prob=curr_fail_prob)

        if not self.mem:
            self.kmer_db.commit()
            # sleep(0.001)

    def check_maker_contingents(self, seq, target, homology):
        # Sequence Legality
        seq_legal, seq_illegal_chars = makerchecks.is_seq_constr_legal(seq)
        if not seq_legal:
            print('\n [ERROR]    Sequence Constraint is not legal due to chars: {}'.format(seq_illegal_chars))
            print(' [SOLUTION] Try correcting Sequence Constraint')
            return False
        # Sequence Sufficiency -- soft check
        seq_sufficient, constrained_motif_locs = makerchecks.is_seq_constr_sufficient(seq, homology, target)
        if not seq_sufficient:
            print('\n [WARNING]  Target size of {} is unreachable from given Sequence Constraint and Lmax of {}'.format(target, homology-1))
            print(' [WARNING]  >> Overly constrained motifs at locations: {}'.format(constrained_motif_locs))
            print(' [WARNING]  Fewer parts will be generated')
        return True

    def non_coding_maker(self,
        homology,
        constraints,
        seq_list,
        formed_seq,
        target_list,
        background_list=None,
        synth_opt=True,
        local_model_fn=None,
        global_model_fn=None,
        jump_count=10,
        fail_count=1000,
        output_file=None,
        verbose=True,
        abortion=True,
        allow_internal_repeat=False):

        print('\n[Null Sequence Generation]')
        
        forseq=forbid_seq.get_forbid_seq()
        #print(forseq)
        constraints['forbid_seq']=forseq

        # Check Maker Inputs
        build_parts = True
        for i, (seq,target) in enumerate(zip(seq_list,target_list)):
            print('\n[Checking Constraints {}]:'.format(i+1))
            print(' Sequence Constraint : {}'.format(seq))
            print('Length of Sequence:{}'.format(len(seq_list[0])))
            print(" GC content:{}".format(constraints['cons'][0]['GC'][0]),"- {}".format(constraints['cons'][0]['GC'][1]))
            print(" Melting Temperature:{}".format(constraints['cons'][0]['tm'][0]),"- {}".format(constraints['cons'][0]['tm'][1]))
            print(' Target Size         : {}'.format(target))
            print(' Lmax                : {}'.format(homology-1))
            check_status = self.check_maker_contingents(
                seq,
                target,
                homology)
            if check_status == False:
                print('\n Check Status : FAIL')
                build_parts = False
            else:
                print('\n Check Status : PASS')
        if not build_parts:
            print('\nPlease Modify Constraints to Eliminate [ERROR]s Found Above.')
            return {} # Empty Dict ... no parts were generated
        print
        # All Checks Completed

        # Execute Maker
        adder.current_uuid     = self.kmer_db_id
        self.mem, self.kmer_db = adder.init_back_db(background_list, homology, seq_list, target_list, verbose=verbose)
        memory_exhausted       = False

        total_nrp_count = 0

        if output_file is None:
            adder.setup_proj_dir()
            output_file = './{}/seq_list.txt'.format(self.kmer_db_id)

        with open(output_file, 'w') as out_file:
            for i, (seq, target) in enumerate(zip(seq_list,target_list)):
                seq    = seq.upper().replace('U', 'T')
                current_nrp_count = 0

                if verbose:
                    print('\nGenerating Inactive Sequences for Constraints {}:'.format(i+1))

                for non_coding_nrp in self._get_non_coding_nrps(i+1,
                    homology,
                    constraints,
                    seq,
                    formed_seq,
                    target,
                    synth_opt,
                    local_model_fn,
                    global_model_fn,
                    jump_count,
                    fail_count,
                    verbose,
                    abortion,
                    allow_internal_repeat):

                    if non_coding_nrp:
                        out_file.write('>constr {} | non-repetitive part {}\n'.format(i+1, current_nrp_count+1))
                        non_coding_nrp = '\n'.join(textwrap.wrap(non_coding_nrp, 80))
                        out_file.write('{}\n'.format(non_coding_nrp))
                        current_nrp_count += 1
                        total_nrp_count += 1
                    else:
                        if non_coding_nrp is None:
                            print('Failure Limits Exceeded or k-mers Exhausted. Cannot Generate more sequences.')
                        else:
                            print('Memory Capacity at Full. Cannot Build More Parts.')
                            memory_exhausted = True
                print('Construction Complete.\n')
                # Memory no longer available ... stop
                if memory_exhausted:
                    break

        # Pack output in dictionary
        parts_dict = {}
        for i,line in enumerate(utils.stream_fasta_seq_list(output_file)):
            line = line.strip()
            parts_dict[i] = line
       
        # Cleanups and return
        adder.remove_proj_dir()
        print('Number of Sequence Generated: {}'.format(current_nrp_count))
        return parts_dict

def main():
    homology = 16
    seq         = 'TAATT TCNNN NNNNN NNNNN NNNNN NNNNN NNNNN NNNNN NNNNN NNNNN' # Nielsen Like
    seq         = ''.join(seq.split(' '))
    output_file = None #'riboz.fa'
    background_list =utils.get_fasta_seq_list(fasta_filename='D:\maker\GCF_000005845.2_ASM584v2_genomic.fna')
    sm_obj = SeqMaker()
    seq_list=[]
    gc_low=0.5
    gc_high=0.5
    tm_low=70
    tm_high=72
    temp=37
    intra_MFE=-21.4
    inter_MFE=-11.0
    SW_score=90
    constraints={'cons':[{'position':0,'range':len(seq),'tm':[tm_low,tm_high],'GC':[gc_low,gc_high]}],'fcons':{'intra_MFE':intra_MFE,'inter_MFE':inter_MFE,'temp':temp},'max_SWscore':SW_score}   

    # Hammerhead Nielsen Paper

    t0 = time()
    sm_obj = SeqMaker(seed=7)
    toolbox = sm_obj.non_coding_maker(homology,constraints, [seq], seq_list,[10],
        background_list=background_list,
        jump_count=100,
        fail_count=1048576,
        synth_opt=True,
        verbose=True,
        abortion=True,
        allow_internal_repeat=False,
        output_file=output_file)
    print(toolbox)

    print('Wall Time {} sec'.format(time()-t0))
    print(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
if __name__ == '__main__':
    main()

