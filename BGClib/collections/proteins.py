#!/usr/bin/env python3

"""
BGC library: ProteinCollection

Handles collections of individual proteins
"""

from ..protein import BGCProtein
from ..annotations.domain import BGCDomain
from ..annotations.hmm import hmmDB
from ..data.utilities import hmmdbs_without_tc
from pathlib import Path
from subprocess import PIPE, Popen
from multiprocessing import Pool, cpu_count
from Bio import SearchIO
from io import StringIO
from pyhmmer.plan7 import HMMFile

class ProteinCollection:
    """Intended for mass-prediction of domains
    """
    
    def __init__(self, fastafile=None):
        self.name = ""              # Name of the collection
        self.proteins = {}          # key = identifier
        self.proteins_by_pid = {}   # key = protein id. Because we can't ensure
                                    #  uniqueness here, only non-empty, 
                                    #  non-repeated keys will be used, so this
                                    #  will be a subset of self.proteins
        
        if fastafile is not None:
            self.fasta_load(fastafile)
    
    
    def __len__(self):
        return len(self.proteins)

            
    def fasta_load(self, fastafile):
        """
        Tries to load a fasta file
        """
        
        # initial string
        random_string = "l5aRx]3DahVr9bnIGVVaMHwpFhOu6iML#DWDbY7DJxTEhc9F4Vw_ZoFRdcA4"
        
        with open(fastafile) as f:
            header = random_string
            sequence = ""
            for line in f:
                line = line.strip()

                if line == "":
                    continue
                
                if line[0] == ">":
                    if header == "":
                        print("Warning: skipping headerless sequence")
                        print(f" ({sequence[:50]}...)")
                        sequence = ""
                        header = line[1:]
                        continue
                    
                    if header != random_string:
                        if header in self.proteins:
                            print(f"Warning: {header} already in Protein "
                                + "Collection. Skipping...")
                            print(f" ({sequence[:50]}...)")
                        else:
                            protein = BGCProtein()
                            protein.identifier = header
                            protein.sequence = sequence
                            self.proteins[header] = protein
                    sequence = ""
                    header = line[1:]
                else:
                    sequence += line
                    
        if header == "":
            print("Warning: skipping headerless sequence")
            print(f" ({sequence[:50]}...)")
        else:
            if header in self.proteins:
                print(f"Warning: {header} already in Protein Collection. "
                    + "Skipping...")
                print(f" ({sequence[:50]})")
            else:
                protein = BGCProtein()
                protein.identifier = header
                protein.sequence = sequence
                self.proteins[header] = protein
        
        
    # TODO: break work on sets of 3 cpus
    # TODO: evaluate whether hmmsearch is better than hmmscan
    # TODO: TEST!
    # TODO: change to pyhmmer
    # TODO: update data structures of parent clusters
    # def predict_domains(self, hmmdb, domtblout_path="", cpus=1, tc=True, filterdoms=True):
    def predict_domains(self, hmmdb, cpus=1, tc=True, filterdoms=True):
        """
        Uses hmmscan to search the protein sequences for hmm models specified
        
        input:
            hmmdb: contains a list to all hmm databases
            domtblout_path: where to store HMMER's output file (optional)
            cpus: number of cpus to pass to HMMER. Remember it uses cpu + 1
            tc: if true, use the model's Trusted Cutoff score (lower score of all
                true positives)
        """
        domtblout_path = "" # TODO: remove this variable?
        # if domtblout_path != "":
        #     try:
        #         assert not domtblout_path.is_file()
        #     except AssertionError:
        #         exit("BGClib.ProteinCollection.predict_domains: domtblout_path should not be a file")
            
        protein_list = []
        for protein_id in self.proteins:
            protein = self.proteins[protein_id]
            protein_list.append(f">{protein.identifier}\n{protein.sequence}")
                    
        if len(protein_list) == 0:
            return
        
        # TODO: 
        #  * merge all hmmdbs and remove this loop
        #  * Divide the protein collection into groups where each group is cpus // 3
        #  * do "with Pool(groups) as pool"
        for db in hmmdb.db_path_list:
            #command = ['hmmscan', '--cpu', str(cpus), '--noali', '--notextw']
            command = ['hmmscan', '--cpu', str(cpus), '--notextw', '--qformat', 'fasta']
            if tc and db.stem not in hmmdbs_without_tc:
                command.append('--cut_tc')
            if domtblout_path != "":
                path = str(domtblout_path / ("output_" + db.stem + ".domtable"))
                command.extend(['--domtblout', path ])
            dbpath = str(db)
            command.extend([ dbpath, '-'])
            print(" ".join(command))

            proc_hmmscan = Popen(command, shell=False, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            out, err = proc_hmmscan.communicate(input="\n".join(protein_list).encode("utf-8"))
            
            # "SearchIO parses a search output file's contents into a hierarchy of four 
            # nested objects: QueryResult, Hit, HSP, and HSPFragment"
            # https://biopython.org/DIST/docs/api/Bio.SearchIO-module.html
            # https://biopython.org/DIST/docs/api/Bio.SearchIO.HmmerIO-module.html
            results = SearchIO.parse(StringIO(out.decode()), 'hmmer3-text')
            #TODO: parallelize this loop
            for qresult in results:
                # if the original header of the fasta sequence contains a space,
                # this thing takes everything after the space as the "description"
                # and the first part as the "id" but probably this isn't what 
                # we want
                if qresult.description != "<unknown description>":
                    seq_identifier = f"{qresult.id} {qresult.description}"
                else:
                    seq_identifier = qresult.id
                                
                for hit in qresult:
                    #print(hit.description)
                    for hsp in hit:
                        hspf = hsp[0] # access to HSPFragment

                        #seq_identifier = qresult.id
                        hmm_id = hit.id
                        try:
                            AC = hmmdb.ID_to_AC[hmm_id]
                        except KeyError:
                            AC = ""

                        try:
                            DE = hmmdb.ID_to_DE[hmm_id]
                        except KeyError:
                            DE = ""
                        
                        try:
                            alias = hmmdb.alias[hmm_id]
                        except KeyError:
                            alias = ""
                        
                        # NOTE: these are 1-based indices but biopython takes
                        # care of this and shifts automatically (e.g. hmmscan
                        # position 1 transforms to array index 0)
                        ali_from = hspf.query_start
                        ali_to = hspf.query_end
                        hmm_from = hspf.hit_start
                        hmm_to = hspf.hit_end
                        #env_from = hsp.env_start
                        #env_to = hsp.env_end
                        
                        # NOTE hspf.query is a SeqRecord object
                        algn_seq = hspf.query.seq
                        
                        # NOTE: The alignment of the query sequence vs the hmm
                        # model can be obtained directly from the hmmscan output.
                        # In order to get same-sized alignments (allowed amino 
                        # acids + deletions in the query), we need to remove all
                        # insertions (of the query, with the model as reference. 
                        # marked as lowercase symbols in the query).
                        # For this we could use the '#=GC RF' line from a proper
                        # Stockholm alignment (e.g. using hmmalign) as a mask. 
                        # hmmscan does not report that line as part of its output
                        # but we can use the "hit" sequence to find the positions
                        # of the insertions (from the perspective of the hmm with
                        # the query as reference, these will be viewed as deletions
                        # which are marked with periods in the hmm sequence).
                        # A second option would be to filter the lowercase
                        # symbols in the sequence
                        #print(hspf.query.seq)
                        #print(hspf.hit.seq)
                        #print(hspf.aln_span, len(hspf.hit.seq))
                        
                        hit_seq = hspf.hit.seq
                        ## this is probably super expensive...
                        #hit_seq_split = hit_seq.split(".")
                        
                        #x = "xxx.x.xxxx.x...xx".split(".")
                        #print([len(a) for a in x])
                        #sys.exit()
                        
                        
                        hmm_size = len(hit_seq)
                        
                        Evalue = hsp.evalue
                        score = hsp.bitscore
                        
                        domain = BGCDomain(self.proteins[seq_identifier], \
                                           hmm_id, AC, DE, alias, \
                                            ali_from, ali_to, 
                                           hmm_from, hmm_to, hmm_size, \
                                           score, Evalue, algn_seq)
                        
                        self.proteins[seq_identifier].domain_list.append(domain)
        
        with Pool(cpus) as pool:
            for p in self.proteins:
                protein = self.proteins[p]
                if filterdoms:
                    pool.apply_async(protein.filter_domains())
                protein.attempted_domain_prediction = True
            pool.close()
            pool.join()
    
        return
    

    def predict_domains_(self, hmmdb):
        """Searches protein sequences against databases of hmm models
        """

        # Prepare databses:
        hmms = []
        for db in hmmdb.db_path_list:
            with HMMFile(db) as hmm_file:
                hmms.append(hmm_file.read())


    def classify_proteins(self, cpus=1):
        """
        Once domains have been predicted, classify proteins from the whole 
        collection
        
        Note that directly classifying a protein will not update BGC.CBPcontent,
         BGC.CBPtypes and BGC.CBPtypes_set of its parent BGC
        """
        
        with Pool(cpus) as pool:
            for protein in self.proteins.values():
                pool.apply_async(protein.classify_sequence())
            pool.close()
            pool.join()
        return
    
    
    def get_fasta(self, sort=True):
        """
        Returns a string with the fasta sequences of all proteins in the
        collection. Header will be the protein accession, if present. Otherwise
        the protein's identifier will be used. 
        """
        sequences = []
        if sort:
            seq_list = sorted(self.proteins.keys())
        else:
            seq_list = self.proteins.keys()
        for p_id in seq_list:
            p = self.proteins[p_id]
            header = f"{p.identifier} ProteinId:{p.protein_id} GeneId:{p.gene}"

            sequences.append(f">{header}\n{p.sequence80()}")

        return "".join(sequences)
