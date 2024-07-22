# this just loads proteins into the database ahead of time.
# this can be done optionally to speed up profiling
# proteins aren't loaded by project so this is particularly useful if one is
# running many projects or many files in a single install

import os
import argparse
import regex as re
from Bio import SeqIO

from results.models import (
    FastaProtein,
    Proteome
)

from .run_command import write_debug, settings

def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_type', choices=['profile', 'full'])  
    args2 = parser.parse_args(args)
    
    fasta_type = args2.fasta_type

    load_fastaproteins(fasta_type)

def load_fastaproteins(fasta_type): 
    print("Loading %s proteins into the database. This may take a while." % (fasta_type))
   
    def load_proteins(fasta):
        proteins_to_add = []
        p1 = re.compile(">?[^|]+\|(?P<accession>[^|]+)\|(?P<description>.+)\sOS=(?P<os>.+)\sOX=(?P<ox>[^\s]+)\s(GN=(?P<gn>[^\s]+)\s)?.+UPId=(?P<upid>[^\s]+)\sPPId=(?P<ppid>.+)$")
        p2 = re.compile(">?[^|]+\|(?P<accession>[^|]+)\|?$")
    
        for record in SeqIO.parse(fasta, "fasta"):
            m1 = p1.search(record.description)
            m2 = p2.search(record.description)
            if m1:
                # reverses are effectively a duplicate
                if "_REVERSED" in m1.group('accession'):
                    continue            
                accession = m1.group('accession')
                description = m1.group('description')
                if m1.group('gn'):
                    gene = m1.group('gn')
                else:
                    gene = "unknown"
                upid = m1.group('upid')
                ppid_m = Proteome.objects.get(proteome=m1.group('ppid'))
            elif m2:
                if "_REVERSED" in m2.group('accession'):
                    continue            
                accession = m2.group('accession')
                description = "CRAP"
                upid = "0"
                gene = "CRAP"
                ppid_m = Proteome.objects.get(proteome=0)
            else:
                print("no regexp match for: %s" % record)
                
            fastaprotein = FastaProtein(accession = accession,
                                        description = description,
                                        gene = gene,
                                        ppid = ppid_m,
                                        length = len(record.seq)
                                    )
                                    
            proteins_to_add.append(fastaprotein)
            if len(proteins_to_add) > 5000:
                FastaProtein.objects.bulk_create(proteins_to_add, 
                                                ignore_conflicts=True
                                                )
                proteins_to_add = []
           
        FastaProtein.objects.bulk_create(proteins_to_add, ignore_conflicts=True)
    
    print("Loading bacterial proteins.")
    if fasta_type == 'profile':
        load_proteins(os.path.join(settings.install_folder, "fasta", "profile.fasta"))
    elif fasta_type == 'full':
        load_proteins(os.path.join(settings.install_folder, "fasta", "full.fasta"))
    
    print("Loading CRAP proteins.")
    load_proteins(os.path.join(settings.install_folder, "fasta", "crap.fasta"))
    
    print("Loading human proteins.")
    load_proteins(os.path.join(settings.install_folder, "fasta", "human.fasta"))