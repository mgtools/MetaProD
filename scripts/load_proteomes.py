from projects.models import Queue, Project, Setting
from results.models import Proteome, FastaProtein

from .run_command import write_debug, settings

import os
import urllib.request
import csv
import regex as re
from Bio import SeqIO
import argparse
import pandas as pd

def run(*args):
    parser = argparse.ArgumentParser()
    args2 = parser.parse_args(args)
    try:
        load_proteomes()
    except:
        return

def load_proteomes():
    print("Loading proteomes into database (this may take some time).")
    if  not os.path.exists(
            os.path.join(settings.install_folder, 
            "fasta", 
            "proteomes.tsv")):
        print("proteomes.tsv does not exist. Run generate_fasta first.")
        return

    proteomes_to_add = []
    reference_proteomes = []

    proteomes = pd.read_csv(os.path.join(settings.install_folder, "fasta", 
                                         "proteomes.tsv"), sep='\t')
                                         
    proteomes['full size'] = proteomes['full size'].fillna(0)                                         
    proteomes['profile size'] = proteomes['profile size'].fillna(0)
    
    for index, row in proteomes.iterrows():
        if row['PPID'] not in reference_proteomes:
            reference_proteomes.append(row['PPID'])
            proteome = Proteome(proteome=row['PPID'],
                                organism=row['OS'],
                                full_size=row['full size'],
                                profile_size=row['profile size'])           
            proteomes_to_add.append(proteome)                        
            if len(proteomes_to_add) > 5000:
                Proteome.objects.bulk_create(proteomes_to_add, 
                                             ignore_conflicts=True)
                proteomes_to_add = []  
 
    proteome = Proteome(proteome="UP000005640",
                        organism="Homo sapiens"
                       )
    proteomes_to_add.append(proteome)

    proteome = Proteome(proteome="0", 
                        organism="CRAP"
                       )
    proteomes_to_add.append(proteome)

    Proteome.objects.bulk_create(proteomes_to_add, ignore_conflicts=True)
    
    print("Finished loading proteomes.")
    
def load_proteins(project_name, job, fasta, accession_list):
    ''' load the proteins into the database so their info can be accessed '''
    # this is kind of slow because it needs to crawl through a FASTA file
    # an option may be to check to see if the proteins are already loaded so we
    # can skip this step. all we need to do is look for 1 protein that isn't
    # loaded
    
    # other ideas
    # trim accession list to only be ones that are loaded
    # pop accession list as we load them and if it's empty, we're done
    
    write_debug("Loading FASTA proteins into database (this may take some time).", 
        job, project_name)
    
    # check for existing proteins because we can save some time sometimes
    # even though this step is time consuming
    fp = list(FastaProtein.objects.filter(accession__in=accession_list).values_list('accession', flat=True))
    
    accessions_needed = []
    for accession in accession_list:
        if accession not in fp:
            accessions_needed.append(accession)
    
    write_debug("%s proteins need to be added." % (len(accessions_needed)), job, project_name)
    
    try:
        if len(accessions_needed) == 0:
            write_debug("Finished loading proteins.", job, project_name)
            return
            
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
                
            if accession not in accessions_needed:
                continue 
       
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
            
            accessions_needed.remove(accession)
            
            if len(accessions_needed) == 0:
                break
                
        FastaProtein.objects.bulk_create(proteins_to_add, ignore_conflicts=True)
    except Exception as e:
        write_debug("Error loading proteins: %s" % (e), job, project_name)
    else:                
        write_debug("Finished loading proteins.", job, project_name)