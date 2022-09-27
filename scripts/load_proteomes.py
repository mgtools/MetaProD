from projects.models import Queue, Project, Setting
from results.models import Proteome, FastaProtein

from .run_command import write_debug, settings

import os
import urllib.request
import csv
import re
from Bio import SeqIO
import argparse

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
            "ref_proteomes_list.tsv")):
        print("ref_proteomes_list.tsv does not exist. Run generate_fasta first.")
        return

    proteomes_to_add = []
    reference_proteomes = []

    with open(os.path.join(
            settings.install_folder, 
            "fasta", 
            "ref_proteomes_list.tsv"), 'r') as file:
        result = csv.reader(file, delimiter='\t')
        header = next(result)
        # uniprot has been inconsistent with how they capitalize their columns
        header = [x.lower() for x in header]
        for row in result:
            if row[header.index('proteome id')] not in reference_proteomes:
                reference_proteomes.append(row[header.index('proteome id')])
                proteome = Proteome(proteome=row[header.index('proteome id')],
                                    organism=row[header.index('organism')])
                proteomes_to_add.append(proteome)
                if len(proteomes_to_add) > 5000:
                    Proteome.objects.bulk_create(proteomes_to_add, 
                                                 ignore_conflicts=True)
                    proteomes_to_add = []
                    
    # add human proteome
    proteome = Proteome(proteome="UP000005640",
                        organism="Homo sapiens"
                       )
    proteomes_to_add.append(proteome)
    #proteome.save()
   
    # this is for CRAP proteins
    proteome = Proteome(proteome="0", 
                        organism="CRAP"
                       )
    proteomes_to_add.append(proteome)
    #proteome.save()

    Proteome.objects.bulk_create(proteomes_to_add, ignore_conflicts=True)
    
    print("Finished loading proteomes.")

def load_proteins(project_name, job, fasta, accession_list):
    write_debug("Loading proteins into database (this may take some time).", 
        job, project_name)
    
    try:
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

            if accession not in accession_list:
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
            
        FastaProtein.objects.bulk_create(proteins_to_add, ignore_conflicts=True)
    except Exception as e:
        write_debug("Error loading proteins: %s" % (e), job, project_name)
    else:                
        write_debug("Finished loading proteins.", job, project_name)