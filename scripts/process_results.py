# generate some initial dataframes with output results

import os
import argparse
import pandas as pd
import math
import time
import warnings
from decimal import Decimal

from django.db.models import Q, Sum, Count
from django.core.exceptions import ObjectDoesNotExist

from results.models import (
    Protein, 
    Protein, 
    Peptide, 
    FastaProtein, 
    Proteome,
    SpeciesSummary, 
    SpeciesFileSummary
)
from projects.models import Queue, SearchSetting, RunTime
from projects.models import Project, MultiplexLabel

from .run_command import write_debug, settings
from .load_proteomes import load_proteins

def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('queue_id', type=str)
    parser.add_argument('fasta_type', choices=['profile', 'proteome', 'custom'])
    args2 = parser.parse_args(args)
    try:
        queue_id = args2.queue_id
        fasta_type = args2.fasta_type
    except:
        return
        
    process_results(queue_id, fasta_type)

def process_results(queue_id, fasta_type):
    try:
        queue = Queue.objects.get(id=queue_id)
    except ObjectDoesNotExist:
        print("Reporter missing queue_id: %s" % queue_id)
        return False
    
    project = queue.project.name
    job = queue.job
    filename = queue.filename

    delete = Protein.objects.filter(queue=queue).filter(type=fasta_type).delete()
    
    try:
        searchsetting=SearchSetting.objects.get(project=project)
    except ObjectDoesNotExist:
        print("Missing searchsetting for project: %s." % project)
        return False
        
    start = time.time()
    
    if (fasta_type != "profile" and fasta_type != "proteome" and fasta_type != "custom"):
        write_debug("process_results() called with invalid type. Valid types are custon, profile, and proteome." % job, project)
        return False
        
    if (fasta_type == "profile" and queue.status < Queue.Status.PROCESS_RESULTS_PROF):
        write_debug("Project: %s, filename %s, type %s is not ready to process results. Correct to continue." % (project, filename, fasta_type), job, project)
        return False
    elif (fasta_type == "proteome" and queue.status < Queue.Status.PROCESS_RESULTS_PROT):
            write_debug("Project: %s, filename %s, type %s is not ready to process results. Correct to continue." % (project, filename, fasta_type), job, project)
            return False
    elif (fasta_type == "custom" and queue.status < Queue.Status.PROCESS_RESULTS_PROF):
        write_debug("Project: %s, filename %s, type %s is not ready to process results. Correct to continue." % (project, filename, fasta_type), job, project)
        return False

    if queue.skip == True:
        write_debug("Project: %s, filename %s, type %s is set to be skipped. Unskip to continue." % (project, filename, fasta_type), job, project)
        return False
        
    if queue.error >= 1 + settings.max_retries:
        write_debug("Project: %s, filename %s, type %s has an error status of %s or higher. Reset error to continue." % (project, filename, fasta_type, (1 + settings.max_retries)), job, project)
        return False
        
    # select the peptides for the file
    query = (Peptide.objects.filter(queue=queue, type=fasta_type))

    # if we have no results, we can't profile, so we need to skip the
    #  file if the profile method is file 
    if not query:
        write_debug("No peptide results in database for project: %s, filename %s, type %s. Skipping process_results.py." % (project, filename, fasta_type), job, project)
        return True
    
    # build the list of potential protein assignments
    accession_list = {}
    for entry in query:
        accessions = entry.accessions.split(',')
        for accession in accessions:
            if accession in accession_list:
                accession_list[accession] += 1
            else:
                accession_list[accession] = 1

    # we have to load the proteins here to get the species assignments
    if fasta_type == "profile":
        try:
            load_proteins(project, job, "%s%s%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "profile"), os.sep, project, "profile"), accession_list.keys())
        except Exception as e:
            write_debug("Failed load_proteins for %s: %s" % (filename, e), job, project)
            return False
    else:
        try:
            load_proteins(project, job, "%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", filename), os.sep, project, filename, "proteome"), accession_list.keys())
        except Exception as e:
            write_debug("Failed load_proteins for %s: %s" % (filename, e), job, project)
            return False
    
    write_debug("Updating inferences.", job, project)
    
    # now we have the species list, which we'll use to break ties, if possible
    # this is how many times a species could be linked to a peptide
    species_list = {}
    for entry in query:
        accessions = entry.accessions.split(',')
        for accession in accessions:
            fp = FastaProtein.objects.get(accession=accession)
            ppid = fp.ppid.proteome
            if ppid not in species_list:
                species_list[ppid] = 1
            else:
                species_list[ppid] += 1
                
    # now we can loop through and assign them
    rows_list = []
    ratio_labels = []
    for entry in query:
        dict1 = {}
        accessions = entry.accessions.split(',')
        top_count = 0
        top_protein = ""
        top_species_count = 0
        top_species = ""
        for accession in accessions:
            fp = FastaProtein.objects.get(accession=accession)
            ppid = fp.ppid.proteome
            if accession_list[accession] > top_count:
                top_protein = accession
                top_count = accession_list[accession]
                top_species_count = species_list[ppid]
                top_species = ppid
            # we have a tie, so try species
            elif accession_list[accession] == top_count:
                if species_list[ppid] > top_species_count:
                    top_species = ppid
                    top_species_count = species_list[ppid]
                    top_protein = accession
                    # top_count doesn't change
        
        # now build the dataframe
        dict1.update({'sequence':entry.sequence,
                      'accession': top_protein,
                      'val_num_psm':entry.val_num_psm,
                      'validation': entry.validation,
                      'type':entry.type,
                      'peptide_id':entry.id,
                      'peak_area':entry.peak_area,
                      'peak_area_psm':entry.peak_area_psm})
     
        rows_list.append(dict1)
    
    # we may have no data, so just move on
    if len(rows_list) == 0:
        write_debug("No results for queue_id: %s" % queue_id, job, project)
        return   

    peptides = pd.DataFrame(rows_list)
    
    # filter out the ones with no valid psms because they probably have all 0s
    # this shouldn't happen but just in case
    peptides = peptides[peptides.val_num_psm > 0]

    #peptides.to_csv('peptides_list.csv')

    proteins = peptides.groupby(peptides['accession']).aggregate({'val_num_psm':'sum',
                                                                  'peak_area':'sum',
                                                                  'peak_area_psm':'sum'})
    # count how many get merged for peptides
    proteins['val_num_peptide'] = peptides.groupby(peptides['accession']).size()
    
    #proteins.to_csv('proteins_list.csv')
    
    # calculate saf, which is #psm / protein_length
    # then create an entry for each protein
    write_debug("Adding proteins to database.", job, project)
    proteins_to_add = []
    for index, row in proteins.iterrows():
        accession = index
        fp = FastaProtein.objects.get(accession=accession)
        saf = Decimal(row['val_num_psm'] / fp.length)
        if searchsetting.mzmine_run_mzmine == False:
            peak_area = None
            peak_area_psm = None
        else:
            peak_area = Decimal(row['peak_area'])
            peak_area_psm = row['peak_area_psm']
        
        protein = Protein(queue=queue,
                          fp=fp,
                          val_num_psm=float(row['val_num_psm']),
                          val_num_peptide=row['val_num_peptide'],
                          saf=saf,
                          type=fasta_type,
                          peak_area=peak_area,
                          peak_area_psm=peak_area_psm)
        proteins_to_add.append(protein)
        #protein.save()
        if len(proteins_to_add) > 1000:
            Protein.objects.bulk_create(proteins_to_add, ignore_conflicts=True)
            proteins_to_add = []

    Protein.objects.bulk_create(proteins_to_add, ignore_conflicts=True)
    proteins_to_add = []
    
    # update the proteininference link for the peptide now that we know what it is
    write_debug("Updating peptides with protein inference.", job, project)
    for index, row in peptides.iterrows():
        peptide = Peptide.objects.get(id=row['peptide_id'])
        fp = FastaProtein.objects.get(accession=row['accession'])
        protein = Protein.objects.get(fp=fp,
                                      type=fasta_type,
                                      queue=queue)
        peptide.protein = protein
        peptide.save()        
    #proteins.to_csv('proteins.csv')

    end = time.time()
    runtime = end - start
    runtimex = RunTime.objects.get(queue=queue)
    if type == 'profile':
        runtimex.process_results_profile = runtime
    elif type == 'proteome':
        runtimex.process_results_proteome = runtime
    runtimex.save()
    
    return True
    
def calculate_nsaf(project, type):
    print("Updating NSAF values for %s %s" % (project, 
                                                 fasta_type))
    
    if type == "profile":
        status = Queue.Status.FINISHED_PROF
    else:
        status = Queue.Status.FINISHED_PROT
        
    queryq = (
        Queue.objects.filter(project__name=project, 
                             status__gte=status)
                     .exclude(error__gte=(1 + settings.max_retries))
                     .exclude(skip=True)
    )
    
    # now that safs are done, need to calculate nsaf
    for entryq in queryq:
        # here we select collect the sum of all the saf values
        query = (Protein.objects.filter(type=fasta_type)
                                .filter(queue=entryq)
                                .aggregate(sum=Sum('saf'))
                )
        # we can end up in a situation where we have no results for a file
        try:
            saf_sum = Decimal(query['sum'])
        except:
            print("Warning: No results for queue_id: %s" % entryq.id)
            saf_sum = 1
    
        query = (Protein.objects.filter(type=fasta_type)
                                .filter(queue=entryq)
                )

        for entry in query:
            nsaf = Decimal(entry.saf / saf_sum)
            entry.nsaf = nsaf
            entry.save()
            
    print("Finished updating NSAF values for %s %s" % (project, fasta_type))

def calculate_species_summary(project, fasta_type):
    print("Calculating species summary for %s %s" % (project, fasta_type))
    delete = SpeciesSummary.objects.filter(project__name=project,
                                           type=fasta_type).delete()
                            
    # start with the protein info
    query = (Protein.objects.filter(queue__project__name=project)
                            .filter(type=fasta_type)
                            .exclude(queue__skip=True)
                            .exclude(queue__error__gte=(1 + settings.max_retries))
                            .values('fp__ppid__proteome')
                            .annotate(total_psm=Sum('val_num_psm'), 
                                      total_pep=Sum('val_num_peptide'), 
                                      total_pro=Count('type'), 
                                      nsaf=Sum('nsaf'),
                                      peak_area=Sum('peak_area'),
                                      peak_area_psm=Sum('peak_area_psm'))
                            .order_by('-total_psm'))
    for entry in query:
        project = Project.objects.get(name=project)
        proteome = Proteome.objects.get(proteome=entry['fp__ppid__proteome'])
        species = SpeciesSummary(project=project,
                                 type=fasta_type,
                                 ppid=proteome,
                                 nsaf = Decimal(entry['nsaf']),
                                 val_num_protein=entry['total_pro'],
                                 val_num_peptide=entry['total_pep'],
                                 val_num_psm=entry['total_psm'],
                                 peak_area=entry['peak_area'],
                                 peak_area_psm=entry['peak_area_psm'])
        species.save()

def calculate_species_file_summary(q_id, fasta_type):
    queue = Queue.objects.get(id=q_id)
    delete = SpeciesFileSummary.objects.filter(queue=queue,
                                               type=fasta_type).delete()
                                               
    # start with the protein info
    query = (Protein.objects.filter(queue=queue)
                            .filter(type=fasta_type)
                            .values('fp__ppid__proteome')
                            .annotate(total_psm=Sum('val_num_psm'), 
                                      total_pep=Sum('val_num_peptide'), 
                                      total_pro=Count('type'), 
                                      nsaf=Sum('nsaf'),
                                      peak_area=Sum('peak_area'),
                                      peak_area_psm=Sum('peak_area_psm'))
                            .order_by('-total_psm'))
    for entry in query:
        proteome = Proteome.objects.get(proteome=entry['fp__ppid__proteome'])
        species = SpeciesFileSummary(queue=queue,
                                     type=fasta_type,
                                     ppid=proteome,
                                     nsaf = Decimal(entry['nsaf']),
                                     val_num_protein=entry['total_pro'],
                                     val_num_peptide=entry['total_pep'],
                                     val_num_psm=entry['total_psm'],
                                     peak_area=entry['peak_area'],
                                     peak_area_psm=entry['peak_area_psm'])
        species.save()    
