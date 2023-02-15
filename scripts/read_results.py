# does not load ratios except protein ratios for now

# relies on proteins, peptides, and psms exporting unique numbers of columns i.e. 11 for proteins, 10 for peptides, 9 for psms
# this is a bad way to do it and may be changed later, such as to look at the numbering in the first column 1 vs 1.1 vs 1.1.1

import os
import csv
import re
import sys
from decimal import Decimal
import argparse
import time
import pandas as pd
import statistics
import math
import numpy as np
import warnings

from django.core.exceptions import ObjectDoesNotExist
from django.db.models import Sum, Count

csv.field_size_limit(sys.maxsize)

from projects.models import Queue, Project, Setting, SearchSetting, RunTime
from results.models import (
    Protein,
    Peptide,
    Psm,
    PsmRatio,
    SpeciesSummary,
    SpeciesFileSummary
)

from .run_command import write_debug, settings

def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('queue_id', type=int)
    parser.add_argument('type', choices=['profile', 'proteome'])  
    args2 = parser.parse_args(args)
    
    queue_id = args2.queue_id
    type = args2.type

    read_results(queue_id, type)

# loop through each queue entry where status is peptideshaker, find result
# files, add to results, then update queue status to finished
def read_results(queue_id, type):
    try:
        queue = Queue.objects.get(id=queue_id)
    except ObjectDoesNotExist:
        print("Reporter missing queue_id: %s" % queue_id)
        return False

    project = queue.project.name
    filename = queue.filename
    job = queue.job
    
    try:
        searchsetting=SearchSetting.objects.get(project=project)
    except ObjectDoesNotExist:
        print("Missing searchsetting for project: %s." % project)
        return False

    start = time.time()

    # update this depending on reporter filenames
    if searchsetting.multiplex == False:
        data_file_psm = r"%s%sps_%s_Default_PSM_Report.txt" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, project)
    else:
        data_file_psm = r"%s%sr_%s_Default_PSM_Report.txt" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, project)
        
    if not os.path.exists(data_file_psm):
        write_debug("Missing PeptideShaker or Reporter output.", job, project)
        return False
    
    delete = Protein.objects.filter(queue=queue).filter(type=type).delete()
    delete = Peptide.objects.filter(queue=queue).filter(type=type).delete()
    delete = Psm.objects.filter(queue=queue).filter(type=type).delete()
    delete = PsmRatio.objects.filter(psm__queue=queue).filter(psm__type=type).delete()
    delete = SpeciesSummary.objects.filter(project__name=project, type=type).delete()
    delete = SpeciesFileSummary.objects.filter(queue=queue, type=type).delete()
    
    write_debug("Reading output.", job, project)
    data_psm  = pd.read_csv(data_file_psm, header=0, sep='\t', low_memory=False)
    # peak area
    if searchsetting.mzmine_run_mzmine == True:
        data_psm_pa = pd.read_csv("%s%s%s_mzexport.csv" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, filename), sep=',', low_memory=False)
    
    # no psms, so we should just skip the file as a way to give a warning
    if data_psm.shape[0] <= 1:
        write_debug("%s has no PSMs. Skipping read_results.py." % queue.filename, job, project)
        return True
    
    # load psms first
    if searchsetting.multiplex == True:
        write_debug("Determining multiplexed headers and labels.", job, project)
        # for multiplexed, reporter uses 2 column headers
        # we need to find the position of the labels we want and then
        # the position of the next set of labels so we can count how many labels
        # there are and the positions (i.e. because we can have tmt9, tmt11, tmt10, etc)
        header1 = 'Deisotoped Intensity'
        header2 = 'Ratios'

        # normal column names
        headers1 = data_psm.columns.tolist()
        # label names
        headers2 = data_psm.iloc[0].tolist()

        pos1 = headers1.index(header1)
        pos2 = headers1.index(header2)
        num_labels = pos2 - pos1
        final_pos = pos2 + num_labels
        # get rid of the row with the multiplex labels
        data_psm = data_psm.drop([0])
  
        # replace 0s with None
        # if we have a none, we ultimately will exclude that value from calculations
        # because we don't know if it's not expressed or not detected
        for i in range(pos2, final_pos):
            data_psm.iloc[:, i] = data_psm.iloc[:, i].astype(float)
            data_psm.iloc[:, i] = data_psm.iloc[:, i].replace(0, np.nan)

    psmratio_list = []
    psm_list = []
    write_debug("Reading and saving PSMs (this may take some time).", job, project)
    for index, row in data_psm.iterrows():
        try:
            peak_area = math.log2(data_psm_pa[data_psm_pa['compound_db_identity:compound_name'] == row['Spectrum Title']]['area'].values[0])
        except:
            # set to 0 if mzmine didn't find a peak
            peak_area = None

        psm = Psm(queue=queue,
                  accessions=row['Protein(s)'],
                  sequence=row['Sequence'],
                  mod_sequence=row['Modified Sequence'],
                  variable_ptm=row['Variable Modifications'],
                  fixed_ptm=row['Fixed Modifications'],
                  rt=row['RT'],
                  mz=row['m/z'],
                  error=row['Precursor m/z Error [ppm]'],
                  charge=row['Measured Charge'],
                  validation=row['Validation'],
                  confidence=Decimal(row['Confidence [%]']),
                  title=row['Spectrum Title'],
                  type=type,
                  peak_area = peak_area)
        # move to bulk_create and just loop again for ratios
        psm_list.append(psm)
        if len(psm_list) > 1000:
            Psm.objects.bulk_create(psm_list, ignore_conflicts=True)
            psm_list = []
            
    Psm.objects.bulk_create(psm_list, ignore_conflicts=True)
    psm_list = []
    
    # we have to loop again because we can't determine the psm foreign key
    # until it has actually been inserted
    if searchsetting.multiplex == True:   
        write_debug("Updating PSM ratios (this may take some time).", job, project)
        for index, row in data_psm.iterrows():
            # queue, title, and type is probably enough but we'll add a few more 
            # things just to make sure it finds the right one
            psm = Psm.objects.get(queue=queue,
                                  sequence=row['Sequence'],
                                  mod_sequence=row['Modified Sequence'],
                                  title=row['Spectrum Title'],
                                  type=type)
            # insert the reporter ratio, which is the deisotoped intensity then normalized vs reference channel
            # (entry for that psm divided by entry for reference for that psm)
    
            for i in range(pos2, final_pos):
                ratio = row[i]
                if math.isnan(ratio):
                    ratio = None
                psmratio = PsmRatio(psm=psm, ratio=ratio, label=headers2[i])
                psmratio_list.append(psmratio)
                #psmratio.save()
    
        PsmRatio.objects.bulk_create(psmratio_list, 5000)
        psmratio_list = []
    
    write_debug("Determing peptides from PSM list.", job, project)
    peptides = data_psm.groupby('Modified Sequence')
    peptide_list = data_psm['Modified Sequence'].unique()
    med_ratios = {}

    rows_list1 = []
    peptide_list_add = []
    peak_area_query = (Psm.objects.filter(queue=queue)
                                .filter(peak_area__gt=0)
                                .filter(type=type).values('mod_sequence')
                                .annotate(peak_area_sum=Sum('peak_area'))
                                .annotate(peak_area_count=Count('peak_area')))
    for p in peptide_list:
        dict1 = {}
        # retrieve the group
        # other than the shape, we can retrieve the peptide info from the 
        # first entry
        group = peptides.get_group(p)
        val_num_psm = group.shape[0]
        first_entry = group.iloc[0]
        if searchsetting.mzmine_run_mzmine == True:
            result = next((item for item in peak_area_query if item['mod_sequence'] == first_entry['Modified Sequence']), None)
            try:
                peak_area = result['peak_area_sum']
                peak_area_psm = result['peak_area_count']
            except:
                peak_area = None
                peak_area_psm = None
        else:
            peak_area = None
            peak_area_psm = None
            
        accessions = first_entry['Protein(s)']
        sequence = first_entry['Sequence']
        mod_sequence = first_entry['Modified Sequence']
        dict1.update({'peptide':mod_sequence})
        variable_ptm = first_entry['Variable Modifications']
        fixed_ptm = first_entry['Fixed Modifications']
        # validation is really a peptideshaker measure
        # if there is one confident PSM, validation is Confident, otherwise
        # Doubtful
        exists = 'Confident' in group.Validation.values
        if exists == True:
            validation = 'Confident'
        else:
            validation = 'Doubtful'
        peptide = Peptide(queue=queue,
                          accessions=accessions,
                          sequence=sequence,
                          mod_sequence=mod_sequence,
                          variable_ptm=variable_ptm,
                          fixed_ptm=fixed_ptm,
                          val_num_psm=val_num_psm,
                          validation=validation,
                          type=type,
                          peak_area=peak_area,
                          peak_area_psm=peak_area_psm)
        peptide_list_add.append(peptide)
        if len(peptide_list) > 1000:
            Peptide.objects.bulk_create(peptide_list_add, ignore_conflicts=True)
            peptide_list_add = []
            
        # peptide.save()

    Peptide.objects.bulk_create(peptide_list_add, ignore_conflicts=True)
    peptide_list = []
    # moved ratios to later because we need to know ppid of protein to remove contaminants

    write_debug("Linking the PSMs to a peptide.", job, project)
    # link the psms to a peptide now
    query = (Psm.objects.filter(queue=queue).filter(type=type))
    for entry in query:
        peptide = Peptide.objects.get(queue=queue, mod_sequence=entry.mod_sequence, type=type)
        entry.peptide = peptide
        #entry.save()
    Psm.objects.bulk_update(query, ['peptide'])
    
    end = time.time()
    runtime = end - start
    runtimex = RunTime.objects.get(queue=queue)
    if type == 'profile':
        runtimex.read_results_profile = runtime
    elif type == 'proteome':
        runtimex.read_results_proteome = runtime
    runtimex.save()
    return True
