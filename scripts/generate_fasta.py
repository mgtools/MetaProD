# downloads proteomes and generates profile and full proteome fasta files

import regex as re
import argparse
import os
import shutil
import gzip
import csv
from pathlib import Path
import urllib.request
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor

from django.db.models import Sum
from django.core.exceptions import ObjectDoesNotExist

from projects.models import Queue, SearchSetting
from results.models import Protein, SpeciesSummary, SpeciesFileSummary
from .run_command import run_command, settings, write_log
from .load_proteomes import load_proteomes

def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('project_name', type=str)
    parser.add_argument('fasta_type', choices=['profile', 'proteome'])
    args2 = parser.parse_args(args)
    project = args2.project_name
    fasta_type = args2.fasta_type

    generate_fasta(project, fasta_type)

# some uniprot code to support the new website
# kept together so if it changes, it's easier to update
import requests
from requests.adapters import HTTPAdapter, Retry

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)
# end uniprot code

def generate_fasta(project, fasta_type):
    write_log("Starting generate_fasta for %s %s." % (project, fasta_type), "fasta", project)

    try:
        searchsetting=SearchSetting.objects.get(project=project)
    except ObjectDoesNotExist:
        write_log("Missing searchsetting for project: %s." % (project), "fasta", project)
        return False

    if searchsetting.profile == False:
        write_log("Generating a FASTA is not supported if profiling is turned off.", "fasta", project)
        write_log("This feature not yet implemented so turn profiling on in the searchsettings.", "fasta", project)
        return 0

    if searchsetting.custom_fasta == True:
        fasta_type = "custom"
        write_log("Project is set to use a custom FASTA.", "fasta", project)
        
    samples = {}
    if fasta_type == "proteome":
        write_log("Including proteomes with at least %s proteins." % (searchsetting.profile_include_above), "fasta", project)
        write_log("Excluding proteomes with fewer than %s proteins." % (searchsetting.profile_exclude_below), "fasta", project)   
    
        if not os.path.exists(os.path.join(settings.data_folder, project, "fasta")):
            os.makedirs(os.path.join(settings.data_folder, project, "fasta"))
        if not os.path.exists(os.path.join(settings.data_folder, project, "fasta", "proteome")):
            os.makedirs(os.path.join(settings.data_folder, project, "fasta", "proteome"))
      
        # list of files in project
        query = (Queue.objects.filter(status__gte=Queue.Status.FINISHED_PROF)
                              .exclude(error__gte=2)
                              .exclude(skip=True)
                              .filter(project__name=project))
        
        # can send them individually
        if searchsetting.profile_type == SearchSetting.ProfileType.FILE:
            write_log("File based profiling.", "fasta", project)
            start = time.time()
        
            for q in query:
                generate_proteome_fasta(project, [q.filename], 0)
                    
            end = time.time()
            runtime = end-start
            print(runtime)
        
        # need to check if a sample has been assigned
        # if no sample, that file has to be processed individually
        # if sample, collect other files with the same sample
        elif searchsetting.profile_type == SearchSetting.ProfileType.SAMPLE:
            write_log("Sample based profiling.", "fasta", project)
            samples = {}
            for q in query:
                # if there's no sample, then it has to be filename
                if q.sample is None:
                    generate_proteome_fasta(project, [q.filename], 0)
                else:
                    if q.sample in samples:
                        samples[q.sample].append(q.filename)
                    else:
                        samples[q.sample] = [q.filename]
            
            if len(samples) > 0:
                for sample in samples:
                    generate_proteome_fasta(project, samples[sample], 1)
                
        elif searchsetting.profile_type == SearchSetting.ProfileType.PROJECT:
            write_log("Project based profiling.", "fasta", project)
            files = []
            for q in query:
                files.append(q.filename)
            
            generate_proteome_fasta(project, files, 0)
            # we generate the FASTA for 1 file and copy it to others
            
    elif fasta_type == "profile":
        # move crap to somewhere accessible without project
        if os.path.exists(os.path.join(os.getcwd(), "scripts", "crap.fasta")):
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "crap.fasta")):
                os.remove(os.path.join(settings.install_folder, "fasta", "crap.fasta"))
            if not os.path.exists(os.path.join(settings.install_folder, "fasta")):
                os.makedirs(os.path.join(settings.install_folder, "fasta"))
            shutil.copy(os.path.join(os.getcwd(), "scripts", "crap.fasta"), os.path.join(settings.install_folder, "fasta", "crap.fasta"))
    
        generate_profile_fasta(project)

        fasta_file_concat = "%s%s%s_profile_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "profile"), os.sep, project)
                
        if not os.path.exists(fasta_file_concat):
            write_log("Missing FASTA file.", "fasta", project)
            write_log("Failed to generate profile FASTA for %s." % (project), "fasta", project)
    
    elif fasta_type == "custom":
        # move crap to somewhere accessible without project
        if os.path.exists(os.path.join(os.getcwd(), "scripts", "crap.fasta")):
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "crap.fasta")):
                os.remove(os.path.join(settings.install_folder, "fasta", "crap.fasta"))
            if not os.path.exists(os.path.join(settings.install_folder, "fasta")):
                os.makedirs(os.path.join(settings.install_folder, "fasta"))
            shutil.copy(os.path.join(os.getcwd(), "scripts", "crap.fasta"), os.path.join(settings.install_folder, "fasta", "crap.fasta"))

        generate_custom_fasta(project)
        
        fasta_file_concat = "%s%s%s_custom_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "custom"), os.sep, project)

        if not os.path.exists(fasta_file_concat):
            write_log("Missing FASTA file.", "fasta", project)
            write_log("Failed to generate profile FASTA for %s." % (project), "fasta", project)

def generate_proteome_fasta(project, filenames, have_sample):

    try:
        searchsetting=SearchSetting.objects.get(project__name=project)
    except:
        write_log("Searchsetting does not exist for %s. Make sure the project and searchsetting are added first." % (project), "fasta", project)
        return 0
 
    write_log("Generating FASTA for %s." % (filenames), "fasta", project)
                     
    banned_species = []
    proteomes = []  
    # these are the proteomes to exclude based on too few proteins
    query_s = (SpeciesFileSummary.objects.filter(fasta_type='profile')
                                         .filter(queue__filename__in=filenames)
                                         .filter(queue__project=project)
                                         .exclude(ppid_id=0)
                                         .exclude(ppid_id='UP000005640')
                                         .values('ppid_id')
                                         .annotate(sum=Sum('val_num_protein'))
                                         .filter(sum__lt=searchsetting.profile_exclude_below)
                                         .values_list('ppid_id', flat=True))
    banned_species = list(query_s)
   
    # sum of all bacterial nsaf for the files being queried
    query_t = (Protein.objects.filter(fasta_type="profile")
                              .filter(queue__filename__in=filenames)
                              .filter(queue__project=project)
                              .exclude(fp__ppid__proteome='0')
                              .exclude(fp__ppid__proteome='UP000005640')
                              .exclude(fp__ppid__proteome__in=banned_species)
                              .aggregate(Sum('nsaf')))
    total = query_t['nsaf__sum']

    if not total is None:
        total = float(total) * float(searchsetting.profile_threshold / 100)
        query = (Protein.objects.filter(fasta_type="profile")
                                .filter(queue__filename__in=filenames)
                                .filter(queue__project=project)
                                .exclude(fp__ppid__proteome='0')
                                .exclude(fp__ppid__proteome='UP000005640')
                                .exclude(fp__ppid__proteome__in=banned_species)
                                .values('fp__ppid__proteome')
                                .annotate(sum=Sum('nsaf'))
                                .order_by('-sum'))
                        
        count = 0
        current_psm = 0       
        for entry in query:
            # running count is below the threshold, so we add it
            if count <= total:
                proteomes.append(entry['fp__ppid__proteome'])
                count += entry['sum']              
                current_psm = entry['sum']
            else:
                if entry['sum'] == current_psm:
                    proteomes.append(entry['fp__ppid__proteome'])
                else:
                    break
            
            # don't bother adding CRAP or any others without headers
            if len(entry['fp__ppid__proteome']) < 1 or entry['fp__ppid__proteome'] == '0':
                continue
                
            if entry['fp__ppid__proteome'] not in proteomes:
                proteomes.append(entry['fp__ppid__proteome'])
        
        # so now we have to check again to look for species that weren't in
        # the proteomes yet but that have proteins above the threshold
        query_s = (SpeciesFileSummary.objects.filter(fasta_type='profile')
                                             .filter(queue__filename__in=filenames)
                                             .filter(queue__project=project)
                                             .exclude(ppid_id=0)
                                             .exclude(ppid_id='UP000005640')
                                             .values('ppid_id')
                                             .annotate(sum=Sum('val_num_protein'))
                                             .filter(sum__gte=(searchsetting.profile_include_above))
                                             .values_list('ppid_id', flat=True))
        approved_species = list(query_s)
        for species in approved_species:
            if species not in proteomes:
                proteomes.append(species)
    
    filename = filenames.pop(0)
    fasta_file = "%s%s%s_%s_%s.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", filename), os.sep, project, filename, "proteome")
    fasta_file_concat = "%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", filename), os.sep, project, filename, "proteome")
    
    if os.path.exists(fasta_file):
        os.remove(fasta_file)

    if os.path.exists(fasta_file_concat):
        os.remove(fasta_file_concat)
 
    if os.path.exists(os.path.join(settings.data_folder, project, "fasta", "proteome", filename)):
        shutil.rmtree(os.path.join(settings.data_folder, project, "fasta", "proteome", filename))
    os.makedirs(os.path.join(settings.data_folder, project, "fasta", "proteome", filename))
            
    if len(proteomes) > 0:
        write_log("Including the proteomes: %s" % (proteomes), "fasta", project)
        for proteome in proteomes:
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta.gz" % proteome)):
                location = "pan"
            elif os.path.exists(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta.gz" % proteome)):
                location = "ref"
            else:
                write_log("Missing FASTA file for %s." % (proteome), "fasta", project)
                return

            with gzip.open(os.path.join(settings.install_folder, "fasta", location, "%s.fasta.gz" % proteome), "rb") as f_in:
                with open(fasta_file, "ab") as f_out:
                    shutil.copyfileobj(f_in, f_out)             

    if searchsetting.use_human == True:
        if not os.path.exists("%s%shuman.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep)):
            write_log("human.fasta does not exist. Remove full.fasta and run generate_fasta again.", "fasta", project)
            return
        with open(fasta_file, "a") as f_out:
            with open(os.path.join(settings.install_folder, "fasta", "human.fasta"), "r") as f_in:
                shutil.copyfileobj(f_in, f_out)
            f_out.write(os.linesep)

    if searchsetting.use_crap == True:
        # copy CRAP to destination file
        with open(fasta_file, "a") as f_out:
            with open(os.path.join(settings.install_folder, "fasta", "crap.fasta"), "r") as f_in:
                shutil.copyfileobj(f_in, f_out)
            f_out.write(os.linesep)
    
    generate_decoy(project, fasta_file)
    
    copy_fasta(project, filename, filenames)

def copy_fasta(project, filename, filenames):
    ''' this just copies the fasta of filename for each of the entries in filename '''
    
    if os.path.exists("%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", filename), os.sep, project, filename, "proteome")):
        for file in filenames:
            write_log("Copying fasta for %s from %s." % (file, filename), "fasta", project)

            if os.path.exists(os.path.join(settings.data_folder, project, "fasta", "proteome", file)):
                shutil.rmtree(os.path.join(settings.data_folder, project, "fasta", "proteome", file))
            os.makedirs(os.path.join(settings.data_folder, project, "fasta", "proteome", file))
    
            if os.path.exists("%s%s%s_%s_%s.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", file), os.sep, project, file, "proteome")):
                os.remove("%s%s%s_%s_%s.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", file), os.sep, project, file, "proteome"))
                
            shutil.copy("%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", filename), os.sep, project, filename, "proteome"),
                       "%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", file), os.sep, project, file, "proteome")) 
                       
def generate_profile_fasta(project):
    try:
        searchsetting=SearchSetting.objects.get(project__name=project)
    except:
        write_log("Searchsetting does not exist for %s. Make sure the project and searchsetting are added first." % (project), "fasta", project)
        return 0
        
    # make sure the fasta dir exists
    if not os.path.exists(os.path.join(settings.data_folder, project, "fasta")):
        os.makedirs(os.path.join(settings.data_folder, project, "fasta"))

    if not os.path.exists(os.path.join(settings.data_folder, project, "fasta", "profile")):
        os.makedirs(os.path.join(settings.data_folder, project, "fasta", "profile"))
        
    if os.path.exists("%s%s%s.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep, "profile")):
        write_log("Profile FASTA already exists. Using existing FASTA.", "fasta", project)
    else:
        write_log("Profile FASTA does not exist. Downloading a new FASTA.", "fasta", project)
        write_log("Note: Uniprot sometimes changes their website. Let the authors know if this does not work anymore.", "fasta", project)
        download_profile_fasta(project)
        create_full_fasta(project)

    fasta_file = "%s%s%s_%s.fasta" % (os.path.join(settings.data_folder, project, "fasta", "profile"), os.sep, project, "profile")
    fasta_file_concat = "%s%s%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "profile"), os.sep, project, "profile")
    
    if os.path.exists(fasta_file):
        os.remove(fasta_file)
 
    if os.path.exists(fasta_file_concat):
        os.remove(fasta_file_concat)
        
    # now copy the base files
    shutil.copy("%s%s%s.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep, "profile"), 
                "%s%s%s_%s.fasta" % (os.path.join(settings.data_folder, project, "fasta", "profile"), os.sep, project, "profile"))

    if searchsetting.use_human == True:
        write_log("Appending human proteome to FASTA file.", "fasta", project)
        if not os.path.exists("%s%shuman.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep)):
            write_log("human.fasta does not exist. Remove profile.fasta and run generate_fasta again.", "fasta", project)
            return
        with open(fasta_file, "a") as f_out:
            with open(os.path.join(settings.install_folder, "fasta", "human.fasta"), "r") as f_in:
                shutil.copyfileobj(f_in, f_out)
            f_out.write(os.linesep)

    if searchsetting.use_crap == True:
        write_log("Appending CRAP database to FASTA file.", "fasta", project)
        # copy CRAP to destination file
        with open(fasta_file, "a") as f_out:
            with open(os.path.join(settings.install_folder, "fasta", "crap.fasta"), "r") as f_in:
                shutil.copyfileobj(f_in, f_out)
            f_out.write(os.linesep)

    success = generate_decoy(project, fasta_file)
    
    if success == True:
        load_proteomes()

        write_log("Done generating profile FASTA.", "fasta", project)
    
def generate_custom_fasta(project):
    try:
        searchsetting=SearchSetting.objects.get(project__name=project)
    except:
        write_log("Searchsetting does not exist for %s. Make sure the project and searchsetting are added first." % (project), "fasta", project)
        return 0
        
    # make sure the fasta dir exists
    if not os.path.exists(os.path.join(settings.data_folder, project, "fasta")):
        os.makedirs(os.path.join(settings.data_folder, project, "fasta"))

    if not os.path.exists(os.path.join(settings.data_folder, project, "fasta", "custom")):
        os.makedirs(os.path.join(settings.data_folder, project, "fasta", "custom"))
        
    if not os.path.exists("%s%s%s.fasta" % (os.path.join(settings.data_folder, project, "fasta"), os.sep, project)):
        write_log("Custom FASTA does not exist.", "fasta", project)
        write_log("Place the custom fasta in: %s%s%s.fasta" % (os.path.join(settings.data_folder, project, "fasta"), os.sep, project), "fasta", project)
        return

    fasta_file = "%s%s%s_%s.fasta" % (os.path.join(settings.data_folder, project, "fasta", "custom"), os.sep, project, "custom")
    fasta_file_concat = "%s%s%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "custom"), os.sep, project, "custom")
    
    if os.path.exists(fasta_file):
        os.remove(fasta_file)
 
    if os.path.exists(fasta_file_concat):
        os.remove(fasta_file_concat)
        
    # now copy the base files
    shutil.copy("%s%s%s.fasta" % (os.path.join(settings.data_folder, project, "fasta"), os.sep, project), 
                "%s%s%s_%s.fasta" % (os.path.join(settings.data_folder, project, "fasta", "custom"), os.sep, project, "custom"))

    if searchsetting.use_human == True:
        write_log("Appending human proteome to FASTA file.", "fasta", project)
        if not os.path.exists("%s%shuman.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep)):
            write_log("human.fasta does not exist. Remove full.fasta and run generate_fasta again.", "fasta", project)
            return
        with open(fasta_file, "a") as f_out:
            with open(os.path.join(settings.install_folder, "fasta", "human.fasta"), "r") as f_in:
                shutil.copyfileobj(f_in, f_out)
            f_out.write(os.linesep)

    if searchsetting.use_crap == True:
        write_log("Appending CRAP database to FASTA file.", "fasta", project)
        # copy CRAP to destination file
        with open(fasta_file, "a") as f_out:
            with open(os.path.join(settings.install_folder, "fasta", "crap.fasta"), "r") as f_in:
                shutil.copyfileobj(f_in, f_out)
            f_out.write(os.linesep)

    success = generate_decoy(project, fasta_file)
    
    # this may not work with custom FASTA files so we may need to add a dummy proteome or something
    #if success == True:
    #    load_proteomes()

    write_log("Done generating custom FASTA.", "fasta", project)
    
# should only need the filename for this as the decoy will end up in the same dir as the filename
def generate_decoy(project, fasta_file):
    ''' generate decoy sequences for FASTA file '''
    write_log("Generating decoy sequences.", "fasta", project)

    if os.path.exists(os.path.join(settings.data_folder, project, "fasta", "temp")):
        shutil.rmtree(os.path.join(settings.data_folder, project, "fasta", "temp"))

    # remake the temp folder
    os.makedirs(os.path.join(settings.data_folder, project, "fasta", "temp", "software"))
    
    if not os.path.exists(os.path.join(settings.install_folder, "software", "SearchGUI-%s" % settings.searchgui_ver, "SearchGUI-%s.jar" % settings.searchgui_ver)):
        write_log("Missing SearchGUI install.", "fasta", project)
        return False
        
    # copy searchgui into temp
    shutil.copytree(os.path.join(settings.install_folder, "software", "SearchGUI-%s" % settings.searchgui_ver), 
                    os.path.join(settings.data_folder, project, "fasta", "temp", "software", "SearchGUI-%s" % settings.searchgui_ver))
    
    run_command(["java", 
                "-cp", os.path.join(settings.data_folder, project, "fasta", "temp", "software", "SearchGUI-%s" % settings.searchgui_ver, "SearchGUI-%s.jar" % settings.searchgui_ver), 
                "eu.isas.searchgui.cmd.FastaCLI",
                "-in", "%s" % fasta_file,
                "-decoy"
                ], "fasta", project)
    
    fasta_file_concat = "%s_concatenated_target_decoy.fasta" % (fasta_file[0:-6])

    if os.path.exists(os.path.join(fasta_file_concat)):
        write_log("Finished generating decoy sequences.", "fasta", project)
        shutil.rmtree(os.path.join(settings.data_folder, project, "fasta", "temp"))
        return True
    else:
        write_log("Failed to generate decoy sequences.", "fasta", project)
        return False

def download_profile_fasta(project):
    ''' downloads zip files '''
    # download ppmembership.txt
    if os.path.exists(os.path.join(settings.install_folder, "fasta", "PPMembership.txt")):
        write_log("PPMembership.txt already exists. Remove to update.", "fasta", project)
    else:
        write_log("Downloading PPMembership.txt from uniprot.", "fasta", project)
        # sometimes need to use different url
        #pp_membership_url = "https://proteininformationresource.org/rps/data/new/PPSeqCurrent/PPMembership.txt"
        pp_membership_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/pan_proteomes/PPMembership.txt"
        urllib.request.urlretrieve(pp_membership_url, os.path.join(settings.install_folder, "fasta", "PPMembership.txt"))
            
    # download reference proteomes
    if os.path.exists(os.path.join(settings.install_folder, "fasta", "ref_proteomes_list.tsv")):
        write_log("ref_proteomes_list.tsv already exists. Remove to update.", "fasta", project)
    else:
        write_log("Downloading list of reference proteomes from uniprot.", "fasta", project)
        #https://www.uniprot.org/proteomes/?query=taxonomy:%22Archaea%20[2157]%22%20OR%20taxonomy:%22Bacteria%20[2]%22&fil=reference%3Ayes&format=list
        # format=tab&columns=id,name
        #ref_proteome_url = "https://www.uniprot.org/proteomes/?query=taxonomy:%22Archaea%20[2157]%22%20OR%20taxonomy:%22Bacteria%20[2]%22&fil=reference%3Ayes&format=tab"
        #urllib.request.urlretrieve(ref_proteome_url, os.path.join(settings.install_folder, "fasta", "ref_proteomes_list.tsv"))
        # uniprot changed the website, so here is the new system
        url = 'https://rest.uniprot.org/proteomes/search?compressed=false&fields=upid%2Corganism%2Corganism_id%2Cprotein_count%2Cbusco%2Ccpd&format=tsv&query=%28%28%28taxonomy_id%3A2%29%20OR%20%28taxonomy_id%3A2157%29%29%20AND%20%28proteome_type%3A1%29%29&size=500'
        progress = 0
        with open(os.path.join(settings.install_folder, "fasta", "ref_proteomes_list.tsv"), 'w') as f:
            for batch, total in get_batch(url):
                lines = batch.text.splitlines()
                if not progress:
                    print(lines[0], file=f)
                for line in lines[1:]:
                    print(line, file=f)
                progress += len(lines[1:])
                print(f'{progress} / {total}')
            
    if not os.path.exists(os.path.join(settings.install_folder, "fasta","pan")):
        os.makedirs(os.path.join(settings.install_folder, "fasta","pan"))
     
    write_log("Generating list of proteomes to download.", "fasta", project)

    member_proteomes = []
    pan_proteomes = []
    ref_proteomes_dl = []
    
    write_log("Including all pan proteomes.", "fasta", project)
    # download pan proteomes first as we include all of them
    with open(os.path.join(settings.install_folder, "fasta", "PPMembership.txt"), 'r') as file:
        i = 1
        result = csv.reader(file, delimiter='\t')
        next(result)
        for row in result:
            # member hasn't been added yet
            if row[1] not in member_proteomes:
                member_proteomes.append(row[1])
                
            # pan proteome hasn't been added yet
            if row[0] not in pan_proteomes:
                pan_proteomes.append(row[0])
    
    write_log("Checking for missing reference proteomes.", "fasta", project)
    ref_count = 0
    # so now we download the pan proteomes but we need to make sure the reference proteomes
    # on uniprot are included as sometimes they won't be a member of a pan proteome if they are
    # the only member of a pan proteome
    with open(os.path.join(settings.install_folder, "fasta", "ref_proteomes_list.tsv"), 'r') as file:
        result = csv.reader(file, delimiter='\t')
        header = next(result)
        header = [x.lower() for x in header]
        for row in result:
            # so this reference is not a member of a pan proteome, so add it as a pan proteome
            if row[header.index('proteome id')] not in member_proteomes:
                ref_count += 1
                member_proteomes.append(row[header.index('proteome id')])
                ref_proteomes_dl.append(row[header.index('proteome id')])
    
    write_log("%s reference proteomes were not a member of a pan-proteome." % (ref_count), "fasta", project)
    # human
    ref_proteomes_dl.append("UP000005640")
        
    # download the pan proteomes
    # https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/pan_proteomes/UP000000212.fasta.gz
    write_log("Downloading needed pan proteomes.", "fasta", project)
    lenp = len(pan_proteomes)
    i = 1
    failed_proteomes = []
    for proteome in pan_proteomes:
        attempts = 10
        for attempt in range(attempts):
            try:
                # if it already exists, don't download it again because uniprot is slow
                # but there could be a problem if the filesize is wrong so might have to check that down the road
                if not os.path.exists(os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta.gz" % proteome)) and not os.path.exists(os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta" % proteome)):
                    write_log("Downloading %s pan-proteome of %s. (%s)" % (i, lenp, proteome), "fasta", project)
                    # fix until uniprot is fixed
                    #proteome_url = "https://proteininformationresource.org/rps/data/new/PPSeqCurrent/" + proteome + ".fasta.gz"
                    proteome_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/pan_proteomes/" + proteome + ".fasta.gz"
                    urllib.request.urlretrieve(proteome_url, os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta.gz" % proteome))
            except Exception as e:
                if attempt == attempts - 1:
                    write_log('Error downloading %s.' % (proteome), "fasta", project)
                    failed_proteoms.append(proteome)
            else:
                break
        i += 1
    
    # now download the reference proteomes remaining
    # put this in another dir so that we know how to change the header later
    # https://www.uniprot.org/uniprot/?query=proteome:UP000027395&format=fasta&compress=yes
    write_log("Downloading needed reference proteomes.", "fasta", project)
    if not os.path.exists(os.path.join(settings.install_folder, "fasta", "ref")):
        os.makedirs(os.path.join(settings.install_folder, "fasta", "ref"))
    i = 1
    lenp = len(ref_proteomes_dl)
    for proteome in ref_proteomes_dl:
        attempts = 10
        for attempt in range(attempts):
            try:
                if not os.path.exists(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta.gz" % proteome)) and not os.path.exists(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta" % proteome)):
                    write_log("Downloading %s reference proteome of %s. (%s)" % (i, lenp, proteome), "fasta", project)
                    proteome_url = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=proteome:" + proteome
                    #proteome_url = "https://www.uniprot.org/uniprot/?query=proteome:" + proteome + "&format=fasta&compress=yes"
                    urllib.request.urlretrieve(proteome_url, os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta.gz" % proteome))
            except Exception as e:
                if attempt == attempts - 1:
                    write_log('Error downloading %s.' % (proteome), "fasta", project)
                    failed_proteoms.append(proteome)
            else:
                break
        i += 1
    
    write_log("The following proteomes failed to download: %s" % (failed_proteomes), "fasta", project)

# we can generate profile.fasta here too
def create_full_fasta(project):
    ''' produces full.fasta with all sequences and no filtering'''
    seq_count = 0

    write_log("Formatting FASTA files and filtering for high profile proteins.", "fasta", project)
    write_log("This can take a long time due to the need to reformat the reference proteomes to match the pan proteome format.", "fasta", project)
 
    if os.path.exists(os.path.join(settings.install_folder, "fasta", "human.fasta")):
        os.remove(os.path.join(settings.install_folder, "fasta", "human.fasta"))
    
    if os.path.exists(os.path.join(settings.install_folder, "fasta", "ref", "temp.fasta")):
        os.remove(os.path.join(settings.install_folder, "fasta", "ref", "temp.fasta"))
    
    if os.path.exists(os.path.join(settings.install_folder, "fasta", "profile.fasta")):
        write_log("profile.fasta already exists. Delete to recreate.", "fasta", project)
        return
        
    accessions = set()

    proteomes = {}
    # pan proteome files
    # unzip pan, look for high profile, then delete unzipped
    write_log("Extracting pan proteomes and filtering for high profile.", "fasta", project)
    p = re.compile("^(?P<start>(?P<db>[^\|]+)\|(?P<accession>[^\|]+)\|(?P<middle>.+)\s)(?P<OS>OS=.+)\s(?P<OX>OX=.+)\s(?P<UPId>UPId=[^\s]+)\s(?P<PPId>PPId=[^\s]+)$")
    for file in sorted(os.listdir(os.path.join(settings.install_folder, "fasta", "pan"))):
        if file.endswith(".fasta.gz"):
            filename = Path(file).stem
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "pan", filename)):
                os.remove(os.path.join(settings.install_folder, "fasta", "pan", filename))
            with gzip.open(os.path.join(settings.install_folder, "fasta", "pan", file), "rb") as f_in:
                with open(os.path.join(settings.install_folder, "fasta", "pan", filename), "ab") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    
            for record in SeqIO.parse(os.path.join(settings.install_folder, "fasta", "pan", filename), "fasta"):           
                m1 = p.search(record.description)         
                if not m1:
                    continue
                    
                accession = m1.group('accession')
                ppid = m1.group('PPId').replace('PPId=','')
                species = m1.group('OS').replace('OS=','')
                
                if accession in accessions:
                    continue
                accessions.add(accession)

                if ppid not in proteomes:
                    proteomes[ppid] = {'OS':species, 'Full Size':1, 'Profile Size':0}
                else:
                    proteomes[ppid]['Full Size'] += 1

                with open(os.path.join(settings.install_folder, "fasta", "pan", "temp.fasta"), "a") as fasta_out:
                    SeqIO.write(record, fasta_out, "fasta")
                        
                m2 = re.search("ribosomal|elongation|chaperon", record.description)              
                if not m2:
                    continue
                    
                proteomes[ppid]['Profile Size'] += 1
                
                with open(os.path.join(settings.install_folder, "fasta", "profile.fasta"), "a") as fasta_out:
                    SeqIO.write(record, fasta_out, "fasta")   
                    
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "pan", filename)):
                os.remove(os.path.join(settings.install_folder, "fasta", "pan", filename))

            with open(os.path.join(settings.install_folder, "fasta", "pan", "temp.fasta"), "rb") as f_in:
                with gzip.open(os.path.join(settings.install_folder, "fasta", "pan", "%s" % file), "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out) 
            
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "pan", "temp.fasta")):
                os.remove(os.path.join(settings.install_folder, "fasta", "pan", "temp.fasta"))
                
    # unzip ref, reformat header, look for high profile
    # once done, zip them all
    write_log("Extracting reference proteomes, reformatting, and filtering for high profile.", "fasta", project)
    p = re.compile("^(?P<start>(?P<db>[^\|]+)\|(?P<accession>[^\|]+)\|(?P<middle>.+))\s(?P<OS>OS=.+)\s(?P<OX>OX=.+)")
    for file in sorted(os.listdir(os.path.join(settings.install_folder, "fasta", "ref"))):
        # exclude human
        if file.endswith(".fasta.gz"):
            filename = Path(file).stem
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "ref", filename)):
                os.remove(os.path.join(settings.install_folder, "fasta", "ref", filename))
            with gzip.open(os.path.join(settings.install_folder, "fasta", "ref", file), "rb") as f_in:
                with open(os.path.join(settings.install_folder, "fasta", "ref", filename), "ab") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                
            for record in SeqIO.parse(os.path.join(settings.install_folder, "fasta", "ref", filename), "fasta"):
                m1 = p.search(record.description)
                if not m1:
                    continue

                accession = m1.group('accession')
                ppid = Path(filename).stem
                species = m1.group('OS').replace('OS=','')                    

                # account for the case where we already ran this on the reference proteomes before
                # we don't want to add part again because it'll already be in m1.group('OX')
                if 'UPId' in record.description or 'PPId' in m1.group('OX'):
                    description = m1.group('start') + " " + m1.group('OS') + " " + m1.group('OX')
                else:
                    description = m1.group('start') + " " + m1.group('OS') + " " + m1.group('OX') + " " + "UPId=" + Path(filename).stem + " " + "PPId=" + Path(filename).stem
                    
                # now write the new fasta header + sequence to a file
                record.description = description
                    
                if file != "UP000005640.fasta.gz":
                    if accession in accessions:
                        continue
                        
                    accessions.add(accession)      
                    
                    if ppid not in proteomes:
                        proteomes[ppid] = {'OS':species, 'Full Size':1, 'Profile Size':0}
                    else:
                        proteomes[ppid]['Full Size'] += 1

                    with open(os.path.join(settings.install_folder, "fasta", "ref", "temp.fasta"), "a") as fasta_out:
                        SeqIO.write(record, fasta_out, "fasta")
                        
                    m2 = re.search("ribosomal|elongation|chaperon", record.description)
                    if not m2:
                        continue
                        
                    proteomes[ppid]['Profile Size'] += 1
                    
                    with open(os.path.join(settings.install_folder, "fasta", "profile.fasta"), "a") as fasta_out:
                        SeqIO.write(record, fasta_out, "fasta")
                                            
                # write human to a different file so we can access it easily later
                else:
                    with open(os.path.join(settings.install_folder, "fasta", "human.fasta"), "a") as fasta_out:
                        SeqIO.write(record, fasta_out, "fasta")
                       
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "ref", filename)):
                os.remove(os.path.join(settings.install_folder, "fasta", "ref", filename))
            
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "ref", "temp.fasta")):
                with open(os.path.join(settings.install_folder, "fasta", "ref", "temp.fasta"), "rb") as f_in:
                    with gzip.open(os.path.join(settings.install_folder, "fasta", "ref", "%s" % file), "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)            
 
                os.remove(os.path.join(settings.install_folder, "fasta", "ref", "temp.fasta"))
                
    proteome_df = pd.DataFrame.from_dict(proteomes, orient='index').reset_index(names=['PPID'])
    proteome_df.to_csv(os.path.join(settings.install_folder, "fasta", "proteomes.tsv"), sep='\t', index=False)
    
    write_log("Finished generating profile.fasta", "fasta", project)