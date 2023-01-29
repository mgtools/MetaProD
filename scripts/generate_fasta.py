# downloads proteomes and generates profile and full proteome fasta files

import re
import argparse
import os
import shutil
import gzip
import csv
from pathlib import Path
import urllib.request
from Bio import SeqIO

from django.db.models import Sum
from django.core.exceptions import ObjectDoesNotExist

from projects.models import Queue, SearchSetting
from results.models import Protein
from .run_command import run_command, settings
from .load_proteomes import load_proteomes

def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('project_name', type=str)
    parser.add_argument('type', choices=['profile', 'proteome'])
    args2 = parser.parse_args(args)
    project = args2.project_name
    type = args2.type

    generate_fasta(project, type)

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

def generate_fasta(project, type):
    print("Starting generate_fasta for %s %s." % (project, type))

    try:
        searchsetting=SearchSetting.objects.get(project=project)
    except ObjectDoesNotExist:
        print("Missing searchsetting for project: %s." % project)
        return False

    if searchsetting.profile == False:
        print("Generating a FASTA is not supported if profiling is turned off.")
        print("This feature not yet implemented so turn profiling on in the searchsettings.")
        return 0

    samples = {}
    if type == "proteome":
        if not os.path.exists(os.path.join(settings.data_folder, project, "fasta")):
            os.makedirs(os.path.join(settings.data_folder, project, "fasta"))
        if not os.path.exists(os.path.join(settings.data_folder, project, "fasta", "proteome")):
            os.makedirs(os.path.join(settings.data_folder, project, "fasta", "proteome"))

        query = (Queue.objects.filter(status__gte=Queue.Status.FINISHED_PROF)
                              .exclude(error__gte=2)
                              .exclude(skip=True)
                              .filter(project__name=project))

        for entry in query:
            # if the sample isn't set, we have to use the filename
            if entry.sample is None:
                sample = entry.filename
                have_sample = 0
            else:
                sample = entry.sample.name
                have_sample = 1

            if os.path.exists(os.path.join(settings.data_folder, project, "fasta", "proteome", entry.filename)):
                shutil.rmtree(os.path.join(settings.data_folder, project, "fasta", "proteome", entry.filename))
            os.makedirs(os.path.join(settings.data_folder, project, "fasta", "proteome", entry.filename))
                    
            if searchsetting.profile_type == SearchSetting.ProfileType.FILE:
                print("File based profiling. Generating proteome FASTA for %s" % entry.filename)
                generate_proteome_fasta(project, entry.filename, "", have_sample)
            elif searchsetting.profile_type == SearchSetting.ProfileType.SAMPLE:
                # haven't made it for the sample yet
                if sample not in samples:
                    print("Sample based profiling. Generating proteome FASTA for %s." % (sample))
                    generate_proteome_fasta(project, entry.filename, sample, have_sample)
                    samples[sample] = entry.filename
                else:
                    print("FASTA already exists. Copying.")
                    previous = samples[sample]
                    shutil.copy("%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", previous), os.sep, project, previous, "proteome"),
                                "%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", entry.filename), os.sep, project, entry.filename, "proteome"))
            else:
                # the first pooled hasn't been made yet
                if len(samples) < 1:
                    print("Project based profiling. Generating proteome FASTA for %s" % (project))
                    generate_proteome_fasta(project, entry.filename, "", have_sample)

                    samples[sample] = entry.filename
                # it's already been made, so just copy the first one everywhere else
                else:
                    print("FASTA already exists. Copying.")
                    previous = list(samples.values())[0]
                    shutil.copy("%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", previous), os.sep, project, previous, "proteome"),
                                "%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", entry.filename), os.sep, project, entry.filename, "proteome"))
          
                
            fasta_file_concat = "%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", entry.filename), os.sep, project, entry.filename, "proteome")
                
            if not os.path.exists(fasta_file_concat):
                print("Missing FASTA file:", project, "fasta.")
                entry.error += 1
                entry.save()
                    
                print("Failed to generate proteome_fasta for %s." % (entry.filename))
            
    elif type == "profile":
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
            print("Missing FASTA file.", project, "fasta")
            print("Failed to generate profile FASTA for %s." % (project))
            
def generate_proteome_fasta(project, filename, sample, have_sample):
    try:
        searchsetting=SearchSetting.objects.get(project__name=project)
    except:
        print("Searchsetting does not exist for %s. Make sure the project and searchsetting are added first." % (project))
        return 0
   
    fasta_file = "%s%s%s_%s_%s.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", filename), os.sep, project, filename, "proteome")
    
    print("Calculating species to include in FASTA.")

    proteomes = []
    count = 0
    current_psm = 0    
    if searchsetting.profile_type == SearchSetting.ProfileType.PROJECT:
        query_t = (
            Protein.objects
            .filter(type="profile")
            .filter(queue__status__gte=Queue.Status.FINISHED_PROF)
            .filter(queue__project__name=project)
            .exclude(queue__skip=True)
            .exclude(fp__ppid__proteome='0')
            .exclude(fp__ppid__proteome='UP000005640')                    
            .aggregate(Sum('nsaf'))
        )
        total = query_t['nsaf__sum']
        if not total is None:
            total = float(total) * float(searchsetting.profile_threshold / 100)
            query = (
                Protein.objects
                .filter(type="profile")
                .filter(queue__status__gte=Queue.Status.FINISHED_PROF)
                .filter(queue__project__name=project)
                .exclude(queue__skip=True)
                .exclude(fp__ppid__proteome='0')
                .exclude(fp__ppid__proteome='UP000005640')                        
                .values('fp__ppid__proteome')
                .annotate(sum=Sum('nsaf'))
                .order_by('-sum')
            )

    # not pooled, so a filter also needs to specify the filename
    elif searchsetting.profile_type == SearchSetting.ProfileType.SAMPLE and have_sample == 1:
        query_t = (
            Protein.objects
            .filter(type="profile")
            .filter(queue__project__name=project)
            .filter(queue__sample__name=sample)
            .filter(queue__sample__project=project)
            .filter(queue__status__gte=Queue.Status.FINISHED_PROF)
            .exclude(queue__error__gte=2)
            .exclude(queue__skip=True)
            .exclude(fp__ppid__proteome='0')
            .exclude(fp__ppid__proteome='UP000005640')
            .aggregate(Sum('nsaf'))
        )
        total = query_t['nsaf__sum']
        if not total is None:
            total = float(total) * float(searchsetting.profile_threshold / 100)
    
            query = (
                Protein.objects
                .filter(type="profile")
                .filter(queue__project__name=project)
                .filter(queue__sample__name=sample)
                .filter(queue__sample__project=project)
                .filter(queue__status__gte=Queue.Status.FINISHED_PROF)
                .exclude(queue__error__gte=2)
                .exclude(queue__skip=True)
                .exclude(fp__ppid__proteome='0')
                .exclude(fp__ppid__proteome='UP000005640')
                .values('fp__ppid__proteome')
                .annotate(sum=Sum('nsaf'))
                .order_by('-sum')
            )

    # either we are doing filename profiling or the file didn't have a sample 
    #  set so the sample is effectively the filename
    else:
        query_t = (
            Protein.objects
            .filter(type="profile")
            .filter(queue__status__gte=Queue.Status.FINISHED_PROF)
            .filter(queue__project__name=project)
            .filter(queue__filename=filename)
            .exclude(queue__error__gte=2)
            .exclude(queue__skip=True)
            .exclude(fp__ppid__proteome='0')
            .exclude(fp__ppid__proteome='UP000005640')
            .aggregate(Sum('nsaf'))
        )
        total = query_t['nsaf__sum']
        if not total is None:
            total = float(total) * float(searchsetting.profile_threshold / 100)

            query = (
                Protein.objects
                .filter(type="profile")
                .filter(queue__status__gte=Queue.Status.FINISHED_PROF)
                .filter(queue__project__name=project)
                .filter(queue__filename=filename)
                .exclude(queue__error__gte=2)
                .exclude(queue__skip=True)
                .exclude(fp__ppid__proteome='0')
                .exclude(fp__ppid__proteome='UP000005640')
                .values('fp__ppid__proteome')
                .annotate(sum=Sum('nsaf'))
                .order_by('-sum')
            )

    if not total is None:
        print("%s%% of bacterial NSAF: %s" % (searchsetting.profile_threshold, total))
        for entry in query:
            # running count is below the threshold, so we add it
            if count <= total:
                proteomes.append(entry['fp__ppid__proteome'])
                count += entry['sum']
                print("Including %s: %s" % (entry['fp__ppid__proteome'], count))                
                current_psm = entry['sum']
            else:
                if entry['sum'] == current_psm:
                    proteomes.append(entry['fp__ppid__proteome'])
                    print("Including %s due to same NSAF as last entry." % (entry['fp__ppid__proteome']))
                else:
                    break
            
            # don't bother adding CRAP or any others without headers
            if len(entry['fp__ppid__proteome']) < 1 or entry['fp__ppid__proteome'] == '0':
                continue
                
            if entry['fp__ppid__proteome'] not in proteomes:
                proteomes.append(entry['fp__ppid__proteome'])

    fasta_file = "%s%s%s_%s_%s.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", filename), os.sep, project, filename, "proteome")
    fasta_file_concat = "%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", filename), os.sep, project, filename, "proteome")
    
    if os.path.exists(fasta_file):
        os.remove(fasta_file)

    if os.path.exists(fasta_file_concat):
        os.remove(fasta_file_concat)
        
    print("Generating full proteome FASTA.")

    accessions = set()
    for proteome in proteomes:
        # check to figure out what dir it is in:
        if os.path.exists(os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta.gz" % proteome)):
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta" % proteome)):
                os.remove(os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta" % proteome))
            with gzip.open(os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta.gz" % proteome), "rb") as f_in:
                with open(os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta" % proteome), "ab") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            # save the file with the corrected taxonomy
            for record in SeqIO.parse(os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta" % proteome), "fasta"):
                p = re.compile("^(?P<start>(?P<db>[^\|]+)\|(?P<accession>[^\|]+)\|(?P<middle>.+)\s)(?P<OS>OS=.+)\s(?P<OX>OX=.+)\s(?P<UPId>UPId=[^\s]+)\s(?P<PPId>PPId=[^\s]+)$")
                m1 = p.search(record.description)
                if m1:
                    accession = m1.group('accession')
                    tax = m1.group('OS')
                    # strip out the =
                    upid = m1.group('UPId').replace("=","")
                    ppid = m1.group('PPId').replace("=","")
                    # add them to OS
                    tax = tax + " " + upid + " " + ppid
                    description = m1.group('start') + " " + tax + " " + m1.group('OX') + " " + m1.group('UPId')  + " " + m1.group('PPId')
                    
                    # now write the new fasta header + sequence to a file
                    record.description = description
                    if accession not in accessions:
                        accessions.add(accession)
                        with open("%s%s%s_%s_%s.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", filename), os.sep, project, filename, "proteome"), "a") as fasta_out:
                            SeqIO.write(record, fasta_out, "fasta")
                else:
                    print("no match for %s" % record.description)
                    
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta" % proteome)):
                os.remove(os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta" % proteome))
        elif os.path.exists(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta.gz" % proteome)):
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta" % proteome)):
                os.remove(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta" % proteome))
            with gzip.open(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta.gz" % proteome), "rb") as f_in:
                with open(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta" % proteome), "ab") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            filename_in = Path(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta.gz" % proteome)).stem
            for record in SeqIO.parse(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta" % proteome), "fasta"):
                p = re.compile("^(?P<start>(?P<db>[^\|]+)\|(?P<accession>[^\|]+)\|(?P<middle>.+)\s)(?P<OS>OS=.+)\s(?P<OX>OX=.+)")
                m1 = p.search(record.description)
                if m1:
                    accession = m1.group('accession')
                    tax = m1.group('OS')
                    # strip out the =
                    upid = 'UPId' + Path(filename_in).stem
                    ppid = 'PPId' + Path(filename_in).stem
                    # add them to OS
                    tax = tax + " " + upid + " " + ppid
                    description = m1.group('start') + " " + tax + " " + m1.group('OX') + " " + "UPId=" + Path(filename_in).stem + " " + "PPId=" + Path(filename_in).stem
                    # now write the new fasta header + sequence to a file
                    record.description = description
                    if accession not in accessions:
                        accessions.add(accession)
                        with open("%s%s%s_%s_%s.fasta" % (os.path.join(settings.data_folder, project, "fasta", "proteome", filename), os.sep, project, filename, "proteome"), "a") as fasta_out:
                            SeqIO.write(record, fasta_out, "fasta")
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta" % proteome)):
                os.remove(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta" % proteome))
        else:
            print("Error. Missing proteome download for %s" % proteome)
            return

    print("Done generating proteome FASTA.")
    
    if searchsetting.use_human == True:
        print("Appending human proteome to FASTA file.")
        with open(fasta_file, "a") as f_out:
            with open(os.path.join(settings.install_folder, "fasta", "human.fasta"), "r") as f_in:
                shutil.copyfileobj(f_in, f_out)
            f_out.write(os.linesep)

    if searchsetting.use_crap == True:
        print("Appending CRAP database to FASTA file.")
        # copy CRAP to destination file
        with open(fasta_file, "a") as f_out:
            with open(os.path.join(settings.install_folder, "fasta", "crap.fasta"), "r") as f_in:
                shutil.copyfileobj(f_in, f_out)
            f_out.write(os.linesep)

    generate_decoy(project, fasta_file)
      
def generate_profile_fasta(project):
    try:
        searchsetting=SearchSetting.objects.get(project__name=project)
    except:
        print("Searchsetting does not exist for %s. Make sure the project and searchsetting are added first." % (project))
        return 0
        
    # make sure the fasta dir exists
    if not os.path.exists(os.path.join(settings.data_folder, project, "fasta")):
        os.makedirs(os.path.join(settings.data_folder, project, "fasta"))

    if not os.path.exists(os.path.join(settings.data_folder, project, "fasta", "profile")):
        os.makedirs(os.path.join(settings.data_folder, project, "fasta", "profile"))
        
    if os.path.exists("%s%s%s.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep, "profile")):
        print("Profile FASTA already exists. Using existing FASTA.")
    else:
        print("Profile FASTA does not exist. Downloading a new FASTA.")
        print("Note: Uniprot sometimes changes their website. Let the authors know if this does not work anymore.")
        download_profile_fasta()
        create_full_fasta()
        filter_high_profile()

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
        print("Appending human proteome to FASTA file.")
        with open(fasta_file, "a") as f_out:
            with open(os.path.join(settings.install_folder, "fasta", "human.fasta"), "r") as f_in:
                shutil.copyfileobj(f_in, f_out)
            f_out.write(os.linesep)

    if searchsetting.use_crap == True:
        print("Appending CRAP database to FASTA file.")
        # copy CRAP to destination file
        with open(fasta_file, "a") as f_out:
            with open(os.path.join(settings.install_folder, "fasta", "crap.fasta"), "r") as f_in:
                shutil.copyfileobj(f_in, f_out)
            f_out.write(os.linesep)

    success = generate_decoy(project, fasta_file)
    
    if success == True:
        load_proteomes()

        print("Done generating profile FASTA.")
    
    
# should only need the filename for this as the decoy will end up in the same dir as the filename
def generate_decoy(project, fasta_file):
    ''' generate decoy sequences for FASTA file '''
    print("Generating decoy sequences.")

    if os.path.exists(os.path.join(settings.data_folder, project, "fasta", "temp")):
        shutil.rmtree(os.path.join(settings.data_folder, project, "fasta", "temp"))

    # remake the temp folder
    os.makedirs(os.path.join(settings.data_folder, project, "fasta", "temp", "software"))
    
    if not os.path.exists(os.path.join(settings.install_folder, "software", "SearchGUI-%s" % settings.searchgui_ver, "SearchGUI-%s.jar" % settings.searchgui_ver)):
        print("Missing SearchGUI install.")
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
    
    fasta_file_concat = "%s%s%s_profile_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", "profile"), os.sep, project)
    
    if os.path.exists(os.path.join(fasta_file_concat)):
        print("Finished generating decoy sequences.")    
        shutil.rmtree(os.path.join(settings.data_folder, project, "fasta", "temp"))
        return True
    else:
        print("Failed to generate decoy sequences.")
        return False

def download_profile_fasta():
    ''' downloads zip files '''
    # download ppmembership.txt
    if os.path.exists(os.path.join(settings.install_folder, "fasta", "PPMembership.txt")):
        print("PPMembership.txt already exists. Remove to update.")
    else:
        print("Downloading PPMembership.txt from uniprot.")
        # need to use different url until uniprot is fixed
        #pp_membership_url = "https://proteininformationresource.org/rps/data/new/PPSeqCurrent/PPMembership.txt"
        pp_membership_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/pan_proteomes/PPMembership.txt"
        urllib.request.urlretrieve(pp_membership_url, os.path.join(settings.install_folder, "fasta", "PPMembership.txt"))
            
    # download reference proteomes
    if os.path.exists(os.path.join(settings.install_folder, "fasta", "ref_proteomes_list.tsv")):
        print("ref_proteomes_list.tsv already exists. Remove to update.")
    else:
        print("Downloading list of reference proteomes from uniprot.")
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
    
    
    print("Generating list of proteomes to download.")

    reference_proteomes = []
    with open(os.path.join(settings.install_folder, "fasta", "ref_proteomes_list.tsv"), 'r') as file:
        result = csv.reader(file, delimiter='\t')
        header = next(result)
        header = [x.lower() for x in header]
        for row in result:
            if row[header.index('proteome id')] not in reference_proteomes:
                reference_proteomes.append(row[header.index('proteome id')])
    
    pan_proteomes = []
    ref_proteomes_added = []

    with open(os.path.join(settings.install_folder, "fasta", "PPMembership.txt"), 'r') as file:
        i = 1
        result = csv.reader(file, delimiter='\t')
        next(result)
        for row in result:
            # so we haven't added this pan proteome yet
            if row[0] not in pan_proteomes:
                # is one of the members in the reference proteome?
                # if so, add the pan proteome
                if row[1] in reference_proteomes:
                    pan_proteomes.append(row[0])
                    ref_proteomes_added.append(row[1])
            else:
                if row[1] in reference_proteomes:
                    ref_proteomes_added.append(row[1])

    ref_proteomes_dl = []
    for proteome in reference_proteomes:
        if proteome not in ref_proteomes_added:
            ref_proteomes_dl.append(proteome)

    # human
    ref_proteomes_dl.append("UP000005640")
        
    # download the pan proteomes
    # https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/pan_proteomes/UP000000212.fasta.gz
    print("Downloading needed pan proteomes.")
    lenp = len(pan_proteomes)
    i = 1
    for proteome in pan_proteomes:
        attempts = 10
        for attempt in range(attempts):
            try:
                # if it already exists, don't download it again because uniprot is slow
                # but there could be a problem if the filesize is wrong so might have to check that down the road
                if not os.path.exists(os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta.gz" % proteome)):
                    print("Downloading %s pan-proteome of %s. (%s)" % (i, lenp, proteome))
                    # fix until uniprot is fixed
                    proteome_url = "https://proteininformationresource.org/rps/data/new/PPSeqCurrent/" + proteome + ".fasta.gz"
                    #proteome_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/pan_proteomes/" + proteome + ".fasta.gz"
                    urllib.request.urlretrieve(proteome_url, os.path.join(settings.install_folder, "fasta", "pan", "%s.fasta.gz" % proteome))
            except Exception as e:
                if attempt == 2:
                    print('Error downloading %s: ' % proteome, e)
            else:
                break
        i += 1
    
    # now download the reference proteomes remaining
    # put this in another dir so that we know how to change the header later
    # https://www.uniprot.org/uniprot/?query=proteome:UP000027395&format=fasta&compress=yes
    print("Downloading needed reference proteomes.")
    if not os.path.exists(os.path.join(settings.install_folder, "fasta", "ref")):
        os.makedirs(os.path.join(settings.install_folder, "fasta", "ref"))
    i = 1
    lenp = len(ref_proteomes_dl)
    for proteome in ref_proteomes_dl:
        attempts = 10
        for attempt in range(attempts):
            try:
                if not os.path.exists(os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta.gz" % proteome)):
                    print("Downloading %s reference proteome of %s. (%s)" % (i, lenp, proteome))
                    proteome_url = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=proteome:" + proteome
                    #proteome_url = "https://www.uniprot.org/uniprot/?query=proteome:" + proteome + "&format=fasta&compress=yes"
                    urllib.request.urlretrieve(proteome_url, os.path.join(settings.install_folder, "fasta", "ref", "%s.fasta.gz" % proteome))
            except Exception as e:
                if attempt == 2:
                    print('Error downloading %s: ' % proteome, e)
            else:
                break
        i += 1

def create_full_fasta():
    ''' produces full.fasta with all sequences and no filtering'''
    seq_count = 0

    if os.path.exists("%s%sfull.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep)):
        print("full.fasta exists already. Delete to recreate.")
        return
   
    if os.path.exists("%s%sfull.fasta.gz" % (os.path.join(settings.install_folder, "fasta"), os.sep)):
        print("full.fasta.gz exists already. Delete to recreate.")
        return
        
    # pan proteome files
    print("Extracting proteins for pan proteomes.")
    for file in os.listdir(os.path.join(settings.install_folder, "fasta", "pan")):
        if file.endswith(".fasta.gz"):
            filename = Path(file).stem
            with gzip.open(os.path.join(settings.install_folder, "fasta", "pan", file), "rb") as f_in:
                with open("%s%sfull.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep), "ab") as f_out:
                    shutil.copyfileobj(f_in, f_out)

    # reference proteomes do not include the proteome in the header, so we 
    #  need to add it
    print("Extracting proteins for reference proteomes (this may take a while).")
    for file in os.listdir(os.path.join(settings.install_folder, "fasta", "ref")):
        # exclude human
        if file.endswith(".fasta.gz"):
            filename = Path(file).stem
            if os.path.exists(os.path.join(settings.install_folder, "fasta", "ref", filename)):
                os.remove(os.path.join(settings.install_folder, "fasta", "ref", filename))
            with gzip.open(os.path.join(settings.install_folder, "fasta", "ref", file), "rb") as f_in:
                with open(os.path.join(settings.install_folder, "fasta", "ref", filename), "ab") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            for record in SeqIO.parse(os.path.join(settings.install_folder, "fasta", "ref", filename), "fasta"):
                p = re.compile("^(?P<start>(?P<db>[^\|]+)\|(?P<accession>[^\|]+)\|(?P<middle>.+)\s)(?P<OS>OS=.+)\s(?P<OX>OX=.+)")
                m1 = p.search(record.description)
                if m1:
                    tax = m1.group('OS')
                    description = m1.group('start') + " " + tax + " " + m1.group('OX') + " " + "UPId=" + Path(filename).stem + " " + "PPId=" + Path(filename).stem
                    # now write the new fasta header + sequence to a file
                    record.description = description
                    if file != "UP000005640.fasta.gz":
                        with open("%s%sfull.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep), "a") as fasta_out:
                            SeqIO.write(record, fasta_out, "fasta")
                        seq_count += 1
                    # write human to a different file so we can access it later
                    else:
                        with open("%s%shuman.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep), "a") as fasta_out:
                            SeqIO.write(record, fasta_out, "fasta")
                else:
                    print("no regexp match for: %s" % record)
                    
            if os.path.exists("%s" % (os.path.join(settings.install_folder, "fasta", "ref", filename))):
                os.remove("%s" % (os.path.join(settings.install_folder, "fasta", "ref", filename)))

    print("Finished generating the full FASTA file.")
    
def filter_high_profile():
    '''filters high profile out of the full fasta'''
    
    if os.path.exists("%s%sfull.fasta.gz" % (os.path.join(settings.install_folder, "fasta"), os.sep)):
        if os.path.exists("%s%sfull.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep)):
            os.remove("%s%sfull.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep))
        
        print("Unzipping full.fasta.gz")
        with gzip.open("%s%sfull.fasta.gz" % (os.path.join(settings.install_folder, "fasta"), os.sep), "rb") as f_in:
            with open("%s%sfull.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep), "ab") as f_out:
                shutil.copyfileobj(f_in, f_out)
        
    print("Filtering HAPs to generate profile FASTA.")
    if os.path.exists(os.path.join(settings.install_folder, "fasta", "profile.fasta")):
        os.remove(os.path.join(settings.install_folder, "fasta", "profile.fasta"))
                
    # need to keep track of accessions. if it exists already, just skip it.
    accessions = set()
    seq_count = 0
    for record in SeqIO.parse("%s%sfull.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep), "fasta"):
        m = re.search("ribosomal", record.description)
        if m:
            # find OS, PPId, UPId
            p = re.compile("^(?P<start>(?P<db>[^\|]+)\|(?P<accession>[^\|]+)\|(?P<middle>.+)\s)(?P<OS>OS=.+)\s(?P<OX>OX=.+)\s(?P<UPId>UPId=[^\s]+)\s(?P<PPId>PPId=[^\s]+)$")
            m1 = p.search(record.description)
            if m1:
                accession = m1.group('accession')
                tax = m1.group('OS')
                description = m1.group('start') + " " + tax + " " + m1.group('OX') + " " + m1.group('UPId')  + " " + m1.group('PPId')
                
                # now write the new fasta header + sequence to a file
                record.description = description
                if accession not in accessions:
                    accessions.add(accession)
                    with open("%s%s%s.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep, "profile"), "a") as fasta_out:
                        SeqIO.write(record, fasta_out, "fasta")
                    seq_count += 1
                else:
                    print("Duplicate accession: %s" % accession)
    
            else:
                print("no regexp match for: %s" % record)

    if os.path.exists("%s%sfull.fasta.gz" % (os.path.join(settings.install_folder, "fasta"), os.sep)):
        os.remove("%s%sfull.fasta.gz" % (os.path.join(settings.install_folder, "fasta"), os.sep))
    
    print("Zipping full.fasta")
    with open("%s%sfull.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep), "rb") as f_in:
        with gzip.open("%s%sfull.fasta.gz" % (os.path.join(settings.install_folder, "fasta"), os.sep), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    if os.path.exists("%s%sfull.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep)):
        os.remove("%s%sfull.fasta" % (os.path.join(settings.install_folder, "fasta"), os.sep))
        
    print("Finished generating profile FASTA with %s sequences." % seq_count)