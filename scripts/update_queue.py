import time
import os
import shutil
from pathlib import Path
import argparse

from django.db.models import Q
from django.core.exceptions import ObjectDoesNotExist

from projects.models import Queue, Project, Setting, SearchSetting
from results.models import SpeciesSummary, SpeciesFileSummary
from .run_command import write_debug, settings
from .generate_fasta import generate_fasta
from .run_queue import cleanup
from .process_results import calculate_nsaf, calculate_species_summary
from .process_results import calculate_species_file_summary

# fasta flag species to download/generate fasta otherwise
# it will use the existing fasta
def run(*args): 
    parser = argparse.ArgumentParser()
    parser.add_argument('project', type=str)
    parser.add_argument('jobs', type=int, default=1, nargs='?')
    args2 = parser.parse_args(args)
    project = args2.project
    jobs = args2.jobs

    generate_file_queue(project, jobs)
    
# generate the queue, as in look for all files that need to be added to the queue
# this is used on script launch
# this looks at raw files matching the format only and not .finished
def generate_file_queue(project, jobs):

    try:
        p = Project.objects.get(name=project)
    except ObjectDoesNotExist:
        print("Project %s does not yet exist. Use create_project project_name to create it." % project)
        return False

    try:
        searchsetting=SearchSetting.objects.get(project=project)
    except ObjectDoesNotExist:
        print("Missing searchsetting for project: %s." % project)
        return False
        
    jobs = int(jobs)
    job = 0

    install_folder = settings.install_folder
    
    # first, we check for any new files. if there are new files, we know it's profile
    files_added = 0

    files_list = list(Queue.objects.filter(project__name=project).values_list('filename', flat=True))

    # nothing added yet, so wipe the log dir
    if len(files_list) == 0:
        if os.path.exists(os.path.join(install_folder, "log", project)):
            shutil.rmtree(os.path.join(install_folder, "log", project))
        if os.path.exists(os.path.join(install_folder, "temp", project)):
            shutil.rmtree(os.path.join(install_folder, "temp", project))
            
    # look for files
    file_count = 0
    for root, dirs, files in os.walk(os.path.join(settings.data_folder, project, "raw")):
        for file in files:
            if file.endswith('.raw') or file.endswith('.mzML'):
                file_count += 1
                filename = Path(file).stem
                if filename not in files_list:
                    try:                
                        add_file_to_queue(os.path.join(settings.data_folder, 
                                            project, "raw", file), 
                                            project, Queue.Status.FILE_ADDED, job)
                        if job < jobs - 1:
                            job += 1
                        else:
                            job = 0
                    # failed so move it back to the raw incoming dir
                    except Exception as e:
                        print("error adding to queue: ", e)
                        write_debug(e, job, project)
                    else:                 
                        files_added += 1
    if files_added > 0:
        print("Files added to queue. Be sure to run: ./generate_fasta %s" % project)
        return True

    if file_count == 0:
        print("No files found in %s." % os.path.join(settings.data_folder, project, "raw"))
        return False
        
    # so we didn't add anything, now check for status 13
    # if this is true, everything is 13 or higher so run update_queue final
    if not Queue.objects.filter(project__name=project, skip=False).exclude(status=Queue.Status.FINISHED_PROT).exists():
        print("All files have finished the proteome step.")
        queue = (Queue.objects.filter(project__name=project))
        for c in queue:
            if c.skip is True:
                print("Warning: %s is set to be skipped." % (c.filename)) 
            
        calculate_nsaf(project, "proteome")
        calculate_species_summary(project, "proteome")
       
        for q in queue:
            if q.skip == True:
                continue
            write_debug("Calculating species file summary for %s." % (q.filename), q.job, q.project.name)
            calculate_species_file_summary(q.id, "proteome")
            q.error = 0
            q.status = Queue.Status.FILE_FINISHED
            q.save()
            
        cleanup(project)
        
        print("Finished update_queue final.")
        
    # if this is true, everything is 7 or higher but not 13 or higher
    elif not Queue.objects.filter(project__name=project, skip=False).exclude(status=Queue.Status.FINISHED_PROF).exists():
        print("All files have finished the profile step.")
        queue = (Queue.objects.filter(project__name=project))
        for c in queue:
            if c.skip is True:
                print("Warning: %s is set to be skipped." % (c.filename)) 
                
        calculate_nsaf(project, "profile")
        calculate_species_summary(project, "profile")
        generate_fasta(project, "proteome")
        
        for q in queue:
            if q.skip == True:
                continue
            write_debug("Calculating species file summary for %s." % (q.filename), q.job, q.project.name)
            calculate_species_file_summary(q.id, "profile")
            q.error = 0
            q.status = Queue.Status.SEARCHGUI_PROT
            q.save()

        cleanup(project)
        
        print("Finished update_queue proteome.")

    else:
        # maybe we should print the statuses here
        print("The project has no new files and has not completed the profile or proteome steps.")
        queue = (Queue.objects.filter(project__name=project))
        for c in queue:
            if c.error >= 2:
                print("Warning: %s has error status 2 and cannot be processed." % (c.filename))
        
def add_file_to_queue(filename, project, status, job):  
    # we only want the basename without path or raw
    filename = Path(filename).stem
    
    try:
        p = Project.objects.get(name=project)
        queue = Queue(filename=filename, project=p, 
                      status=status, job=job)
        queue.save()
        
        write_debug("Adding file to queue: Filename: %s Project: %s Queue ID: %s Job: %s" % (filename, 
                project, queue.id, job), job, project)
    except Exception as e:
        write_debug("adding to models failed for: %s %s" % (filename, e), job, project)
        pass

    return 1