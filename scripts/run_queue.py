import os
import shutil
import sys
import argparse
from django.utils import timezone
from django.core.exceptions import ObjectDoesNotExist

from projects.models import Queue, Project, Setting, SearchSetting
from projects.models import RunTime
from results.models import Protein, Peptide, Psm
from results.models import PsmRatio, SpeciesSummary, SpeciesFileSummary

from .run_command import run_command, write_debug, settings
from .run_msconvert import run_msconvert
from .run_searchgui import run_searchgui
from .run_peptideshaker import run_peptideshaker
from .read_results import read_results
from .generate_fasta import generate_fasta
from .run_reporter import run_reporter
from .process_results import process_results
from .run_mzmine import run_mzmine

def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('project', type=str)
    parser.add_argument('job', type=int, nargs='?', default=0)
    args2 = parser.parse_args(args)
    project = args2.project
    job = args2.job
    
    run_queue(project, job)

def run_queue(project, job):
    ''' main queue processing '''
    try:
        p = Project.objects.get(name=project)
    except ObjectDoesNotExist:
        print("Unable to find project %s." % project)
        return False
    
    try:
        searchsetting=SearchSetting.objects.get(project=project)
    except ObjectDoesNotExist:
        write_debug("Missing searchsetting for project: %s." % project, job, project)
        return False
    
    write_debug("Starting queue processing for project %s and job %s." % 
        (project, job), job, project
    )
    rerun = 1
    while rerun == 1:
        queue = (Queue.objects.filter(project__name=project, job=job)
                              .exclude(error__gte=(1 + settings.max_retries))
                              .exclude(status=Queue.Status.FINISHED_PROF)
                              .exclude(status=Queue.Status.FINISHED_PROT)
                              .exclude(status=Queue.Status.FILE_FINISHED)
                              .exclude(skip=True)
                              .order_by('-status')
                )
        if not queue:
            write_debug("No remaining entries left in the queue for project %s and job %s." % (project, job), job, project)
            return

        queue = queue[:1][0]
        filename = queue.filename
        project = queue.project.name
        install_folder = settings.install_folder
        
        # done with proteome step, so hold here
        if queue.status == Queue.Status.FINISHED_PROT:
            return
            
        elif queue.status == Queue.Status.PROCESS_RESULTS_PROT:
            write_debug("Starting process_results for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
            if process_results(queue.id, "proteome") == True:
                queue.status = Queue.Status.FINISHED_PROT
                queue.date_finished_proteome = timezone.now()
                # calculate the runtime now
                write_debug("Calculating total runtime for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
                runtimex = RunTime.objects.get(queue=queue)
                runtime = (runtimex.msconvert
                           + runtimex.searchgui_profile
                           + runtimex.peptideshaker_profile
                           + runtimex.reporter_profile
                           + runtimex.mzmine_profile
                           + runtimex.read_results_profile
                           + runtimex.process_results_profile
                           + runtimex.searchgui_proteome
                           + runtimex.peptideshaker_proteome
                           + runtimex.reporter_proteome
                           + runtimex.mzmine_proteome
                           + runtimex.read_results_proteome
                           + runtimex.process_results_proteome
                           )
                queue.total_runtime = runtime
                queue.save()
                write_debug("Finished process_results for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
            else:
                write_debug("Failed process_results for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
                queue.error += 1
                queue.save()
              
        elif queue.status == Queue.Status.READ_RESULTS_PROT:
            write_debug("Starting read_results for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
            if read_results(queue.id, "proteome") == True:
                queue.status = Queue.Status.PROCESS_RESULTS_PROT
                queue.error = 0
                queue.save()
                write_debug("Finished read_results for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
            else:
                write_debug("Failed read_results for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
                queue.error += 1 
                queue.save()

        elif queue.status == Queue.Status.MZMINE_PROT:
            if searchsetting.mzmine_run_mzmine == True:
                write_debug("Starting run_mzmine for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
                if run_mzmine(queue.id) == True:
                    queue.status = Queue.Status.READ_RESULTS_PROT
                    queue.error = 0
                    queue.save()
                    write_debug("Finished run_mzmine for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
                else:
                    write_debug("Failed run_mzmine for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
                    queue.error += 1
                    queue.save()
            else:
                queue.status = Queue.Status.READ_RESULTS_PROT
                queue.save()
                
        elif queue.status == Queue.Status.REPORTER_PROT:
            # reporter is only needed for multiplexed
            if searchsetting.multiplex == True:
                write_debug("Starting run_reporter for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
                if run_reporter(queue.id) == True:
                    queue.status = Queue.Status.MZMINE_PROT
                    queue.error = 0
                    queue.save()
                    write_debug("Finished run_reporter for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
                else:
                    write_debug("Failed run_reporter for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
                    queue.error += 1
                    queue.save()
            else:
                queue.status = Queue.Status.MZMINE_PROT
                queue.save()

        elif queue.status == Queue.Status.PEPTIDESHAKER_PROT:
            # peptideshaker needs to be run for reporter
            write_debug("Starting run_peptideshaker for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
            if run_peptideshaker(queue.id) == True:
                queue.status = Queue.Status.REPORTER_PROT
                queue.error = 0
                queue.save()
                write_debug("Finished run_peptideshaker for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
            else:
                write_debug("Failed run_peptideshaker for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
                queue.error += 1
                queue.save()

        elif queue.status == Queue.Status.SEARCHGUI_PROT:
            write_debug("Starting run_searchgui for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
            if run_searchgui(queue.id) == True:
                queue.status = Queue.Status.PEPTIDESHAKER_PROT
                queue.error = 0
                queue.save()
                write_debug("Finished run_searchgui for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
            else:
                write_debug("Failed run_searchgui for project: %s, filename: %s, type: proteome." % (project, filename), job, project)
                queue.error += 1
                queue.save()

        # all this status indicates is that the file is done, so nothing should be done
        # this shouldn't be selected though but is here for reference
        elif queue.status == Queue.Status.FINISHED_PROF:
            return
            
        elif queue.status == Queue.Status.PROCESS_RESULTS_PROF:
            write_debug("Starting process_results for project: %s, filename: %s" % (project, filename), job, project)
            if process_results(queue.id, "profile") == True:
                queue.status = Queue.Status.FINISHED_PROF
                queue.date_finished_profile = timezone.now()
                # calculate the runtime now
                runtimex = RunTime.objects.get(queue=queue)
                runtime = (runtimex.msconvert
                           + runtimex.searchgui_profile
                           + runtimex.peptideshaker_profile
                           + runtimex.reporter_profile
                           + runtimex.mzmine_profile
                           + runtimex.read_results_profile
                           + runtimex.process_results_profile
                           )
                queue.total_runtime = runtime
                queue.save()
                queue.date_finished_profile = timezone.now()
                write_debug("Finished process_results for project: %s, filename: %s, type: profile." % (project, filename), job, project)
            else:
                write_debug("Failed process_results for project: %s, filename: %s, type: profile." % (project, filename), job, project)
                queue.error += 1
                queue.save()
                
        elif queue.status == Queue.Status.READ_RESULTS_PROF:
            if searchsetting.custom_fasta == True:
                fasta_type = "custom"
            else:
                fasta_type = "profile"        
            write_debug("Starting read_results for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
            if read_results(queue.id, fasta_type) == True:
                queue.status = Queue.Status.PROCESS_RESULTS_PROF
                queue.error = 0
                queue.save()
                write_debug("Finished read_results for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
            else:
                write_debug("Failed read_results for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
                queue.error += 1 
                queue.save()

        elif queue.status == Queue.Status.MZMINE_PROF:
            if searchsetting.custom_fasta == True:
                fasta_type = "custom"
            else:
                fasta_type = "profile"
            if searchsetting.mzmine_run_mzmine == True:
                write_debug("Starting run_mzmine for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
                if run_mzmine(queue.id) == True:
                    queue.status = Queue.Status.READ_RESULTS_PROF
                    queue.error = 0
                    queue.save()
                    write_debug("Finished run_mzmine for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
                else:
                    write_debug("Failed run_mzmine for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
                    queue.error += 1
                    queue.save()
            else:
                queue.status = Queue.Status.READ_RESULTS_PROF
                queue.save()
                
        elif queue.status == Queue.Status.REPORTER_PROF:
            # reporter is only needed for multiplexed
            if searchsetting.multiplex == True:
                write_debug("Starting run_reporter for project: %s, filename: %s" % (project, filename), job, project)
                if run_reporter(queue.id) == True:
                    queue.status = Queue.Status.MZMINE_PROF
                    queue.error = 0
                    queue.save()
                    write_debug("Starting run_reporter for project: %s, filename: %s" % (project, filename), job, project)
                else:
                    write_debug("Finished run_reporter for project: %s, filename: %s, type: profile." % (project, filename), job, project)
                    queue.error += 1
                    queue.save()
            else:
                queue.status = Queue.Status.MZMINE_PROF
                queue.save()
                
        elif queue.status == Queue.Status.PEPTIDESHAKER_PROF:
            if searchsetting.custom_fasta == True:
                fasta_type = "custom"
            else:
                fasta_type = "profile"        
            # peptideshaker needs to be run for reporter
            write_debug("Starting run_peptideshaker for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
            if run_peptideshaker(queue.id) == True:
                queue.status = Queue.Status.REPORTER_PROF
                queue.error = 0
                queue.save()
                write_debug("Finished run_peptideshaker for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
            else:
                write_debug("Failed run_peptideshaker for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
                queue.error += 1
                queue.save()
                
        elif queue.status == Queue.Status.SEARCHGUI_PROF:
            if searchsetting.custom_fasta == True:
                fasta_type = "custom"
            else:
                fasta_type = "profile"         

            write_debug("Starting run_searchgui for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
            if run_searchgui(queue.id) == True:
                queue.status = Queue.Status.PEPTIDESHAKER_PROF
                queue.error = 0
                queue.save()
                write_debug("Finished run_searchgui for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
            else:
                write_debug("Failed run_searchgui for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
                queue.error += 1
                queue.save()
                
        # this only needs to be done once per file
        elif queue.status == Queue.Status.THERMO:
            if searchsetting.custom_fasta == True:
                fasta_type = "custom"
            else:
                fasta_type = "profile"         
                
            write_debug("Starting run_msconvert for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
            if run_msconvert(queue.id) == True:
                if searchsetting.profile == True:
                    queue.status = Queue.Status.SEARCHGUI_PROF
                else:
                    queue.status = Queue.Status.SEARCHGUI_PROT
                queue.error = 0
                queue.save()
                write_debug("Finished run_msconvert for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)                
            else:
                write_debug("Failed run_msconvert for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
                queue.error += 1
                queue.save()
                
        # clean up/set up
        # wipe existing outputs and make sure the proper dirs exist
        # also wipe up the database
        elif queue.status == Queue.Status.FILE_ADDED:
            if searchsetting.custom_fasta == True:
                fasta_type = "custom"
            else:
                fasta_type = "profile"        

            write_debug("Cleaning up existing entries for project: %s, filename: %s, type: %s." % (project, filename, fasta_type), job, project)
            delete = Protein.objects.filter(queue=queue).delete()
            delete = Peptide.objects.filter(queue=queue).delete()
            delete = Psm.objects.filter(queue=queue).delete()
            delete = PsmRatio.objects.filter(psm__queue=queue).delete()
            delete = RunTime.objects.filter(queue=queue).delete()
            delete = SpeciesFileSummary.objects.filter(queue=queue).delete()
            # create the runtimex table
            runtimex = RunTime(queue=queue)
            runtimex.save()
            
            if os.path.exists(os.path.join(settings.data_folder, project, "out", filename)):
                shutil.rmtree(os.path.join(settings.data_folder, project, "out", filename))
                
            os.makedirs(os.path.join(settings.data_folder, project, "out", filename))
            
            if searchsetting.custom_fasta == True:
                os.makedirs(os.path.join(settings.data_folder, project, "out", filename, "custom"))
            else:
                os.makedirs(os.path.join(settings.data_folder, project, "out", filename, "profile"))
            
                os.makedirs(os.path.join(settings.data_folder, project, "out", filename, "proteome"))

            # clean up the temp dir
            if os.path.exists(os.path.join(install_folder, "temp", project, str(job))):
                shutil.rmtree(os.path.join(install_folder, "temp", project, str(job)))
                     
            queue.date_finished_profile = None
            queue.date_finished_proteome = None
            queue.status = Queue.Status.THERMO
            queue.error = 0
            queue.total_runtime = 0
            queue.save()
            
        else:
            write_debug("Invalid status for project: %s, filename: %s" % (project, filename), job, project)
            queue.error = 1 + settings.max_retries
            queue.save()

# once we run update_queue final, wipe the old files
# we no longer remove files needed for other steps so they can be re-run
def cleanup(project):
    print("Done processing. Cleaning up %s." % (project))
    # remove temp folder
    if os.path.exists(os.path.join(settings.install_folder, "temp", project)):
        shutil.rmtree(os.path.join(settings.install_folder, "temp", project))