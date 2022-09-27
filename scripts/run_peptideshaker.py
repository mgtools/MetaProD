# runs peptideshaker
# python3 manage.py runscript run_peptideshaker --script-args 1

import shutil
import os
import psutil
import argparse
import time

from django.core.exceptions import ObjectDoesNotExist

from projects.models import Setting, Queue, RunTime

from .run_command import run_command, write_debug, settings

def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('queue_id', type=int)
    args2 = parser.parse_args(args)
    
    queue_id = args2.queue_id
        
    run_peptideshaker(queue_id)

def run_peptideshaker(queue_id):
    try:
        queue = Queue.objects.get(id=queue_id)
    except ObjectDoesNotExist:
        print("peptideshaker missing queue_id: %s" % queue_id)
        return False

    filename = queue.filename
    project = queue.project.name
 
    install_folder = settings.install_folder
    job = queue.job
    
    start = time.time()
    
    if queue.status == Queue.Status.PEPTIDESHAKER_PROF:
        type = "profile"
        fasta_file = "%s%s%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", type), os.sep, project, type)
    elif queue.status == Queue.Status.PEPTIDESHAKER_PROT:
        type = "proteome"
        fasta_file = "%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", type, filename), os.sep, project, filename, type)

    if not os.path.exists(os.path.join(install_folder, "software", "PeptideShaker-%s" % settings.peptideshaker_ver, "PeptideShaker-%s.jar" % settings.peptideshaker_ver)):
        write_debug("Missing PeptideShaker install.", job, project)
        return False
        
    if not os.path.exists(fasta_file):
        write_debug("Missing %s FASTA file for %s." % (type, filename), job, project)
        return False
        
    # remake the temp folder
    if not os.path.exists(os.path.join(install_folder, "temp", project, str(job), "software")):
        os.makedirs(os.path.join(install_folder, "temp", project, str(job), "software"))

    if os.path.exists(os.path.join(install_folder, "temp", project, str(job), "temp", "PeptideShaker")):
        shutil.rmtree(os.path.join(install_folder, "temp", project, str(job), "temp", "PeptideShaker"))
    
    # copy peptideshaker into temp
    if os.path.exists(os.path.join(install_folder, "temp", project, str(job), "software", "PeptideShaker-%s" % settings.peptideshaker_ver)):
        shutil.rmtree(os.path.join(install_folder, "temp", project, str(job), "software", "PeptideShaker-%s" % settings.peptideshaker_ver))
        
    shutil.copytree(os.path.join(settings.install_folder, "software", "PeptideShaker-%s" % settings.peptideshaker_ver), 
                    os.path.join(install_folder, "temp", project, str(job), "software", "PeptideShaker-%s" % settings.peptideshaker_ver))

    # remove the old output if it exists
    if os.path.exists(os.path.join(settings.data_folder, project, "out", filename, type, "%s.psdb" % filename)):
        os.remove(os.path.join(settings.data_folder, project, "out", filename, type, "%s.psdb" % filename))
        
    if os.path.exists(os.path.join(settings.data_folder, project, "out", filename, type, "ps_%s_Default_PSM_Report.txt" % project)):
        os.remove(os.path.join(settings.data_folder, project, "out", filename, type, "ps_%s_Default_PSM_Report.txt" % project))
        
    write_debug("Starting PathSettingsCLI: %s" % (os.path.join(settings.data_folder, project)), job, project)
    success = run_command(["timeout", "86400", 
                    "java", "-Xms%s" % settings.memory, "-Xmx%s" % settings.memory, 
                    "-cp", os.path.join(install_folder, "temp", project, str(job), "software", "PeptideShaker-%s" % settings.peptideshaker_ver, "PeptideShaker-%s.jar" % settings.peptideshaker_ver), 
                    "eu.isas.peptideshaker.cmd.PathSettingsCLI",
                    "-temp_folder", "%s" % os.path.join(install_folder, "temp", project, str(job), "temp", "PeptideShaker"),
                    "-identification_parameters", "%s" % os.path.join(install_folder, "temp", project, str(job), "temp", "PeptideShaker"),
                    ], job, project) 

    if (success == 0):
            write_debug("peptideShakerPathSettingsCLI failed", job, project)
            return False
            
    if settings.threads == -1:
        threads = psutil.cpu_count()
    else:
        threads = settings.threads
     
    # timeout after 8 hours
    write_debug("Starting PeptideShakerCLI: %s" % (os.path.join(settings.data_folder, project, "out", filename, type, "searchgui_out.zip")), job, project)
    success = run_command(["timeout", "86400", 
                            "java", "-Xms%s" % settings.memory, "-Xmx%s" % settings.memory, 
                            "-cp", os.path.join(install_folder, "temp", project, str(job), "software", "PeptideShaker-%s" % settings.peptideshaker_ver, "PeptideShaker-%s.jar" % settings.peptideshaker_ver),
                            "eu.isas.peptideshaker.cmd.PeptideShakerCLI",
                            "-out", "%s%s%s.psdb" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, filename),
                            "-reference", project,
                            "-identification_files", "%s" % (os.path.join(settings.data_folder, project, "out", filename, type, "searchgui_out.zip")),
                            "-threads", "%s" % threads,
                            "eu.isas.peptideshaker.cmd.ReportCLI",
                            "-reports", "3", # this has to be generated in the gui other than defaults
                            "-out_reports", "%s" % (os.path.join(settings.data_folder, project, "out", filename, type)),
                            "-report_prefix", "ps_"
                          ], job, project)
        
    if success == 0 or not os.path.exists(r"%s%s%s.psdb" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, filename)):
        write_debug("Missing PeptideShakerCLI output: %s%s%s.psdb" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, filename), job, project)
        return False
        
    if (not os.path.exists(r"%s%sps_%s_Default_PSM_Report.txt" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, project))):
        write_debug("PeptideShaker output missing %s" % (os.path.join(settings.data_folder, project, "out", filename, type)), job, project)
        return False
    else:
        end = time.time()
        runtime = end - start
        runtimex = RunTime.objects.get(queue=queue)
        if type == 'profile':
            runtimex.peptideshaker_profile = runtime
        elif type == 'proteome':
            runtimex.peptideshaker_proteome = runtime
        runtimex.save()
        
        return True
