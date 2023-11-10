import shutil
import os
import argparse
import time
import zipfile

from django.core.exceptions import ObjectDoesNotExist

from projects.models import (
    Setting, 
    Queue, 
    SearchSetting, 
    RunTime,
    LabelChoice,
)

from .run_command import run_command, write_debug, settings, write_error

def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('queue_id', type=int)
    parser.add_argument('job', type=int)
    args2 = parser.parse_args(args)
    
    queue_id = args2.queue_id
        
    run_reporter(queue_id)

def run_reporter(queue_id):
    ''' runs reporter for multiplexed data to get the psm ratios '''
    try:
        queue = Queue.objects.get(id=queue_id)
    except ObjectDoesNotExist:
        print("Reporter missing queue_id: %s" % queue_id)
        return False
    
    filename = queue.filename
    project = queue.project.name
    
    install_folder = settings.install_folder
    job = queue.job
    
    start = time.time()
    
    try:
        searchsetting=SearchSetting.objects.get(project=project)
    except ObjectDoesNotExist:
        print("Missing searchsetting for project: %s." % project)
        return False

    if queue.status == Queue.Status.REPORTER_PROT:
        type = "proteome"
    else:
        type = "profile"

    if not os.path.exists(os.path.join(install_folder, "software", "Reporter-%s" % settings.reporter_ver, "Reporter-%s.jar" % settings.reporter_ver)):
        write_debug("Missing Reporter install.", job, project)
        return False
        
    # remake the temp folder
    if not os.path.exists(os.path.join(install_folder, "temp", project, str(job), "software")):
        os.makedirs(os.path.join(install_folder, "temp", project, str(job), "software"))

    if os.path.exists(os.path.join(install_folder, "temp", project, str(job), "temp", "Reporter")):
        shutil.rmtree(os.path.join(install_folder, "temp", project, str(job), "temp", "Reporter"))
        
    # copy reporter into temp
    if os.path.exists(os.path.join(install_folder, "temp", project, str(job), "software", "Reporter-%s" % settings.reporter_ver)):
        shutil.rmtree(os.path.join(install_folder, "temp", project, str(job), "software", "Reporter-%s" % settings.reporter_ver))
        
    shutil.copytree(os.path.join(install_folder, "software", "Reporter-%s" % settings.reporter_ver), 
                    os.path.join(install_folder, "temp", project, str(job), "software", 
                    "Reporter-%s" % settings.reporter_ver))

    # remove the old output if it exists
    if os.path.exists(os.path.join(settings.data_folder, project, "out", filename, type, "%s_reporter.psdb" % filename)):
        os.remove(os.path.join(settings.data_folder, project, "out", filename, type, "%s_reporter.psdb" % filename))    
    
    if os.path.exists(os.path.join(settings.data_folder, project, "out", filename, type, "r_%s_Default_PSM_Report.txt" % project)):
        os.remove(os.path.join(settings.data_folder, project, "out", filename, type, "r_%s_Default_PSM_Report.txt" % project))
        
    # determine the reference positions
    sample = queue.sample
    labels = LabelChoice.objects.filter(multiplexlabel__sample=sample)
    i = 1
    ref_sample_list = []
    for label in labels:
        if label.tag.t_type == 'Reference':
            ref_sample_list.append(i)
        i += 1
    
    if len(ref_sample_list) == 0:
        write_debug("Warning: no reference sample set so ratios will not be normalized.", job, project)
        
    ref_samples = ",".join([str(ref) for ref in ref_sample_list])

    ############ temporary reporter bug workaround ##############
    # we need to unzip the searchgui_out.zip folder
    write_debug("Extracting data folder from searchgui_out.zip", job, project)
    archive = zipfile.ZipFile(os.path.join(settings.data_folder, project, "out", filename, type, "searchgui_out.zip"))

    for file in archive.namelist():
        if file.startswith('data/'):
            archive.extract(file, os.path.join(settings.data_folder, project, "out", filename, type))
    
    # we need to figure out what the reference samples are
    # we will pull it from the multiplex label information so we need
    # this to be set ahead of time. error out if it's not set for a file
    

        
    success = run_command(["timeout", "600", 
                            "java", "-Xms%s" % settings.memory, "-Xmx%s" % settings.memory, 
                            "-cp", os.path.join(install_folder, "temp", project, str(job), "software", "Reporter-%s" % settings.reporter_ver, "Reporter-%s.jar" % settings.reporter_ver),
                            "eu.isas.reporter.cli.PathSettingsCLI",
                            "-temp_folder", "%s" % os.path.join(install_folder, "temp", project, str(job), "temp", "Reporter"),
                        ], job, project)
    
    if len(ref_samples) >= 1:
        # note that some settings, such as isotope correction have the ability
        # to be changed but aren't supported right now. 
        success = run_command(["timeout", "3600", 
                                "java", "-Xms%s" % settings.memory, "-Xmx%s" % settings.memory, 
                                "-cp", os.path.join(install_folder, "temp", project, str(job), "software", "Reporter-%s" % settings.reporter_ver, "Reporter-%s.jar" % settings.reporter_ver), 
                                "eu.isas.reporter.cli.ReporterCLI",
                                "-id", "%s.psdb" % (os.path.join(settings.data_folder, project, "out", filename, type, filename)),
                                "-out", "%s%s%s_reporter.psdb" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, filename),
                                "-reports", "3",
                                "-report_prefix", "r_",
                                "-ref_samples", "%s" % ref_samples,
                            ], job, project)
    else:
        success = run_command(["timeout", "3600", 
                                "java", "-Xms%s" % settings.memory, "-Xmx%s" % settings.memory, 
                                "-cp", os.path.join(install_folder, "temp", project, str(job), "software", "Reporter-%s" % settings.reporter_ver, "Reporter-%s.jar" % settings.reporter_ver), 
                                "eu.isas.reporter.cli.ReporterCLI",
                                "-id", "%s.psdb" % (os.path.join(settings.data_folder, project, "out", filename, type, filename)),
                                "-out", "%s%s%s_reporter.psdb" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, filename),
                                "-reports", "3",
                                "-report_prefix", "r_",
                            ], job, project)
                            
    if success == 0 or not os.path.exists("%s%s%s_reporter.psdb" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, filename)):
        write_debug("Reporter ReportCLI failed", job, project)
        return False
           
    if (not os.path.exists(r"%s%sr_%s_Default_PSM_Report.txt" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, project))):
        write_debug("Reporter output missing %s" % (os.path.join(settings.data_folder, project, "out", filename, type)), job, project)
        return False
    else:
        end = time.time()
        runtime = end - start
        runtimex = RunTime.objects.get(queue=queue)
        if type == 'profile':
            runtimex.reporter_profile = runtime
        elif type == 'proteome':
            runtimex.reporter_proteome = runtime
        runtimex.save()
        shutil.rmtree(os.path.join(settings.data_folder, project, "out", filename, type, "data"))
        return True