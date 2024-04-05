# runs msconvert (thermorawfileparser)
# python3 manage.py runscript run_msconvert --script-args 1

import shutil
import os
import argparse
import time

from django.core.exceptions import ObjectDoesNotExist

from projects.models import Setting, Queue, RunTime, SearchSetting

from .run_command import run_command, write_debug, settings

def run(*args):
    # here if we know queue_id, we don't care about project
    parser = argparse.ArgumentParser()
    parser.add_argument('queue_id', type=int)
    args2 = parser.parse_args(args)
    
    queue_id = args2.queue_id
        
    run_msconvert(queue_id)
      
# args: filename
# this will only get called when processing the queue
def run_msconvert(queue_id):
    try:
        queue = Queue.objects.get(id=queue_id)
    except ObjectDoesNotExist:
        print("run_msconvert missing queue entry for ID: %s." % queue_id)
        return False
    
    success = 0
    filename = queue.filename
    project = queue.project.name
    install_folder = settings.install_folder
    job = queue.job

    try:
        searchsetting=SearchSetting.objects.get(project=project)
    except ObjectDoesNotExist:
        write_debug("Missing searchsetting for project: %s." % project, job, project)
        return False
    
    if not os.path.exists(os.path.join(install_folder, "software", "SearchGUI-%s" % settings.searchgui_ver, "SearchGUI-%s.jar" % settings.searchgui_ver)):
        write_debug("Missing SearchGUI install.", job, project)
        return False
        
    start = time.time()

    if os.path.exists(os.path.join(install_folder, "temp", project, str(job))):
        shutil.rmtree(os.path.join(install_folder, "temp", project, str(job)))
    os.makedirs(os.path.join(install_folder, "temp", project, str(job)))
  
    if os.path.exists(os.path.join(settings.data_folder, project, "raw", "%s.raw" % filename)):
        # copy searchgui into temp
        shutil.copytree(os.path.join(install_folder, "software", "SearchGUI-%s" % settings.searchgui_ver), 
                        os.path.join(install_folder, "temp", project, str(job), "software", "SearchGUI-%s" % settings.searchgui_ver))
                        
        write_debug("Starting msconvert: %s" % (filename), job, project)
         
        # remove old file
        if os.path.exists(r"%s.mzML" % (os.path.join(settings.data_folder, project, "out", filename, filename))):
            os.remove(r"%s.mzML" % (os.path.join(settings.data_folder, project, "out", filename, filename)))
        
        success = run_command(["timeout", "1800", 
                                "mono", os.path.join(install_folder, "temp", project, str(job), "software", "SearchGUI-%s" % settings.searchgui_ver, "resources", "ThermoRawFileParser", "ThermoRawFileParser.exe"),
                                "-i=%s" % os.path.join(settings.data_folder, project, "raw", "%s.raw" % filename),
                                "-o=%s" % os.path.join(settings.data_folder, project, "out", filename),
                                #"-e", # ignore instrument errors
                                #"-b=%s" % OUTPUT FILE,
                                "-f=2", # format: 0 MGF, 1 mzml, 2 imzml
                                #"-g", # gzip
                                #"-m=0", # meta data format
                                #"-c=%s" % META_FILE,
                                #"-z", # no zlib compression
                              ], job, project)

        if os.path.exists(os.path.join(install_folder, "temp", project, str(job))):
            shutil.rmtree(os.path.join(install_folder, "temp", project, str(job)))                              
    # if ends with .mzML, assume success and copy file
    elif os.path.exists(os.path.join(settings.data_folder, project, "raw", "%s.mzML" % filename)):
        write_debug("Found an mzML file already converted. Using it instead.", job, project)
        try:
            shutil.copyfile(os.path.join(settings.data_folder, project, "raw", "%s.mzML" % filename), 
                            os.path.join(settings.data_folder, project, "out", filename, "%s.mzML" % filename))
        except:
            write_debug("Failed to copy the mzML file.", job, project)
            return False
        else:
            success = 1
    
    else:
        write_debug("Queue file does not exist.", job, project)
      
    # need database update here to change the attempt count and status
    if success == 0 or not os.path.exists(r"%s.mzML" % (os.path.join(settings.data_folder, project, "out", filename, filename))):
        write_debug("Missing msconvert output: %s.mzML" % (os.path.join(settings.data_folder, project, "out", filename, filename)), job, project)
        queue.error += 1
        queue.save()
        return False
    else:
        end = time.time()
        runtime = end - start
        queue.error = 0
        queue.save()
        runtimex = RunTime.objects.get(queue=queue)
        runtimex.msconvert = runtime
        runtimex.save()
        return True
