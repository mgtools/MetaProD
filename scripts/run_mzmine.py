import shutil
import os
import argparse
import time
import zipfile
import pandas as pd
import xml.etree.ElementTree as ET

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
        
    run_mzmine(queue_id)

def run_mzmine(queue_id):
    ''' runs mzmine for peak areas for identified psms '''
    try:
        queue = Queue.objects.get(id=queue_id)
    except ObjectDoesNotExist:
        print("MZmine missing queue_id: %s" % queue_id)
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

    if queue.status == Queue.Status.MZMINE_PROT:
        type = "proteome"
    else:
        type = "profile"

    if not os.path.exists(os.path.join(install_folder, "software", 
                          "MZmine-%s" % settings.mzmine_ver, "lib", "app", 
                          "mzmine3-%s.jar" % settings.mzmine_ver)):
        write_debug("Missing MZmine install.", job, project)
        return False
        
    # remake the temp folder
    if not os.path.exists(os.path.join(install_folder, "temp", project, str(job), "software")):
        os.makedirs(os.path.join(install_folder, "temp", project, str(job), "software"))

    if os.path.exists(os.path.join(install_folder, "temp", project, str(job), "temp", "MZmine")):
        shutil.rmtree(os.path.join(install_folder, "temp", project, str(job), "temp", "MZmine"))
        
    # copy reporter into temp
    if os.path.exists(os.path.join(install_folder, "temp", project, str(job), "software", 
                                   "MZmine-%s" % settings.mzmine_ver)):
        shutil.rmtree(os.path.join(install_folder, "temp", project, str(job), "software", 
                                   "MZmine-%s" % settings.mzmine_ver))
        
    shutil.copytree(os.path.join(install_folder, "software", "MZmine-%s" % settings.mzmine_ver), 
                    os.path.join(install_folder, "temp", project, str(job), "software", 
                    "MZmine-%s" % settings.mzmine_ver))
                    
    # remove the old output if it exists
    if os.path.exists(os.path.join(settings.data_folder, project, "out", filename, type, "%s_mzexport.csv" % filename)):
        os.remove(os.path.join(settings.data_folder, project, "out", filename, type, "%s_mzexport.csv" % filename))    
    
    # also remove the csv generated from the psm results
    if os.path.exists(os.path.join(settings.data_folder, project, "out", filename, type, "%s_mzmine_tpd.csv" % filename)):
        os.remove(os.path.join(settings.data_folder, project, "out", filename, type, "%s_mzmine_tpd.csv" % filename))

    # remove the batch file
    if os.path.exists(os.path.join(settings.data_folder, project, "out", filename, type, "%s_mzmine_batch.xml" % filename)):
        os.remove(os.path.join(settings.data_folder, project, "out", filename, type, "%s_mzmine_batch.xml" % filename))
        
    # now generate the csv to import into mzmine
    # results aren't in the database yet so we need to use the ps output
    ps_output = pd.read_csv(r"%s%sps_%s_Default_PSM_Report.txt" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, project), sep='\t')
    mz_output = pd.DataFrame()
    mz_output['mz'] = round(ps_output['m/z'], 4)
    mz_output['rt'] = ps_output['RT']
    mz_output['rt'] = round(mz_output['rt'] / 60, 2)
    mz_output['name'] = ps_output['Spectrum Title']
    # mzmine bug workaround
    mz_output.loc[len(mz_output.index)] = ['0','0','None']
    mz_output.to_csv(r"%s%s%s_mzmine_tpd.csv" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, filename), sep=',')

    # then copy the xml to run mzmine
    if not os.path.exists(os.path.join(os.getcwd(), "scripts", "mzmine_batch.xml")):
        write_debug("Missing mzmine_batch.xml", job, project)
        return False
    else:
        shutil.copy(os.path.join(os.getcwd(), "scripts", "mzmine_batch.xml"), 
                    os.path.join(settings.data_folder, project, "out", filename, type, "%s_mzmine_batch.xml" % filename))

    # now edit the xml to add the filename and eventually the settings    
    tree = ET.parse(os.path.join(settings.data_folder, project, "out", filename, type, "%s_mzmine_batch.xml" % filename))
    root = tree.getroot()
    # update the filename for importing
    for batchstep in root.findall('batchstep'):
        if batchstep.attrib['method'] == 'io.github.mzmine.modules.io.import_rawdata_mzml.MSDKmzMLImportModule':
            for parameter in batchstep.findall('parameter'):
                if parameter.attrib['name'] == 'File names':
                    file = parameter.find('file')
                    file.text = os.path.join(settings.data_folder, project, "out", filename, "%s.mzML" % filename)
                    file.set('updated', 'yes')
                    
        # here we should also update various settings but we'll add that later
        # also need to update the import of the tpd.csv
        elif batchstep.attrib['method'] == 'io.github.mzmine.modules.dataprocessing.featdet_targeted.TargetedFeatureDetectionModule':
            for parameter in batchstep.findall('parameter'):
                if parameter.attrib['name'] == 'Database file':
                    file = parameter.find('current_file')
                    file.text = os.path.join(settings.data_folder, project, "out", filename, type, "%s_mzmine_tpd.csv" % filename)
                    file.set('updated', 'yes')
                    # should also update options here for tolerances, etc
                elif parameter.attrib['name'] == 'Intensity tolerance':
                    parameter.text = "%s" % searchsetting.mzmine_tpd_intensity
                elif parameter.attrib['name'] == 'm/z tolerance':
                    mztol = parameter.find('absolutetolerance')
                    mztol.text = "%s" % searchsetting.mzmine_tpd_mztolerance
                    ppmtol = parameter.find('ppmtolerance')
                    ppmtol.text = "%s" % searchsetting.mzmine_tpd_ppmtolerance
        # finally update the final export dir
        elif batchstep.attrib['method'] == 'io.github.mzmine.modules.io.export_features_csv.CSVExportModularModule':
            for parameter in batchstep.findall('parameter'):
                if parameter.attrib['name'] == 'Filename':
                    file = parameter.find('current_file')
                    file.text = os.path.join(settings.data_folder, project, "out", filename, type, "%s_mzexport.csv" % filename)
                    file.set('updated', 'yes')
                    
    tree.write(os.path.join(settings.data_folder, project, "out", filename, type, "%s_mzmine_batch.xml" % filename))
    
    # no mzmine tpd so mzmine will fail
    # make an empty file
    success = 0
    if len(mz_output.index) <= 1:
        mzexport_csv = pd.DataFrame(columns=['1', '2', '3'])
        mzexport_csv.loc[len(mzexport_csv.index)] = ['','','']
        mzexport_csv.to_csv("%s%s%s_mzexport.csv" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, filename), sep=',')
        success = 1
    else:
        # now run the command
        # java --enable-preview -cp "/home/jamie/metaprod_projects/software/MZmine-3.3.0/lib/app/*" io.github.mzmine.main.MZmineCore -b 
        success = run_command(["timeout", "86400", 
                                "java", "-Xms%s" % settings.memory, "-Xmx%s" % settings.memory,
                                "--enable-preview",
                                "-cp", os.path.join(install_folder, "temp", project, str(job), "software", "MZmine-%s" % settings.mzmine_ver, "lib", "app", "*"),
                                "io.github.mzmine.main.MZmineCore",
                                "-b", "%s" % os.path.join(settings.data_folder, project, "out", filename, type, "%s_mzmine_batch.xml" % filename),
                            ], job, project)
                                
    if success == 0 or not os.path.exists("%s%s%s_mzexport.csv" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, filename)):
        write_debug("MZmine failed", job, project)
        return False
    else:
        end = time.time()
        runtime = end - start
        runtimex = RunTime.objects.get(queue=queue)
        if type == 'profile':
            runtimex.mzmine_profile = runtime
        elif type == 'proteome':
            runtimex.mzmine_proteome = runtime
        runtimex.save()
        return True
