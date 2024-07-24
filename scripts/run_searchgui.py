import shutil
import os
import psutil
import argparse
import time
import zipfile

from django.core.exceptions import ObjectDoesNotExist

from projects.models import Setting, Queue, SearchSetting, ModChoice, EnzymeChoice, RunTime

from .run_command import run_command, write_debug, settings

def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('queue_id', type=int)
    parser.add_argument('job', type=int)
    args2 = parser.parse_args(args)
    
    queue_id = args2.queue_id
        
    run_searchgui(queue_id)

def run_searchgui(queue_id):
    try:
        queue = Queue.objects.get(id=queue_id)
    except ObjectDoesNotExist:
        print("searchgui missing queue_id: %s" % queue_id)
        return False
    
    job = queue.job
    filename = queue.filename
    project = queue.project.name
    
    install_folder = settings.install_folder
    
    start = time.time()

    try:
        searchsetting=SearchSetting.objects.get(project=project)
    except ObjectDoesNotExist:
        write_debug("Missing searchsetting for project: %s." % project, job, project)
        return False

    if queue.status == Queue.Status.SEARCHGUI_PROF:
        if searchsetting.custom_fasta == True:
            type = "custom"
            fasta_file = "%s%s%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", type), os.sep, project, type)     
        else:
            type = "profile"
            fasta_file = "%s%s%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", type), os.sep, project, type)
    elif queue.status == Queue.Status.SEARCHGUI_PROT:
        type = "proteome"
        fasta_file = "%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", type, filename), os.sep, project, filename, type)
    


    if not os.path.exists(os.path.join(install_folder, "software", "SearchGUI-%s" % settings.searchgui_ver, "SearchGUI-%s.jar" % settings.searchgui_ver)):
        write_debug("Missing SearchGUI install.", job, project)
        return False
        
    if not os.path.exists(fasta_file):
        write_debug("Missing %s FASTA file for %s." % (type, filename), job, project)
        return False
        
    mods = searchsetting.mods.all()
    enzymes = searchsetting.enzymes.all()
        
    modchoice = ModChoice.objects.filter(searchsetting=searchsetting)
    mod_list_f = ",".join([mod.mod.name for mod in modchoice if mod.modtype=='Fixed'])
    mod_list_v = ",".join([mod.mod.name for mod in modchoice if mod.modtype=='Variable'])

    enzymechoice = EnzymeChoice.objects.filter(searchsetting=searchsetting)
    enzyme_list_name = ",".join([enzyme.enzyme.name for enzyme in enzymechoice])
    enzyme_list_specificity = ",".join([str(enzyme.specificity) for enzyme in enzymechoice])
    enzyme_list_mc = ",".join([str(enzyme.mc) for enzyme in enzymechoice])

    # remove the temp dir to prepare for running the programs
    if os.path.exists(os.path.join(install_folder, "temp", project, str(job))):
        shutil.rmtree(os.path.join(install_folder, "temp", project, str(job)))
    os.makedirs(os.path.join(install_folder, "temp", project, str(job)))

    # remake the temp software folder
    if not os.path.exists(os.path.join(install_folder, "temp", project, str(job), "software")):
        os.makedirs(os.path.join(install_folder, "temp", project, str(job), "software"))
        
    # copy searchgui into temp
    shutil.copytree(os.path.join(install_folder, "software", "SearchGUI-%s" % settings.searchgui_ver), 
                    os.path.join(install_folder, "temp", project, str(job), "software", "SearchGUI-%s" % settings.searchgui_ver))
    
    # remove the parameter file
    if os.path.exists("%s%s%s_%s.par" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, project, type)):
        os.remove("%s%s%s_%s.par" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, project, type))

    # remove the folder we check results in
    if os.path.exists(os.path.join(settings.data_folder, project, "out", filename, type, "temp")):
        shutil.rmtree(os.path.join(settings.data_folder, project, "out", filename, type, "temp"))
        
    # remove the old output if it exists
    if os.path.exists(os.path.join(settings.data_folder, project, "out", filename, type, "searchgui_out.zip")):
        os.remove(os.path.join(settings.data_folder, project, "out", filename, type, "searchgui_out.zip"))

    write_debug("Starting SearchGUI PathSettingsCLI: %s" % (os.path.join(settings.data_folder, project)), job, project)
    success = run_command(["timeout", "86400", 
                    "java", "-Xms%s" % settings.memory, "-Xmx%s" % settings.memory, 
                    "-cp", os.path.join(install_folder, "temp", project, str(job), "software", "SearchGUI-%s" % settings.searchgui_ver, "SearchGUI-%s.jar" % settings.searchgui_ver), 
                    "eu.isas.searchgui.cmd.PathSettingsCLI",
                    "-temp_folder", "%s" % os.path.join(install_folder, "temp", project, str(job), "temp", "SearchGUI"),
                    "-identification_parameters", "%s" % os.path.join(install_folder, "temp", project, str(job), "temp", "SearchGUI"),
                    "-gene_mapping", "%s" % os.path.join(install_folder, "temp", project, str(job), "temp", "SearchGUI"),
                    "-pride_annotation", "%s" % os.path.join(install_folder, "temp", project, str(job), "temp", "SearchGUI"),
                    "-use_log_folder", "0"
                    ], job, project) 

    if (success == 0):
            write_debug("SearchGUI PathSettingsCLI failed", job, project)
            return False
            
    # put this in the temp folder from now on
    write_debug("Generating SearchGUI PAR file.", job, project)
    command = ["timeout", "600",
               "java", "-Xms%s" % settings.memory, "-Xmx%s" % settings.memory, 
               "-cp", os.path.join(install_folder, "temp", project, str(job), "software", "SearchGUI-%s" % settings.searchgui_ver, "SearchGUI-%s.jar" % settings.searchgui_ver),
               "eu.isas.searchgui.cmd.IdentificationParametersCLI",
               "-out", "%s%s%s_%s.par" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, project, type),
               "-min_charge", "%s" % searchsetting.min_charge,
               "-max_charge", "%s" % searchsetting.max_charge,
               "-import_peptide_length_min", "%s" % searchsetting.min_peptide_length, # min peptide length
               "-import_peptide_length_max", "%s" % searchsetting.max_peptide_length, #max peptide length
               "-psm_fdr", "%s" % searchsetting.psm_fdr, # fdr precent
               "-peptide_fdr", "%s" % searchsetting.peptide_fdr,
               "-protein_fdr", "%s" % searchsetting.protein_fdr,
               "-useGeneMapping", "0",
               "-updateGeneMapping", "0",
               "-prec_tol", "%s" % searchsetting.prec_tol,
               "-prec_ppm", "%s" % searchsetting.prec_ppm,
               "-frag_tol", "%s" % searchsetting.frag_tol,
               "-frag_ppm", "%s" % searchsetting.frag_ppm,
               "-min_isotope", "%s" % searchsetting.isotope_min,
               "-max_isotope", "%s" % searchsetting.isotope_max,
               "-simplify_groups", "0",               
               "-msgf_instrument", "%s" % searchsetting.instrument,
               "-msgf_fragmentation", "%s" % searchsetting.fragmentation,
               "-msgf_num_tasks", "%s" % searchsetting.msgf_num_tasks,
               "-msgf_min_pep_length", "%s" % searchsetting.min_peptide_length,
               "-msgf_max_pep_length", "%s" % searchsetting.max_peptide_length,               
               "-myrimatch_min_pep_length", "%s" % searchsetting.min_peptide_length,
               "-myrimatch_max_pep_length", "%s" % searchsetting.max_peptide_length,
               "-omssa_min_pep_length", "%s" % searchsetting.min_peptide_length,
               "-omssa_max_pep_length", "%s" % searchsetting.max_peptide_length,
               "-comet_min_pep_length", "%s" % searchsetting.min_peptide_length,
               "-comet_max_pep_length", "%s" % searchsetting.max_peptide_length,
               "-comet_batch_size", "%s" % searchsetting.comet_batch_size,               
               "-meta_morpheus_min_pep_length", "%s" % searchsetting.min_peptide_length,
               "-meta_morpheus_max_pep_length", "%s" % searchsetting.max_peptide_length,
               "-sage_min_pep_length", "%s" % searchsetting.min_peptide_length,
               "-sage_max_pep_length", "%s" % searchsetting.max_peptide_length,
               "-digestion", "%s" % searchsetting.digestion
               ]
   
    if len(mod_list_f) > 0:
        command.extend(["-fixed_mods", mod_list_f])
        
    if len(mod_list_v) > 0:
        command.extend(["-variable_mods", mod_list_v])
        
    if len(enzyme_list_name) > 0:
        command.extend(["-digestion", "0"])
        command.extend(["-enzyme", "%s" % enzyme_list_name])
        command.extend(["-specificity", "%s" % enzyme_list_specificity])
        command.extend(["-mc", "%s" % enzyme_list_mc])
        
        
    # timeout after 10 mins
    success = run_command(command, job, project)
    
    if success == 0 or not os.path.exists("%s%s%s_%s.par" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, project, type)):
        write_debug("Missing searchgui PAR file.", job, project)
        return False
    # this is a workaround for the latest searchgui
    else:
        if not os.path.exists(os.path.join(install_folder, "temp", project, str(job), "identification_parameters_4")):
            os.makedirs(os.path.join(install_folder, "temp", project, str(job), "identification_parameters_4"))
        if not os.path.exists("%s%s%s_%s.par" % (os.path.join(settings.install_folder, "temp", project, str(job), "identification_parameters_4"), os.sep, project, type)):
            shutil.copy("%s%s%s_%s.par" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, project, type),
                        "%s%s%s_%s.par" % (os.path.join(settings.install_folder, "temp", project, str(job), "identification_parameters_4"), os.sep, project, type))
    if settings.threads == -1:
        threads = psutil.cpu_count()
    else:
        threads = settings.threads
        
    write_debug("Starting searchgui: %s" % (filename), job, project)            
    if type == "profile" or type == "custom":
        xtandem = int(searchsetting.xtandem_profile)
        msgf = int(searchsetting.msgf_profile)
        comet = int(searchsetting.comet_profile)
        omssa = int(searchsetting.omssa_profile)
        metamorpheus = int(searchsetting.metamorpheus_profile)
        myrimatch = int(searchsetting.myrimatch_profile)
        sage = int(searchsetting.sage_profile)
    else:
        xtandem = int(searchsetting.xtandem_proteome)
        msgf = int(searchsetting.msgf_proteome)
        comet = int(searchsetting.comet_proteome)
        omssa = int(searchsetting.omssa_proteome)
        metamorpheus = int(searchsetting.metamorpheus_proteome)
        myrimatch = int(searchsetting.myrimatch_proteome)
        sage = int(searchsetting.sage_proteome)
    
    # this is a workaround for myrimatch
    if myrimatch == 1:
        #success = run_command(["export LC_ALL=C"], job, project)
        os.environ['LC_ALL'] = 'C'
    
    # issues with getting these two to run so disabled
    ms_amanda = 0
    tide = 0
    
    if (xtandem == 0) and (msgf == 0) and (comet == 0) and (omssa == 0) and (metamorpheus == 0) and (myrimatch == 0) and (sage == 0):
        write_debug("No search engines are selected. Make sure at least one search engine is selected for both profile and proteome steps.", job, project)
        return False
        
    success = run_command(["timeout", "172800", 
                            "java", "-Xms%s" % settings.memory, "-Xmx%s" % settings.memory, 
                            "-cp", os.path.join(install_folder, "temp", project, str(job), "software", "SearchGUI-%s" % settings.searchgui_ver, "SearchGUI-%s.jar" % settings.searchgui_ver), 
                            "eu.isas.searchgui.cmd.SearchCLI",
                            "-spectrum_files", "%s.mzML" % (os.path.join(settings.data_folder, project, "out", filename, filename)),
                            "-output_folder", "%s" % (os.path.join(settings.data_folder, project, "out", filename, type)),
                            "-id_params", "%s%s%s_%s.par" % (os.path.join(settings.data_folder, project, "out", filename, type), os.sep, project, type),
                            "-fasta_file", "%s" % fasta_file,
                            "-xtandem", "%s" % xtandem,
                            "-msgf", "%s" % msgf,
                            "-omssa", "%s" % omssa,
                            "-comet", "%s" % comet,
                            "-meta_morpheus", "%s" % metamorpheus,
                            "-myrimatch", "%s" % myrimatch,
                            "-ms_amanda", "%s" % ms_amanda,
                            "-sage", "%s" % sage,
                            "-tide", "%s" % tide,
                            "-output_option", "0",
                            "-output_data", "1",
                            "-output_date", "0",
                            "-threads", "%s" % threads
                          ], job, project)

    # need database update here to change the attempt count and status
    if success == 0 or not os.path.exists(r"%s" % (os.path.join(settings.data_folder, project, "out", filename, type, "searchgui_out.zip"))):
        write_debug("missing searchgui output: %s" % (os.path.join(settings.data_folder, project, "out", filename, type, "searchgui_out.zip")), job, project)
        return False
    
    # to verify that the search engines work, we'll have to unzip the searchgui file and look
    write_debug("Verifying correct search engine outputs.", job, project)
    if not os.path.exists(os.path.join(settings.data_folder, project, "out", filename, type, "temp")):
        os.makedirs(os.path.join(settings.data_folder, project, "out", filename, type, "temp"))
        
    with (zipfile.ZipFile(os.path.join(settings.data_folder, project, "out", 
                          filename, type, "searchgui_out.zip"),"r")) as zip_ref:
        zip_ref.extractall(os.path.join(settings.data_folder, project, "out", 
                                        filename, type, "temp"))

    if (msgf == 1):
        if not os.path.exists(
            os.path.join(settings.data_folder, project, "out", filename, 
                type, "temp", "%s.msgf.mzid.gz" % filename)
        ):
            write_debug("Missing MSGF output from run_searchgui.", 
                job, project)
            return False
    if (comet == 1):
        if not os.path.exists(
            os.path.join(settings.data_folder, project, "out", filename, 
            type, "temp", "%s.comet.pep.xml.gz" % filename)
        ):
            write_debug("Missing Comet output from run_searchgui.", 
                job, project)
            return False
    if (omssa == 1):
        if not os.path.exists(
            os.path.join(settings.data_folder, project, "out", filename, 
            type, "temp", "%s.omx.gz" % filename)
        ):
            write_debug("Missing OMSSA output from run_searchgui.", 
                job, project)
            return False
    if (xtandem == 1):
        if not os.path.exists(
            os.path.join(settings.data_folder, project, "out", filename, 
            type, "temp", "%s.t.xml.gz" % filename)
        ):
            write_debug("Missing XTandem output from run_searchgui.", 
                job, project)
            return False
    if (myrimatch == 1):
        if not os.path.exists(
            os.path.join(settings.data_folder, project, "out", filename, 
            type, "temp", "%s.myrimatch.mzid.gz" % filename)
        ):
            write_debug("Missing MyriMatch output from run_searchgui.", 
                job, project)
            return False
    if (metamorpheus == 1):
        if not os.path.exists(
            os.path.join(settings.data_folder, project, "out", filename, 
            type, "temp", "%s.mzID.gz" % filename)
        ):
            write_debug("Missing MetaMorpheus output from run_searchgui.", 
                job, project)
            return False
    if (sage == 1):
        if not os.path.exists(
            os.path.join(settings.data_folder, project, "out", filename,
            type, "temp", "%s.sage.tsv.gz" % filename)
        ):
            write_debug("Missing Sage output from run_searchgui.",
                job, project)
            return False
            
    if os.path.exists(os.path.join(settings.data_folder, project, "out", filename, type, "temp")):
        shutil.rmtree(os.path.join(settings.data_folder, project, "out", filename, type, "temp"))
    
    end = time.time()
    runtime = end-start
    runtimex = RunTime.objects.get(queue=queue)
    if type == 'profile' or type == 'custom':
        runtimex.searchgui_profile = runtime
    elif type == 'proteome':
        runtimex.searchgui_proteome = runtime
    runtimex.save()
    
    return True
