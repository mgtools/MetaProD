import shutil
import os
import psutil
import argparse
import time
import zipfile

from django.core.exceptions import ObjectDoesNotExist

from projects.models import Setting, Queue, SearchSetting, ModChoice, EnzymeChoice, RunTime, EngineStatus

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
            fasta_type = "custom"
            fasta_file = "%s%s%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", fasta_type), os.sep, project, fasta_type)     
        else:
            fasta_type = "profile"
            fasta_file = "%s%s%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", fasta_type), os.sep, project, fasta_type)
    elif queue.status == Queue.Status.SEARCHGUI_PROT:
        fasta_type = "proteome"
        fasta_file = "%s%s%s_%s_%s_concatenated_target_decoy.fasta" % (os.path.join(settings.data_folder, project, "fasta", fasta_type, filename), os.sep, project, filename, fasta_type)
    


    if not os.path.exists(os.path.join(install_folder, "software", "SearchGUI-%s" % settings.searchgui_ver, "SearchGUI-%s.jar" % settings.searchgui_ver)):
        write_debug("Missing SearchGUI install.", job, project)
        return False
        
    if not os.path.exists(fasta_file):
        write_debug("Missing %s FASTA file for %s." % (fasta_type, filename), job, project)
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
    if os.path.exists("%s%s%s_%s.par" % (os.path.join(settings.data_folder, project, "out", filename, fasta_type), os.sep, project, fasta_type)):
        os.remove("%s%s%s_%s.par" % (os.path.join(settings.data_folder, project, "out", filename, fasta_type), os.sep, project, fasta_type))

    # remove the folder we check results in
    if os.path.exists(os.path.join(settings.data_folder, project, "out", filename, fasta_type, "temp")):
        shutil.rmtree(os.path.join(settings.data_folder, project, "out", filename, fasta_type, "temp"))
        
    # remove the old output if it exists
    if os.path.exists(os.path.join(settings.data_folder, project, "out", filename, fasta_type, "searchgui_out.zip")):
        os.remove(os.path.join(settings.data_folder, project, "out", filename, fasta_type, "searchgui_out.zip"))

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
               "-out", "%s%s%s_%s.par" % (os.path.join(settings.data_folder, project, "out", filename, fasta_type), os.sep, project, fasta_type),
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
    
    if success == 0 or not os.path.exists("%s%s%s_%s.par" % (os.path.join(settings.data_folder, project, "out", filename, fasta_type), os.sep, project, fasta_type)):
        write_debug("Missing searchgui PAR file.", job, project)
        return False
    # this is a workaround for the latest searchgui
    else:
        if not os.path.exists(os.path.join(install_folder, "temp", project, str(job), "identification_parameters_4")):
            os.makedirs(os.path.join(install_folder, "temp", project, str(job), "identification_parameters_4"))
        if not os.path.exists("%s%s%s_%s.par" % (os.path.join(settings.install_folder, "temp", project, str(job), "identification_parameters_4"), os.sep, project, fasta_type)):
            shutil.copy("%s%s%s_%s.par" % (os.path.join(settings.data_folder, project, "out", filename, fasta_type), os.sep, project, fasta_type),
                        "%s%s%s_%s.par" % (os.path.join(settings.install_folder, "temp", project, str(job), "identification_parameters_4"), os.sep, project, fasta_type))
    if settings.threads == -1:
        threads = psutil.cpu_count()
    else:
        threads = settings.threads
        
    write_debug("Starting searchgui: %s" % (filename), job, project)            
    if fasta_type == "profile" or fasta_type == "custom":
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

    # run the search
    # check for success
    # if failure, try that engine again until max retries
    # update status table with results
    
    enginestatus = EngineStatus.objects.get(queue=queue)
    enginestatus.comet_tries = 0
    enginestatus.xtandem_tries = 0
    enginestatus.omssa_tries = 0
    enginestatus.msgf_tries = 0
    enginestatus.myrimatch_tries = 0
    enginestatus.metamorpheus_tries = 0
    enginestatus.sage_tries = 0
    if fasta_type == "profile":
        enginestatus.comet_profile = False
        enginestatus.msgf_profile = False
        enginestatus.xtandem_profile = False
        enginestatus.xtandem_profile = False
        enginestatus.myrimatch_profile = False
        enginestatus.metamorpheus_profile = False
        enginestatus.sage_profile = False
    elif fasta_type == "proteome":
        enginestatus.comet_proteome = False
        enginestatus.msgf_proteome = False
        enginestatus.xtandem_proteome = False
        enginestatus.xtandem_proteome = False
        enginestatus.myrimatch_proteome = False
        enginestatus.metamorpheus_proteome = False
        enginestatus.sage_proteome = False    
    
    enginestatus.save()
    
    def run_search(**engine):
        success = run_command(["timeout", "172800", 
                                "java", "-Xms%s" % settings.memory, "-Xmx%s" % settings.memory, 
                                "-cp", os.path.join(install_folder, "temp", project, str(job), "software", "SearchGUI-%s" % settings.searchgui_ver, "SearchGUI-%s.jar" % settings.searchgui_ver), 
                                "eu.isas.searchgui.cmd.SearchCLI",
                                "-spectrum_files", "%s.mzML" % (os.path.join(settings.data_folder, project, "out", filename, filename)),
                                "-output_folder", "%s" % (os.path.join(settings.data_folder, project, "out", filename, fasta_type)),
                                "-id_params", "%s%s%s_%s.par" % (os.path.join(settings.data_folder, project, "out", filename, fasta_type), os.sep, project, fasta_type),
                                "-fasta_file", "%s" % fasta_file,
                                "-xtandem", "%s" % engine['xtandem'],
                                "-msgf", "%s" % engine['msgf'],
                                "-omssa", "%s" % engine['omssa'],
                                "-comet", "%s" % engine['comet'],
                                "-meta_morpheus", "%s" % engine['metamorpheus'],
                                "-myrimatch", "%s" % engine['myrimatch'],
                                "-ms_amanda", "%s" % engine['ms_amanda'],
                                "-sage", "%s" % engine['sage'],
                                "-tide", "%s" % engine['tide'],
                                "-output_option", "3",
                                "-output_data", "0",
                                "-output_date", "0",
                                "-threads", "%s" % threads
                            ], job, project)
                            
        if (engine['comet'] == 1):
            if not os.path.exists(
                os.path.join(settings.data_folder, project, "out", filename, 
                    fasta_type, "%s.comet.pep.xml.gz" % filename)
            ):
                write_debug("Missing Comet output from run_searchgui.", 
                    job, project)
                enginestatus.comet_tries += 1
                enginestatus.save()
                
                if enginestatus.comet_tries < settings.max_retries:
                    run_search(xtandem=0, msgf=0, comet=1, omssa=0, metamorpheus=0, myrimatch=0, sage=0, ms_amanda=0, tide=0)
                    
            else:
                enginestatus.comet_tries = 0
                if fasta_type == "profile":                
                    enginestatus.comet_profile = True
                elif fasta_type == "proteome":
                    enginestatus.comet_proteome = True
                enginestatus.save()

        if (engine['msgf'] == 1):
            if not os.path.exists(
                os.path.join(settings.data_folder, project, "out", filename, 
                    fasta_type, "%s.msgf.mzid.gz" % filename)
            ):
                write_debug("Missing MSGF+ output from run_searchgui.", 
                    job, project)
                enginestatus.msgf_tries += 1
                enginestatus.save()
                
                if enginestatus.msgf_tries < settings.max_retries:
                    run_search(xtandem=0, msgf=1, comet=0, omssa=0, metamorpheus=0, myrimatch=0, sage=0, ms_amanda=0, tide=0)
                    
            else:
                enginestatus.msgf_tries = 0
                if fasta_type == "profile":                
                    enginestatus.msgf_profile = True
                elif fasta_type == "proteome":
                    enginestatus.msgf_proteome = True                
                enginestatus.save()
                
 
        if (engine['xtandem'] == 1):
            if not os.path.exists(
                os.path.join(settings.data_folder, project, "out", filename, 
                    fasta_type, "%s.t.xml.gz" % filename)
            ):
                write_debug("Missing XTandem output from run_searchgui.", 
                    job, project)
                enginestatus.xtandem_tries += 1
                enginestatus.save()
                
                if enginestatus.xtandem_tries < settings.max_retries:
                    run_search(xtandem=1, msgf=0, comet=0, omssa=0, metamorpheus=0, myrimatch=0, sage=0, ms_amanda=0, tide=0)
                    
            else:
                enginestatus.xtandem_tries = 0
                if fasta_type == "profile":                
                    enginestatus.xtandem_profile = True
                elif fasta_type == "proteome":
                    enginestatus.xtandem_proteome = True
                enginestatus.save()
 
        if (engine['omssa'] == 1):
            if not os.path.exists(
                os.path.join(settings.data_folder, project, "out", filename, 
                    fasta_type, "%s.omx.gz" % filename)
            ):
                write_debug("Missing OMSSA output from run_searchgui.", 
                    job, project)
                enginestatus.omssa_tries += 1
                enginestatus.save()
                
                if enginestatus.omssa_tries < settings.max_retries:
                    run_search(xtandem=0, msgf=0, comet=0, omssa=1, metamorpheus=0, myrimatch=0, sage=0, ms_amanda=0, tide=0)
                    
            else:
                enginestatus.omssa_tries = 0
                if fasta_type == "profile":                
                    enginestatus.omssa_profile = True
                elif fasta_type == "proteome":
                    enginestatus.omssa_proteome = True
                enginestatus.save()

        if (engine['myrimatch'] == 1):
            if not os.path.exists(
                os.path.join(settings.data_folder, project, "out", filename, 
                    fasta_type, "%s.myrimatch.mzid.gz" % filename)
            ):
                write_debug("Missing MyriMatch output from run_searchgui.", 
                    job, project)
                enginestatus.myrimatch_tries += 1
                enginestatus.save()
                
                if enginestatus.myrimatch_tries < settings.max_retries:
                    run_search(xtandem=0, msgf=0, comet=0, omssa=0, metamorpheus=0, myrimatch=1, sage=0, ms_amanda=0, tide=0)
                    
            else:
                enginestatus.myrimatch_tries = 0
                if fasta_type == "profile":
                    enginestatus.myrimatch_profile = True
                elif fasta_type == "proteome":
                    enginestatus.myrimatch_proteome = True
                enginestatus.save()

        if (engine['metamorpheus'] == 1):
            if not os.path.exists(
                os.path.join(settings.data_folder, project, "out", filename, 
                    fasta_type, "%s.mzID.gz" % filename)
            ):
                write_debug("Missing Metamorpheus output from run_searchgui.", 
                    job, project)
                enginestatus.metamorpheus_tries += 1
                enginestatus.save()
                
                if enginestatus.metamorpheus_tries < settings.max_retries:
                    run_search(xtandem=0, msgf=0, comet=0, omssa=0, metamorpheus=1, myrimatch=0, sage=0, ms_amanda=0, tide=0)
                    
            else:
                enginestatus.metamorpheus_tries = 0
                if fasta_type == "profile":                
                    enginestatus.metamorpheus_profile = True
                elif fasta_type == "proteome":
                    enginestatus.metamorpheus_proteome = True
                enginestatus.save()

        if (engine['sage'] == 1):
            if not os.path.exists(
                os.path.join(settings.data_folder, project, "out", filename, 
                    fasta_type, "%s.sage.tsv.gz" % filename)
            ):
                write_debug("Missing Sage output from run_searchgui.", 
                    job, project)
                enginestatus.sage_tries += 1
                enginestatus.save()
                
                if enginestatus.sage_tries < settings.max_retries:
                    run_search(xtandem=0, msgf=0, comet=0, omssa=0, metamorpheus=0, myrimatch=0, sage=1, ms_amanda=0, tide=0)
                    
            else:
                enginestatus.sage_tries = 0
                enginestatus.sage_profile=True
                enginestatus.save()
                
    if xtandem == 1:
        run_search(xtandem=1, msgf=0, comet=0, omssa=0, metamorpheus=0, myrimatch=0, sage=0, ms_amanda=0, tide=0)
        
    if msgf == 1:
        run_search(xtandem=0, msgf=1, comet=0, omssa=0, metamorpheus=0, myrimatch=0, sage=0, ms_amanda=0, tide=0)
        
    if comet == 1:
        run_search(xtandem=0, msgf=0, comet=1, omssa=0, metamorpheus=0, myrimatch=0, sage=0, ms_amanda=0, tide=0)
        
    if omssa == 1:
        run_search(xtandem=0, msgf=0, comet=0, omssa=1, metamorpheus=0, myrimatch=0, sage=0, ms_amanda=0, tide=0)
        
    if metamorpheus == 1:
        run_search(xtandem=0, msgf=0, comet=0, omssa=0, metamorpheus=1, myrimatch=0, sage=0, ms_amanda=0, tide=0)
        
    if myrimatch == 1:
        run_search(xtandem=0, msgf=0, comet=0, omssa=0, metamorpheus=0, myrimatch=1, sage=0, ms_amanda=0, tide=0)
        
    if sage == 1:
        run_search(xtandem=0, msgf=0, comet=0, omssa=0, metamorpheus=0, myrimatch=0, sage=1, ms_amanda=0, tide=0)

    end = time.time()
    runtime = end-start
    runtimex = RunTime.objects.get(queue=queue)
    if fasta_type == 'profile' or fasta_type == 'custom':
        runtimex.searchgui_profile = runtime
    elif fasta_type == 'proteome':
        runtimex.searchgui_proteome = runtime
    runtimex.save()
    
    return True
