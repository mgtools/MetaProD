import os
import argparse
import pandas as pd
import statistics
import math
import numpy as np
from decimal import Decimal
import warnings

from django.core.exceptions import ObjectDoesNotExist
from django.db.models import Sum

from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

base = importr('base')
utils = importr('utils')
    
from results.models import (
    Protein, 
    Protein, 
    Peptide, 
    FastaProtein, 
    PsmRatio,
    Psm,
    DiffProtein
)

from projects.models import (
    Queue, 
    SearchSetting, 
    MultiplexLabel, 
    Project, 
    Sample,
    LabelChoice,
    Tag
)

from .run_command import write_debug, settings

# we have protein inferences in the proteininference table

# in the future, support nsaf quantification for non-multiplexed, but this is much easier
def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('project_name', type=str)
    args2 = parser.parse_args(args)
    project = args2.project_name

    analyze_results(project)

def analyze_results(project):
    print("Starting analyze_results for %s" % (project))

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
        
    # we only store the diff proteins for full-proteome
    delete = DiffProtein.objects.filter(project=project).delete()
    
    if searchsetting.multiplex == True:
        print("Updating peptide ratios for multiplexed data (this may take some time).")
    
        psm_ratio_list, columns, labelchoices = load_ratios(project)
    
        peptides_n = generate_initial_peptides(psm_ratio_list, columns)
    
        peptides_initial = peptides_to_phenotypes(peptides_n, columns, labelchoices)
    else:
        print("Updating peptide ratios for label-free data (this may take some time).")
        peptides_initial = generate_lf_peptides(project)
        
    if not os.path.exists(os.path.join(settings.data_folder, project, 'results')):
        os.makedirs(os.path.join(settings.data_folder, project, 'results'))
    peptides_initial.to_csv(os.path.join(settings.data_folder, project, 'results', '%s_peptides_initial.tsv' % project), index=False, sep='\t')
    
    def run_deqms_lf(phenotype):
        print("Calculating differentially expressed proteins for %s." % phenotype)
        
        query = (Tag.objects.filter(project=project)
                            .filter(t_type='Control')
                            .values('name'))
        if len(query) == 0:
            print("There must be at least 1 control tag.")
            return
            
        count_t = (Queue.objects.filter(project=project)
                                .filter(tag__name=phenotype)
                                .count())
        if count_t == 0:
            print("There are no samples with the %s phenotype. Skipping." % phenotype)
            return
            
        ro.r('TMT_columns2_control <- data.frame(matrix(ncol=0, nrow=nrow(df.prot2)))')
        ro.r('count_columns2_control <- data.frame(matrix(ncol=0, nrow=nrow(df.prot2)))')
        for q in query:
            ro.r('TMT_columns2_control = cbind(TMT_columns2_control, df.prot2[, grep("Peak.Area.*.%s$", colnames(df.prot2))])' % q['name'])
            ro.r('count_columns2_control = cbind(count_columns2_control, df.prot2[, grep("psm.*.%s$", colnames(df.prot2), ignore.case=TRUE)])' % q['name'])
        ro.r('TMT_columns2_treatment = df.prot2[, grep("Peak.Area.*.%s$", colnames(df.prot2))]' % phenotype)
        ro.r('dat2 = cbind(TMT_columns2_control, TMT_columns2_treatment)')
        ro.r('rownames(dat2) = df.prot2$accession')
        ro.r('count_columns2_treatment = df.prot2[, grep("psm.*.%s$", colnames(df.prot2), ignore.case=TRUE)]' % phenotype)
        ro.r('count_columns2 = cbind(count_columns2_control, count_columns2_treatment)')
        ro.r('psm.count.table2 = data.frame(count = rowMins(as.matrix(count_columns2)), row.names =  df.prot2$accession)')
        ro.r('control2 = rep("control", each=length(TMT_columns2_control))')
        ro.r('treatment2 = rep("treatment", each=length(TMT_columns2_treatment))')
        ro.r('cond2 = as.factor(c(control2, treatment2))')
        ro.r('design2 <- model.matrix(~0+cond2)')
        ro.r('colnames(design2) = gsub("cond2","",colnames(design2))')
        ro.r('fit12 <- lmFit(dat2, design2)')
        ro.r('x2 <- c("treatment-control")')
        ro.r('contrast2 = makeContrasts(contrasts=x2, levels=design2)')
        ro.r('fit22 <- contrasts.fit(fit12, contrasts = contrast2)')
        ro.r('fit32 <- eBayes(fit22)')
        ro.r('fit32$count = psm.count.table2[rownames(fit32$coefficients),"count"]')
        ro.r('fit42 = spectraCounteBayes(fit32)')
        ro.r('DEqMS.results2 = outputResult(fit42,coef_col = 1)')
        ro.r('fit42$p.value = fit42$sca.p')
        ro.r('prots2 = rownames(DEqMS.results2)')
        ro.r('protein_info2 = protein_list[protein_list$accession %in% prots2, ]')
        ro.r('protein_info2 <- protein_info2[order(protein_info2$accession),]')
        ro.r('DEqMS.results2["accession"] = rownames(DEqMS.results2)')
        ro.r('DEqMS.results2 <- DEqMS.results2[order(DEqMS.results2$accession),]')
        ro.r('protein_info2$accession <- NULL')
        ro.r('DEqMS.results2.final = cbind(protein_info2, DEqMS.results2)')
        ro.r('DEqMS.results2.final <- DEqMS.results2.final[, c("accession", names(DEqMS.results2.final)[names(DEqMS.results2.final) != "accession"])]')
        ro.r('DEqMS.results2.final["gene.1"] <- NULL')
        DEqMS_results2_final_r = ro.r('DEqMS.results2.final')
        # write.table(DEqMS.results2.final, file="x:/DEqMS_results2_final.tsv", sep="\t", row.names=FALSE, quote=FALSE)
        with localconverter(ro.default_converter + pandas2ri.converter):
            DEqMS_results2_final = ro.conversion.rpy2py(DEqMS_results2_final_r)
        DEqMS_results2_final.to_csv(os.path.join(settings.data_folder, project, 'results', '%s_DEqMS_results_final_%s.tsv' % (project, phenotype)), index=False, sep='\t')
        # load results to database  
        # in this case, we're just going to link the files and not store it in a table
        # however, we can store the diff protein results because the columns are the same each time  
        diffprotein_list = []
        project_ = Project.objects.get(name=project)
        for index, row in DEqMS_results2_final.iterrows():
            fp = FastaProtein.objects.get(accession=row['accession'])
            diffprotein = DiffProtein(fp=fp, project=project_, logfc=row['logFC'],
                                    p_value=row['adj.P.Val'], d_p_value=row['sca.adj.pval'])
            diffprotein_list.append(diffprotein)
        DiffProtein.objects.bulk_create(diffprotein_list, 5000)

    def run_deqms_mp(phenotype):
        # df.prot2 = read.table('Z:/data/colon_projects/co1/results/co1_proteins_final.tsv', sep="\t", quote="", header=TRUE)
        print("Calculating differentially expressed proteins for %s." % phenotype)
        
        query = (Tag.objects.filter(project=project)
                            .filter(t_type='Control')
                            .values('name'))
        if len(query) == 0:
            print("There must be at least 1 control tag.")
            return
            
        count_t = (LabelChoice.objects.filter(multiplexlabel__project__name=project)
                                      .filter(tag__name=phenotype)
                                      .count())
        if count_t == 0:
            print("There are no samples with the %s phenotype. Skipping." % phenotype)
            return
            
        ro.r('TMT_columns2_control <- data.frame(matrix(ncol=0, nrow=nrow(df.prot2)))')
        ro.r('count_columns2_control <- data.frame(matrix(ncol=0, nrow=nrow(df.prot2)))')
        for q in query:
            ro.r('TMT_columns2_control = cbind(TMT_columns2_control, df.prot2[, grep(".ratio.%s$", colnames(df.prot2))])' % q['name'])
            ro.r('count_columns2_control = cbind(count_columns2_control, df.prot2[, grep("psm.*.%s$", colnames(df.prot2), ignore.case=TRUE)])' % q['name'])
        ro.r('TMT_columns2_treatment = df.prot2[, grep(".ratio.%s$", colnames(df.prot2))]' % phenotype)
        ro.r('dat2 = cbind(TMT_columns2_control, TMT_columns2_treatment)')
        ro.r('rownames(dat2) = df.prot2$accession')
        ro.r('count_columns2_treatment = df.prot2[, grep("psm.*.%s$", colnames(df.prot2), ignore.case=TRUE)]' % phenotype)
        ro.r('count_columns2 = cbind(count_columns2_control, count_columns2_treatment)')
        ro.r('psm.count.table2 = data.frame(count = rowMins(as.matrix(count_columns2)), row.names =  df.prot2$accession)')        
        #ro.r('psm.count.table2 = data.frame(count = df.prot2$psm, row.names =  df.prot2$accession)')
        # normalize
        ro.r('dat2 = equalMedianNormalization(dat2)')
        # drop na although with pemm, this shouldn't happen
        ro.r('dat2 = na.omit(dat2)')
        ro.r('control2 = rep("control", each=length(TMT_columns2_control))')
        ro.r('treatment2 = rep("treatment", each=length(TMT_columns2_treatment))')
        ro.r('cond2 = as.factor(c(control2, treatment2))')
        ro.r('design2 <- model.matrix(~0+cond2)')
        ro.r('colnames(design2) = gsub("cond2","",colnames(design2))')
        ro.r('x2 = "treatment-control"')
        ro.r('contrast2 = makeContrasts(contrasts=x2, levels=design2)')
        ro.r('fit12 <- lmFit(dat2, design2)')
        ro.r('fit22 <- contrasts.fit(fit12, contrasts = contrast2)')
        ro.r('fit32 <- eBayes(fit22)')
        ro.r('fit32$count = psm.count.table2[rownames(fit32$coefficients),"count"]')
        ro.r('fit42 = spectraCounteBayes(fit32)')
        ro.r('DEqMS.results2 = outputResult(fit42,coef_col = 1)')
        ro.r('fit42$p.value = fit42$sca.p')
        ro.r('prots2 = rownames(DEqMS.results2)')
        ro.r('protein_info2 = protein_list[protein_list$accession %in% prots2, ]')
        ro.r('protein_info2 <- protein_info2[order(protein_info2$accession),]')
        ro.r('DEqMS.results2["accession"] = rownames(DEqMS.results2)')
        ro.r('DEqMS.results2 <- DEqMS.results2[order(DEqMS.results2$accession),]')
        ro.r('protein_info2$accession <- NULL')
        ro.r('DEqMS.results2.final = cbind(protein_info2, DEqMS.results2)')
        # move accession to the front
        ro.r('DEqMS.results2.final <- DEqMS.results2.final[, c("accession", names(DEqMS.results2.final)[names(DEqMS.results2.final) != "accession"])]')
        ro.r('DEqMS.results2.final["gene.1"] <- NULL')
        DEqMS_results2_final_r = ro.r('DEqMS.results2.final')
        # write.table(DEqMS.results2.final, file="x:/DEqMS_results2_final.tsv", sep="\t", row.names=FALSE, quote=FALSE)
        with localconverter(ro.default_converter + pandas2ri.converter):
            DEqMS_results2_final = ro.conversion.rpy2py(DEqMS_results2_final_r)
        DEqMS_results2_final.to_csv(os.path.join(settings.data_folder, project, 'results', '%s_DEqMS_results_final_%s.tsv' % (project, phenotype)), index=False, sep='\t')
        # load results to database
        # in this case, we're just going to link the files and not store it in a table
        # however, we can store the diff protein results because the columns are the same each time  
        diffprotein_list = []
        project_ = Project.objects.get(name=project)
        for index, row in DEqMS_results2_final.iterrows():
            fp = FastaProtein.objects.get(accession=row['accession'])
            diffprotein = DiffProtein(fp=fp, project=project_, logfc=row['logFC'],
                                    p_value=row['adj.P.Val'], d_p_value=row['sca.adj.pval'])
            diffprotein_list.append(diffprotein)
        DiffProtein.objects.bulk_create(diffprotein_list, 5000)        
        
    print("Note: There may be some R warrnings or messages that can be ignored.")
    print("Running PEMM (this may take some time).")
    load_pemm()
    # df.pep = read.table('z:cd2_peptides_initial.tsv', sep="\t", quote="", header=TRUE)
    with localconverter(ro.default_converter + pandas2ri.converter):
        peptides_initial_r = ro.conversion.py2rpy(peptides_initial)
    # assign the dataframe to df.pep
    ro.r.assign('df.pep', peptides_initial_r)

    ### start of R part
    if searchsetting.multiplex == True:
        query = (Tag.objects.filter(project=project)
                            .filter(t_type='Reference')
                            .values('name'))
        # we don't want the reference columns in analysis but store them for later
        ro.r('reference <- data.frame(matrix(ncol=0, nrow=nrow(df.pep)))')
        for q in query:
            ro.r('reference = cbind(reference, df.pep[, grep("%s", colnames(df.pep))])' % q['name'])
            # we only want the log of the Reference ratio columns and not PSM
            ro.r('reference[grep("ratio.%s", colnames(reference))] = log2(reference[grep("ratio.%s", colnames(reference))])')
            ro.r('df.pep = df.pep[, -grep("%s", colnames(df.pep))]' % q['name'])
        # find the ratio columns    
        ro.r('ratio_columns = grep(".ratio.", colnames(df.pep), ignore.case=TRUE)')
        ro.r('psm_columns = grep(".psm.", colnames(df.pep), ignore.case=TRUE)')
        # find psm columns
    else:
        ro.r('ratio_columns = grep("Peak.Area.", colnames(df.pep))')
        ro.r('psm_columns = grep("psm.", colnames(df.pep), ignore.case=TRUE)')
    # select the ratio columns from the data
    ro.r('dat.pep.ratio=df.pep[ratio_columns]')        
    # count non ratio columns
    ro.r('dat.pep.ratio["count"] = (1 - rowSums(is.na(dat.pep.ratio)) / ncol(dat.pep.ratio))')
    # select the rows meating the non-empty criteria (0.5 = 50%)
    threshold = searchsetting.imput_threshold / 100
    ro.r('dat.filtered = dat.pep.ratio[dat.pep.ratio["count"] >= %s,]' % threshold)
    # remove the count column
    ro.r('dat.filtered["count"] <- NULL')
    # label-free is already transformed
    ro.r('dat.filtered.log2 = log2(dat.filtered)')
    ro.r('dat.filtered.log2.normalized = sweep(dat.filtered.log2, 2, colMedians(as.matrix(dat.filtered.log2), na.rm=TRUE))')
    # run PEMM on the log2 transformed data
    ro.r('PEM.result = PEMM_fun(as.matrix(dat.filtered.log2.normalized), phi=0)')
    print('Generating final peptides.')
    # PEM.result$Xhat contains the results
    ro.r('PEM.final = as.data.frame(PEM.result$Xhat)')
    ro.r('rows = row.names(PEM.final)')
    # gather the filtered psm data
    ro.r('psm_data = df.pep[rows, psm_columns]')
    # add 1 to all PSM counts
    # therefore, imputed peptides will count as having 1 PSM
    ro.r('psm_data[is.na(psm_data)] <- 0')
    ro.r('psm_data = psm_data + 1')
    # gather peptide info
    ro.r('data_columns = df.pep[rows, 1:6]')
    ro.r('peptide_final = cbind(data_columns, PEM.final, psm_data)')
    # re-add the reference data in case somebody wants to use it
    if searchsetting.multiplex == True:
        ro.r('reference = reference[rows,]')
        peptide_final = ro.r('peptide_final = cbind(peptide_final, reference)')
    else:
        peptide_final = ro.r('peptide_final')
    with localconverter(ro.default_converter + pandas2ri.converter):
        peptides_final = ro.conversion.rpy2py(peptide_final)
    peptides_final.to_csv(os.path.join(settings.data_folder, project, 'results', '%s_peptides_final.tsv' % project), index=False, sep='\t')

    print('Generating proteins.')
    #ro.r('write.table(peptide_final, file="x:/peptides_final.tsv", sep="\t", row.names=FALSE)')
    ro.r('PEM.final["accession"] = df.pep[rows, "accession"]')
    ro.r('psm_data["accession"] = df.pep[rows, "accession"]')
    # aggregate all the data
    # these should result in things being in the same order so we should be able
    # to merge the two
    ro.r('df.prot2.ratio = aggregate(PEM.final[,1:ncol(PEM.final)-1], by=list(PEM.final$accession), FUN=median)')

    ro.r('df.prot2.psm = aggregate(psm_data[,1:ncol(psm_data)-1], by=list(psm_data$accession), FUN=sum)')
    # now we need to pull the protein information from a previous table
    ro.r('protein_list = df.pep[!duplicated(df.pep$accession),][2:6]')
    ro.r('prots = df.prot2.ratio["Group.1"]')
    ro.r('protein_info = protein_list[protein_list$accession %in% prots$Group.1, ]')
    # make sure they're ordered the same
    ro.r('df.prot2.ratio <- df.prot2.ratio[order(df.prot2.ratio$Group.1),]')
    ro.r('df.prot2.ratio$Group.1 <- NULL')
    ro.r('df.prot2.psm <- df.prot2.psm[order(df.prot2.psm$Group.1),]')
    ro.r('df.prot2.ratio.normalized = sweep(df.prot2.ratio, 2, colMedians(as.matrix(df.prot2.ratio), na.rm=TRUE))')
    ro.r('protein_info <- protein_info[order(protein_info$accession),]')
    # merge it
    ro.r('df.prot2 = cbind(protein_info, df.prot2.ratio.normalized, df.prot2.psm)')
    # drop group.1
    ro.r('df.prot2$Group.1 <- NULL')
    df_prot2 = ro.r('df.prot2')
    # save the results and also bring the dataframe into python
    with localconverter(ro.default_converter + pandas2ri.converter):
        proteins_final = ro.conversion.rpy2py(df_prot2)
    proteins_final.to_csv(os.path.join(settings.data_folder, project, 'results', '%s_proteins_final.tsv' % project), index=False, sep='\t')

    # this is useful if one wants to manually analyze the data
    #if searchsetting.run_deqms == False:
    #    return

    print("Running DEqMS.")
    load_deqms()
    # deqms part
    # select the ratio and treatment columns
    # if there are more phenotypes, one would select those too
    ro.r('library(matrixStats)')

    phenotypes = (Tag.objects.filter(project=project)
                             .filter(t_type='Treatment')
                             .values('name'))
    if len(phenotypes) == 0:
        print("There must be at least 1 treatment tag.")
        return
     
    query = (Tag.objects.filter(project=project)
                            .filter(t_type='Control')
                            .values('name'))
    if len(query) == 0:
        print("There must be at least 1 control tag.")
        return
            
    if searchsetting.multiplex == True:
        for pt in phenotypes:
            run_deqms_mp(pt['name'])
    else:
        for pt in phenotypes:
            run_deqms_lf(pt['name'])

def load_deqms():
    utils.chooseCRANmirror(ind=1)
    packnames = ('BiocManager', 'matrixStats')
    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]

    if len(names_to_install) > 0:
        utils.install_packages(StrVector(names_to_install))
        
    if not rpackages.isinstalled('DEqMS'):
        biocmanager = importr("BiocManager")
        biocmanager.install("DEqMS")
    
    global deqms 
    deqms = importr('DEqMS')

# note that PEMM isn't on CRAN anymore so it has to be installed manually
def load_pemm():
    pemm_dir = os.path.join(settings.install_folder, "software", "PEMM_1.0.tar.gz")
    utils.install_packages(pemm_dir, repos = ro.r("NULL"), type = "source")
    global pemm
    pemm = importr('PEMM')
    
    utils.chooseCRANmirror(ind=1)
    packnames = ('matrixStats')
    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
    if len(names_to_install) > 0:
        utils.install_packages(StrVector(names_to_install))
        
    global matrixstats
    matrixstats = importr('matrixStats')

def load_ratios(project):
    print("Loading ratios for multiplexed data.")
    query = (PsmRatio.objects.filter(psm__queue__project__name=project)
                             .filter(psm__type="proteome")
                             .exclude(psm__peptide__protein__fp__ppid='0')
                             .exclude(psm__queue__skip=True)
                             .exclude(psm__queue__error__gte=(1 + settings.max_retries))
                             .values('psm_id', 'psm__peptide__id', 'ratio', 'label',
                                     'psm__mod_sequence',                             
                                     'psm__peptide__protein__fp__accession',
                                     'psm__peptide__protein__fp__gene',
                                     'psm__peptide__protein__fp__description',                                     
                                     'psm__peptide__protein__fp__ppid',
                                     'psm__peptide__protein__fp__ppid__organism',
                                     'psm__queue__sample__name'))

    # convert to dataframe
    psm_ratio_list = pd.DataFrame(list(query))

    del query
    
    psm_ratio_list = psm_ratio_list.rename(columns={'psm_id': 'psm_id',
                                                    'psm__peptide__id': 'peptide_id',
                                                    'psm__mod_sequence': 'sequence',                                                    
                                                    'psm__peptide__protein__fp__accession': 'accession',
                                                    'psm__peptide__protein__fp__gene': 'gene',                                                    
                                                    'psm__peptide__protein__fp__description': 'description',
                                                    'psm__peptide__protein__fp__ppid': 'ppid',                                                    
                                                    'psm__peptide__protein__fp__ppid__organism': 'organism',
                                                    'psm__queue__sample__name': 'sample'})
                    
    # turn ratio/label columns into separate columns per label
    psm_ratio_list = psm_ratio_list.pivot(index=['psm_id', 
                                                 'peptide_id',
                                                 'sequence',                                                 
                                                 'accession',
                                                 'gene',
                                                 'description',                                                 
                                                 'ppid',
                                                 'organism',
                                                 'sample'], 
                                            columns='label', values='ratio')
    psm_ratio_list = psm_ratio_list.reset_index()

    query = (LabelChoice.objects.filter(
                multiplexlabel__project__name=project,
            ).values('multiplexlabel__sample__name', 'label__name', 'identifier', 'tag__t_type', 'tag__name'))
    labelchoices = pd.DataFrame(list(query))
    del query
    columns = [x for x in psm_ratio_list.columns[9:]]
    # drop anything PSMs missing all ratios
    psm_ratio_list.dropna(subset=columns, inplace=True, how='all')
    # also drop any PSMs missing the reference because it can't be normalized
    # have to lookup the reference for a given sample
    references = list(labelchoices[labelchoices['tag__t_type'] == 'Reference'][['multiplexlabel__sample__name', 'label__name']].index)
    psm_ratio_list.drop(index=references, inplace=True)

    for col in columns:
        psm_ratio_list[col] = psm_ratio_list[col].astype(float)
        
    return(psm_ratio_list, columns, labelchoices)
    
def generate_initial_peptides(psm_ratio_list, columns):
    print("Generating initial peptides.")
    # find the median of the psms for the raw peptide ratio
    peptides = psm_ratio_list.groupby(['peptide_id',
                                       'sample',
                                       'sequence',
                                       'accession',
                                       'gene',
                                       'description',                                       
                                       'ppid',
                                       'organism',
                                      ])
    # median of PSM is peptide
    peptides_r = peptides[columns].median()
    # peptide PSM count is total PSMs
    peptides_r['psm'] = peptides.size()
    peptides_r = peptides_r.reset_index()
    # normalize peptide ratio by median of the channel (sets the channel median to 1)
    peptides_n = peptides_r.groupby('sample')[columns].transform(lambda x: x/x.median())

    # pull the peptide id and sample from other dataframe
    peptides_n = pd.merge(peptides_r[['peptide_id', 
                                      'sample',
                                      'sequence',
                                      'accession',
                                      'gene',
                                      'description',
                                      'ppid',
                                      'organism',
                                      'psm',
                                    ]], peptides_n, left_index=True, right_index=True)
                                    
    return(peptides_n)
    
def peptides_to_phenotypes(peptides_n, columns, labelchoices):
    print("Converting peptides into phenotypes.")
    print("1")
    sample_list = list(peptides_n['sample'].unique())
    # map peptides to phenotypes
    series_dict = {}
    series_dict_p = {}
    for sample in sample_list:
        for col in columns:
            mapping = labelchoices[(labelchoices['multiplexlabel__sample__name'] == sample)
                                    & (labelchoices['label__name'] == col)]

            phenotype_t = "%s ratio %s" % (mapping.iloc[0]['identifier'], mapping.iloc[0]['tag__name'])
            phenotype_p = "%s psm %s" % (mapping.iloc[0]['identifier'], mapping.iloc[0]['tag__name'])
            new_series = peptides_n[(peptides_n['sample'] == sample)][col]
            new_series.name = phenotype_t
            new_series_p = peptides_n[(peptides_n['sample'] == sample)]['psm']
            new_series.name = phenotype_p
            if phenotype_t not in series_dict:
                series_dict[phenotype_t] = new_series
            else:
                series_dict[phenotype_t] = pd.concat(objs=[series_dict[phenotype_t], new_series])
            if phenotype_p not in series_dict_p:
                series_dict_p[phenotype_p] = new_series_p
            else:
                series_dict_p[phenotype_p] = pd.concat(objs=[series_dict_p[phenotype_p], new_series_p])
    print("2")
    # in some cases, the phenotype will be the same so merge the numbers
    for s in series_dict:
        series_dict[s] = series_dict[s].groupby(series_dict[s].index).median()

    for s in series_dict_p:
        series_dict_p[s] = series_dict_p[s].groupby(series_dict_p[s].index).median()
    
    series_df = pd.DataFrame(series_dict)
    del series_dict
    series_df_p = pd.DataFrame(series_dict_p)
    del series_dict_p
    print("3")    
    peptide_samples = pd.merge(peptides_n, series_df, left_index=True, right_index=True)
    del series_df
    peptide_samples = pd.merge(peptide_samples, series_df_p, left_index=True, right_index=True)
    del series_df_p
    peptide_samples = peptide_samples.reset_index()
    print("4")
    # 'ColonRef Not Reported ratio Reference'
    peptide_samples = peptide_samples.drop(columns=columns)
    data_cols = list(sorted(peptide_samples.columns[10:]))
    peptide_samples = peptide_samples.groupby(['sequence',
                                               'accession',
                                               'gene',
                                               'description',
                                               'ppid',
                                               'organism']).agg({'psm': 'sum', **{e: 'median' for e in data_cols}})
    peptide_samples = peptide_samples.reset_index()
 
    return(peptide_samples)
    
def generate_lf_peptides(project):
    query = (Peptide.objects.filter(queue__project__name=project)
                            .filter(type="proteome")
                            .filter(peak_area_psm__gt=0)
                            .exclude(protein__fp__ppid='0')
                            .exclude(queue__skip=True)
                            .exclude(queue__error__gte=(1 + 2))
                            .values('mod_sequence',
                                    'protein__fp__accession',
                                    'protein__fp__gene',
                                    'protein__fp__description',
                                    'protein__fp__ppid',
                                    'protein__fp__ppid__organism',
                                    'peak_area_psm',
                                    'queue__sample__name',
                                    'queue__filename',
                                    'queue__tag__name',
                                    'peak_area'))
            
    peptide_list = pd.DataFrame(list(query))

    peptide_list = peptide_list.rename(columns={'mod_sequence': 'sequence',    
                                                'protein__fp__accession': 'accession',
                                                'protein__fp__gene': 'gene',
                                                'protein__fp__description': 'description',
                                                'protein__fp__ppid': 'ppid',
                                                'protein__fp__ppid__organism': 'organism',
                                                'peak_area_psm': 'psm',
                                                'queue__sample__name': 'sample',
                                                'queue__filename': 'filename',
                                                'queue__tag__name': 'tag',
                                                })
    
    rows_list = []
    for index, row in peptide_list.iterrows():
        dict1 = {}
        # maybe try to do this as a dict/list to build a new dataframe then concat them
        if row['sample'] is not None:
            new_column_1 = "Peak Area %s %s" % (str(row['sample']), str(row['tag']))
            new_column_2 = "PSM %s %s" % (str(row['sample']), str(row['tag']))
        else:
            new_column_1 = "Peak Area %s %s" % (str(row['filename']), str(row['tag']))
            new_column_2 = "PSM %s %s" % (str(row['filename']), str(row['tag']))
            
        dict1.update({new_column_1:float(row['peak_area']), new_column_2:float(row['psm'])})
        rows_list.append(dict1)

    peak_area_df = pd.DataFrame(rows_list)
    peptide_list = pd.concat([peptide_list, peak_area_df], axis='columns')
    peptide_list = peptide_list.drop(columns=['sample', 'tag', 'peak_area', 'filename', 'psm'])
         
    data_cols = list(sorted(peptide_list.columns[6:]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)  
        peptide_list = peptide_list.groupby(['sequence',
                                            'accession',
                                            'gene',
                                            'description',
                                            'ppid',
                                            'organism']).agg({**{e: 'median' for e in data_cols}})
    peptide_list = peptide_list.reset_index()

    return(peptide_list)