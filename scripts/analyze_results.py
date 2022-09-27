import os
import argparse
import pandas as pd
import statistics
import math
import numpy as np
from decimal import Decimal

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
    LabelChoice
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
        
    if searchsetting.multiplex == False:
        print("Non-multiplex not supported yet.")
        return

    # we only store the diff proteins for full-proteome
    delete = DiffProtein.objects.filter(project=project).delete()
    
    print("Updating peptide ratios (this may take some time).")

    psm_ratio_list, columns, labelchoices = load_ratios(project)

    peptides_n = generate_initial_peptides(psm_ratio_list, columns)

    peptide_samples = peptides_to_phenotypes(peptides_n, columns, labelchoices)

    if not os.path.exists(os.path.join(settings.data_folder, project, 'results')):
        os.makedirs(os.path.join(settings.data_folder, project, 'results'))
    peptide_samples.to_csv(os.path.join(settings.data_folder, project, 'results', 'peptide_samples.tsv'), index=False, sep='\t')
    
    if searchsetting.run_deqms == False:
        return

    print("Note: There may be some R warrnings or messages that can be ignored.")
    print("Running PEMM (this may take some time).")
    load_pemm()
    with localconverter(ro.default_converter + pandas2ri.converter):
        peptide_samples_r = ro.conversion.py2rpy(peptide_samples)
    # assign the dataframe to df.prot2
    ro.r.assign('df.pep', peptide_samples_r)

    ### start of R part
    # load peptide ratios
    #df.pep = read.table('x:/co_3_c4_10ppm/peptide_samples.csv', sep='\t', quote="", header=TRUE)
    # df.pep = read.table('z:/data/arsenic_c4_10ppm/results/peptide_samples_modified.tsv', sep='\t', quote="", header=TRUE)
    # we don't want the reference columns in analysis
    ro.r('reference = df.pep[, grep("Reference", colnames(df.pep))]')
    ro.r('reference[colnames(reference)[2]] = log2(reference[colnames(reference)[2]])')
    ro.r('df.pep = df.pep[, -grep("Reference", colnames(df.pep))]')
    # find the ratio columns    
    ro.r('ratio_columns = grep(".ratio.", colnames(df.pep))')
    # find psm columns
    ro.r('psm_columns = grep(".psm.", colnames(df.pep))')
    # select the ratio columns from the data
    ro.r('dat.pep.ratio=df.pep[ratio_columns]')
    # count non ratio columns
    ro.r('dat.pep.ratio["count"] = (1 - rowSums(is.na(dat.pep.ratio)) / ncol(dat.pep.ratio))')
    # select the rows meating the non-empty criteria (0.5 = 50%)
    threshold = searchsetting.imput_threshold / 100
    ro.r('dat.filtered = dat.pep.ratio[dat.pep.ratio["count"] >= %s,]' % threshold)
    # remove the count column
    ro.r('dat.filtered["count"] <- NULL')
    # log2 transform
    ro.r('dat.filtered.log2 = log2(dat.filtered)')
    # run PEMM on the log2 transformed data
    ro.r('PEM.result = PEMM_fun(as.matrix(dat.filtered.log2), phi=0)')
    # PEM.result$Xhat contains the results
    ro.r('PEM.final = as.data.frame(PEM.result$Xhat)')
    ro.r('rows = row.names(PEM.final)')
    # gather the filtered psm data
    ro.r('psm_data = df.pep[rows, psm_columns]')
    # replace NA values with 1  
    # therefore, imputed peptides will count as having 1 PSM, but we don't actually
    # know what it will be. could also try median or mean value for that peptide across all data
    ro.r('psm_data[is.na(psm_data)] <- 1')
    # gather peptide info
    ro.r('data_columns = df.pep[rows, 1:5]')
    ro.r('peptide_final = cbind(data_columns, PEM.final, psm_data)')
    # re-add the reference data in case somebody wants to use it
    ro.r('reference = reference[rows,]')
    peptide_final = ro.r('peptide_final = cbind(peptide_final, reference)')
    with localconverter(ro.default_converter + pandas2ri.converter):
        peptides_final = ro.conversion.rpy2py(peptide_final)
    peptides_final.to_csv(os.path.join(settings.data_folder, project, 'results', 'peptides_final.tsv'), index=False, sep='\t')
    #ro.r('write.table(peptide_final, file="x:/peptides_final.tsv", sep="\t", row.names=FALSE)')
    ro.r('PEM.final["accession"] = df.pep[rows, "accession"]')
    ro.r('psm_data["accession"] = df.pep[rows, "accession"]')
    # aggregate all the data
    # these should result in things being in the same order so we should be able
    # to merge the two
    ro.r('df.prot2.ratio = aggregate(PEM.final[,1:ncol(PEM.final)-1], by=list(PEM.final$accession), FUN=median)')
    ro.r('df.prot2.psm = aggregate(psm_data[,1:ncol(psm_data)-1], by=list(psm_data$accession), FUN=sum)')
    # now we need to pull the protein information from a previous table
    ro.r('protein_list = df.pep[!duplicated(df.pep$accession),][2:5]')
    ro.r('prots = df.prot2.ratio["Group.1"]')
    ro.r('protein_info = protein_list[protein_list$accession %in% prots$Group.1, ]')
    # make sure they're ordered the same
    ro.r('df.prot2.ratio <- df.prot2.ratio[order(df.prot2.ratio$Group.1),]')
    ro.r('df.prot2.psm <- df.prot2.psm[order(df.prot2.psm$Group.1),]')
    ro.r('protein_info <- protein_info[order(protein_info$accession),]')
    # merge it
    ro.r('df.prot2 = cbind(protein_info, df.prot2.ratio, df.prot2.psm)')
    # drop group.1
    ro.r('df.prot2$Group.1 <- NULL')
    ro.r('df.prot2$Group.1 <- NULL')
    df_prot2 = ro.r('df.prot2')
    # save the results and also bring the dataframe into python
    with localconverter(ro.default_converter + pandas2ri.converter):
        proteins_final = ro.conversion.rpy2py(df_prot2)
    proteins_final.to_csv(os.path.join(settings.data_folder, project, 'results', 'proteins_final.tsv'), index=False, sep='\t')

    print("Running DEqMS.")
    load_deqms()        
    # deqms part
    # select the ratio and treatment columns
    # if there are more phenotypes, one would select those too
    ro.r('library(matrixStats)')
    ro.r('TMT_columns2_treatment = grep("ratio.Treatment$", colnames(df.prot2))')
    ro.r('TMT_columns2_control = grep("ratio.Control$", colnames(df.prot2))')
    ro.r('TMT_columns2 = c(TMT_columns2_control, TMT_columns2_treatment)')
    ro.r('dat2=df.prot2[TMT_columns2]')
    ro.r('rownames(dat2) = df.prot2$accession')
    # normalize
    ro.r('dat2 = equalMedianNormalization(dat2)')
    # drop na although with pemm, this shouldn't happen
    ro.r('dat2 = na.omit(dat2)')
    # for more phenotypes, repeat this for each name above
    ro.r('control2 = rep("control", each=length(TMT_columns2_control))')
    ro.r('treatment2 = rep("treatment", each=length(TMT_columns2_treatment))')
    # add factors for more phenotypes
    ro.r('cond2 = as.factor(c(control2, treatment2))')
    ro.r('design2 <- model.matrix(~0+cond2)')
    ro.r('colnames(design2) = gsub("cond2","",colnames(design2))')
    # for additional phenotypes, add them here, eg. x2 <- c("p1-ctrl", "p2-ctrl")
    ro.r('x2 = "treatment-control"')
    ro.r('contrast2 = makeContrasts(contrasts=x2, levels=design2)')
    ro.r('fit12 <- lmFit(dat2, design2)')
    ro.r('fit22 <- contrasts.fit(fit12, contrasts = contrast2)')
    ro.r('fit32 <- eBayes(fit22)')
    # add psm counts for more phenotypes
    ro.r('count_columns2_control = grep("psm.Control$", colnames(df.prot2))')
    ro.r('count_columns2_treatment = grep("psm.Treatment$", colnames(df.prot2))')
    ro.r('count_columns2 = c(count_columns2_control, count_columns2_treatment)')
    ro.r('psm.count.table2 = data.frame(count = rowMins(as.matrix(df.prot2[,count_columns2])), row.names =  df.prot2$accession)')
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
    DEqMS_results2_final.to_csv(os.path.join(settings.data_folder, project, 'results', 'DEqMS_results_final.tsv'), index=False, sep='\t')

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

def load_ratios(project):
    query = (PsmRatio.objects.filter(psm__queue__project__name=project)
                             .filter(psm__type="proteome")
                             .exclude(psm__peptide__protein__fp__ppid='0')
                             .values('psm_id', 'psm__peptide__id', 'ratio', 'label',
                                     'psm__peptide__protein__fp__accession',
                                     'psm__peptide__protein__fp__ppid',
                                     'psm__peptide__protein__fp__description',
                                     'psm__peptide__protein__fp__gene',
                                     'psm__mod_sequence',
                                     'psm__queue__sample__name'))

    # convert to dataframe
    psm_ratio_list = pd.DataFrame(list(query))

    psm_ratio_list = psm_ratio_list.rename(columns={'psm_id': 'psm_id',
                                                    'psm__peptide__id': 'peptide_id',
                                                    'psm__peptide__protein__fp__accession': 'accession',
                                                    'psm__peptide__protein__fp__ppid': 'ppid',
                                                    'psm__peptide__protein__fp__gene': 'gene',
                                                    'psm__mod_sequence': 'sequence',
                                                    'psm__peptide__protein__fp__description': 'description',
                                                    'psm__queue__sample__name': 'sample'})
                    
    # turn ratio/label columns into separate columns per label
    psm_ratio_list = psm_ratio_list.pivot(index=['psm_id', 
                                                 'peptide_id',
                                                 'accession',
                                                 'ppid',
                                                 'gene',
                                                 'sequence',
                                                 'description',
                                                 'sample'], 
                                            columns='label', values='ratio')
    psm_ratio_list = psm_ratio_list.reset_index()

    query = (LabelChoice.objects.filter(
                multiplexlabel__project__name=project,
            ).values('multiplexlabel__sample__name', 'label__name', 'identifier', 'phenotype'))
    labelchoices = pd.DataFrame(list(query))
    columns = [x for x in psm_ratio_list.columns[8:]]
    # drop anything PSMs missing all ratios
    psm_ratio_list.dropna(subset=columns, inplace=True, how='all')
    # also drop any PSMs missing the reference because it can't be normalized
    # have to lookup the reference for a given sample
    references = list(labelchoices[labelchoices['phenotype'] == 'Reference'][['multiplexlabel__sample__name', 'label__name']].index)
    psm_ratio_list.drop(index=references, inplace=True)

    for col in columns:
        psm_ratio_list[col] = psm_ratio_list[col].astype(float)
        
    return(psm_ratio_list, columns, labelchoices)
    
def generate_initial_peptides(psm_ratio_list, columns):
    # find the median of the psms for the raw peptide ratio
    peptides = psm_ratio_list.groupby(['peptide_id',
                                       'sample',
                                       'sequence',
                                       'accession',
                                       'gene',
                                       'ppid',
                                       'description',
                                      ])
    peptides_r = peptides[columns].median()
    peptides_r['psm'] = peptides.size()
    peptides_r = peptides_r.reset_index()
    peptides_n = peptides_r.groupby('sample')[columns].transform(lambda x: x/x.median())

    # pull the peptide id and sample from other dataframe
    peptides_n = pd.merge(peptides_r[['peptide_id', 
                                      'sample',
                                      'sequence',
                                      'accession',
                                      'gene',
                                      'description',
                                      'ppid',
                                      'psm',
                                    ]], peptides_n, left_index=True, right_index=True)
                                    
    return(peptides_n)
    
def peptides_to_phenotypes(peptides_n, columns, labelchoices):
    sample_list = list(peptides_n['sample'].unique())
    # map peptides to phenotypes
    series_dict = {}
    series_dict_p = {}
    for sample in sample_list:
        for col in columns:
            mapping = labelchoices[(labelchoices['multiplexlabel__sample__name'] == sample)
                                    & (labelchoices['label__name'] == col)]

            phenotype_t = "%s ratio %s" % (mapping.iloc[0]['identifier'], mapping.iloc[0]['phenotype'])
            phenotype_p = "%s psm %s" % (mapping.iloc[0]['identifier'], mapping.iloc[0]['phenotype'])
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

    # in some cases, the phenotype will be the same so merge the numbers
    for s in series_dict:
        series_dict[s] = series_dict[s].groupby(series_dict[s].index).median()

    for s in series_dict_p:
        series_dict_p[s] = series_dict_p[s].groupby(series_dict_p[s].index).median()
    
    series_df = pd.DataFrame(series_dict)
    series_df_p = pd.DataFrame(series_dict_p)
    
    peptide_samples = pd.merge(peptides_n, series_df, left_index=True, right_index=True)
    peptide_samples = pd.merge(peptide_samples, series_df_p, left_index=True, right_index=True)
    peptide_samples = peptide_samples.reset_index()

    # 'ColonRef Not Reported ratio Reference'
    peptide_samples = peptide_samples.drop(columns=columns)
    data_cols = list(sorted(peptide_samples.columns[9:]))
    peptide_samples = peptide_samples.groupby(['sequence',
                                               'accession',
                                               'gene',
                                               'description',
                                               'ppid']).agg({'psm': 'sum', **{e: 'median' for e in data_cols}})
    peptide_samples = peptide_samples.reset_index()
 
    return(peptide_samples)