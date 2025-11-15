from django.shortcuts import render
from django.http import HttpResponse
from django.views import generic
from django.http import Http404
from django.shortcuts import get_list_or_404
from django.db.models import Sum, Count
from django_tables2.export.views import ExportMixin
from django_tables2.export.export import TableExport

from django_tables2.views import SingleTableMixin, MultiTableMixin
from django_tables2 import SingleTableView
from django_filters.views import FilterView

from projects.models import Project, Queue, RunTime, SearchSetting
from results.models import (
    Protein, 
    Peptide, 
    Psm, 
    SpeciesSummary, 
    SpeciesFileSummary,
    DiffProtein,
)
from website.tables import (
    ProjectListTable, 
    ProteinListTable, 
    PeptideListTable, 
    PsmListTable, 
    FileListTable, 
    SpeciesListTable,
    FileSummaryTable,
    DiffProteinListTable
)
from website.filters import (
    ProteinListFilter, 
    PeptideListFilter, 
    PsmListFilter, 
    FileListFilter, 
    SpeciesListFilter,
    FileSummaryFilter,
    DiffProteinListFilter,
)
from website.forms import (
    ProteinListFilterFormHelper, 
    PeptideListFilterFormHelper, 
    PsmListFilterFormHelper, 
    FileListFilterFormHelper, 
    SpeciesListFilterFormHelper,
    FileSummaryFilterFormHelper,
    DiffProteinListFilterFormHelper,
)

# Create your views here.

class FilteredSingleTableView(SingleTableMixin, FilterView):
    formhelper_class = None

    def get_filterset(self, filterset_class):
        kwargs = self.get_filterset_kwargs(filterset_class)
        filterset = filterset_class(**kwargs)
        filterset.form.helper = self.formhelper_class()
        return filterset

class FilteredMultiTableView(MultiTableMixin, FilterView):
    formhelper_class = None

    def get_filterset(self, filterset_class):
        kwargs = self.get_filterset_kwargs(filterset_class)
        filterset = filterset_class(**kwargs)
        filterset.form.helper = self.formhelper_class()
        return filterset
        
def index(request):
    return render(request, 'website/index.html')

def help(request):
    return render(request, 'website/help.html')
    
class ProjectListView(SingleTableView):
    model = Project
    template_name = 'website/project_list.html'
    table_class = ProjectListTable
    paginate_by = 100

def project_list(request):
    projects = (Project.objects.all().values('name', 'description'))
    project_dict = {}
    for project in projects:
        project_dict[project['name']] = {}
        project_dict[project['name']]['description'] = project['description']
        try:
            multiplexed = (SearchSetting.objects.get(project__name=project['name']))
            project_dict[project['name']]['multiplexed'] = multiplexed.multiplex
        except:
            project_dict[project['name']]['multiplexed'] = False
    return render(request, 
        'website/project_list.html', 
        {'project_dict': project_dict}
    )
    
class ProteinListView(ExportMixin, FilteredSingleTableView):
    model = Protein
    template_name = 'website/protein_list.html'
    table_class = ProteinListTable
    paginate_by = 250
    filterset_class = ProteinListFilter
    formhelper_class = ProteinListFilterFormHelper
    export_formats = ['csv']

class DiffProteinListView(ExportMixin, FilteredSingleTableView):
    model = DiffProtein
    template_name = 'website/diffprotein_list.html'
    table_class = DiffProteinListTable
    paginate_by = 250
    filterset_class = DiffProteinListFilter
    formhelper_class = DiffProteinListFilterFormHelper
    export_formats = ['csv']
    
class PeptideListView(ExportMixin, FilteredSingleTableView):
    model = Peptide
    template_name = 'website/peptide_list.html'
    table_class = PeptideListTable
    paginate_by = 250
    filterset_class = PeptideListFilter
    formhelper_class = PeptideListFilterFormHelper
    export_formats = ['csv']
    
class PsmListView(ExportMixin, FilteredSingleTableView):
    model = Psm
    template_name = 'website/psm_list.html'
    table_class = PsmListTable
    paginate_by = 250
    filterset_class = PsmListFilter
    formhelper_class = PsmListFilterFormHelper
    export_formats = ['csv']

class FileListView(FilteredSingleTableView):
    model = Queue
    template_name = 'website/file_list.html'
    table_class = FileListTable
    paginate_by = 100
    filterset_class = FileListFilter
    formhelper_class = FileListFilterFormHelper

class SpeciesListView(ExportMixin, FilteredSingleTableView):
    model = SpeciesSummary
    template_name = 'website/species_list.html'
    table_class = SpeciesListTable
    paginate_by = 100
    filterset_class = SpeciesListFilter
    formhelper_class = SpeciesListFilterFormHelper

class FileSummaryView(ExportMixin, FilteredSingleTableView):
    model = SpeciesFileSummary
    template_name = 'website/file_summary.html'
    table_class = FileSummaryTable
    paginate_by = 100
    filterset_class = FileSummaryFilter
    formhelper_class = FileSummaryFilterFormHelper

def summary(request, project):
    files = (Queue.objects.filter(project__name=project,
                                  skip=False)
                          .values('id', 'filename'))
    searchsetting=SearchSetting.objects.get(project=project)
    if searchsetting.custom_fasta == True:
        fasta_type = "custom"
    else:
        fasta_type = "not_custom"
    custom_dict = {}
    profile_dict = {}
    proteome_dict = {}
    trfp = 0
    sgui = 0
    peps = 0
    repo = 0
    rere = 0
    prre = 0
    mzre = 0
    total_psm_cust = 0
    total_psm_prof = 0
    total_psm_prot = 0
    total_pep_cust = 0
    total_pep_prof = 0
    total_pep_prot = 0
    total_pro_cust = 0
    total_pro_prof = 0
    total_pro_prot = 0
    for file in files:
        custom_dict[file['id']] = {}
        custom_dict[file['id']]['project'] = project
        custom_dict[file['id']]['filename'] = file['filename']        
        proteins_cust = (
            Protein.objects.filter(
                queue__project__name=project, 
                fasta_type='custom', 
                queue__id=file['id']
            ).count()
        )
        custom_dict[file['id']]['protein'] = proteins_cust
        peptides_cust = (
            Peptide.objects.filter(
                queue__project__name=project, 
                fasta_type='custom', 
                queue__id=file['id']
            ).count()
        )
        custom_dict[file['id']]['peptide'] = peptides_cust
        psms_cust = (
            Psm.objects.filter(
                queue__project__name=project, 
                fasta_type='custom', 
                queue__id=file['id']
            ).count()
        )    
        custom_dict[file['id']]['psm'] = psms_cust

        profile_dict[file['id']] = {}
        profile_dict[file['id']]['project'] = project
        profile_dict[file['id']]['filename'] = file['filename']    
        proteins_prof = (
            Protein.objects.filter(
                queue__project__name=project, 
                fasta_type='profile', 
                queue__id=file['id']
            ).count()
        )
        profile_dict[file['id']]['protein'] = proteins_prof
        peptides_prof = (
            Peptide.objects.filter(
                queue__project__name=project, 
                fasta_type='profile', 
                queue__id=file['id']
            ).count()
        )
        profile_dict[file['id']]['peptide'] = peptides_prof
        psms_prof = (
            Psm.objects.filter(
                queue__project__name=project, 
                fasta_type='profile', 
                queue__id=file['id']
            ).count()
        )
        profile_dict[file['id']]['psm'] = psms_prof
        
        proteome_dict[file['id']] = {}
        proteome_dict[file['id']]['project'] = project
        proteome_dict[file['id']]['filename'] = file['filename']
        proteins_prot = (
            Protein.objects.filter(
                queue__project__name=project, 
                fasta_type='proteome', 
                queue__id=file['id']
            ).count()
        )
        proteome_dict[file['id']]['protein'] = proteins_prot
        peptides_prot = (
            Peptide.objects.filter(
                queue__project__name=project, 
                fasta_type='proteome', 
                queue__id=file['id']
            ).count()
        )
        proteome_dict[file['id']]['peptide'] = peptides_prot
        psms_prot = (
            Psm.objects.filter( 
                queue__project__name=project, 
                fasta_type='proteome', 
                queue__id=file['id']
            ).count()
        )
        proteome_dict[file['id']]['psm'] = psms_prot
        
        try:
            runtime = RunTime.objects.get(queue__id=file['id'])
            trfp += runtime.msconvert
            sgui += runtime.searchgui_profile+runtime.searchgui_proteome
            peps += runtime.peptideshaker_profile+runtime.peptideshaker_proteome
            repo += runtime.reporter_profile+runtime.reporter_proteome
            rere += runtime.read_results_profile+runtime.read_results_proteome
            prre += runtime.process_results_profile+runtime.process_results_proteome
            mzre += runtime.mzmine_profile+runtime.mzmine_proteome
        except:
            pass
        
        total_psm_cust += psms_cust
        total_psm_prof += psms_prof
        total_psm_prot += psms_prot
        
        total_pep_cust += peptides_cust
        total_pep_prof += peptides_prof
        total_pep_prot += peptides_prot
        
        total_pro_cust += proteins_cust
        total_pro_prof += proteins_prof
        total_pro_prot += proteins_prot
        
    runtime_dict = {'trfp':trfp, 'sgui':sgui, 'peps':peps, 'repo':repo,
                    'rere':rere, 'prre':prre, 'mzre':mzre}
    totals_dict = {'pep_cust': total_pep_cust, 'pep_prof': total_pep_prof, 'pep_prot': total_pep_prot,
                   'psm_cust': total_psm_cust, 'psm_prof': total_psm_prof, 'psm_prot': total_psm_prot,
                   'pro_cust': total_pro_cust, 'pro_prof': total_pro_prof, 'pro_prot': total_pro_prot,
                  }
    
    # omit this for custom because it's not necessarily possible
    top_x_ppid_psm_profile = (
        Protein.objects.filter(fasta_type='profile')
                       .filter(queue__project__name=project)
                       .values('fp__ppid__proteome', 'fp__ppid__organism')
                       .annotate(total=Sum('val_num_psm'))
                       .order_by('-total')[:10]
    )
    top_x_ppid_psm_proteome = (
        Protein.objects.filter(fasta_type='proteome')
                       .filter(queue__project__name=project)
                       .values('fp__ppid__proteome', 'fp__ppid__organism')
                       .annotate(total=Sum('val_num_psm'))
                       .order_by('-total')[:10]
    )                                             
    top_x_ppid_nsaf_profile = (
        Protein.objects.filter(fasta_type='profile')
                       .filter(queue__project__name=project)
                       .values('fp__ppid__proteome', 'fp__ppid__organism')
                       .annotate(total=Sum('nsaf'))
                       .order_by('-total')[:10]
    )
    top_x_ppid_nsaf_proteome = (
        Protein.objects.filter(fasta_type='proteome')
                       .filter(queue__project__name=project)
                       .values('fp__ppid__proteome', 'fp__ppid__organism')
                       .annotate(total=Sum('nsaf'))
                       .order_by('-total')[:10]
    )
                                             
    return render(request, 
        'website/summary.html',
        {'custom_dict': custom_dict,
         'profile_dict': profile_dict,
         'proteome_dict': proteome_dict,
         'top_x_ppid_psm_profile': top_x_ppid_psm_profile,
         'top_x_ppid_psm_proteome': top_x_ppid_psm_proteome,
         'top_x_ppid_nsaf_profile': top_x_ppid_nsaf_profile,
         'top_x_ppid_nsaf_proteome': top_x_ppid_nsaf_proteome,
         'runtimes': runtime_dict,
         'totals': totals_dict
        }
    )
    
def graphs(request, project):
    # we can load the species summary and use that and generate species graphs by nsaf

    species_dict = {}
    ss = SpeciesSummary.objects.filter(
        project=project,
        fasta_type='proteome'
    ).values(
        'ppid__organism',
        'nsaf'
    )
    for entry in ss:
       species_dict[entry['ppid__organism']] = entry['nsaf']
       
    return render(request, 
        'website/graphs.html',
        {'species_dict': species_dict,
        }
    )    