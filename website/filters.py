import django_filters

from projects.models import Project, Queue, Tag, SearchSetting

from results.models import (
    Protein, 
    Peptide, 
    Psm, 
    SpeciesSummary,
    SpeciesFileSummary,
    DiffProtein,
)

STEP_CHOICES = (
    ('profile', 'Profile'),
    ('proteome', 'Proteome'),
)
    
class SummaryFilter(django_filters.FilterSet):
    project = django_filters.CharFilter(field_name='project__name')
    class Meta:
        model = Project
        fields = ['project',]

def get_files_protein(request):
    if request is None:
        return Queue.objects.none()
    
    project = request.GET.get('project','')
    
    return Queue.objects.filter(project__name=project)

def get_tags_protein(request):
    if request is None:
        return Tag.objects.none()
        
    project = request.GET.get('project','')
    
    return Tag.objects.filter(project__name=project)


class ProteinListFilter(django_filters.FilterSet):
    project = django_filters.CharFilter(field_name='queue__project__name')
    gene = django_filters.CharFilter(field_name='fp__gene')
    proteome = django_filters.CharFilter(field_name='fp__ppid__proteome')
    organism = django_filters.CharFilter(field_name='fp__ppid__organism')
    file = django_filters.ModelChoiceFilter(field_name='queue__filename', 
                                            queryset=get_files_protein,
                                            label='Filename')
    fasta_type = django_filters.AllValuesFilter(lookup_expr='iexact')
    tag = django_filters.ModelChoiceFilter(field_name='queue__tag__name', 
                                           queryset=get_tags_protein,
                                           label='Tag')
    class Meta:
        model = Protein
        fields = ['project', 'gene', 'proteome', 'organism', 'file', 'fasta_type']

class DiffProteinListFilter(django_filters.FilterSet):
    project = django_filters.CharFilter(field_name='project__name')
    gene = django_filters.CharFilter(field_name='fp__gene')
    proteome = django_filters.CharFilter(field_name='fp__ppid__proteome')
    organism = django_filters.CharFilter(field_name='fp__ppid__organism')
    logfc__gte = django_filters.NumberFilter(field_name='logfc', lookup_expr='gte')
    logfc__lte = django_filters.NumberFilter(field_name='logfc', lookup_expr='lte')
    d_p_value__lt = django_filters.NumberFilter(field_name='d_p_value', lookup_expr='lt')
    class Meta:
        model = DiffProtein
        fields = ['project', 'gene', 'proteome', 'organism', 'logfc', 'd_p_value__lt']
        
class PeptideListFilter(django_filters.FilterSet):
    project = django_filters.CharFilter(field_name='queue__project__name')
    protein__id = django_filters.CharFilter(field_name='protein__id')
    organism = django_filters.CharFilter(field_name='protein__fp__ppid__organism')
    proteome = django_filters.CharFilter(field_name='protein__fp__ppid__proteome')
    file = django_filters.ModelChoiceFilter(field_name='queue__filename',
                                            queryset=get_files_protein,
                                            label='Filename')
    fasta_type = django_filters.AllValuesFilter(lookup_expr='iexact')
    tag = django_filters.ModelChoiceFilter(field_name='queue__tag__name', 
                                           queryset=get_tags_protein,
                                           label='Tag')
    class Meta:
        model = Peptide
        fields = ['project', 'protein__id','sequence','fasta_type', 'organism', 'file']
        
class PsmListFilter(django_filters.FilterSet):
    accession = django_filters.CharFilter(field_name='peptide__protein__fp__accession')
    project = django_filters.CharFilter(field_name='queue__project__name')
    organism = django_filters.CharFilter(field_name='peptide__protein__fp__ppid__organism')
    proteome = django_filters.CharFilter(field_name='peptide__protein__fp__ppid__proteome')
    file = django_filters.ModelChoiceFilter(field_name='queue__filename',
                                            queryset=get_files_protein,
                                            label='Filename')
    fasta_type = django_filters.AllValuesFilter(lookup_expr='iexact')
    tag = django_filters.ModelChoiceFilter(field_name='queue__tag__name', 
                                           queryset=get_tags_protein,
                                           label='Tag')
    class Meta:
        model = Psm
        fields = ['project', 'sequence', 'mod_sequence', 'fasta_type', 'accession', 'organism', 'proteome', 'file']
        
class FileListFilter(django_filters.FilterSet):
    project = django_filters.CharFilter(field_name='project__name')
    class Meta:
        model = Queue
        fields = ['project', 'status', 'error']

class SpeciesListFilter(django_filters.FilterSet):
    project = django_filters.CharFilter(field_name='project__name')
    fasta_type = django_filters.AllValuesFilter(lookup_expr='iexact')
    class Meta:
        model = SpeciesSummary
        fields = ['project', 'fasta_type']

class FileSummaryFilter(django_filters.FilterSet):    
    project = django_filters.CharFilter(field_name='queue__project__name')
    file = django_filters.ModelChoiceFilter(field_name='queue__filename', queryset=get_files_protein)
    fasta_type = django_filters.AllValuesFilter(lookup_expr='iexact')
    class Meta:
        model = SpeciesFileSummary
        fields = ['project', 'file', 'fasta_type']