import os

from django.utils.html import format_html
from django_tables2 import Column, Table, A, LinkColumn, TemplateColumn

from projects.models import Project, Queue
from results.models import (
    Peptide, 
    Psm, 
    Protein, 
    SpeciesSummary,
    SpeciesFileSummary,
    DiffProtein
)

class ProjectListTable(Table):
    summary = TemplateColumn(
        orderable=False, 
        verbose_name='Summary', 
        template_code=
            '''<a href="{% url 'summary' record.name %}">Link</a>'''
    )
    files = TemplateColumn(
        orderable=False, 
        verbose_name='Files', 
        template_code=
            '''<a href="{% url 'file_list' %}?project={{ record.name }}">Link</a>'''
    )
    species = TemplateColumn(
        orderable=False, 
        verbose_name='Species', 
        template_code=
            '''<a href="{% url 'species_list' %}?project={{ record.name }}">Link</a>'''
    )        
    proteins = TemplateColumn(
        orderable=False, 
        verbose_name='Proteins', 
        template_code=
            '''<a href="{% url 'protein_list' %}?project={{ record.name }}">Link</a>'''
    )
    peptides = TemplateColumn(
        orderable=False, 
        verbose_name='Peptides', 
        template_code=
            '''<a href="{% url 'peptide_list' %}?project={{ record.name }}">Link</a>'''
    )
    psms = TemplateColumn(
        orderable=False, 
        verbose_name='PSMs', 
        template_code=
            '''<a href="{% url 'psm_list' %}?project={{ record.name }}">Link</a>'''
    )
    diffproteins = TemplateColumn(
        orderable=False, 
        verbose_name='Diff. Proteins', 
        template_code=
            '''<a href="{% url 'diffprotein_list' %}?project={{ record.name }}&d_p_value__lt=0.05">Link</a>'''
    )
    class Meta:
        model = Project
        fields = ('name', 'description', 'summary', 'files', 'species', 
                  'proteins', 'peptides', 'psms', 'diffproteins')
        attrs = {"class": "table table-hover"}
        order_by = 'name'

# https://www.ncbi.nlm.nih.gov/gene?term=(tns1%5BGene%20Name%5D)%20AND%20homo%20sapiens%5BOrganism%5D        
class ProteinListTable(Table):
    fp__accession = TemplateColumn(
        verbose_name='Accession', 
        template_code=
            '''<a href="https://www.uniprot.org/uniprot/{{ record.fp.accession }}">{{ record.fp.accession }}</a>'''
    )
    val_num_peptide = TemplateColumn(
        verbose_name='# Peptides', 
        template_code=
            '''<a href="{% load spurl %}{% spurl base="{\% url 'peptide_list' %\}" add_query="protein__id={{ record.id }}" add_query_from=request.get_full_path %}">{{ record.val_num_peptide }}</a>'''
    )
    fp__gene = TemplateColumn(
        verbose_name='Gene', 
        template_code=
            '''<a href="{% load spurl %}{% spurl base="{\% url 'protein_list' %\}" add_query="gene={{ record.fp.gene }}" add_query_from=request.get_full_path %}">{{ record.fp.gene }}</a>''')
    fp__ppid__proteome = TemplateColumn(
        verbose_name='Proteome', 
        template_code=
            '''<a href="{% load spurl %}{% spurl base="{\% url 'protein_list' %\}" add_query="proteome={{ record.fp.ppid.proteome }}" add_query_from=request.get_full_path %}">{{ record.fp.ppid.proteome }}</a>''')
    fp__ppid__organism = TemplateColumn(
        verbose_name='Organism', 
        template_code=
            '''<a href="{% load spurl %}{% spurl base="{\% url 'protein_list' %\}" add_query="organism={{ record.fp.ppid.organism }}" add_query_from=request.get_full_path %}">{{ record.fp.ppid.organism }}</a>'''
    )
    val_num_psm = TemplateColumn(
        verbose_name='# PSM', 
        template_code=
            '''<a href="{% load spurl %}{% spurl base="{\% url 'psm_list' %\}" add_query="accession={{ record.fp.accession }}" add_query_from=request.get_full_path %}">{{ record.val_num_psm }}</a>'''
    )
    nsaf = Column(verbose_name='NSAF')
    class Meta:
        model = Protein
        fields = ('fp__accession', 'fp__gene', 'fp__ppid__proteome', 'fp__ppid__organism', 
                  'fp__description', 'val_num_peptide', 'val_num_psm', 'nsaf')
        attrs = {"class": "table table-hover"}
        order_by = '-nsaf'
        export_formats = ['csv']

class DiffProteinListTable(Table):
    fp__accession = TemplateColumn(
        verbose_name='Accession', 
        template_code=
            '''<a href="https://www.uniprot.org/uniprot/{{ record.fp.accession }}">{{ record.fp.accession }}</a>'''
    )
    fp__gene = TemplateColumn(
        verbose_name='Gene', 
        template_code=
            '''<a href="{% load spurl %}{% spurl base="{\% url 'diffprotein_list' %\}" add_query="gene={{ record.fp.gene }}" add_query_from=request.get_full_path %}">{{ record.fp.gene }}</a>''')
    fp__ppid__proteome = TemplateColumn(
        verbose_name='Proteome', 
        template_code=
            '''<a href="{% load spurl %}{% spurl base="{\% url 'diffprotein_list' %\}" add_query="proteome={{ record.fp.ppid.proteome }}" add_query_from=request.get_full_path %}">{{ record.fp.ppid.proteome }}</a>''')
    fp__ppid__organism = TemplateColumn(
        verbose_name='Organism', 
        template_code=
            '''<a href="{% load spurl %}{% spurl base="{\% url 'diffprotein_list' %\}" add_query="organism={{ record.fp.ppid.organism }}" add_query_from=request.get_full_path %}">{{ record.fp.ppid.organism }}</a>'''
    )
    d_p_value = Column(verbose_name='P-Value (DEqMS)')
    logfc = Column(verbose_name='Log(2) FC')
    class Meta:
        model = DiffProtein
        fields = ('fp__accession', 'fp__gene', 'fp__ppid__proteome', 'fp__ppid__organism', 
                  'fp__description', 'logfc', 'd_p_value')
        attrs = {"class": "table table-hover"}
        order_by = 'accession'
        export_formats = ['csv']
        
class PeptideListTable(Table):
    protein__fp__accession = TemplateColumn(
        verbose_name='Accession', 
        template_code=
            '''<a href="https://www.uniprot.org/uniprot/{{ record.protein.fp.accession }}">{{ record.protein.fp.accession }}</a>'''
    )
    val_num_psm = TemplateColumn(
        verbose_name='# PSMs', 
        template_code=
            '''<a href="{% load spurl %}{% spurl base="{\% url 'psm_list' %\}" add_query="mod_sequence={{ record.mod_sequence }}" add_query_from=request.get_full_path %}">{{ record.val_num_psm }}</a>'''
    )
    organism = Column(accessor='protein__fp__ppid__organism', 
                      verbose_name='Organism', visible=False)
    proteome = Column(accessor='protein__fp__ppid__proteome', 
                      verbose_name='Proteome', visible=False)
    class Meta:
        model = Peptide
        fields = ('proteome', 'organism', 'protein__fp__accession', 
                  'sequence', 'mod_sequence', 'val_num_psm')
        attrs = {"class": "table table-hover"}
        order_by = '-val_num_psm'
        
class PsmListTable(Table):
    peptide__protein__fp__accession = TemplateColumn(
        verbose_name='Accession', 
        template_code=
            '''<a href="https://www.uniprot.org/uniprot/{{ record.peptide.protein.fp.accession }}">{{ record.peptide.protein.fp.accession }}</a>'''
    )
    rt = Column(verbose_name='RT')
    fixed_ptm = Column(verbose_name='Fixed PTMs')
    variable_ptm = Column(verbose_name='Variable PTMs')
    mod_sequence = Column(verbose_name='Modified Sequence')
    organism = Column(accessor='peptide__protein__fp__ppid__organism', 
                      verbose_name='Organism', 
                      visible=False)
    proteome = Column(accessor='peptide__protein__fp__ppid__proteome', 
                      verbose_name='Proteome', 
                      visible=False)
    class Meta:
        model = Psm
        fields = ('proteome', 'organism', 'peptide__protein__fp__accession', 
                  'sequence', 'mod_sequence', 'fixed_ptm', 'variable_ptm', 
                  'rt', 'charge')
        attrs = {"class": "table table-hover"}
        order_by = 'sequence'
        
class FileListTable(Table):
    summary = TemplateColumn(
        orderable=False, 
        verbose_name='File Summary', 
        template_code=
            '''<a href="{% url 'file_summary' %}?project={{ record.project.name }}&file={{ record.id }}">Link</a>'''
    )  
    proteins = TemplateColumn(
        orderable=False, 
        verbose_name='Proteins', 
        template_code=
            '''<a href="{% url 'protein_list' %}?project={{ record.project.name }}&file={{ record.id }}">Link</a>'''
    )
    peptides = TemplateColumn(
        orderable=False, 
        verbose_name='Peptides', 
        template_code=
            '''<a href="{% url 'peptide_list' %}?project={{ record.project.name }}&file={{ record.id }}">Link</a>'''
    )
    psms = TemplateColumn(
        orderable=False, 
        verbose_name='PSMs', 
        template_code=
            '''<a href="{% url 'psm_list' %}?project={{ record.project.name }}&file={{ record.id }}">Link</a>'''
    )      
    class Meta:
        model = Queue
        fields = ('filename', 'project', 'total_runtime', 'summary', 'proteins', 'peptides', 'psms')
        attrs = {"class": "table table-hover"}
        order_by = 'id'
  
class FileSummaryTable(Table):
    proteome = TemplateColumn(
        verbose_name="Proteome",
        template_code=
            '''<a href="https://www.uniprot.org/proteomes/{{ record.ppid.proteome }}">{{ record.ppid.proteome }}</a>'''
    )
    organism = TemplateColumn(
        verbose_name="Organism",
        template_code=
            '''<a href="https://www.uniprot.org/proteomes/{{ record.ppid.proteome }}">{{ record.ppid.organism }}</a>'''
    )    
    proteins = TemplateColumn(
        orderable=False, 
        verbose_name='Proteins', 
        template_code=
            '''<a href="{% url 'protein_list' %}?project={{ record.queue.project.name }}&file={{ record.queue.id }}&proteome={{ record.ppid.proteome }}&type={{ record.type|capfirst }}">{{ record.val_num_protein }}</a>'''
    )
    peptides = TemplateColumn(
        orderable=False, 
        verbose_name='Peptides', 
        template_code=
            '''<a href="{% url 'peptide_list' %}?project={{ record.queue.project.name }}&file={{ record.queue.id }}&proteome={{ record.ppid.proteome }}&type={{ record.type|capfirst }}">{{ record.val_num_peptide }}</a>'''
    )
    psms = TemplateColumn(
        orderable=False, 
        verbose_name='PSMs', 
        template_code=
            '''<a href="{% url 'psm_list' %}?project={{ record.queue.project.name }}&file={{ record.queue.id }}&proteome={{ record.ppid.proteome }}&type={{ record.type|capfirst }}">{{ record.val_num_psm }}</a>'''
    )
    nsaf = Column(verbose_name='NSAF') 
    class Meta:
        model = SpeciesFileSummary
        fields = ('proteome', 'organism', 'proteins', 'peptides', 'psms', 'nsaf')
        attrs = {"class": "table table-hover"}
        order_by = '-nsaf'
        
class SpeciesListTable(Table):
    proteome = TemplateColumn(
        verbose_name="Proteome",
        template_code=
            '''<a href="https://www.uniprot.org/proteomes/{{ record.ppid.proteome }}">{{ record.ppid.proteome }}</a>'''
    )
    organism = TemplateColumn(
        verbose_name="Organism",
        template_code=
            '''<a href="https://www.uniprot.org/proteomes/{{ record.ppid.proteome }}">{{ record.ppid.organism }}</a>'''
    )
    val_num_protein = TemplateColumn(
        verbose_name='# Proteins', 
        template_code=
            '''<a href="{% load spurl %}{% spurl base="{\% url 'protein_list' %\}" add_query="proteome={{ record.ppid.proteome }}" add_query_from=request.get_full_path %}">{{ record.val_num_protein }}</a>'''
    )
    val_num_peptide = TemplateColumn(
        verbose_name='# Peptides', 
        template_code=
            '''<a href="{% load spurl %}{% spurl base="{\% url 'peptide_list' %\}" add_query="proteome={{ record.ppid.proteome }}" add_query_from=request.get_full_path %}">{{ record.val_num_peptide }}</a>'''
    )
    val_num_psm = TemplateColumn(verbose_name='# PSMs', template_code=
        '''<a href="{% load spurl %}{% spurl base="{\% url 'psm_list' %\}" add_query="proteome={{ record.ppid.proteome }}" add_query_from=request.get_full_path %}">{{ record.val_num_psm }}</a>''')
    nsaf = Column(verbose_name="NSAF")
    class Meta:
        model = SpeciesSummary
        fields = ('proteome', 'organism', 'val_num_protein', 
                  'val_num_peptide', 'val_num_psm', 'nsaf')
        attrs = {"class": "table table-hover"}
        order_by = '-nsaf'