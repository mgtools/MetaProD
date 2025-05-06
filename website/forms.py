from crispy_forms.helper import FormHelper
from crispy_forms.layout import (
    Layout, 
    Submit, 
    Row, 
    Column, 
    Reset, 
    Div, 
    Field, 
    ButtonHolder,
    Fieldset
)
from crispy_forms.bootstrap import FormActions
from django import forms

class ProteinListFilterFormHelper(FormHelper):
    form_method = 'GET'
    layout = Layout(
        Div(
            Field('organism', wrapper_class='col-md-1', type="hidden"),
            Field('project', wrapper_class='col-md-1', type="hidden"),
            Field('file', title='Filename', wrapper_class='col-md-3'),
            Field('fasta_type', wrapper_class='col-md-1'),
            Field('tag', wrapper_class='col-md-1'),
            Submit('submit', 'Apply Filter'),
            css_class='form-row',
        ),
    )

class DiffProteinListFilterFormHelper(FormHelper):
    form_method = 'GET'
    layout = Layout(
        Div(
            Field('organism', wrapper_class='col-md-1', type="hidden"),
            Field('project', wrapper_class='col-md-1', type="hidden"),
            Field('logfc__gte', title='Log(2) FC', wrapper_class='col-md-2'),
            Field('logfc__lte', title='Log(2) FC', wrapper_class='col-md-2'),
            Field('d_p_value__lt', title='P-Value', wrapper_class='col-md-2'),
            Submit('submit', 'Apply Filter'),
            css_class='form-row',
        ),
    )
    
class PeptideListFilterFormHelper(FormHelper):
    form_method = 'GET'
    layout = Layout(
        Div(
            Field('protein__id', wrapper_class='col-md-1', type="hidden"),
            Field('organism', wrapper_class='col-md-1', type="hidden"),
            Field('project', wrapper_class='col-md-1', type="hidden"),
            Field('file', title='Filename', wrapper_class='col-md-3'),
            Field('fasta_type', wrapper_class='col-md-1'),
            Field('tag', wrapper_class='col-md-1'),
            Submit('submit', 'Apply Filter'),
            css_class='form-row',
        ),
    )
    
class PsmListFilterFormHelper(FormHelper):
    form_method = 'GET'
    layout = Layout(
        Div(
            Field('accession', wrapper_class='col-md-1', type="hidden"),
            Field('project', wrapper_class='col-md-1', type="hidden"),
            Field('file', title='Filename', wrapper_class='col-md-3'),
            Field('fasta_type', wrapper_class='col-md-1'),
            Field('tag', wrapper_class='col-md-1'),
            Submit('submit', 'Apply Filter'),
            css_class='form-row',
        ),
    )
    
class FileListFilterFormHelper(FormHelper):
    form_method = 'GET'
    layout = Layout(
        Div(
            Submit('submit', 'Apply Filter'),
            css_class='form-row',
        ),
    )

class SpeciesListFilterFormHelper(FormHelper):
    form_method = 'GET'
    layout = Layout(
        Div(
            Field('project', wrapper_class='col-md-1', type="hidden"),
            Field('fasta_type', wrapper_class='col-md-1'),
            Submit('submit', 'Apply Filter'),
            css_class='form-row',
        ),
    )

class FileSummaryFilterFormHelper(FormHelper):
    form_method = 'GET'
    layout = Layout(
        Div(
            Field('project', wrapper_class='col-md-1', type="hidden"),
            Field('file', title='Filename', wrapper_class='col-md-3'),
            Field('fasta_type', wrapper_class='col-md-1'),
            Submit('submit', 'Apply Filter'),
            css_class='form-row',
        ),
    )