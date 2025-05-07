from django.contrib import admin
from django.utils.translation import gettext_lazy as _
from django import forms
from django.forms import ModelForm
from dal import autocomplete
from django.contrib import messages

from .models import (
    Queue, 
    Project,
    Setting, 
    RunTime, 
    Sample,
    SearchSetting, 
    ModList, 
    ModChoice, 
    EnzymeList, 
    EnzymeChoice,
    MultiplexLabel,
    Label,
    LabelChoice,
    Reagent,
    Tag,
    MetaData,
    MetaDataChoice,
    EngineStatus
)

from .forms import MultiplexLabelForm, MultiplexLabelInlineForm, QueueForm, MetaDataChoiceForm

admin.site.site_header = 'MetaProD Admin'
admin.site.site_title = "MetaProD Admin"

class SettingAdmin(admin.ModelAdmin):
    list_display = ('server', 'memory', 'threads', 'default')
    list_display_links = ('server', 'memory', 'threads', 'default')

class RunTimeInline(admin.TabularInline):
    model = RunTime
    extra = 0
    max_num = 0
    can_delete = False
    readonly_fields = ('msconvert','searchgui_profile','peptideshaker_profile',
        'reporter_profile', 'mzmine_profile', 'read_results_profile', 'process_results_profile',
        'searchgui_proteome', 'peptideshaker_proteome', 'reporter_proteome', 'mzmine_proteome',
        'read_results_proteome', 'process_results_proteome')

class EngineStatusInline(admin.TabularInline):
    model = EngineStatus
    extra = 0
    max_num = 0
    can_delete = False
    readonly_fields = ('comet_profile', 'xtandem_profile', 'msgf_profile', 'omssa_profile', 'metamorpheus_profile', 'myrimatch_profile', 'sage_profile',
                       'comet_proteome', 'xtandem_proteome', 'msgf_proteome', 'omssa_proteome', 'metamorpheus_proteome', 'myrimatch_proteome', 'sage_proteome')
    exclude = ('comet_tries', 'xtandem_tries', 'msgf_tries', 'omssa_tries', 'metamorpheus_tries', 'myrimatch_tries', 'sage_tries')
    
class ProjectListFilter(admin.SimpleListFilter):
    title = _('Project')
    parameter_name = 'project__name'

    def lookups(self, request, model_admin):
        projects = Project.objects.all()
        return [(project.name, project.name) for project in projects]

    def queryset(self, request, queryset):
        value = self.value()
        if value is not None:
            return queryset.filter(project__name=self.value())
        return queryset

class QueueProjectListFilter(admin.SimpleListFilter):
    title = _('Project')
    parameter_name = 'queue__project__name'

    def lookups(self, request, model_admin):
        projects = Project.objects.all()
        return [(project.name, project.name) for project in projects]

    def queryset(self, request, queryset):
        value = self.value()
        if value is not None:
            return queryset.filter(queue__project__name=self.value())
        return queryset
 
class LabelReagentListFilter(admin.SimpleListFilter):
    title = _('Reagent')
    parameter_name = 'reagent__name'

    def lookups(self, request, model_admin):
        reagents = Reagent.objects.all()
        return [(reagent.name, reagent.name) for reagent in reagents]

    def queryset(self, request, queryset):
        value = self.value()
        if value is not None:
            return queryset.filter(reagent__name=self.value())
        return queryset

class MetaDataChoiceInline(admin.TabularInline):
    model = MetaDataChoice
    form = MetaDataChoiceForm
    extra = 0

    def metadata(self, obj):
        return obj.metadata.name
        
    def project(self, obj):
        return obj.queue.project
        
class QueueAdmin(admin.ModelAdmin):

    inlines = (RunTimeInline, EngineStatusInline, MetaDataChoiceInline)
    list_filter = (ProjectListFilter,)
    fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': (
                ('project', 'job'),
                'filename',
                'sample', 
                'status', 
                'error',
                'tag',
                'description',
                'skip',
                ('date_added', 'date_finished_profile', 'date_finished_proteome'),
                'total_runtime',
            )
        }),
    )
    form = QueueForm
    list_display = ('id', 'project', 'filename', 'get_sample', 'status', 'tag', 
                    'error', 'skip', 'job')
    list_display_links = ('id', 'project', 'filename', 'get_sample', 'status', 
                          'tag', 'error', 'skip', 'job')
    readonly_fields = ('total_runtime', 'filename', 'date_added', 'date_finished_profile', 'date_finished_proteome')

    @admin.display(ordering='sample__name', description='Sample')
    def get_sample(self, obj):
        if obj.sample is not None:
            return obj.sample.name
        else:
            return obj.sample

class EnzymeChoiceInline(admin.TabularInline):
    model = EnzymeChoice
    extra = 0

class ModChoiceInline(admin.TabularInline):
    model = ModChoice
    extra = 0
    
class SearchSettingInline(admin.StackedInline):
    model = SearchSetting
    extra = 0
    max_num = 1
    fieldsets = (
        ('MetaProD', {
            'classes': ('collapse in'),
            'fields': (('project'),
                ('use_crap'), ('use_human', 'human_fasta'),
                ('custom_fasta'),
                ('perform_second_step'),
                ('profile_type', 'profile_threshold', 'profile_method'),
                ('multiplex'), ('run_deqms'), ('mzmine_run_mzmine'),
                ('imput_threshold')
            ),
            'description': 'MetaProD specific options',
        }),    
        ('SearchGUI/PeptideShaker', {
            'classes': ('collapse', 'extrapretty'),
            'fields': (('min_peptide_length', 'max_peptide_length'),
                ('min_charge', 'max_charge'),
                ('prec_tol', 'prec_ppm'),
                ('frag_tol', 'frag_ppm'),
                ('isotope_min', 'isotope_max'),
                ('instrument', 'fragmentation'),
                ('psm_fdr', 'peptide_fdr', 'protein_fdr'),
                ('digestion')
            ),
            'description': 'SearchGUI/PeptideShaker specific options',
        }), 
        ('Profile/custom search engines', {
            'classes': ('collapse', 'extrapretty'),
            'fields': (('comet_profile'), ('metamorpheus_profile'), ('msgf_profile'), ('myrimatch_profile'), ('omssa_profile'), ('xtandem_profile'), ('sage_profile'),
            ),
        }),
        ('Proteome search engines', {
            'classes': ('collapse', 'extrapretty'),
            'fields': (('comet_proteome'), ('metamorpheus_proteome'), ('msgf_proteome'), ('myrimatch_proteome'), ('omssa_proteome'), ('xtandem_proteome'), ('sage_proteome'),
            ),
        }),
        ('MZmine', {
            'classes': ('collapse', 'extrapretty'),
            'fields': (
                'mzmine_tpd_intensity',
                'mzmine_tpd_mztolerance',
                'mzmine_tpd_ppmtolerance'
            )
        }),
    )
    readonly_fields = ('profile_method',)
        
class ProjectAdmin(admin.ModelAdmin):
    list_display = ('name', 'description')
    list_display_links = ('name', 'description')
    #inlines = (SearchSettingInline,)

    # this will end up also showing the success message but that's a minor issue
    def save_model(self, request, obj, form, change):
        disallowed = ['log', 'software', 'fasta', 'temp', 'metaprod']
        if obj.name in disallowed:
            messages.error(request, 'Log, software, fasta, and metaprod cannot be the name of a project.')
        else:
            super(ProjectAdmin, self).save_model(request, obj, form, change)
        
class SearchSettingAdmin(admin.ModelAdmin):
    inlines = (EnzymeChoiceInline, ModChoiceInline)
    list_display = ('project', 'instrument', 'fragmentation', 'modification_list', 'enzyme_list', 'multiplex',
                    'mzmine_run_mzmine')
    list_display_links = ('project',)
    
    def modification_list(self, obj):
        return ", ".join([a.name for a in obj.mods.all()])

    def enzyme_list(self, obj):
        return ", ".join([a.name for a in obj.enzymes.all()])

    fieldsets = (
        ('MetaProD', {
            'classes': ('collapse in'),
            'fields': (('project'),
                ('use_crap'), ('use_human', 'human_fasta'),
                ('custom_fasta'),
                ('perform_second_step'),
                ('profile_type', 'profile_method'),
                ('profile_threshold', 'profile_exclude_below', 'profile_include_above'),
                ('multiplex'), ('run_deqms'), ('mzmine_run_mzmine'),
                ('imput_threshold')
            ),
            'description': 'MetaProD specific options',
        }),    
        ('SearchGUI/PeptideShaker', {
            'classes': ('collapse', 'extrapretty'),
            'fields': (('min_peptide_length', 'max_peptide_length'),
                ('min_charge', 'max_charge'),
                ('prec_tol', 'prec_ppm'),
                ('frag_tol', 'frag_ppm'),
                ('isotope_min', 'isotope_max'),
                ('instrument', 'fragmentation'),
                ('psm_fdr', 'peptide_fdr', 'protein_fdr'),
                ('digestion')
            ),
            'description': 'SearchGUI/PeptideShaker specific options',
        }), 
        ('Profile/custom search engines', {
            'classes': ('collapse', 'extrapretty'),
            'fields': (('comet_profile'), ('metamorpheus_profile'), ('msgf_profile'), ('myrimatch_profile'), ('omssa_profile'), ('xtandem_profile'), ('sage_profile'),
            ),
        }),
        ('Proteome search engines', {
            'classes': ('collapse', 'extrapretty'),
            'fields': (('comet_proteome'), ('metamorpheus_proteome'), ('msgf_proteome'), ('myrimatch_proteome'), ('omssa_proteome'), ('xtandem_proteome'), ('sage_proteome'),
            ),
        }),
        ('MZmine', {
            'classes': ('collapse', 'extrapretty'),
            'fields': (
                'mzmine_tpd_intensity',
                'mzmine_tpd_mztolerance',
                'mzmine_tpd_ppmtolerance'
            )
        }),
        ('Search engine specific settings', {
            'classes': ('collapse', 'extrapretty'),
            'fields': (
                ('comet_batch_size'),
                ('msgf_num_tasks'),
            )
        }),        
    )
    readonly_fields = ('profile_method',)

class ModListAdmin(admin.ModelAdmin):
    list_display = ('name','description')
    list_display_links = ('name',)

class ModChoiceAdmin(admin.ModelAdmin):
    list_filter = (ProjectListFilter,)
    list_display = ('project', 'mod', 'modtype',)
    list_display_links = ('project', 'mod', 'modtype')

    def mod(self, obj):
        return obj.mod.name
        
    def project(self, obj):
        return obj.searchsetting.project

class EnzymeListAdmin(admin.ModelAdmin):
    list_display = ('name','description')
    list_display_links = ('name',)

class EnzymeChoiceAdmin(admin.ModelAdmin):
    list_display = ('project', 'enzyme', 'specificity', 'mc',)
    list_display_links = ('project', 'enzyme', 'specificity', 'mc')

    def enzyme(self, obj):
        return obj.enzyme.name
        
    def project(self, obj):
        return obj.searchsetting.project

class RunTimeAdmin(admin.ModelAdmin):
    list_filter = (QueueProjectListFilter,)
    list_display = ('project', 'queue', 'msconvert', 
                    'searchgui_profile', 'peptideshaker_profile', 
                    'reporter_profile', 'mzmine_profile', 
                    'read_results_profile', 'process_results_profile', 
                    'searchgui_proteome', 'peptideshaker_proteome', 
                    'reporter_proteome', 'mzmine_proteome',
                    'read_results_proteome', 'process_results_proteome')
    list_display_links = ('project', 'queue', 'msconvert', 
                          'searchgui_profile', 'peptideshaker_profile', 
                          'reporter_profile', 'mzmine_profile',
                          'read_results_profile', 'process_results_profile', 
                          'searchgui_proteome', 'peptideshaker_proteome',
                          'reporter_proteome', 'mzmine_proteome',
                          'read_results_proteome', 'process_results_proteome')
    
    def project(self, obj):
        return obj.queue.project.name
        
class SampleAdmin(admin.ModelAdmin):
    list_filter = (ProjectListFilter,)
    list_display = ('project', 'name')
    list_display_links = ('project', 'name')
    ordering = ['name']

class LabelChoiceInline(admin.TabularInline, autocomplete.Select2QuerySetView):
    model = LabelChoice
    form = MultiplexLabelInlineForm
    extra = 11
    max_num = 11
        
class LabelInline(admin.TabularInline):
    model = Label
    extra = 11
    max_num = 11

class MultiplexLabelAdmin(admin.ModelAdmin):
    list_filter = (ProjectListFilter,)
    form = MultiplexLabelForm
    #fields = ('project', 'sample', 'reagent')
    list_display = ('project', 'sample', 'reagent')
    list_display_links = ('project', 'sample', 'reagent')
    inlines = (LabelChoiceInline,)

class LabelAdmin(admin.ModelAdmin):
    list_filter = (LabelReagentListFilter,)
    list_display = ('reagent', 'name', 'description')
    list_display_links = ('reagent', 'name', 'description')

class PhenotypeAdmin(admin.ModelAdmin):
    list_display = ('phenotype', 'description')
    list_display_links = ('phenotype', 'description')

class ReagentAdmin(admin.ModelAdmin):
    list_display = ('name', 'description')
    list_display_links = ('name', 'description')
    inlines = (LabelInline,)

class TagAdmin(admin.ModelAdmin):
    list_filter = (ProjectListFilter,)
    list_display = ('project', 'name', 't_type')
    list_display_links = ('project', 'name', 't_type')

class MetaDataAdmin(admin.ModelAdmin):
    list_filter = (ProjectListFilter,)
    list_display = ('project', 'name')
    list_display_links = ('project', 'name')

class MetaDataChoiceAdmin(admin.ModelAdmin):
    list_display = ('metadata', 'value')
    list_display_links = ('metadata', 'value')

    form = MetaDataChoiceForm
    
    def metadata(self, obj):
        return obj.metadata.name
        
    def project(self, obj):
        return obj.queue.project
        
admin.site.register(Queue, QueueAdmin)
admin.site.register(Project, ProjectAdmin)
admin.site.register(Setting, SettingAdmin)
admin.site.register(SearchSetting, SearchSettingAdmin)
admin.site.register(ModList, ModListAdmin)
admin.site.register(EnzymeList, EnzymeListAdmin)
admin.site.register(MultiplexLabel, MultiplexLabelAdmin)
admin.site.register(Sample, SampleAdmin)
admin.site.register(Label, LabelAdmin)
admin.site.register(Reagent, ReagentAdmin)
admin.site.register(Tag, TagAdmin)
admin.site.register(MetaData, MetaDataAdmin)