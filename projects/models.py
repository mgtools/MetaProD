from django.db import models
from django.utils.translation import gettext_lazy as _
from django.contrib import messages
from django.core.validators import MaxValueValidator, MinValueValidator

class QueueManager(models.Manager):
    def get_by_natural_key(self, project, filename):
        return self.get(project=project, filename=filename)
        
# queue of files to run
class Queue(models.Model):

    objects = QueueManager()
    
    class Status(models.TextChoices):
        FILE_ADDED = 'FILE_ADDED', _('File added')
        THERMO = 'THERMO', _('Ready for thermorawfileparser')
        SEARCHGUI_PROF = 'SEARCHGUI_PROF', _('Ready for SearchGUI profile')
        PEPTIDESHAKER_PROF = 'PEPTIDESHAKER_PROF', _('Ready for PeptideShaker profile')
        REPORTER_PROF = 'REPORTER_PROF', ('Ready for Reporter profile')
        MZMINE_PROF = 'MZMINE_PROF', ('Ready for MZmine profile')
        READ_RESULTS_PROF = 'READ_RESULTS_PROF', _('Ready for read_results profile')
        PROCESS_RESULTS_PROF = 'PROCESS_RESULTS_PROF', _('Ready for process_results profile')
        FINISHED_PROF = 'FINISHED_PROF', _('File is finished profile step')
        SEARCHGUI_PROT = 'SEARCHGUI_PROT', _('Ready for SearchGUI proteome')
        PEPTIDESHAKER_PROT = 'PEPTIDESHAKER_PROT', _('Ready for PeptideShaker proteome')
        REPORTER_PROT = 'REPORTER_PROT', _('Ready for Reporter proteome')
        MZMINE_PROT = 'MZMINE_PROT', _('Ready for MZmine proteome')
        READ_RESULTS_PROT = 'READ_RESULTS_PROT', _('Ready for read_results proteome')
        PROCESS_RESULTS_PROT = 'PROCESS_RESULTS_PROT', _('Ready for process_results proteome')
        FINISHED_PROT = 'FINISHED_PROT', _('File is finished proteome step')
        FILE_FINISHED = 'FILE_FINISHED', _('File is finished and cleaned up')
    
    class Error(models.IntegerChoices):
        READY = 0, _('No Error')
        RETRY = 1, _('Retry')
        FAILED = 2, _('Failed')
        
    # full filename including any fraction data
    filename = models.CharField(
        max_length=255, 
        help_text="Filename without extension."
    )
    # sample grouping if there are fractions
    # this is unused for profiling but could be used for grouping of results
    # we'll set null if the sample is deleted because there could be reasons 
    # to delete the sample but want to keep the file
    sample = models.ForeignKey(
        'Sample', 
        on_delete=models.SET_NULL, 
        blank=True, 
        null=True, 
        help_text="Sample this file is associated with (for fractions)."
    )
    project = models.ForeignKey(
        'Project', 
        on_delete=models.CASCADE, 
        help_text="Project this file is associated with."
    )
    date_added = models.DateTimeField(
        auto_now_add=True, 
        help_text="Date file was added for processing."
    )
    date_finished_profile = models.DateTimeField(
        null=True, 
        blank=True, 
        help_text="Date file completed processing."
    )    
    date_finished_proteome = models.DateTimeField(
        null=True, 
        blank=True, 
        help_text="Date file completed processing."
    )
    status = models.TextField(
        choices=Status.choices, 
        help_text="File progress status."
    )
    error = models.IntegerField(
        choices=Error.choices, 
        default=0, 
        help_text="Error status."
    )
    skip = models.BooleanField(
        default=False, 
        help_text="Skip file from both queue and results."
    )
    job = models.IntegerField(
        default=0, 
        help_text="File job for HPC situations."
    )
    # this is just a calculated sum of the info in the runtimex
    total_runtime = models.IntegerField(
        default=0, 
        help_text="Total runtime including processing steps."
    )
    tag = models.ForeignKey(
        'Tag', 
        on_delete=models.SET_NULL, 
        blank=True, 
        null=True, 
        help_text="Optional tag for this file (e.g., phenotype)."
    )
    description = models.CharField(
        max_length=100,
        null=True,
        blank=True,
        help_text="Optional description for the file."
    )
    class Meta:
        unique_together = ('project', 'filename')
        
    def __str__(self):
        return self.filename
        
    def natural_key(self):
        return (self.project.name, self.filename)
        
    natural_key.dependencies = ['projects.project']

class SearchSetting(models.Model):
    # these are mostly for MSGF+
    class Fragmentation(models.IntegerChoices):
        AUTO = 0, _('AUTO')
        CID = 1, _('CID')
        ETD = 2, _('ETD')
        HCD = 3, _('HCD')
        UVPD = 4, _('UVPD')

    class Instrument(models.IntegerChoices):
        LTQ = 0, _('LTQ')
        ORBITRAP = 1, _('Orbitrap')
        TOF = 2, _('TOF')
        QEXACTIVE = 3, _('Q-Exactive')
    
    class Digestion(models.IntegerChoices):
        ENZYME = 0, ('Enzyme')
        UNSPECIFIC = 1, _('Unspecific')
        WHOLE_PROTEIN = 2, _('Whole Protein')
     
    class FragError(models.IntegerChoices):
        DA = 0, _('Da')
        PPM = 1, _('PPM')
    
    class ProfileType(models.IntegerChoices):
        FILE = 0, _('File')
        SAMPLE = 1, _('Sample')
        PROJECT = 2, _('Project')
    
    class ProfileMethod(models.IntegerChoices):
        NSAF = 0, _('NSAF')
        PSM = 1, _('PSM')
        PEAK_AREA = 2, _('Peak Area')
        
    project = models.OneToOneField(
        'Project', 
        on_delete=models.CASCADE, 
        help_text="Project name",
        primary_key=True
    )
    digestion = models.IntegerField(choices=Digestion.choices, default=0)
    min_peptide_length = models.IntegerField(
        default=8, 
        help_text="Minimum peptide length to include"
    )
    max_peptide_length = models.IntegerField(
        default=30, 
        help_text="Maximum peptide length to include"
    )
    min_charge = models.IntegerField(
        default=2, 
        help_text="Minimum precursor charge to include"
    )
    max_charge = models.IntegerField(
        default=4, 
        help_text="Maximum precursor charge to include"
        )
    psm_fdr = models.IntegerField(
        "PSM FDR", default=1, 
        help_text="False discovery rate at a PSM level"
    )
    peptide_fdr = models.IntegerField(
        "Peptide FDR", 
        default=1, 
        help_text="False discovery rate at a peptide level"
    )
    protein_fdr = models.IntegerField(
        "Protein FDR", 
        default=1, 
        help_text="False discovery rate at a protein level"
    )
    prec_tol = models.DecimalField(
        "Precursor Tolerance", 
        max_digits=5, 
        decimal_places=2, 
        default=10, 
        help_text="Precursor tolerance"
    )
    prec_ppm = models.IntegerField(
        "Precursor Unit", 
        choices=FragError.choices, 
        default=1, 
        help_text="Use daltons or PPM for precursor tolerance"
    )
    frag_tol = models.DecimalField(
        "Fragment Tolerance", 
        max_digits=5, 
        decimal_places=2, 
        default=10, 
        help_text="Fragment tolerance"
    )
    frag_ppm = models.IntegerField(
        "Fragment Unit", 
        choices=FragError.choices, 
        default=1, 
        help_text="Use daltons or PPM for fragment tolerance"
    )
    isotope_min = models.IntegerField(
        "Minimum precursor isotope",
        default=0,
        help_text = "Minumum precursor isotope."
    )
    isotope_max = models.IntegerField(
        "Maximum precursor isotope",
        default=1,
        help_text = "Maximum precursor isotope."
    )
    # frag/inst only affects msgf+
    # 0: default/cid, 1: cid, 2: etd, 3: hcd, 4: uvpd
    fragmentation = models.IntegerField(
        choices=Fragmentation.choices, 
        default=3, 
        help_text="Fragmentation method used by instrument"
    )
    # 0 ltq/ltq, 1: orbitrap, 2: tof, 3: qe
    instrument = models.IntegerField(
        choices=Instrument.choices, 
        default=3, 
        help_text="Type of instrument"
    )
    # use search engines
    # note other engines are supported in searchgui but may not work in linux
    xtandem_proteome = models.BooleanField(
        "Xtandem", 
        default=True, 
        help_text="Use Xtandem Search Algorithm"
    )
    msgf_proteome = models.BooleanField(
        "MSGF+", 
        default=True, 
        help_text="Use MSGF+ Search Algorithm"
    )
    omssa_proteome = models.BooleanField(
        "OMSSA", 
        default=False, 
        help_text="Use OMSSA Search Algorithm"
    )
    comet_proteome = models.BooleanField(
        "Comet", 
        default=True, 
        help_text="Use Comet Search Algorithm"
    )
    metamorpheus_proteome = models.BooleanField(
        "MetaMorpheus", 
        default=True, 
        help_text="Use MetaMorpheus Search Algorithm"
    )
    myrimatch_proteome = models.BooleanField(
        "MyriMatch", 
        default=False, 
        help_text="Use Myrimatch Search Algorithm"
    )              
    xtandem_profile = models.BooleanField(
        "Xtandem", 
        default=False, 
        help_text="Use Xtandem Search Algorithm"
    )
    msgf_profile = models.BooleanField(
        "MSGF+", 
        default=False, 
        help_text="Use MSGF+ Search Algorithm"
    )
    omssa_profile = models.BooleanField(
        "OMSSA", 
        default=False, 
        help_text="Use OMSSA Search Algorithm"
    )
    comet_profile = models.BooleanField(
        "Comet", 
        default=True, 
        help_text="Use Comet Search Algorithm"
    )
    metamorpheus_profile = models.BooleanField(
        "MetaMorpheus", 
        default=False, 
        help_text="Use MetaMorpheus Search Algorithm"
    )
    myrimatch_profile = models.BooleanField(
        "MyriMatch", 
        default=False, 
        help_text="Use Myrimatch Search Algorithm"
    )          
    mods = models.ManyToManyField(
        'ModList', 
        related_name="SearchSettingMods",
        db_index=True,
        through='ModChoice'
    )                                
    enzymes = models.ManyToManyField(
        'EnzymeList', 
        related_name="SearchSettingEnzymes",
        db_index=True, 
        through="EnzymeChoice"
    )
    # this dataset is multiplexed
    multiplex = models.BooleanField(
        default=False, 
        help_text="Data is multiplexed"
    )
    # use the homo sapien reference proteome in the FASTA
    use_human = models.BooleanField(
        default=True, 
        help_text="Include human proteome in FASTA file"
    )
    # pool results from project to profile versus profile on individual samples
    profile_type = models.IntegerField(
        default=1, 
        choices=ProfileType.choices, 
        help_text="Profiling method."
    )
    # use contaminants database in FASTA
    use_crap = models.BooleanField(
        default=True, 
        help_text="Include CRAP database in FASTA file"
    )
    # profile first? if False, then we do a single-step search
    profile = models.BooleanField(
        default=True, 
        help_text="Include a profile step"
    )
    profile_threshold = models.IntegerField(
        default=90, 
        help_text="Top percent of NSAF to include when profiling"
    )
    run_deqms = models.BooleanField(
        "Run DEqMS",
        default=True,
        help_text = "Run DEqMS on multiplexed data?"
    )
    imput_threshold = models.IntegerField(
        default=50,
        help_text="Percent of channels to require when filtering results before PEMM."
    )
    profile_method = models.IntegerField(
        choices=ProfileMethod.choices,
        default=0,
        help_text="LFQ method to use for profiling."
    )
    mzmine_run_mzmine = models.BooleanField(
        "Run MZmine",
        default=False,
        help_text = "Run MZmine to calculate peak areas."
    )
    mzmine_tpd_intensity = models.DecimalField(max_digits=3, decimal_places=2,
        default=0.50,
        validators=[
            MaxValueValidator(1),
            MinValueValidator(0)
        ],
        help_text="MZmine targeted peak detection intensity tolerance (from 0 to 1)."
    )
    mzmine_tpd_mztolerance = models.DecimalField(max_digits=4, decimal_places=3,
        default=0.001,
        validators=[
            MaxValueValidator(1),
            MinValueValidator(0.0001)
        ],
        help_text = "MZmine targeted peak detection absolute m/z tolerance."
    )
    mzmine_tpd_ppmtolerance = models.DecimalField(max_digits=3, decimal_places=1,
        default=10.0,
        validators=[
            MaxValueValidator(99.9),
            MinValueValidator(0.1)
        ],
        help_text = "MZmine targeted peak detection m/z PPM tolerance."
    )
        
    def __str__(self):
        return self.project.name
    
# keep a table of runtimes for each step
class RunTime(models.Model):
    queue = models.OneToOneField(
        'projects.Queue', 
        on_delete=models.CASCADE, 
    )
    msconvert = models.IntegerField(default=0)
    searchgui_profile = models.IntegerField(default=0)
    peptideshaker_profile = models.IntegerField(default=0)
    reporter_profile = models.IntegerField(default=0)
    mzmine_profile = models.IntegerField(default=0)
    read_results_profile = models.IntegerField(default=0)
    process_results_profile = models.IntegerField(default=0)
    searchgui_proteome = models.IntegerField(default=0)
    peptideshaker_proteome = models.IntegerField(default=0)
    reporter_proteome = models.IntegerField(default=0)
    mzmine_proteome = models.IntegerField(default=0)
    read_results_proteome = models.IntegerField(default=0)
    process_results_proteome = models.IntegerField(default=0)
    
    def natural_key(self):
        return (self.queue.natural_key(),)
        
class Project(models.Model):
    name = models.CharField(max_length=100, primary_key=True)
    description = models.CharField(max_length=250, blank=True, null=True)

    def __str__(self):
        return self.name
        
    def save(self, *args, **kwargs):
        disallowed = ['log', 'software', 'fasta', 'temp', 'metaprod']
        if self.name not in disallowed:
            super().save(*args, **kwargs)
        else:
            return
        
class Setting(models.Model):
    server = models.CharField(max_length=100, primary_key=True)
    install_folder = models.CharField(
        max_length=100, 
        blank=False, 
        null=False, 
        default='/home/metaprod/projects',
        help_text='Directory main files are installed in.'
    )
    memory = models.CharField(max_length=100, 
        blank=False, 
        null=False, 
        default='25G',
        help_text='How much memory to use?'
    )
    searchgui_ver = models.CharField(
        max_length=100, 
        blank=False, 
        null=False, 
        default='4.1.24',
        help_text = "SearchGUI version"
    )
    peptideshaker_ver = models.CharField(
        max_length=100, 
        blank=False, 
        null=False, 
        default='2.2.17',
        help_text = "PeptideShaker version"
    )
    reporter_ver = models.CharField(
        max_length=100, 
        blank=False, 
        null=False, 
        default='0.9.14',
        help_text="Reporter version"
    )
    mzmine_ver = models.CharField(
        max_length=100,
        blank=False,
        null=False,
        default='3.3.0',
        help_text="MZmine version"
    )
    threads = models.IntegerField(
        default=-1,
        help_text="Number of threads to use. -1 for all available."
    )
    debug_mode = models.BooleanField(default=1)
    data_folder = models.CharField(
        max_length=100, 
        blank=False, 
        null=False, 
        default='/home/metaprod/data'
    )
    max_retries = models.IntegerField(default=1)
    default = models.BooleanField(default=0)

class EnzymeList(models.Model):
    name = models.CharField(
        max_length=100, 
        primary_key=True
    )
    description = models.CharField(max_length=100, null=True, blank=True)
    
    def __str__(self):
        return self.name

       
class EnzymeChoice(models.Model):
    class Specificity(models.IntegerChoices):
        SPECIFIC = 0, _('Specific')
        SEMISPECIFIC = 1, _('Semi-Specific')
        NTERM = 2, _('N-Term Specific')
        CTERM = 3, _('C-Term Specific')
        
    enzyme = models.ForeignKey('EnzymeList', on_delete=models.CASCADE)
    searchsetting = models.ForeignKey(
        'SearchSetting', 
        on_delete=models.CASCADE
    )
    specificity = models.IntegerField(choices=Specificity.choices, default=0)
    mc = models.IntegerField(default=2, help_text='Missed Cleavages Allowed')
    
    def __str__(self):
        return self.enzyme.name
        
    class Meta:
        unique_together = ('searchsetting', 'enzyme')
    
    def natural_key(self):
        return(self.enzyme,) + self.searchsetting.natural_key()
    natural_key.dependencies = ['projects.project', 'projects.searchsetting']        
        
class ModList(models.Model):
    name = models.CharField(
        max_length=100,  
        primary_key=True
    )
    description = models.CharField(max_length=100, null=True, blank=True)

    def __str__(self):
        return self.name

class ModChoice(models.Model):
    class ModType(models.TextChoices):
        FIXED = 'Fixed'
        VARIABLE = 'Variable'
    mod = models.ForeignKey('ModList', on_delete=models.CASCADE)
    searchsetting = models.ForeignKey('SearchSetting', on_delete=models.CASCADE)
    modtype = models.CharField(choices=ModType.choices, max_length=10)

    def __str__(self):
        return self.mod.name
        
    class Meta:
        unique_together = ('searchsetting', 'mod')

    def natural_key(self):
        return(self.mod,) + self.searchsetting.natural_key()
    natural_key.dependencies = ['projects.searchsetting']
        

class SampleManager(models.Manager):
    def get_by_natural_key(self, project, name):
        return self.get(project=project, name=name)
        
# table of samples, i.e. for fractions
class Sample(models.Model):

    objects = SampleManager()
    
    name = models.CharField(
        max_length=255, 
        help_text="Name of sample (fractions should have same sample name)."
    )
    project = models.ForeignKey(
        'Project', 
        on_delete=models.CASCADE, 
        help_text="Project name"
    )
    description = models.CharField(max_length=100, blank=True, null=True)

    def __str__(self):
        return self.name
        
    class Meta:
        unique_together = ('name', 'project')
        
    def natural_key(self):
        return(self.project.name, self.name,)
    natural_key.dependencies = ['projects.project']

class TagManager(models.Manager):
    def get_by_natural_key(self, project, name):
        return self.get(project=project, name=name)
        
class Tag(models.Model):

    objects = TagManager()
        
    name = models.CharField(
        max_length=50,
        help_text="Name of tag."
    )
    project = models.ForeignKey(
        'Project', 
        on_delete=models.CASCADE, 
        help_text="Project name"
    )
    description = models.CharField(max_length=100, blank=True, null=True)
    
    def __str__(self):
        return self.name
        
    class Meta:
        unique_together = ('name', 'project')

    def natural_key(self):
        return(self.project.name, self.name,)
    natural_key.dependencies = ['projects.project']

class MetaDataManager(models.Manager):
    def get_by_natural_key(self, project, name):
        return self.get(project=project, name=name)
        
class MetaData(models.Model):

    objects = MetaDataManager()
    
    name = models.CharField(
        max_length=50,
        help_text="MetaData name."
    )
    project = models.ForeignKey(
        'Project',
        on_delete=models.CASCADE,
        help_text="Project name"
    )

    def __str__(self):
        return self.name
        
    class Meta:
        unique_together = ('name', 'project')

    def natural_key(self):
        return(self.project.name, self.name,)
    natural_key.dependencies = ['projects.project']

class MetaDataChoice(models.Model):

    metadata = models.ForeignKey('MetaData', on_delete=models.CASCADE)
    queue = models.ForeignKey('Queue', on_delete=models.CASCADE)
    value = models.CharField(max_length=25)

    def __str__(self):
        return self.metadata.name
        
    class Meta:
        unique_together = ('metadata', 'queue')

    def natural_key(self):
        return(self.project.name,) + self.metadata.natural_key() + self.queue.natural_key()
        
    natural_key.dependencies = ['projects.metadata', 'projects.queue']
    
# table of phenotypes, i.e. for multiplexed samples
class Phenotype(models.Model):
    phenotype = models.CharField(
        max_length=255, 
        help_text="Phenotype (i.e. normal, cancer, reference)",
        primary_key=True
    )
    description = models.CharField(
        max_length=100, 
        blank=True, 
        null=True, 
        help_text="Phenotype description."
    )

    def __str__(self):
        return self.phenotype
        
class MultiplexLabel(models.Model):
    project = models.ForeignKey('projects.Project', on_delete=models.CASCADE)
    # filename: 01CPTAC_COprospective_Proteome_PNNL_20170123
    sample = models.ForeignKey(
        'projects.Sample', 
        on_delete=models.SET_NULL, 
        null=True
    )
    reagent = models.ForeignKey('projects.reagent', null=True, 
                                on_delete=models.SET_NULL)
        
class LabelChoice(models.Model):
    multiplexlabel = models.ForeignKey('MultiplexLabel', on_delete=models.CASCADE)
    label = models.ForeignKey('Label', on_delete=models.CASCADE)
    # phenotype: colon, normal, reference, etc.
    phenotype = models.ForeignKey(
        'projects.Phenotype', 
        on_delete=models.SET_NULL, 
        null=True
    )
    identifier = models.CharField(
        max_length=255, 
        help_text=("Unique patient identifier. Must be consistent so that the "
            "same patient has the same ID for each sample they are in."
        )
    )

class Reagent(models.Model):
    name = models.CharField(max_length=255, primary_key=True)
    description = models.CharField(max_length=250, null=True, blank=True)

    def __str__(self):
        return self.name
        
class Label(models.Model):
    reagent = models.ForeignKey('projects.reagent', on_delete=models.CASCADE)
    name = models.CharField(max_length=100)
    description = models.CharField(max_length=250, null=True, blank=True)

    class Meta:
        unique_together = ('reagent', 'name')
        
    def natural_key(self):
        return(self.reagent, self.name,)
        
    def __str__(self):
        return self.name or ''