from django.db import models

TYPE_CHOICES = (
    ('Profile', 'Profile'),
    ('Proteome', 'Proteome'),
)

class ProteinManager(models.Manager):
    def get_by_natural_key(self, accession, type, queue_project, queue_filename):
        return self.get(fp__accession=accession, type=type, queue__project=queue_project, queue__filename=queue_filename)
        
# table of calculated inferences (as opposed to using peptideshaker inferences)
class Protein(models.Model):
    objects = ProteinManager()
    queue = models.ForeignKey('projects.Queue', on_delete=models.CASCADE)
    fp = models.ForeignKey('results.FastaProtein', on_delete=models.CASCADE)
    # this needs to be recalculated if we split PSMs for ambiguous inferences
    val_num_psm = models.DecimalField(max_digits=15, decimal_places=3)
    val_num_peptide = models.IntegerField()
    saf = models.DecimalField(null=True, max_digits=15, decimal_places=6)
    nsaf = models.DecimalField(null=True, max_digits=15, decimal_places=6)
    validation = models.CharField(max_length=20)
    type = models.CharField(max_length=20, choices=TYPE_CHOICES)
    peak_area = models.DecimalField(null=True, max_digits=19, decimal_places=3)
    peak_area_psm = models.IntegerField(null=True)
    class Meta:
        unique_together = ('queue', 'fp', 'type')
        
    def natural_key(self):
        return (self.fp.accession, self.type) + self.queue.natural_key()

    natural_key.dependencies = ['projects.queue', 'results.fastaprotein']

class PeptideManager(models.Manager):
    def get_by_natural_key(self, mod_sequence, type, queue_project, queue_filename):
        return self.get(mod_sequence=mod_sequence, type=type, queue__project=queue_project, queue__filename=queue_filename)
        
# set_null on delete of protein as we may be re-running process_results
class Peptide(models.Model):
    objects = PeptideManager()
    queue = models.ForeignKey('projects.Queue', on_delete=models.CASCADE)
    # parsimony inference
    protein = models.ForeignKey(
        'results.Protein', 
        on_delete=models.SET_NULL, 
        null=True
    )
    # possible proteins as decided by peptideshaker
    accessions = models.TextField()
    # sequence without mods (may not be unique)
    sequence = models.CharField(max_length=255)
    # sequence with mods (should be unique)
    mod_sequence = models.CharField(max_length=255)
    variable_ptm = models.CharField(max_length=255, blank=True)
    fixed_ptm = models.CharField(max_length=255, blank=True)
    val_num_psm = models.IntegerField()
    validation = models.CharField(max_length=20)
    type = models.CharField(max_length=20, choices=TYPE_CHOICES)
    peak_area = models.DecimalField(null=True, max_digits=19, decimal_places=3)
    # how many PSMs had a peak area as it may not be the same as val_num_psm
    peak_area_psm = models.IntegerField(null=True)
    class Meta:
        unique_together = ('queue', 'mod_sequence', 'type')
        
    def natural_key(self):
        return (self.mod_sequence, self.type) + self.queue.natural_key()

    natural_key.dependencies = ['projects.queue']

class PsmManager(models.Manager):
    def get_by_natural_key(self, type, title, queue_project, queue_filename):
        return self.get(type=type, title=title, queue__project=queue_project, queue__filename=queue_filename)
        
class Psm(models.Model):
    objects = PsmManager()
    queue = models.ForeignKey('projects.Queue', on_delete=models.CASCADE)
    peptide = models.ForeignKey(
        'results.Peptide', 
        on_delete=models.CASCADE, 
        null=True
    )
    accessions = models.TextField(max_length=255)
    sequence = models.CharField(max_length=255)
    mod_sequence = models.CharField(max_length=255)
    variable_ptm = models.CharField(max_length=255, blank=True)
    fixed_ptm = models.CharField(max_length=255, blank=True)
    rt = models.DecimalField(max_digits=15, decimal_places=7)
    mz = models.DecimalField(max_digits=22, decimal_places=15)
    error = models.DecimalField(max_digits=25, decimal_places=15)
    charge = models.CharField(max_length=255)
    validation = models.CharField(max_length=20)
    confidence = models.DecimalField(max_digits=15, decimal_places=6)
    title = models.CharField(max_length=255)
    type = models.CharField(max_length=20, choices=TYPE_CHOICES)
    peak_area = models.DecimalField(null=True, max_digits=19, decimal_places=3)

    def natural_key(self):
        return (self.type, self.title) + self.queue.natural_key()
        
class Proteome(models.Model):
    proteome = models.CharField(max_length=20, primary_key=True)
    organism = models.TextField()
    # number of unique proteins in proteome
    full_size = models.IntegerField(default=0)
    # number in the profile.fastaprotein
    # note that some proteomes will actually have 0
    profile_size = models.IntegerField(default=0)

# table of proteins loaded from the FASTA file
class FastaProtein(models.Model):
    accession = models.CharField(max_length=100, primary_key=True)
    description = models.TextField()
    # proteome
    # pan proteome (sometimes the same as upid, sometimes not)
    ppid = models.ForeignKey(
        'results.Proteome', 
        on_delete=models.CASCADE, 
        related_name='ppid_fasta'
    )
    length = models.IntegerField(default=1)
    gene = models.CharField(max_length=50)

class PsmRatio(models.Model):
    psm = models.ForeignKey('results.Psm', on_delete=models.CASCADE)
    # the number
    ratio = models.DecimalField(max_digits=25, decimal_places=15, null=True)
    # the sample ID, e.g. TMT-127C
    label = models.CharField(max_length=255)
 
    def natural_key(self):
        return(self.psm, self.label)
        
# to make it easier later, we will calculate these values and save them
# rather than calculate them on demand
# this lets us more easily combine table + filter + export + sort
# and we can use it both for the species page and the summary page
class SpeciesSummary(models.Model):
    project = models.ForeignKey('projects.Project', on_delete=models.CASCADE)
    type = models.CharField(max_length=20, choices=TYPE_CHOICES)
    ppid = models.ForeignKey(
        'results.Proteome', 
        on_delete=models.CASCADE, 
        related_name='ppid_species'
    )
    val_num_protein = models.IntegerField()
    val_num_psm = models.DecimalField(max_digits=15, decimal_places=3)
    val_num_peptide = models.IntegerField()
    nsaf = models.DecimalField(max_digits=15, decimal_places=6)
    peak_area = models.DecimalField(null=True, max_digits=19, decimal_places=3)
    peak_area_psm = models.IntegerField(null=True)

    def natural_key(self):
        return(self.project, self.ppid, self.type)
        
class SpeciesFileSummary(models.Model):
    queue = models.ForeignKey('projects.Queue', on_delete=models.CASCADE)
    type = models.CharField(max_length=20, choices=TYPE_CHOICES)
    ppid = models.ForeignKey(
        'results.Proteome', 
        on_delete=models.CASCADE, 
        related_name='ppid_species_file_summary'
    )
    val_num_protein = models.IntegerField()
    val_num_psm = models.DecimalField(max_digits=15, decimal_places=3)
    val_num_peptide = models.IntegerField()
    nsaf = models.DecimalField(max_digits=15, decimal_places=6)
    peak_area = models.DecimalField(null=True, max_digits=19, decimal_places=3)
    peak_area_psm = models.IntegerField(null=True)
 
    def natural_key(self):
        return(self.type, self.ppid) + self.queue.natural_key()

# just keep tabs on the list of files we are loading onto the website for a project
class ResultsFiles(models.Model):
    project = models.ForeignKey('projects.Project', on_delete=models.CASCADE)
    file = models.CharField(max_length=100)
    
    def natural_key(self):
        return(self.project, self.file)
        
class DiffProtein(models.Model):
    project = models.ForeignKey('projects.Project', on_delete=models.CASCADE)
    fp = models.ForeignKey('results.FastaProtein', on_delete=models.CASCADE)
    logfc = models.DecimalField(max_digits=15, decimal_places=10)
    # adj.P.Val (LIMMA)
    p_value = models.DecimalField(max_digits=15, decimal_places=10)
    # sca.adj.pval (DEqMS)
    d_p_value = models.DecimalField(max_digits=15, decimal_places=10)

    def natural_key(self):
        return (self.fp.accession, self.project)
        
    natural_key.dependencies = ['results.fastaprotein']