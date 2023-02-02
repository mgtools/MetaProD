from django.urls import path
from django.urls import re_path as url
from dal import autocomplete

from website import views
from django.contrib import admin
from django.conf.urls.static import static
from django.conf import settings
from projects.models import MultiplexLabel, Sample, Label, Queue, Tag, MetaData

class QueueSample(autocomplete.Select2QuerySetView):
    def get_queryset(self):
        qs = Sample.objects.none()
        project = self.forwarded.get('project', None)
        if project:
            qs = Sample.objects.filter(project=project)
            #qs = qs.filter(project=project)

        return qs

class QueueTag(autocomplete.Select2QuerySetView):
    def get_queryset(self):
        qt = Tag.objects.none()
        project = self.forwarded.get('project', None)
        if project:
            qt = Tag.objects.filter(project=project)

        return qt

class QueueMetaData(autocomplete.Select2QuerySetView):
    def get_queryset(self):
        qt = MetaData.objects.none()
        project = self.forwarded.get('project', None)
        if project:
            qt = MetaData.objects.filter(project=project)

        return qt
        
class MultiplexLabelSample(autocomplete.Select2QuerySetView):
    def get_queryset(self):
        qs = Sample.objects.none()
        #qs = super(LinkedDataView, self).get_queryset()
        project = self.forwarded.get('project', None)

        if project:
            qs = Sample.objects.filter(project=project)
            #qs = qs.filter(project=project)

        return qs

class MultiplexLabelReagent(autocomplete.Select2QuerySetView):
  
    def get_queryset(self):
        qs = Label.objects.none()
        reagent = self.forwarded.get('reagent', None)
        l_0 = self.forwarded.get('labelchoice_set-0-label', None)
        l_1 = self.forwarded.get('labelchoice_set-1-label', None)
        l_2 = self.forwarded.get('labelchoice_set-2-label', None)
        l_3 = self.forwarded.get('labelchoice_set-3-label', None)
        l_4 = self.forwarded.get('labelchoice_set-4-label', None)
        l_5 = self.forwarded.get('labelchoice_set-5-label', None)
        l_6 = self.forwarded.get('labelchoice_set-6-label', None)
        l_7 = self.forwarded.get('labelchoice_set-7-label', None)
        l_8 = self.forwarded.get('labelchoice_set-8-label', None)
        l_9 = self.forwarded.get('labelchoice_set-9-label', None)
            
        #print(self.forwarded.get())
        if reagent:
            qs = (Label.objects.filter(reagent=reagent)
                              .exclude(name=l_0)
                              .exclude(name=l_1)
                              .exclude(name=l_2)
                              .exclude(name=l_3)
                              .exclude(name=l_4)
                              .exclude(name=l_5)
                              .exclude(name=l_6)
                              .exclude(name=l_7)
                              .exclude(name=l_8)
                              .exclude(name=l_9))

        return qs

class InitialStateAutocomplete(autocomplete.Select2QuerySetView):
    def get_queryset(self):
        print(self.forwarded.get())
        
urlpatterns = [
    path('', views.index, name='index'),
    #path('projects/', views.ProjectListView.as_view(), name='project_list'),
    path('projects/', views.project_list, name='project_list'),
    path('protein/', views.ProteinListView.as_view(), name='protein_list'),
    path('peptide/', views.PeptideListView.as_view(), name='peptide_list'),
    path('psm/', views.PsmListView.as_view(), name='psm_list'),
    path('diffprotein/', views.DiffProteinListView.as_view(), name='diffprotein_list'),
    path('files/', views.FileListView.as_view(), name='file_list'),
    path('summary/<str:project>/', views.summary, name='summary'),
    path('file_summary/', views.FileSummaryView.as_view(), name='file_summary'),
    path('species/', views.SpeciesListView.as_view(), name='species_list'),
    path('help', views.help, name='help'),
    path('admin/', admin.site.urls, name="admin"),
    url(
        '^multiplexlabelsample/$',
        MultiplexLabelSample.as_view(model=MultiplexLabel),
        name='multiplexlabelsample'
    ),
    url(
        '^multiplexlabelreagent/$',
        MultiplexLabelReagent.as_view(model=MultiplexLabel),
        name='multiplexlabelreagent'
    ),
    url(
        '^initial-state-autocomplete/$',
        InitialStateAutocomplete.as_view(model=MultiplexLabel),
        name='initial_state_autocomplete'
    ),
    url(
        '^queuesample/$',
        QueueSample.as_view(model=Queue),
        name='queuesample'
    ),
    url(
        '^queuetag/$',
        QueueTag.as_view(model=Queue),
        name='queuetag'
    ),
    url(
        '^metadata/$',
        QueueMetaData.as_view(model=Queue),
        name='metadata'
    )
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
