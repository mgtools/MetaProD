from dal import autocomplete
from dal import forward

from django import forms

from .models import MultiplexLabel, LabelChoice, Queue

class MultiplexLabelForm(forms.ModelForm):

    class Meta:
        model = MultiplexLabel
        fields = ('project', 'sample', 'reagent')
        widgets = {
            'sample': autocomplete.ModelSelect2(url='multiplexlabelsample',
                                              forward=('project','reagent'))
        }
    class Media:
        js = (
            'linked_data.js',
        )
        
class MultiplexLabelInlineForm(forms.ModelForm):

    class Meta:
        model = LabelChoice
        fields = ('__all__')
        forwards = ['labelchoice_set-'+str(i)+'-label' for i in range(10)]
        forwards.append('reagent')
        widgets = {
            'label': autocomplete.ModelSelect2(url='multiplexlabelreagent',
                                              forward=forwards),
        }        
        
    class Media:
        js = (
            'linked_data.js',
        )

class QueueForm(forms.ModelForm):

    class Meta:
        model = Queue
        fields = ('__all__')
        widgets = {
            'sample': autocomplete.ModelSelect2(url='queuesample',
                                                forward=('project',)),
            'tag': autocomplete.ModelSelect2(url='queuetag',
                                             forward=('project',))
        }
    class Media:
        js = (
            'linked_data.js',
        )