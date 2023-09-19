from dal import autocomplete
from dal import forward

from django import forms

from .models import MultiplexLabel, LabelChoice, Queue, MetaDataChoice

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
        # this is a workaround for a limitation
        # 20 needs to be increased if we have more labels
        forwards = ['labelchoice_set-'+str(i)+'-label' for i in range(20)]
        forwards.append('reagent')
        forwards.append('project')
        widgets = {
            'label': autocomplete.ModelSelect2(url='multiplexlabelreagent',
                                              forward=forwards),
            'tag': autocomplete.ModelSelect2(url='queuetag',
                                             forward=forwards)
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
        
class MetaDataChoiceForm(forms.ModelForm):

    class Meta:
        model = MetaDataChoice
        fields = ('__all__')
        widgets = {
            'metadata': autocomplete.ModelSelect2(url='metadata',
                                                  forward=('project',)),
        }

    class Media:
        js = (
            'linked_data.js',
        )