### unused at this point
### was used as a way to automatically load from a metadata file but these files
### are often in different formats 

import csv
import os
import argparse

from results.models import MultiplexLabel

from .run_command import write_debug, settings

def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('client_name', type=str)
    parser.add_argument('project_name', type=str)
    args2 = parser.parse_args(args)
    try:
        client = args2.client_name
        project = args2.project_name
    except:
        return
        
    load_multiplex_labels(client, project)
    
def load_multiplex_labels(client, project):
    print("Loading Multiplex Labels into database.")

    if  not os.path.exists(os.path.join(settings.data_folder, client, project, "%s_%s_multiplexlabels.tsv" % (client, project))):
        print("%s_%s_multiplexlabels.tsv does not exist. Generate first." % (client, project))
        return
    
    delete = MultiplexLabel.objects.all().delete()
    
    with open(os.path.join(settings.data_folder, client, project, "%s_%s_multiplexlabels.tsv" % (client, project)), 'r') as file:
        result = csv.reader(file, delimiter='\t')
        header = next(result)
        for row in result:
            try:
                multiplexlabel = MultiplexLabel(sample = row[header.index('filename')],
                                                  label = row[header.index('sample')].lower(),
                                                  phenotype = row[header.index('label')].lower(),
                                                  identifier = row[header.index('identifier')].lower()
                                                 )
                multiplexlabel.save()
            except Exception as e:
                print("Failure to load labels: %s" % e)
                return
