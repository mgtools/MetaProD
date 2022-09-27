import shutil
import os
import argparse
import time
import zipfile

from projects.models import Project, SearchSetting

from .run_command import run_command, write_debug, settings, write_error

def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('project_name', type=str)
    parser.add_argument('description', type=str, nargs='?', default="")
    args2 = parser.parse_args(args)
    
    project_name = args2.project_name
    description = args2.description
        
    create_project(project_name, description)

def create_project(project_name, description):

    data_folder = settings.data_folder
    
    if Project.objects.filter(name=project_name).exists() is True:
        print("Project %s already exists." % project_name)
        return
    else:
        print("Creating project %s." % project_name)
        project = Project(name=project_name, description=description)
        project.save()
        
        print("Creating SearchSetting for project %s with default values." % project_name)
        print("Change settings and add mods and enzymes in the web interface.")
        searchsetting = SearchSetting(project=project)
        searchsetting.save()
        
        print("Finished creating the project.")