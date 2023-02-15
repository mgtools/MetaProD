# should sort the order of the display by the queue order
import argparse
from projects.models import Queue


def run(*args):
    parser = argparse.ArgumentParser()
    parser.add_argument('project', type=str)
    args2 = parser.parse_args(args)
    project = args2.project
    
    check_status(project)

def check_status(project):
    queue = Queue.objects.filter(project__name=project)
    statuses = {
        'File added': 0,
        'Ready for thermorawfileparser': 0,
        'Ready for SearchGUI profile': 0,
        'Ready for PeptideShaker profile': 0,
        'Ready for Reporter profile': 0,
        'Ready for read_results profile': 0,
        'Ready for process_results profile': 0,
        'File is finished profile step': 0,
        'Ready for SearchGUI proteome': 0,
        'Ready for PeptideShaker proteome': 0,
        'Ready for Reporter proteome': 0,
        'Ready for read_results proteome': 0,
        'Ready for process_results proteome': 0,
        'File is finished proteome step': 0,
        'File is finished and cleaned up': 0,
        'Ready for MZmine profile': 0,
        'Ready for MZmine proteome': 0,
        'Skipped': 0
    }
    errors = 0
    jobs = []
    for q in queue:
        if q.skip == True:
            statuses['Skipped'] += 1
        else:
            statuses[q.get_status_display()] += 1
            if q.error == 2:
                errors += 1
        if (q.get_status_display() != 'File is finished proteome step'
            and q.get_status_display() != 'File is finished and cleaned up'
            and q.skip == False):
            if q.job not in jobs:
                jobs.append(q.job)
    print("Queue statuses for %s:" % project)
    for s in statuses:
        if statuses[s] > 0:
            print("%s: %s" % (s, statuses[s]))
    
    print("Queue entries with errors: %s" % errors)
    print("Job IDs remaining to be processed: %s " % (','.join(str(j) for j in sorted(jobs))))