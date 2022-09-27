import subprocess
import time
import os

from django.core.exceptions import ObjectDoesNotExist

from projects.models import Setting

try:
    settings = Setting.objects.get(default=True)
except ObjectDoesNotExist:
    print("Create default settings for this server first.")
    raise SystemExit

# program outputs only
def write_log(msg, job, project):
    ts = time.localtime()
    msg = "%s %s" % (time.strftime("%Y-%m-%d %H:%M:%S", ts), msg)
    
    filename = "%s%s%s_%s.log" % (os.path.join(settings.install_folder, "log", project), os.sep, project, job)

    if not os.path.exists(os.path.join(settings.install_folder, "log", project)):
        os.makedirs(os.path.join(settings.install_folder, "log", project))
        
    f = open(filename, 'a')
    ts = time.localtime()
    f.write(msg)
    f.close()
    #print(msg, end='\r\n', flush=True)
    print(msg.strip("\n"), flush=True)

# debug mode messages
def write_debug(msg, job, project):
    # this turns off all output for now but is intended to do something else
    #   in the future
    if settings.debug_mode == '0':
        return

    ts = time.localtime()
    msg = "%s %s" % (time.strftime("%Y-%m-%d %H:%M:%S", ts), msg)
    
    filename = "%s%s%s_%s.log" % (os.path.join(settings.install_folder, "log", project), os.sep, project, job)

    if not os.path.exists(os.path.join(settings.install_folder, "log", project)):
        os.makedirs(os.path.join(settings.install_folder, "log", project))
        
    f = open(filename, 'a')
    ts = time.localtime()
    f.write("%s\r\n" % msg)
    f.close()
    #print(msg, end='\r\n', flush=True)
    print(msg, flush=True)

# writes an error message relating to entire project rathern specifically one queue entry
# this would focus on things that would break the script from working
def write_error(msg):
    ts = time.localtime()
    msg = "%s %s" % (time.strftime("%Y-%m-%d %H:%M:%S", ts), msg)
    
    filename = "%s" % (os.path.join(settings.install_folder, "log", "error.log"))
    
    f = open(filename, 'a')
    ts = time.localtime()
    f.write("%s\r\n" % msg)
    f.close()
    #print(msg, end='\r\n', flush=True)
    print(msg, flush=True)
    
# execute a command and log each individual line of output as it comes out    
def run_command(cmd, job, project):
    write_debug("runCommand: %s" % (cmd), job, project)

    with subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                            stderr=subprocess.STDOUT, universal_newlines=True, bufsize=1) as p:
        for line in p.stdout:
            write_log(line, job, project)
    if p.returncode != 0:
        print("returned error:" + str(p.returncode) + str(p.args))
        return 0
    else:
        return 1
