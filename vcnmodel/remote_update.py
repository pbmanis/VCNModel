#!/usr/bin/python
from __future__ import print_function
__author__ = 'pbmanis'
import string
import getpass

"""
remote_update.py is for updating the software on remote machines in the lab
for modeling
The program:
1. logs in through ssh
2. changes directories to the right directory (PyNeuronLibrary/Cortex_STDP in this case)
3. uploads the source files in the list sourceFiles
4. changes directories to Desktop/Python, copies: pylibrary, nrnlibrary,
4. submits the job using bsub and requesting an appropriate number of processors
4. (optionally) patiently waits for the results to become available.
5. (optionally) downloads the result file
if run in the non-interactive mode, 4 and 5 are skipped, and you must manually download
once the email arrives...
Note you must be on campus or using a VPN to access the machine this way
"""
import paramiko
import sys
print(len(sys.argv))
print(sys.argv[0])
if len(sys.argv) > 1:
    machine = sys.argv[1]
else:
    machine = 'Killdevil'

# system selection:
model = 'CalyxModel'

# list the source files that should be updated for a run
sourceFiles = ['analyze_run.py', 'channel_decorate.py', 'channel_manager.py', 'generate_run.py',
               'listDirs.py', 'h_reader.py', 'model_run.py', 'remap_hoc.py', 'parallel_model_run.py']

# define system directory information...
sysChoice = {'Killdevil': {'uname': 'pmanis', 'dir': 'PyNeuronLibrary/%s' % model, 'addr': 'killdevil.unc.edu'}, #'152.19.187.84'},
             'Lytle': {'uname': 'pbmanis', 'dir': '/Users/pbmanis/Desktop/Python/PyNeuronLibrary/%s' % model, 'addr': '152.19.86.111'},
             'Tule': {'uname': 'pbmanis', 'dir': '/Users/pbmanis/Desktop/Python/PyNeuronLibrary/%s' % model, 'addr': '152.19.86.116'},
            }

if machine not in sysChoice.keys():
    print('Machine %s not recognized' % machine)
    exit()

mypwd = getpass.getpass("Password for %s: " % machine)

ssh=paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
print('connecting to: ', sysChoice[machine]['addr'])
print('as user: ', sysChoice[machine]['uname'])
conn = ssh.connect(sysChoice[machine]['addr'], username=sysChoice[machine]['uname'], password=mypwd)  # hate putting pwd in file like this...
if conn is False:
    print('Connection failed')
    exit()

cdcmd = 'cd ' + sysChoice[machine]['dir']
# Note tht the exe command is a "single session", so all commands that will be done need to be concatenated
# using semicolons as if you were doing it from a command line on a terminal.
# cd will not "hold"
#stdin, stdout, stderr = ssh.exec_command('cd ~/PyNeuronLibrary/Cortex-STDP; ls -la')
#for l in stdout.readlines():
#    print l,

print('dir : ', sysChoice[machine]['dir'])
print('open ftp:')
ftp = ssh.open_sftp()
ftp.chdir(sysChoice[machine]['dir'])  # however this works for the sftp
print(ftp.getcwd())
for f in sourceFiles:
    ftp.put(f, f)  # update the source files
ftp.close()

# print 'execute ls'
# stdin, stdout, stderr = ssh.exec_command('cd ~/PyNeuronLibrary/Cortex-STDP ; ls -la')
# for l in stdout.readlines():
#      print l,

#

ssh.close()

