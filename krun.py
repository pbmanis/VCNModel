#!/usr/bin/python
__author__ = 'pbmanis'
import string
import getpass
import stdp_test_plot
"""
krun.py is for submitting jobs to killdevil.unc.edu
The program:
1. logs in through ssh
2. changes directories to the right directory (PyNeuronLibrary/Cortex_STDP in this case)
3. uploads the source file (for example, stdp_test_parallel.py)
4. submits the job using bsub and requesting an appropriate number of processors
4. patiently waits for the results to become available.
5. downloads the result file

Note you must be on campus or using a VPN to access the machine this way
"""
import paramiko
import sys
print len(sys.argv)
print sys.argv[0]
if len(sys.argv) > 1:
    machine = sys.argv[1]
else:
    machine = 'Killdevil'

# system selection:
model = 'CalyxModel'
sourceFiles = ['analyze_run.py', 'channel_decorate.py', 'channel_manager.py', 'generate_run.py',
               'listDirs.py', 'model_run.py', 'remap_hoc.py', ]
sysChoice = {'Killdevil': {'uname': 'pmanis', 'dir': 'PyNeuronLibrary/%s' % model, 'addr': 'killdevil.unc.edu'}, #'152.19.187.84'},
             'Lytle': {'uname': 'pbmanis', 'dir': '/Users/pbmanis/Desktop/Python/PyNeuronLibrary/%s', model, 'addr': '152.19.86.111'},
             'Tule': {'uname': 'pbmanis', 'dir': '/Users/pbmanis/Desktop/Python/PyNeuronLibrary/%s', model, 'addr': '152.19.86.116'},
            }

if machine not in sysChoice.keys():
    print 'Machine %s not recognized' % machine
    exit()

mypwd = getpass.getpass("Password for %s: " % machine)

ssh=paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
print 'connecting to: ', sysChoice[machine]['addr']
print 'as user: ', sysChoice[machine]['uname']
conn = ssh.connect(sysChoice[machine]['addr'], username=sysChoice[machine]['uname'], password=mypwd)  # hate putting pwd in file like this...
if conn is False:
    print 'Connection failed'
    exit()

cdcmd = 'cd ' + sysChoice[machine]['dir']
# Note tht the exe command is a "single session", so all commands that will be done need to be concatenated
# using semicolons as if you were doing it from a command line on a terminal.
# cd will not "hold"
#stdin, stdout, stderr = ssh.exec_command('cd ~/PyNeuronLibrary/Cortex-STDP; ls -la')
#for l in stdout.readlines():
#    print l,

print 'dir : ', sysChoice[machine]['dir']
print 'open ftp:'
ftp = ssh.open_sftp()
ftp.chdir(sysChoice[machine]['dir'])  # however this works for the sftp
print ftp.getcwd()
for f in sourceFiles:
    ftp.put(f, f)  # update the source files
ftp.close()

# print 'execute ls'
# stdin, stdout, stderr = ssh.exec_command('cd ~/PyNeuronLibrary/Cortex-STDP ; ls -la')
# for l in stdout.readlines():
#      print l,
run = True
interactive = False
if run is True:
    if machine == 'Killdevil':
        if interactive:
            print 'run the bsub... interactive'
            stdin, stdout, stderr = ssh.exec_command(cdcmd + ' ; bsub -Ip python parallel_model_run.py')
            for l in stdout.readlines():
                print l,
        else:
            print 'run bsub, non-interactive'
            stdin, stdout, stderr = ssh.exec_command(cdcmd + '  ; bsub python parallel_model_run.py')
    else:
        print 'run standard job'
        stdin, stdout, stderr = ssh.exec_command(cdcmd + ' ;  python stdp_test_parallel.py')
        for l in stdout.readlines():
            print l,

if interactive:
    ftp = ssh.open_sftp()
    ftp.chdir(sysChoice[machine]['dir'])
    ftp.get('parallel.p', 'parallel.p')
    ftp.get('parallel_analyzed.p', 'parallel_analyzed.p')
    ftp.close()

ssh.close()

if interactive:
    stdp_test_plot.plot()