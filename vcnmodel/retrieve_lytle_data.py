#!/usr/bin/python

__author__ = "pbmanis"
import string
import getpass
from pathlib import Path
import toml
import subprocess

"""
remote_update.py is for pulling data from lytle for the vcnmodel results

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
    machine = "Lytle"

gradeACells = [2, 5, 6, 9, 10, 11, 13, 17, 24, 29, 30]


class Retiever(object):
    def __init__(self):
        self.setPaths()
        self.sysChoice = {  #'152.19.187.84'},
            "Lytle": {
                "uname": "pbmanis",
                "dir": self.datapaths,
                "addr": "152.19.86.111",
            },
            "Tule": {
                "uname": "pbmanis",
                "dir": self.datapaths,
                "addr": "152.19.86.116",
            },
        }

        print("Machine: ", machine)
        print(self.datapaths["baseDirectory"])
        if machine not in self.sysChoice.keys():
            print("Machine %s not recognized" % machine)
            exit()

        mypwd = getpass.getpass("Password for %s: " % machine)

        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        print("connecting to: ", self.sysChoice[machine]["addr"])
        print("as user: ", self.sysChoice[machine]["uname"])
        conn = self.ssh.connect(
            self.sysChoice[machine]["addr"],
            username=self.sysChoice[machine]["uname"],
            password=mypwd,
        )  # hate putting pwd in file like this...
        if conn is False:
            print("Connection failed")
            exit()
        self.conn = conn

    def setPaths(self, stimtype="AN", cell=11):
        where_is_data = Path("wheres_my_data.toml")
        if where_is_data.is_file():
            self.datapaths = toml.load("wheres_my_data.toml")
        else:
            self.datapaths = {"baseDirectory": Path("../VCN-SBEM-Data", "VCN_Cells",)}
        self.basepath = self.datapaths["baseDirectory"]

    # define system directory information...

    def getCellData(self):
        # gradeACells = [30]
        for ic in gradeACells:
            self.setPaths(cell=ic)
            dpath = str(Path(
                    self.sysChoice[machine]["dir"]["baseDirectory"],
                    f"VCN_c{ic:02d}",
                    "Simulations",
                    "AN",
                ))
            cdcmd = "cd " + dpath
            backupDir = dpath # important: end path slash!!!! 
            print('cdcmd: ', cdcmd)
            # print('backup dir: ', backupDir)

            # Note that the exe command is a "single session", so all commands that will be done need to be concatenated
            # using semicolons as if you were doing it from a command line on a terminal.
            # cd will not "hold"
            # stdin, stdout, stderr = self.ssh.exec_command(f"{cdcmd:s}; ls -la")
            # sourcefiles = []
            # for l in stdout.readlines():
            #     lsp = l.split(' ')
            #     print(lsp[-1])
            #     sourcefiles.append(lsp[-1])
            #
            bwLimit = ''
            excludeStr = ''
            trash = ''
            dryrun = '' # '--dry-run'
            print('dpath:     ', dpath)
            backupDir = str(Path(backupDir).parent) + '/'
            print('backupDir: ', backupDir)
            #continue
            cmd = f"rsync -rtDuv {dryrun:s} --chmod=Du=rwx,Dg=rwx,Do=,Fu=rw,Fg=rw,Fo= {bwLimit:s}"
            cmd += f" --itemize-changes --out-format='%n' -e ssh pbmanis@152.19.86.111:{dpath:s} {backupDir:s}"
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            logstr = ''
            lo = ''
            while True:
                l = p.stdout.read(1)
                if l == b'':
                    break
                lo += l.decode("utf-8")
                # log.write(l)
                # log.flush()
                # logstr += l
            ret = p.wait()
            print(lo)
            print(ret)
            
            # print('open ftp:')
#             ftp = self.ssh.open_sftp()
#             ftp.chdir(dpath)  # however this works for the sftp
#             print(ftp.getcwd())
#             for f in sourceFiles:
#                 ftp.get(f, f)  # update the source files
#             ftp.close()
#
#             print ('execute ls after ftp')
#             stdin, stdout, stderr = ssh.exec_command(f"{cdcmd:s}; ls -la")
#             for l in stdout.readlines():
#                  print(l)



        self.ssh.close()


if __name__ == "__main__":
    R = Retiever()
    R.getCellData()
