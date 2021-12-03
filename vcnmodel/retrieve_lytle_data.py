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
2. changes directories to the right directory
3. uploads the source files in the list sourceFiles
4. submits the job using bsub and requesting an appropriate number of processors
4. (optionally) patiently waits for the results to become available.
5. (optionally) downloads the result file
if run in the non-interactive mode, 4 and 5 are skipped, and you must manually download
once the email arrives...
Note you must be on campus or using a VPN to access the machine this way
"""
from dataclasses import dataclass
from typing import Union
import paramiko
import sys

print(len(sys.argv))
print(sys.argv[0])
if len(sys.argv) > 1:
    machine = sys.argv[1]
else:
    machine = "Lytle"



@dataclass
class SysInfo:
    uname: str='pbmanis'
    directory: Union[str, Path, None]=None
    addr: str=''


class Syncer(object):
    def __init__(self, remotemachine, push_pull="pull"):
        self.push_pull = push_pull
        self.machine = remotemachine
        self.gradeACells = [2, 5, 6, 9, 10, 11, 13, 17, 30]
        self.set_paths()
        if self.machine == 'Lytle':
            self.remote =SysInfo(directory=self.datapaths, addr="192.168.1.195")
        elif self.machine == 'Tamalpais2':
            self.remote = SysInfo(directory=self.datapaths, addr="192.168.1.218")
        else:
            print("Machine {self.machine:s} is not recognized")
            exit()
            
        print("Machine: ", self.machine)
        print(self.datapaths["baseDirectory"])

        # set password to empty string if have passwordless ssh setup already
        mypwd = "" # getpass.getpass("Password for %s: " % self.machine)

        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        print("connecting to: ", self.remote.addr)
        print("as user: ", self.remote.uname)
        conn = self.ssh.connect(
            self.remote.addr,
            username=self.remote.uname,
            password= mypwd,
        )  # hate putting pwd in file like this...
        if conn is False:
            print("Connection failed")
            exit()
        self.conn = conn
    
    def done(self):
        self.ssh.close()

    def set_paths(self, stimtype="AN", cell=11):
        where_is_data = Path("wheres_my_data.toml")
        if where_is_data.is_file():
            self.datapaths = toml.load("wheres_my_data.toml")
        else:
            self.datapaths = {"baseDirectory": Path("../VCN-SBEM-Data", "VCN_Cells",)}
        self.basepath = self.datapaths["baseDirectory"]

    # define system directory information...
    def get_all_cell_data(self):
        for ic in self.gradeACells:
            self.get_cell_data(self, ic)
    
    def transfer_cell_data(self, ic:int, push_pull="push", dryRun:bool=True, daysback:int=3):
    # gradeACells = [30]
        self.set_paths(cell=ic)
        remote_dir = str(Path(
                self.remote.directory["baseDirectory"],
                f"VCN_c{ic:02d}",
                "Simulations",
                "AN",
            ))
        local_dir = str(Path(remote_dir).parent) + '/'
        cdcmd = "cd " + remote_dir
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
        if dryRun:
            dryrun = '--dry-run'
        else:
            dryrun = ''
        datelimit = f"--files-from=<(find {remote_dir:s} -mtime -{daysback:d}"
        print('Remote_dir:     ', remote_dir)
        datelimit +=  "-type f -exec basename {} \;)"

        print('Local_dir: ', local_dir)
        # return
        # note this uses the brew installed rsync, V3.2.3 as of 8/28/2020. The "system" version is 2.6.9
        cmd = f"/usr/local/bin/rsync -rtDuv {dryrun:s} --chmod=Du=rwx,Dg=rwx,Do=,Fu=rw,Fg=rw,Fo= {bwLimit:s}"
        if self.push_pull == "pull":
            cmd += f" --itemize-changes --out-format='%n' -e ssh pbmanis@{self.remote.addr:s}:{remote_dir:s} {local_dir:s}"
        elif self.push_pull == "push":
            cmd += f" --itemize-changes --out-format='%n' -e ssh {local_dir:s} pbmanis@{self.remote.addr:s}:{remote_dir:s}"
        else:
            raise ValueError('push-pull must be either push or pull')
        
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        logstr = ''
        lo = ''
        lo_line = ''
        while True:
            l = p.stdout.read(1)
            if l == b'':
                break
            lo_line += l.decode("utf-8")
            lo += lo_line
            if lo_line[-1] == '\n':
                print(lo_line,)
                lo_line = ''
            # log.write(l)
            # log.flush()
            # logstr += l
        ret = p.wait()
        # print(lo)
        # print(ret)


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






if __name__ == "__main__":
    R = Syncer(remotemachine = machine, push_pull='push')
    R.transfer_cell_data(5, dryRun=True)
    R.done()