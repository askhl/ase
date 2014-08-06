import os
import sys
import shutil

def smart_open(filename=None):
    if filename and filename != '-':
        if os.path.isfile(filename):
            shutil.move(filename, filename+'.bak')
        fh = open(filename, 'w')
    else:
        fh = sys.stdout
    return fh

class free_printer():
    def __init__(self, outp='-'):
	self.fil = smart_open(outp)

    def header(self):
        fil = self.fil
        header='\n\t\tFree Energy Calculator ver. S2014\n\n'
        fil.write(header)

    def sub_header(self, str):
        fil = self.fil
        header='\n\t\t##################################\n\t\tRunning '+str+'\n\t\t##################################\n\n'
        fil.write(header)
    
    def sub_footer(self, str):
        fil = self.fil
        footer='\n\t\t##################################\n\t\tLeaving '+str+'\n\t\t##################################\n\n'
        fil.write(footer)
    
    def f_print(self, str):
        fil = self.fil
        fil.write(str+'\n')
