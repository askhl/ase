"""
This is conditional on the existence of the VASP_COMMAND or VASP_SCRIPT
environment variables.
"""

def vasp_installed():
    import os
    from ase.test import NotAvailable

    vcmd = os.getenv('VASP_COMMAND')
    vscr = os.getenv('VASP_SCRIPT')
    if vcmd == None and vscr == None:
        raise NotAvailable('Neither VASP_COMMAND nor VASP_SCRIPT defined')

if __name__ == '__main__':
    vasp_installed()
