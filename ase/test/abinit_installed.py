def abinit_installed():
    import os
    from ase.test import NotAvailable

    abinitcmd = os.getenv('ABINIT_SCRIPT')
    psppath = os.getenv('ABINIT_PP_PATH')
    if abinitcmd == None or psppath == None:
        raise NotAvailable('ABINIT_SCRIPT or ABINIT_PP_PATH not defined')

if __name__ == '__main__':

    abinit_installed()
