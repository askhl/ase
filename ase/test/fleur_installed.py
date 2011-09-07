def fleur_installed():
    import os
    from ase.test import NotAvailable

    inpgencmd = os.getenv('FLEUR_INPGEN')
    fleurcmd = os.getenv('FLEUR')
    if inpgencmd == None or fleurcmd == None:
        raise NotAvailable('FLEUR_INPGEN or FLEUR not defined')

if __name__ == '__main__':
    elk_installed()
