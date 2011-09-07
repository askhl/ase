def elk_installed():
    import os
    from ase.test import NotAvailable

    elkcmd = os.getenv('ELK')
    speciespath = os.getenv('ELK_SPECIES_PATH')
    if elkcmd == None or speciespath == None:
        raise NotAvailable('ELK or ELK_SPECIES_PATH not defined')

if __name__ == '__main__':
    elk_installed()
