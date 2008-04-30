# -*- coding: utf-8 -*-
import os
from os.path import join
from stat import ST_MTIME
from docutils import nodes, utils
from docutils.parsers.rst.roles import set_classes

def svn_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
    if text[-1] == '>':
        i = text.index('<')
        name = text[:i - 1]
        text = text[i + 1:-1]
    else:
        name = text
    ref = 'http://trac.fysik.dtu.dk/projects/ase/browser/trunk/' + text
    set_classes(options)
    node = nodes.reference(rawtext, name, refuri=ref,
                           **options)
    return [node], []

def epydoc_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
    name = None
    if text[-1] == '>':
        i = text.index('<')
        name = text[:i - 1]
        text = text[i + 1:-1]
        
    components = text.split('.')
    if components[0] != 'ase':
        components.insert(0, 'ase')

    if name is None:
        name = components[-1]
        
    try:
        for n in range(2, len(components) + 1):
            module = __import__('.'.join(components[:n]))
    except ImportError:
        for component in components[1:n]:
            module = getattr(module, component)
            ref = '.'.join(components[:n])
            if isinstance(module, type):
                ref += '-class.html'
            else:
                ref += '-module.html'
        if n < len(components):
            ref += '#' + components[-1]
    else:
        ref = '.'.join(components) + '-module.html'

    ref = 'http://web2.fysik.dtu.dk/ase/epydoc/' + ref
    set_classes(options)
    node = nodes.reference(rawtext, name,
                           refuri=ref,
                           **options)
    return [node], []

def setup(app):
    app.add_role('svn', svn_role)
    app.add_role('epydoc', epydoc_role)
    #import atexit
    #atexit.register(fix_sidebar)
    create_png_files()

def create_png_files():
    for dirpath, dirnames, filenames in os.walk('.'):
        for filename in filenames:
            if filename.endswith('.py'):
                path = join(dirpath, filename)
                line = open(path).readline()
                if line.startswith('# creates:'):
                    t0 = os.stat(path)[ST_MTIME]
                    run = False
                    for file in line.split()[2:]:
                        try:
                            t = os.stat(join(dirpath, file))[ST_MTIME]
                        except OSError:
                            run = True
                            break
                        else:
                            if t < t0:
                                run = True
                                break
                    if run:
                        print 'running:', join(dirpath, filename)
                        os.system('cd %s; python %s' % (dirpath, filename))
        if '.svn' in dirnames:
            dirnames.remove('.svn')
            
"""
def fix_sidebar():
    print 'fixing sidebars...',
    for docname in ['index', 'contents']:
        t = open('.build/%s.html' % docname).read()
        i1 = t.find('<p class="first sidebar-title">')
        if i1 != -1:
            print docname,
            i0 = t.find('<div class="sidebar">')
            i1b = t.find('</p>', i1)
            name = t[i1 + 31:i1b]
            i2 = t.index('</div>', i1)
            i3 = t.index('<div class="sidebarwrapper">', i2)
            i4 = t.index('<h3>This Page</h3>', i3)
            t = (t[:i0] + t[i2 + 6:i3] + '<div class="sidebarwrapper">' +
                 ('<h4>%s</h4>' % name) + 
                 t[i1b + 5:i2] + t[i4:])
            open('.build/%s.html' % docname, 'w').write(t)
    print
"""
    
