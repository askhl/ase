"""
This modules defines the generic chemical element
represented by the object ChemicalElement
"""
# Import the needed common check functins.
from ase.commonchecks import check_string, check_positive_integer, \
                             check_positive_float, check_rgb_tuple
from ase import exceptions                             

class ChemicalElement:
    """
    Class for representing a chemical element.
    In order to construct a chemical element the following parameters are 
    required:

    element: ChemicalElement
         Given if the new ChemicalElement should retain all the
         properties of the given element, except those additional
         specified.
    symbol: string
         The symbol to be used for this element
    name: string
        The name of the element given as a string.
    atomic_number: int
        The atomic number -  It must be given as positive integer.

    Example of usage:
    Iron = ChemicalElement('Fe', 'Iron', 26)

    Another example of usage:
    Carbon = ChemicalElement(symbol='C',
                             name='Carbon',
                             atomic_number=6)
    """
    def __init__(self,
                 element=None,
                 symbol=None,
                 name=None,
                 atomic_number=None):
        # Check if the element is given, otherwise the user is forced to
        # provided all the parameters.
        check_valid_input(element, symbol, name, atomic_number)
        # After this point we know, that if any parameters is none,
        # the user has given a base element.
    
        # Check if we should use the symbol from the base element.
        if symbol is None:
            self._symbol = element.get_symbol()
        else:            
            # Check that the symbol is a valid string, and store if it is.
            self._symbol = check_string(symbol,'symbol')
           
        # Check if we should use the name from the base element.
        if name is None:
            self._name = element.get_name()
        else:
            # Check that the name is a valid string, and store if it is.
            self._name = check_string(name,'name')
       
        # Check if the we should use the atomic number from the base element.
        if atomic_number is None:
            self._atomic_number = element.get_atomic_number()
        else:
            # Check that the atomic number is a positive integer, and
            # if it is stored.
            self._atomic_number = check_positive_integer(atomic_number,
                                                          'atomic_number')
    
    def get_symbol(self):
        """ Query method for getting the symbol of this element. """
        return self._symbol
    symbol = property(get_symbol,
                 doc='The symbol used to represent this ChemicalElement.')

    def get_name(self):
        """ Query method for getting the name of this element. """
        return self._name
    name = property(get_name,
                    doc='The convential name of this ChemicalElement.')

    def get_atomic_number(self):
        """ Query method for getting the atomic number of this element. """
        return self._atomic_number
    atomic_number = property(
            get_atomic_number,
            doc='The atomic number of this ChemicalElement.')  

def check_valid_input(element, symbol, name, atomic_number):
    """    
    Small private helper function for checking that user has
    given the right combination of the parameters.
    """
    # Check if the element is given, otherwise the user is forced to
    # provided all the parameters.
    if element is None:
        if symbol is None:
            raise exceptions.ASEValueError('When no base element is'\
               ' given to ChemicalElement, the symbol must be specified.')
        if name is None:
            raise exceptions.ASEValueError('When no base element is'\
               ' given to ChemicalElement, the name must be specified.')
        if atomic_number is None:
            raise exceptions.ASEValueError('When no base element is'\
               ' given to ChemicalElement, the atomic_number'\
               ' must be specified.')
    # Since the element is given, check it is a ChemicalElement.
    else:
        if not(isinstance(element,ChemicalElement)):
            raise exceptions.ASETypeError('The base element must be '\
                    'a ChemicalElement.')

