"""
Unit tests for the chemicalelement module.
"""
# Import the unittest module.
import unittest
# Import numpy
import numpy

from ase import exceptions
from ase.chemicalelement import ChemicalElement

class ChemicalElementTest(unittest.TestCase):
    """ Test class for the ChemicalElement """
    
    def test_valid_input1(self):
        """ Test construction and query methods on the ChemicalElement"""
        # Create a valid element.            
        my_element = ChemicalElement(None, 'H', 'Hydrogen',6)

        # Check that the member data is correct.
        self.assertEqual(my_element.get_symbol(), 'H')
        self.assertEqual(my_element.get_name(), 'Hydrogen')
        self.assertEqual(my_element.get_atomic_number(), 6)
    
    def test_query_method_return_correct(self):
        """ Test that the query method returns the proper member data """
        # Create a valid element.            
        my_element = ChemicalElement(None, 'He', 'Helimum', 2)
        # Check that the member data and query methods return the same.
        self.assertEqual(my_element.symbol, my_element.get_symbol())
        self.assertEqual(my_element.name, my_element.get_name())
        self.assertEqual(my_element.atomic_number, 
                         my_element.get_atomic_number())
    def test_error_1(self):
        """ Test that the symbol must be string. """
        self.assertRaises(exceptions.ASETypeError,
            lambda: ChemicalElement(None,1,'Hydrogen',1))

    def test_error_2(self):
        """ Test that the name must be string. """
        self.assertRaises(exceptions.ASETypeError,
            lambda: ChemicalElement(None,'H',1,1))

    def test_error_3(self):
        """ Test that the atomic number must be an positive integer """
        self.assertRaises(exceptions.ASETypeError,
            lambda: ChemicalElement(None,'H','Hyd',1.0))
        self.assertRaises(exceptions.ASEValueError,
            lambda: ChemicalElement(None,'H','Hyd',-1,))
    
    def test_base_element_1(self):
        """ Test that base element transfer all values that are not given """
        # Create an chemical element.
        my_element = ChemicalElement(None,'H','Hydrogen',1)

        # Create another chemical element.
        another_element = ChemicalElement(my_element)

        # Check that the properties are correctly transfered.
        self.assertEqual(my_element.symbol,
                         another_element.symbol)
        self.assertEqual(my_element.name,
                         another_element.name)
        self.assertEqual(my_element.atomic_number,
                         another_element.atomic_number)

    def test_base_element_2(self):
        """ Test that base element tranfer values are correct replaced. """
        # Create an chemical element.
        my_element = ChemicalElement(None,'H','Hydrogen',1)

        # Check that the properties are correctly transfered.
        self.assertEqual('H2',
             ChemicalElement(my_element,symbol='H2').symbol)
        self.assertEqual('Hydrogen+',
             ChemicalElement(my_element,name='Hydrogen+').name)
        self.assertEqual(10,
             ChemicalElement(my_element,atomic_number=10).atomic_number)

    def test_base_element_error(self):
        """ Test that the base element must be a ChemicalElement """
        self.assertRaises(exceptions.ASETypeError,
            lambda: ChemicalElement('element','H','Hydrogen',1))

if __name__ == '__main__':
    unittest.main()
