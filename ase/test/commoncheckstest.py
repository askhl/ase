"""
Unit tests for the commonchecks module
"""
# Import the unittest module.
import unittest
# Import numpy
import numpy
# Import the needed functions from commonchecks.
from ase.commonchecks import is_sequence, check_string, \
                             check_positive_integer, check_positive_float, \
                             check_rgb_tuple
from ase import exceptions

class CommonChecksTest(unittest.TestCase):
    """ Test class for the commonchecks """

    def test_is_sequence1(self):
        """ Test that is_sequence return true for various sequences. """
        self.assertTrue( is_sequence(numpy.zeros(3)) )
        self.assertTrue( is_sequence((1,2,3,4)) )
        self.assertTrue( is_sequence([1.2,0.3,2.0]) )               
        self.assertTrue( is_sequence('test'))

    def test_is_sequence2(self):
        """ Test that is_sequence return false for a string. """
        self.assertFalse( is_sequence('array') )

    def test_is_sequence3(self):
        """ Test that is_sequence return false for a single value. """
        self.assertFalse( is_sequence(1.0) )

    def test_check_string1(self):
        "Test that check_string raises an error,if it is not given a string. "
        self.assertRaises(exceptions.ASETypeError,
                          lambda: check_string(1.0,'test') )

    def test_check_string2(self):
        """ Test that check_string return the correct string, if valid """
        self.assertEqual( 'somestring', check_string('somestring','test'))

    def test_check_positive_integer1(self):
        "Test that check_positive_integer return the correct integer, if ok"
        self.assertEqual( 3, check_positive_integer(3,'test') )

    def test_check_positive_integer2(self):
        "Test that check_positive_integer will not accept a negative integer"
        self.assertRaises(exceptions.ASEValueError,
                          lambda: check_positive_integer(-5,'test') )

    def test_check_positive_integer3(self):
        "Test that check_positive_integer will not accept a float value"
        self.assertRaises(exceptions.ASETypeError,
                          lambda: check_positive_integer(1.0,'test') )

    def test_check_positive_float1(self):        
        " Test that check_positive_float will accept both float and integer"
        self.assertEqual( 2.7, check_positive_float(2.7,'test') )
        self.assertEqual( 3, check_positive_float(3,'test') )
    
    def test_check_positive_float2(self):
        """ Test that check_positive_float will not accept not numbers """
        self.assertRaises(exceptions.ASETypeError,
                          lambda: check_positive_float('string','test') )

    def test_check_positive_float3(self):
        " Test that check_positive_float raises error if number is negative "
        self.assertRaises(exceptions.ASEValueError,
                          lambda: check_positive_float(-1.0,'test') )
        # Try with an integer - just to be sure.
        self.assertRaises(exceptions.ASEValueError,
                          lambda: check_positive_float(0,'test') )
    
    def test_check_rgb_tuple1(self):
        " Test that rgb tuple return the correct tuple if given correct input "
        self.assertEqual( (0.2,0.3,0.4), check_rgb_tuple((0.2,0.3,0.4), 'test') )

    def test_check_rgb_tuple2(self):
        " Test that check_rgb_tuple raises an error, if given a non-sequence "
        self.assertRaises(exceptions.ASETypeError,
                          lambda: check_rgb_tuple(1.0,'test') )

    def test_check_rgb_tuple3(self):
        " Test that the check_rgb_tuple only accept sequence with length 3 "
        self.assertRaises(exceptions.ASEValueError,
                          lambda : check_rgb_tuple((1.0,1.0),'test'))
        self.assertRaises(exceptions.ASEValueError,
                          lambda : check_rgb_tuple((1.0,1.0,1.0,1.0),'test'))
    
    def test_check_rgb_tuple4(self):
        " Test that all values in rgb_tuple must be positive integer values "
        self.assertRaises(exceptions.ASEValueError,
                          lambda : check_rgb_tuple((-1.0,1.0,1.0),'test') )
        self.assertRaises(exceptions.ASEValueError,
                          lambda : check_rgb_tuple((1.0,-1.0,1.0),'test') )
        self.assertRaises(exceptions.ASEValueError,
                          lambda : check_rgb_tuple((1.0,1.0,-1.0),'test') )

    def test_check_rgb_tuple5(self):
        " Test that all values in the rgb_tuple must be below 1. "
        self.assertRaises(exceptions.ASEValueError,
                          lambda : check_rgb_tuple((1.1,1.0,1.0),'test') )
        self.assertRaises(exceptions.ASEValueError,
                          lambda : check_rgb_tuple((1.0,1.1,1.0),'test') )
        self.assertRaises(exceptions.ASEValueError,
                          lambda : check_rgb_tuple((1.0,1.0,1.1),'test') )

if __name__ == '__main__':
    unittest.main()

