"""
Module for making the most basic check of the user input
"""
# Import the ase exceptions.
from ase import exceptions

def is_sequence(test_object):
    """
    Function for checking that the test object is a sequence
    @param test_object: The object that should be checked if it is a sequence.
    @return Boolean: True if the object is a sequence.
    """
    return '__getitem__' in dir(test_object) and not isinstance(test_object,str)

def check_string(test_object,variable_name):
    """
    Method for checking if the given test_object is a valid string.
    If the given test_object is a valid string, the function will return the 
    string, otherwise it will raise an error.
    @param test_object : The object to be tested.
    @param variable_name : The name of the variable that is being tested, 
                           helpfull for generating a good error message.
    """                           
    # Check if the given object is a string, if not raise an error.
    if not(isinstance(test_object,str)):
        raise exceptions.ASETypeError, \
                    "The parameter, %s, must be a string."%(variable_name)
    # Since the test_object is a string, return the string.        
    else:       
        return test_object

def check_positive_integer(test_object,variable_name):
    """
    Method for checking if the given test_object is a valid integer,
    and positive.   If the given test_object is such an integer,
    the function will return the integer, otherwise it will raise an error.
    @param test_object : The object to be tested.
    @param variable_name : The name of the variable that is being tested, 
                           helpfull for generating a good error message.
    """                           
    # Check if the given object is an integer, if not raise an error.
    if not(isinstance(test_object,int)):
        raise exceptions.ASETypeError, \
                    "The parameter, %s, must be a positive integer."\
                    %(variable_name)
    # Check if it is an integer, that is positive.
    elif not(test_object > 0):
        raise exceptions.ASEValueError, \
                  "The parameter, %s, must be a positive integer, but is %s." \
                  %(variable_name,test_object)
    # Since the test_object is a a positive integer, return it.
    else:       
        return test_object

def check_positive_float(test_object,variable_name):
    """
    Method for checking if the given test_object is a valid float,
    and positive. It will also accept an integer as valid float.
    If it is invalid, it will raise an error otherwise return the float.
    @param test_object   : The object to be tested.
    @param variable_name : The name of the variable that is being tested, 
                           helpfull for generating a good error message.
    """                           
    # Check if the given object is an integer or float, if not raise an error.
    if not(isinstance(test_object,int) or isinstance(test_object,float)):
        raise exceptions.ASETypeError, \
                    "The parameter, %s, must be a positive float."\
                    %(variable_name)
    # Check if it is an integer or float, that is positive.
    elif not(test_object > 0):
        raise exceptions.ASEValueError, \
                  "The parameter, %s, must be a positive float, but is %s." \
                  %(variable_name,test_object)
    # Since the test_object is a positive float, return it.
    else:       
        return test_object


def check_rgb_tuple(test_object, variable_name):
    """
    Function for checking a RGB-tuple.
    The function will check that the tuple has the correct length,
    and consist of the elements with a float value between 0
    and 1. If this is the case it will return the color tuple,
    otherwise it will raise an error.
    @param test_object   : The object to be tested.
    @param variable_name : The name of the variable that is being tested, 
                           helpfull for generating a good error message.
    """
    # First check it is sequence.
    if not(is_sequence(test_object)):
        raise exceptions.ASETypeError, \
                    "The parameter, %s, must be a sequence of 3 floats " \
                    "with a value between 0 and 1."%(variable_name)

    # Check that the lenght is correct.
    elif not(len(test_object)==3):
        raise exceptions.ASEValueError, \
                    "The parameter, %s, must be a sequence of 3 floats " \
                    "with a value between 0 and 1.\n"%(variable_name) + \
                    " "*30 + "It has %s elements."%(len(test_object))
    
    # Check that the first element is valid RBG value.                    
    elif not( isinstance(test_object[0],float) and 
              test_object[0] >= 0 and test_object[0] <= 1.0):
        raise exceptions.ASEValueError, \
                    "The parameter, %s, must be a sequence of 3 floats " \
                    "with a value between 0 and 1.\n"%(variable_name) + \
                    " "*30 + "The first element has a value of %s."\
                    %(test_object[0])
    
    # Check that the second element is valid RBG value.                    
    elif not( isinstance(test_object[1],float) and 
              test_object[1] >= 0.0 and test_object[1] <= 1.0):
        raise exceptions.ASEValueError, \
                    "The parameter, %s, must be a sequence of 3 floats " \
                    "with a value between 0 and 1.\n"%(variable_name) + \
                    " "*30 + "The second element has a value of %s."\
                    %(test_object[1])

    # Check that the thrid element is valid RBG value.                    
    elif not( isinstance(test_object[2],float) and 
              test_object[2] >= 0.0 and test_object[2] <= 1.0):
        raise exceptions.ASEValueError, \
                    "The parameter, %s, must be a sequence of 3 floats " \
                    "with a value between 0 and 1.\n"%(variable_name) + \
                    " "*30 + "The thrid element has a value of %s."\
                    %(test_object[2])

    #Since everything thing is okay and no error was encountered, return it.
    else:             
        return test_object
