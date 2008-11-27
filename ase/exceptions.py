"""
Module for implementing all the exceptions handling.
"""

class ASEError(Exception):
    """
    Class for the basic exception in Atomic Simulation Environment"
    """

class ASEUnitError(ASEError):
    """
    Class for representing when an error is raised
    due to a incorrect unit.
    """

class ASETypeError(ASEError):
    """
    Class for representing when an error is raised
    due to a incorrect unit.
    """

class ASEValueError(ASEError):
    """
    Class for representing when an error is raised
    due to a incorrect value.
    """


