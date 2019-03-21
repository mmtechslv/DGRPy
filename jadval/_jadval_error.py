#!/usr/bin/env python3

"""
DGRPy is a python package for DGRP microbiome data analysis
The project is under early development and this is my initial release.
NOTE: This is an unstable early release so as soon as I publish stable version I will also add documentation.
Many functions are not at final form at all. For example, JadvalOTU can only load OTU data only in format as OTUs_sample.csv

@author: Farid MUSA(mmtechslv)
"""

__author__ = "Farid MUSA"
__copyright__ = "Copyright (C) 2019, DGRPy Project"
__credits__ = ["Farid MUSA"]
__license__ = "GPLv3"
__version__ = "1.1"
__maintainer__ = "Farid MUSA"
__email__ = "farid.musa.h@gmail.com"
__status__ = "Development"


from functools import wraps as metawrapper
from datetime import datetime
import inspect

from ._jadval_constants import _ERROR_CODES as ERROR_CODES
from ._jadval_constants import _WARNING_CODES as WARNING_CODES

ERROR_LOG = []



def any_error():
    return True if len(ERROR_LOG)>0 else False

def get_error_log():
    return ERROR_LOG if any_error else None

def get_error(error_no):
    if type(error_no)==int:
        if error_no<len(ERROR_LOG):
            return ERROR_LOG[error_no]
    return False

def get_last_error():
    return ERROR_LOG[-1] if any_error else None

def raise_error(code,message=None):
    curframe = inspect.currentframe()
    prevframe = inspect.getouterframes(curframe, 2)
    raise JadvalFreeError(error_code=code,error_message=message,frame_pass=prevframe[1].frame)

def explain_warning_code(warning_code):
    warning_desc = 'N/A'
    if warning_code:
        if warning_code in WARNING_CODES.keys():
            warning_desc = WARNING_CODES[warning_code]
    return warning_desc

def raise_warning(code,message):
    warning_desc = explain_warning_code(code)
    print('Warning! %s [EC: %s]: %s' % (message, code, warning_desc))


def error_handler(error_code=None,error_msj=None):
    def handler_decorator(func_pass):

        @metawrapper(func_pass)
        def func_pass_wrapper(*args, **kwargs):
            try:
                return func_pass(*args, **kwargs)
            except Exception as exception:
                raise JadvalMethodError(source_exception=exception,error_code=error_code,error_message=error_msj,source_pass=func_pass,source_args={'args':args,'kwargs':kwargs})

        return func_pass_wrapper
    return handler_decorator


class JadvalBaseException(Exception):
    """This is base exception class for Jadval Exceptions"""
    def __init__(self,error_code):
        global ERROR_LOG
        self.jadval_error_time = str(datetime.now().strftime('%Y-%m-%d(%H-%M-%S)'))
        self.jadval_error_code = error_code if error_code else False
        self.jadval_error_no = len(ERROR_LOG)
        Exception.__init__(self,self.get_error_message())

    def get_when(self):
        return self.jadval_error_time

    def _log_self_as_error(self):
        global ERROR_LOG
        ERROR_LOG.append(self)

    def get_error_no(self):
        return self.jadval_error_no

    def get_error_message(self):
        error_no = self.jadval_error_no
        if self.jadval_error_code:
            error_code = str(self.jadval_error_code)
            error_desc = str(self.explain_error_code(self.jadval_error_code))
            return '<ErrorNo: %s> [EC: %s]: %s' %(error_no,error_code,error_desc)
        else:
            return 'INCORRECT ERROR MESSAGE'

    def print_report(self):
        report = self.get_error_report()
        spacer = '\n====================' + self.__class__.__name__ + '====================\n'
        print(spacer + ('\n'.join(report)) + '\n')
        return

    def print_message(self):
        message = self.get_error_message()
        print(message)
        return

    def get_code(self):
        return self.jadval_error_code

    @staticmethod
    def explain_error_code(error_code):
        error_desc = 'N/A'
        if error_code:
            if error_code in ERROR_CODES.keys():
                error_desc = ERROR_CODES[error_code]
        return error_desc


class JadvalMethodError(JadvalBaseException):
    """This is main exception class that handles all exceptions that occur in Jadval module via error_handler() decorator."""
    def __init__(self,source_exception=None,error_code=None,error_message=None,source_pass=None,source_args=None):
        JadvalBaseException.__init__(self,error_code)
        self.jadval_additional_message = error_message if error_message else False
        self.jadval_exception = source_exception if source_exception else False
        self.jadval_source_obj = source_pass if source_pass else False
        self.jadval_source_args = source_args['args'] if source_args['args'] else []
        self.jadval_source_kwargs = source_args['kwargs'] if source_args['kwargs'] else []
        self.remake_source_obj()
        self._log_self_as_error()

    def remake_source_obj(self):
        if self.get_source_type():
            self_str = 'self.jadval_source_args[0]'
            method_str = self.jadval_source_obj.__name__
            if self.jadval_source_obj.__name__[0:2] == '__':
                method_str  = '_'+(type(self.jadval_source_args[0]).__name__)+method_str
            obj_str = self_str+'.'+method_str
            self.jadval_source_obj = eval(obj_str)

    def get_source_type(self):
        if 'self' in inspect.getfullargspec(self.jadval_source_obj).args:
            return True
        else:
            return False

    def get_error_report(self):
        report = []
        if self.jadval_error_time:
            report.append('Exception Log Number: ' + str(self.jadval_error_no))
            report.append('Record Time: '  + self.jadval_error_time)
        if self.jadval_exception:
            report.append('Exception Details: ' + str(self.jadval_exception) if isinstance(self.jadval_exception,Exception) else 'Unknown')
        if self.jadval_error_code:
            report.append('Error Code: ' + str(self.jadval_error_code))
            report.append('Error Description: ' + str(self.explain_error_code(self.jadval_error_code)))
        if self.jadval_additional_message:
            report.append('Additional Message: ' + self.jadval_additional_message)
        if self.jadval_source_obj is not None:
            if self.get_source_type():
                report.append('Instance Class Type: ' + str(self.jadval_source_obj.__self__.__class__))
                report.append('Method Name: ' + self.jadval_source_obj.__func__.__name__)
            else:
                report.append('Function Name: ' + self.jadval_source_obj.__name__)
            report.append('Module Name: '+ inspect.getmodule(self.jadval_source_obj).__file__)
        if len(self.jadval_source_args)>0:
            report.append('Args: ['+', '.join([str(arg) for arg in self.jadval_source_args])+']')
        if len(self.jadval_source_kwargs)>0:
            report.append('Kwargs: {' + ', '.join([str(k)+':'+str(v) for k,v in self.jadval_source_kwargs.items()])+'}')
        return report



class JadvalFreeError(JadvalBaseException):
    """This is exception class that handles all exceptions that were raised through raise_error()"""
    def __init__(self,error_code=None,error_message=None,frame_pass=None):
        JadvalBaseException.__init__(self,error_code)
        self.jadval_additional_message = error_message if error_message else False
        self.jadval_source_frame = frame_pass if frame_pass  else False
        self.jadval_source_args = False
        self.make_args()
        self._log_self_as_error()

    def make_args(self):
        if self.jadval_source_frame.f_code.co_argcount>0:
            args = inspect.getargvalues(self.jadval_source_frame).args
            locals = inspect.getargvalues(self.jadval_source_frame).locals
            arg_values = [locals[arg] for arg in args]
            self.jadval_source_args = {'args': args,'values':arg_values}
        return

    def get_source_type(self):
        if 'self' in inspect.getargvalues(self.jadval_source_frame).locals.keys():
            return True
        else:
            return False

    def get_error_report(self):
        report = []
        if self.jadval_error_time:
            report.append('Exception Log Number: ' + str(self.jadval_error_no))
            report.append('Record Time: '  + self.jadval_error_time)
        if self.jadval_error_code:
            report.append('Error Code: ' + str(self.jadval_error_code))
            report.append('Error Description: ' + str(self.explain_error_code(self.jadval_error_code)))
        if self.jadval_additional_message:
            report.append('Additional Message: ' + self.jadval_additional_message)
        if self.jadval_source_frame is not None:
            if self.get_source_type():
                report.append('Instance Class Type: ' + str(self.jadval_source_args['values'][0].__class__.__name__))
                report.append('Method Name: ' + self.jadval_source_frame.f_code.co_name)
            else:
                report.append('Function Name: ' + self.jadval_source_frame.f_code.co_name)
            report.append('Module Name: '+ inspect.getmodulename(inspect.getframeinfo(self.jadval_source_frame).filename))
        if self.jadval_source_args:
            report.append('Args: ['+','.join(self.jadval_source_args['args'])+']')
            report.append('Values: [' + ','.join([str(e) for e in self.jadval_source_args['values']]) + ']')
        return report
