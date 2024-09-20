# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_base import Stack_base
from .PRM import PRM

class Stack_prm(Stack_base):

    def PRM(self, date=None, subswath=None):
        """
        Open a PRM (Parameter) file.

        Parameters
        ----------
        date : str, optional
            The date of the PRM file. If None or equal to self.reference, return the reference PRM file. Default is None.
        multi : bool, optional
            If True, open a multistem PRM file. If False, open a stem PRM file. Default is True.
        singleswath : bool, optional
            If True, open a single-digit subswath PRM file instead of a merged (multi-digit) one. Default is False.

        Returns
        -------
        PRM
            An instance of the PRM class representing the opened PRM file.
        """
        import os

        # check if subswath exists or return a single subswath for None
        # workaround for stack_rep_subswath()
        if subswath is None:
            subswath = self.get_subswath()

        if date is None:
            date == self.reference

        prefix = self.multistem_stem(subswath, date)
        filename = os.path.join(self.basedir, f'{prefix}.PRM')
        #print (filename)
        return PRM.from_file(filename)