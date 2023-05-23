# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2022, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from distributed.diagnostics.progressbar import ProgressBar
from distributed.utils import LoopRunner
from tqdm.auto import tqdm
from tornado.ioloop import IOLoop

class TqdmDaskProgress(ProgressBar):
    """
    Class to create a progress bar for Dask computation using tqdm.

    Inherits from Dask's ProgressBar class.
    """
    __loop = None

    def __init__(
        self,
        keys,
        scheduler=None,
        interval="100ms",
        loop=None,
        complete=True,
        start=True,
        **tqdm_kwargs
    ):
        """
        Initialize the progress bar.

        Parameters
        ----------
        keys : iterable
            Iterable of keys that uniquely identify the tasks for which to track progress.
        scheduler : Scheduler, optional
            Dask scheduler instance. Default is None.
        interval : str, optional
            Time interval for updates. Default is "100ms".
        loop : IOLoop, optional
            Tornado IOLoop to use for updates. Default is None.
        complete : bool, optional
            If True, show completed tasks. Default is True.
        start : bool, optional
            If True, start listening for updates immediately. Default is True.
        **tqdm_kwargs : dict
            Additional keyword arguments for tqdm.
        """
        self._loop_runner = loop_runner = LoopRunner(loop=loop)
        super().__init__(keys, scheduler, interval, complete)
        self.tqdm = tqdm(total=9e6,**tqdm_kwargs)

        if start:
            loop_runner.run_sync(self.listen)
        
    @property
    def loop(self):
        """
        Get the current Tornado IOLoop.
        
        Returns
        -------
        IOLoop
            The current Tornado IOLoop.
        """
        loop = self.__loop
        if loop is None:
            # If the loop is not running when this is called, the LoopRunner.loop
            # property will raise a DeprecationWarning
            # However subsequent calls might occur - eg atexit, where a stopped
            # loop is still acceptable - so we cache access to the loop.
            self.__loop = loop = self._loop_runner.loop
        return loop
        
    @loop.setter
    def loop(self, value):
        """
        Set the current Tornado IOLoop.

        Parameters
        ----------
        value : IOLoop
            Tornado IOLoop instance.
        """
        warnings.warn(
            "setting the loop property is deprecated", DeprecationWarning, stacklevel=2
        )
        self.__loop = value

    def _draw_bar(self, remaining, all, **kwargs):
        """
        Update the progress bar.

        Parameters
        ----------
        remaining : int
            Number of tasks remaining.
        all : int
            Total number of tasks.
        **kwargs : dict
            Additional keyword arguments.
        """
        if self.tqdm.total == 9e6:
            self.tqdm.reset(total=all)
        update_ct = (all - remaining) - self.tqdm.n
        self.tqdm.update(update_ct)

    def _draw_stop(self, **kwargs):
        """
        Close the progress bar.
        
        **kwargs : dict
            Additional keyword arguments.
        """
        self.tqdm.close()

def tqdm_dask(futures, **kwargs):
    """
    Helper function to track progress of Dask computation using tqdm.

    Parameters
    ----------
    futures : list
        List of Dask futures.
    **kwargs : dict
        Additional keyword arguments for TqdmDaskProgress.
    """
    from distributed.client import futures_of
    
    futures = futures_of(futures)
    if not isinstance(futures, (set, list)):
        futures = [futures]
    TqdmDaskProgress(futures, **kwargs)
