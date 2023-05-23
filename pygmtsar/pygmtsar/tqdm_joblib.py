# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2022, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------

class tqdm_joblib:
    """
    Class to provide a context manager for `joblib` that allows progress reporting via `tqdm`.

    This class is primarily used to facilitate the display of progress information during 
    the execution of parallel tasks with `joblib`.
    """
    import contextlib

    @staticmethod
    @contextlib.contextmanager
    def tqdm_joblib(tqdm_object):
        """
        Context manager to patch `joblib` to report into `tqdm` progress bar given as argument.

        Parameters
        ----------
        tqdm_object : tqdm object
            The `tqdm` progress bar object to be updated.
        """
        import joblib

        class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
            """
            Subclass of `joblib.parallel.BatchCompletionCallBack` to update the `tqdm` 
            progress bar upon the completion of a batch of tasks.
            """
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)

            def __call__(self, *args, **kwargs):
                """
                Overridden method from `joblib.parallel.BatchCompletionCallBack` to update 
                the `tqdm` progress bar and then call the super class method.
                """
                tqdm_object.update(n=self.batch_size)
                return super().__call__(*args, **kwargs)

        old_batch_callback = joblib.parallel.BatchCompletionCallBack
        joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
        try:
            yield tqdm_object
        finally:
            joblib.parallel.BatchCompletionCallBack = old_batch_callback
            tqdm_object.close()
