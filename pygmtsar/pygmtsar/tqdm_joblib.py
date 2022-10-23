#!/usr/bin/env python3
# Alexey Pechnikov, Oct, 2022, https://github.com/mobigroup/gmtsar

class tqdm_joblib:

    import contextlib
    @staticmethod
    @contextlib.contextmanager
    def tqdm_joblib(tqdm_object):
        import joblib
        """Context manager to patch joblib to report into tqdm progress bar given as argument"""
        class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)

            def __call__(self, *args, **kwargs):
                tqdm_object.update(n=self.batch_size)
                return super().__call__(*args, **kwargs)

        old_batch_callback = joblib.parallel.BatchCompletionCallBack
        joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
        try:
            yield tqdm_object
        finally:
            joblib.parallel.BatchCompletionCallBack = old_batch_callback
            tqdm_object.close()
