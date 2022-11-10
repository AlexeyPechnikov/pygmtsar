from distributed.diagnostics.progressbar import ProgressBar
from distributed.utils import LoopRunner
from tqdm.auto import tqdm
from tornado.ioloop import IOLoop

class TqdmDaskProgress(ProgressBar):
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
        self._loop_runner = loop_runner = LoopRunner(loop=loop)
        super().__init__(keys, scheduler, interval, complete)
        self.tqdm = tqdm(total=9e6,**tqdm_kwargs)

        if start:
            loop_runner.run_sync(self.listen)
        
    @property
    def loop(self):
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
        warnings.warn(
            "setting the loop property is deprecated", DeprecationWarning, stacklevel=2
        )
        self.__loop = value

    def _draw_bar(self, remaining, all, **kwargs):
        if self.tqdm.total == 9e6:
            self.tqdm.reset(total=all)
        update_ct = (all - remaining) - self.tqdm.n
        self.tqdm.update(update_ct)

    def _draw_stop(self, **kwargs):
        self.tqdm.close()

def tqdm_dask(futures, **kwargs):
    from distributed.client import futures_of
    
    futures = futures_of(futures)
    if not isinstance(futures, (set, list)):
        futures = [futures]
    TqdmDaskProgress(futures, **kwargs)
