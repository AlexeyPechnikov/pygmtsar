# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2024, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
class MultiInstanceManager:
    def __init__(self, *instances):
        self.instances = instances
        self.context_params = {}

    def __getattr__(self, name):
        def method_wrapper(*args, **kwargs):
            results = []
            for instance in self.instances:
                # Adjust arguments with context-specific parameters if applicable
                instance_kwargs = {**kwargs}
                for key, values in self.context_params.items():
                    if len(values) == len(self.instances):
                        instance_kwargs[key] = values[self.instances.index(instance)]
                instance_method = getattr(instance, name)
                results.append(instance_method(*args, **instance_kwargs))
            return results

        return method_wrapper

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.context_params = {}

    def apply(self, **kwargs):
        """Prepare specific attributes or arguments for each instance."""
        for key, value in kwargs.items():
            if isinstance(value, list) and len(value) == len(self.instances):
                self.context_params[key] = value
            else:
                raise ValueError(f"Each key in apply must have a list of values equal to the number of instances")
        return self
