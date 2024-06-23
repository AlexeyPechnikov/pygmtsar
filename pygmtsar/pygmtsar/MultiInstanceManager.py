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

    def set_context(self, **kwargs):
        """Prepare specific attributes or arguments for each instance."""
        # import collections.abc to check for iterable
        import collections.abc
        for key, value in kwargs.items():
            if isinstance(value, collections.abc.Iterable) and not isinstance(value, (str, bytes)):
                value_list = list(value)
                if len(value_list) != len(self.instances):
                    raise ValueError(f"Expected a list of {len(self.instances)} values for '{key}', but got {len(value_list)}")
                self.context_params[key] = value_list
            else:
                raise ValueError(f"Provided value for '{key}' is not an appropriate iterable")
        return self

    def run_callable(self, func):
        """
        Execute a lambda function or any callable across all managed instances, 
        passing each instance and its context-specific arguments to the callable.
    
        This method allows for flexible execution of any function that requires instance-level context. 
        It is particularly useful for operations that need to dynamically interact with instance attributes 
        or methods during execution.
    
        Args:
            func (callable): A function or lambda to execute. The function should
                             accept the instance as its first argument, followed by any
                             number of keyword arguments.
    
        Returns:
            list: The results from executing the function across all instances.
    
        Example:
            # Assuming 'func' is a function defined to operate on an instance 'inst' with additional arguments
            with sbas.set_context(arg1=value1, arg2=value2):
                results = sbas.run_callable(lambda inst, arg1, arg2: inst.custom_method(arg1, arg2))
    
        Note:
            The function passed to this method should be capable of handling the specific attributes
            or the state of the instances as passed. Misalignment between the expected instance state
            and the function's requirements can lead to runtime errors.
        """
        results = []
        for instance in self.instances:
            instance_args = {k: v[self.instances.index(instance)] for k, v in self.context_params.items()}
            results.append(func(instance, **instance_args))
        return results

    def run_method(self, method_name, **method_args):
        """
        Execute a specified method on all managed instances with given arguments,
        considering any context-specific adjustments.

        Args:
            method_name (str): The name of the method to execute.
            **method_args: Arbitrary keyword arguments to pass to the method.

        Returns:
            list: The results from each instance's method execution.

        Example:
            with sbas.apply(da=psf):
                psf = sbas.run_method('conncomp_main', data=xr.where(np.isfinite(da), 1, 0).chunk(-1))

        Raises:
            AttributeError: If the specified method is not found on an instance.
        """
        results = []
        for instance in self.instances:
            if hasattr(instance, method_name):
                method = getattr(instance, method_name)
                instance_args = {**{k: v[self.instances.index(instance)] for k, v in self.context_params.items()}, **method_args}
                results.append(method(**instance_args))
            else:
                raise AttributeError(f"{instance} does not have a method named {method_name}")
        return results
