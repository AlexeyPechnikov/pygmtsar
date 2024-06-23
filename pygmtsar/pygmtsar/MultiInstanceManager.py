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
        # import collections.abc to check for iterable
        import collections.abc
        for key, value in kwargs.items():
            # Check if the value is an iterable (but not a string)
            if isinstance(value, collections.abc.Iterable) and not isinstance(value, (str, bytes)):
                value_list = list(value)
                if len(value_list) != len(self.instances):
                    raise ValueError(f"Each key in apply must have a list of values equal to the number of instances")
                self.context_params[key] = value_list
            else:
                raise ValueError("Non-iterable or incorrect number of elements provided")
        return self

#     def lambda_(self, func):
#         results = []
#         for instance in self.instances:
#             instance_args = {k: v[self.instances.index(instance)] for k, v in self.context_params.items()}
#             # Provide the instance directly to the function
#             results.append(func(instance, **instance_args))
#         return results

    def lambda_(self, func):
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
            with sbas.apply(arg1=value1, arg2=value2):
                results = sbas.lambda_(lambda inst, arg1, arg2: inst.custom_method(arg1, arg2))
    
        Note:
            The function passed to this method should be capable of handling the specific attributes
            or the state of the instances as passed. Misalignment between the expected instance state
            and the function's requirements can lead to runtime errors.
        """
        results = []
        for instance in self.instances:
            instance_args = {k: v[self.instances.index(instance)] for k, v in self.context_params.items()}
            # Provide the instance directly to the function
            results.append(func(instance, **instance_args))
        return results

#     def method(self, method_name, **method_args):
#         """
#         with sbas.apply(da=psf):
#             psf = sbas.method('conncomp_main', data=xr.where(np.isfinite(da), 1, 0).chunk(-1))
#         """
#         results = []
#         for instance in self.instances:
#             # Fetch the method from the instance
#             method = getattr(instance, method_name)
#             # Call the method with the arguments processed for each instance
#             processed_args = {k: v[self.instances.index(instance)] for k, v in self.context_params.items()}
#             results.append(method(**processed_args, **method_args))
#         return results

    def execute(self, method_name, **method_args):
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
                psf = sbas.execute('conncomp_main', data=xr.where(np.isfinite(da), 1, 0).chunk(-1))

        Raises:
            AttributeError: If the specified method is not found on an instance.
        """
        results = []
        for instance in self.instances:
            # Fetch the method from the instance
            if hasattr(instance, method_name):
                method = getattr(instance, method_name)
                # Combine context-specific and direct method arguments
                instance_args = {**{k: v[self.instances.index(instance)] for k, v in self.context_params.items()}, **method_args}
                results.append(method(**instance_args))
            else:
                raise AttributeError(f"{instance} does not have a method named {method_name}")
        return results
