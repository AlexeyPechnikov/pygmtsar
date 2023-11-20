# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------

# the code adopted from my project https://github.com/mobigroup/ParaView-plugins
class NCubeVTK:
    """
    This class provides a static method for mapping a 2D image onto a 3D topography.

    Example usage:

    ```python
    from pygmtsar import NCubeVTK
    import vtk

    vtk_ugrid = NCubeVTK.ImageOnTopography(sbas.get_dem().to_dataset().rename({'lat': 'y', 'lon': 'x'}))
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName("DEM_WGS84.vtk")
    writer.SetInputData(vtk_ugrid)
    writer.Write()
    vtk_ugrid
    ```
    ```python
    from pygmtsar import NCubeVTK
    import pyvista as pv

    vtk_ugrid = NCubeVTK.ImageOnTopography(sbas.get_dem().to_dataset().rename({'lat': 'y', 'lon': 'x'}))
    vtk_ugrid = pv.UnstructuredGrid(vtk_ugrid)
    vtk_ugrid.save('DEM_WGS84.vtk')
    vtk_ugrid
    ```
    """

    @staticmethod
    def ImageOnTopography(dataset, band_mask=None, use_sealevel=False):
        """
        Map a 2D image onto a 3D topography.

        Parameters:
            dataset: A file name or Xarray Dataset to use for mapping the image.
            band_mask (optional): A string representing a variable in the dataset to use as a mask. Default is None.
            use_sealevel (optional): A boolean flag indicating whether to replace negative topography by sea level (z=0). Default is False.

        Returns:
            A vtkUnstructuredGrid representing the image mapped onto the topography.

        Raises:
            Error if the band_mask variable is not found in the dataset, or if the dataset is not an Xarray Dataset.

        Note: fill NODATA by NAN for float variables.
        """
        from vtk import vtkPoints, vtkStructuredGrid, vtkThreshold, vtkDataObject, \
            VTK_FLOAT, VTK_UNSIGNED_CHAR, vtkStringArray, vtkFloatArray, vtkIntArray
        from vtk.util import numpy_support as vn
        import numpy as np
        import xarray as xr
        # https://corteva.github.io/rioxarray/stable/getting_started/getting_started.html
        #import rioxarray as rio

        if band_mask is not None:
            if band_mask not in dataset.data_vars:
                print (f'ERROR: {band_mask} variable not found')
                return
            #if 'band' not in dataset[band_mask].dims:
            #    print (f'ERROR: {band_mask} variable has no band dimension')
            #    return

        if not isinstance(dataset, xr.Dataset):
            print ('ERROR: Unsupported ds argument, should be supported file name or Xarray Dataset')
            return

        # fill NODATA by NAN
        for data_var in dataset.data_vars:
            #if dataset[data_var].values.dtype in [np.dtype('float16'),np.dtype('float32'),np.dtype('float64')]:
            if np.issubdtype(dataset[data_var].dtype, np.floating):
                dataset[data_var].values = dataset[data_var].values.astype('float32')
                if '_FillValue' in dataset[data_var] and not np.isnan(dataset[data_var]._FillValue):
                    dataset[data_var].values[dataset[data_var].values == dataset[data_var]._FillValue] = np.nan
        
        xs = dataset.x.values
        ys = dataset.y.values        
        (yy,xx) = np.meshgrid(ys, xs)
        if 'z' in dataset:
            zs = dataset.z.values
            znanmask = np.isnan(dataset.z)
            # replace negative topography by sea level (z=0)
            if use_sealevel:
                zs[zs <= 0] = 0
        else:
            # set z coordinate to zero
            zs = np.zeros_like(xx)
            znanmask = 0
            
        # create raster mask by geometry and for NaNs
        vtk_points = vtkPoints()
        points = np.column_stack((xx.ravel('F'), yy.ravel('F'), zs.ravel('C')))
        vtk_points.SetData(vn.numpy_to_vtk(points, deep=True))

        sgrid = vtkStructuredGrid()
        sgrid.SetDimensions(len(xs), len(ys), 1)
        sgrid.SetPoints(vtk_points)

        for data_var in dataset.data_vars:
            if 'z' == data_var:
                # variable is already included into coordinates
                continue
            da = dataset[data_var]
            dims = da.dims
            #print (data_var, dims, da.dtype)
            if 'band'in dims:
                bands = da.band.shape[0]
                #print ('bands', bands)
                if bands in [3,4]:
                    # RGB or RGBA, select 3 bands only
                    array = vn.numpy_to_vtk(da[:3].values.reshape(3,-1).T, deep=True, array_type=VTK_UNSIGNED_CHAR)
                    array.SetName(da.name)
                    sgrid.GetPointData().AddArray(array)
                elif bands == 1:
                    array = vn.numpy_to_vtk(da.values.reshape(1,-1).T, deep=True, array_type=VTK_FLOAT)
                    array.SetName(da.name)
                    sgrid.GetPointData().AddArray(array)
                else:
                    print (f'ERROR: Unsupported bands count (should be 1,3 or 4) {bands} for variable {data_var}')
                    return
            else:
                array = vn.numpy_to_vtk(da.values.ravel(), deep=True, array_type=VTK_FLOAT)
                array.SetName(da.name)
                sgrid.GetPointData().AddArray(array)

        if band_mask is not None and 'band' in dataset[band_mask]:
            # mask NaNs for band variable
            nanmask = np.any(np.isnan(dataset[band_mask]),axis=0)
            # mask zero areas for color bands
            # a bit of magic to enhance image (mask single channel zeroes)
            # that's more correct way, only totally black pixels ignored
            #zeromask = (~np.all(image.values==0,axis=0)).astype(float)
            # that's a magic for better borders
            zeromask = np.any(ds[band_mask].values==0,axis=0)
            mask = nanmask | zeromask | znanmask
        elif band_mask is not None:
            # mask NaNs
            nanmask = np.isnan(dataset[band_mask])
            mask = nanmask | znanmask
        #else:
        #    mask = znanmask

        if band_mask is not None:
            array = vn.numpy_to_vtk(mask.values.ravel(), deep=True, array_type=VTK_UNSIGNED_CHAR)
            array.SetName('band_mask')
            sgrid.GetPointData().AddArray(array)

        for coord in dataset.coords:
            #print (coord, dataset[coord].dims, len(dataset[coord].dims))
            if len(dataset[coord].dims) > 0:
                # real coordinate, ignore it
                continue
            #print (coord, dataset[coord].dtype)
            #print ('np.datetime64', np.issubdtype(dataset[coord].dtype, np.datetime64))
            #print ('np.int64', np.issubdtype(dataset[coord].dtype, np.int64))
            #print ('np.float64', np.issubdtype(dataset[coord].dtype, np.float64))
            if np.issubdtype(dataset[coord].dtype, np.datetime64):
                data_array = vtkStringArray()
                data_value = str(dataset[coord].values)
            elif np.issubdtype(dataset[coord].dtype, np.int64):
                data_array = vtkIntArray()
                data_value = dataset[coord].values.astype(np.int64)
            elif np.issubdtype(dataset[coord].dtype, np.float64):
                data_array = vtkFloatArray
                data_value = dataset[coord].values.astype(np.float64)
            else:
                print(f'NOTE: unsupported attribute {coord} datatype {dataset[coord].dtype}, miss it')
                continue
            data_array.SetName(coord)
            data_array.InsertNextValue(data_value)
            sgrid.GetFieldData().AddArray(data_array)

        # required for 3D interactive rendering on Google Colab
        if band_mask is None:
            return sgrid

        thresh = vtkThreshold()
        thresh.SetInputData(sgrid)
        thresh.SetInputArrayToProcess(0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_POINTS, 'band_mask')
        # drop float NANs
        #thresh.SetLowerThreshold(-1e30)
        #thresh.SetUpperThreshold(1e30)
        thresh.SetUpperThreshold(0)
        thresh.Update()
        ugrid = thresh.GetOutput()
        # remove internal mask variable from output
        ugrid.GetPointData().RemoveArray('band_mask')
        return ugrid

#vtk_ugrid = NCubeVTK.ImageOnTopography(ds)
#pv.UnstructuredGrid(vtk_ugrid)
