# ----------------------------------------------------------------------------
# InSAR.dev
# 
# This file is part of the InSAR.dev project: https://InSAR.dev
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .tqdm_joblib import tqdm_joblib
#from .S1 import S1

class ASF(tqdm_joblib):
    import pandas as pd
    from datetime import timedelta

    # check for downloaded burst files
    # like S1_370328_IW1_20150121T134421_VV_DBBE-BURST
    template_burst = '{burstId}/*/{burst}.*'

    def __init__(self, username=None, password=None):
        import asf_search
        import getpass
        if username is None:
            username = getpass.getpass('Please enter your ASF username and press Enter key:')
        if password is None:
            password = getpass.getpass('Please enter your ASF password and press Enter key:')
        self.username = username
        self.password = password

    def _get_asf_session(self):
        import asf_search
        return asf_search.ASFSession().auth_with_creds(self.username, self.password)

    # https://asf.alaska.edu/datasets/data-sets/derived-data-sets/sentinel-1-bursts/
    def download(self, basedir, bursts, session=None, n_jobs=8, joblib_backend='loky', skip_exist=True,
                        retries=30, timeout_second=3, debug=False):
        """
        Downloads the specified bursts extracted from Sentinel-1 SLC scenes.

        Parameters
        ----------
        basedir : str
            The directory where the downloaded bursts will be saved.
        bursts : list of str
            List of burst identifiers to download.
        session : asf_search.ASFSession, optional
            The session object for authentication. If None, a new session is created.
        n_jobs : int, optional
            The number of concurrent download jobs. Default is 8.
        joblib_backend : str, optional
            The backend for parallel processing. Default is 'loky'.
        skip_exist : bool, optional
            If True, skips downloading bursts that already exist. Default is True.
        debug : bool, optional
            If True, prints debugging information. Default is False.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the list of downloaded bursts.
        """
        import rioxarray as rio
        from tifffile import TiffFile
        import xmltodict
        from xml.etree import ElementTree
        import pandas as pd
        import asf_search
        import joblib
        from tqdm.auto import tqdm
        import os
        import glob
        from datetime import datetime, timedelta
        import time
        import warnings
        # supress asf_search 'UserWarning: File already exists, skipping download'
        warnings.filterwarnings("ignore", category=UserWarning)

        def filter_azimuth_time(items, start_utc_dt, stop_utc_dt, delta=3):
            return [item for item in items if
                 datetime.strptime(item['azimuthTime'], '%Y-%m-%dT%H:%M:%S.%f') >= start_utc_dt - timedelta(seconds=delta) and
                 datetime.strptime(item['azimuthTime'], '%Y-%m-%dT%H:%M:%S.%f') <= stop_utc_dt + timedelta(seconds=delta)]

        if isinstance(bursts, str):
            bursts = list(filter(None, map(str.strip, bursts.split('\n'))))

        # create the directory if needed
        os.makedirs(basedir, exist_ok=True)

        # skip existing bursts
        if skip_exist:
            bursts_missed = []
            for burst in bursts:
                #print (burst)
                # orbital path is not included into burst
                fakeBurstId = '_'.join(['*'] + burst.split('_')[1:3])
                template = self.template_burst.format(burstId=fakeBurstId, burst=burst)
                #print ('template', template)
                files = glob.glob(template, root_dir=basedir)
                #print ('files', files)
                exts =[ os.path.splitext(file)[-1] for file in files]
                #print ('exts', exts)
                if '.tiff' in exts and '.xml' in exts and len(exts)==4:
                    #print ('pass')
                    pass
                else:
                    bursts_missed.append(burst)
        else:
            # process all the defined scenes
            bursts_missed = bursts
        #print ('bursts_missed', len(bursts_missed))
        # do not use internet connection, work offline when all the scenes already available
        if len(bursts_missed) == 0:
            return

        def download_burst(result, basedir, session):
            properties = result.geojson()['properties']
            #print ('result properties', properties)
            burst = properties['fileID']
            burstId = properties['burst']['fullBurstID']
            burstIndex = properties['burst']['burstIndex']
            platform = properties['platform'][-2:]
            polarization = properties['polarization']
            #print ('polarization', polarization)
            subswath = properties['burst']['subswath']

            # create the directories if needed
            burst_dir = os.path.join(basedir, burstId)
            tif_dir = os.path.join(burst_dir, 'measurement')
            xml_annot_dir = os.path.join(burst_dir, 'annotation')
            xml_noise_dir = os.path.join(burst_dir, 'noise')
            xml_calib_dir = os.path.join(burst_dir, 'calibration')
            # save annotation using the burst and scene names
            xml_file = os.path.join(xml_annot_dir, f'{burst}.xml')
            xml_noise_file = os.path.join(xml_noise_dir, f'{burst}.xml')
            xml_calib_file = os.path.join(xml_calib_dir, f'{burst}.xml')
            #rint ('xml_file', xml_file)
            tif_file = os.path.join(tif_dir, f'{burst}.tiff')
            #print ('tif_file', tif_file)
            for dirname in [burst_dir, tif_dir, xml_annot_dir, xml_noise_dir, xml_calib_dir]:
                os.makedirs(dirname, exist_ok=True)

            # download tif
            # properties['bytes'] is not an accurate file size but it looks about 40 kB smaller
            if os.path.exists(tif_file) and os.path.getsize(tif_file) >= int(properties['bytes']):
                #print (f'pass {tif_file}')
                pass
            else:
                #print ('YYY', os.path.getsize(tif_file), properties['bytes'])
                # remove potentially incomplete file if needed
                if os.path.exists(tif_file):
                    os.remove(tif_file)
                # check if we can open the downloaded file without errors
                tmp_file = os.path.join(burst_dir, os.path.basename(tif_file))
                # remove potentially incomplete data file if needed
                if os.path.exists(tmp_file):
                    os.remove(tmp_file)
                result.download(burst_dir, filename=os.path.basename(tif_file), session=session)
                if not os.path.exists(tmp_file):
                    raise Exception(f'ERROR: TiFF file is not downloaded: {tmp_file}')
                if os.path.getsize(tmp_file) == 0:
                    raise Exception(f'ERROR: TiFF file is empty: {tmp_file}')
                # check TiFF file validity opening it
                with TiffFile(tmp_file) as tif:
                    # get TiFF file information
                    page = tif.pages[0]
                    tags = page.tags
                    data = page.asarray()
                # attention: rasterio can crash the interpreter on a corrupted TIFF file
                # perform this check as the final step
                with rio.open_rasterio(tmp_file) as raster:
                    raster.load()
                # TiFF file is well loaded
                if not os.path.exists(tmp_file):
                    raise Exception(f'ERROR: TiFF file is missed: {tmp_file}')
                # move to persistent name
                if os.path.exists(tmp_file):
                    os.rename(tmp_file, tif_file)

            # download xml
            if os.path.exists(xml_file) and os.path.getsize(xml_file) > 0 \
                and os.path.exists(xml_noise_file) and os.path.getsize(xml_noise_file) > 0 \
                and os.path.exists(xml_calib_file) and os.path.getsize(xml_calib_file) > 0:
                #print (f'pass {xml_file}')
                pass
            else:
                # get TiFF file information
                with TiffFile(tif_file) as tif:
                    page = tif.pages[0]
                    offset = page.dataoffsets[0]
                #print ('offset', offset)
                # get the file name
                basename = os.path.basename(properties['additionalUrls'][0])
                #print ('basename', '=>', basename)
                manifest_file = os.path.join(burst_dir, basename)
                # remove potentially incomplete manifest file if needed
                if os.path.exists(manifest_file):
                    os.remove(manifest_file)
                asf_search.download_urls(urls=properties['additionalUrls'], path=burst_dir, session=session)
                if not os.path.exists(manifest_file):
                    raise Exception(f'ERROR: manifest file is not downloaded: {manifest_file}')
                if os.path.getsize(manifest_file) == 0:
                    raise Exception(f'ERROR: manifest file is empty: {manifest_file}')
                # check XML file validity parsing it
                with open(manifest_file, 'r') as file:
                    xml_content = file.read()
                    _ = ElementTree.fromstring(xml_content)
                # xml file is well parsed
                if not os.path.exists(manifest_file):
                    raise Exception(f'ERROR: manifest file is missed: {manifest_file}')
                # parse xml
                with open(manifest_file, 'r') as file:
                    xml_content = file.read()
                # remove manifest file
                if os.path.exists(manifest_file):
                    os.remove(manifest_file)

                subswathidx = int(subswath[-1:]) - 1
                content = xmltodict.parse(xml_content)['burst']['metadata']['product'][subswathidx]
                assert polarization == content['polarisation'], 'ERROR: XML polarization differs from burst polarization'
                annotation = content['content']

                annotation_burst = annotation['swathTiming']['burstList']['burst'][burstIndex]
                start_utc = annotation_burst['azimuthTime']
                start_utc_dt = datetime.strptime(start_utc, '%Y-%m-%dT%H:%M:%S.%f')
                #print ('start_utc', start_utc, start_utc_dt)

                length = int(annotation['swathTiming']['linesPerBurst'])
                #print (f'length={length}, burstIndex={burstIndex}')
                azimuth_time_interval = annotation['imageAnnotation']['imageInformation']['azimuthTimeInterval']
                burst_time_interval = timedelta(seconds=(length - 1) * float(azimuth_time_interval))
                stop_utc_dt = start_utc_dt + burst_time_interval
                stop_utc = stop_utc_dt.strftime('%Y-%m-%dT%H:%M:%S.%f')
                #print ('stop_utc', stop_utc, stop_utc_dt)

                # output xml
                product = {}

                adsHeader = annotation['adsHeader']
                adsHeader['startTime'] = start_utc
                adsHeader['stopTime'] = stop_utc
                adsHeader['imageNumber'] = '001'
                product = product   | {'adsHeader': adsHeader}

                qualityInformation = {'productQualityIndex': annotation['qualityInformation']['productQualityIndex']} |\
                                      {'qualityDataList':     annotation['qualityInformation']['qualityDataList']}
                product = product   | {'qualityInformation': qualityInformation}

                generalAnnotation = annotation['generalAnnotation']
                # filter annotation['generalAnnotation']['replicaInformationList'] by azimuthTime
                product = product   | {'generalAnnotation': generalAnnotation}

                imageAnnotation = annotation['imageAnnotation']
                imageAnnotation['imageInformation']['productFirstLineUtcTime'] = start_utc
                imageAnnotation['imageInformation']['productLastLineUtcTime'] = stop_utc
                imageAnnotation['imageInformation']['productComposition'] = 'Assembled'
                imageAnnotation['imageInformation']['sliceNumber'] = '0'
                imageAnnotation['imageInformation']['sliceList'] = {'@count': '0'}
                imageAnnotation['imageInformation']['numberOfLines'] = str(length)
                # imageStatistics and inputDimensionsList are not updated
                product = product   | {'imageAnnotation': imageAnnotation}

                dopplerCentroid = annotation['dopplerCentroid']
                items = filter_azimuth_time(dopplerCentroid['dcEstimateList']['dcEstimate'], start_utc_dt, stop_utc_dt)
                dopplerCentroid['dcEstimateList'] = {'@count': len(items), 'dcEstimate': items}
                product = product   | {'dopplerCentroid': dopplerCentroid}

                antennaPattern = annotation['antennaPattern']
                items = filter_azimuth_time(antennaPattern['antennaPatternList']['antennaPattern'], start_utc_dt, stop_utc_dt)
                antennaPattern['antennaPatternList'] = {'@count': len(items), 'antennaPattern': items}
                product = product   | {'antennaPattern': antennaPattern}

                swathTiming = annotation['swathTiming']
                items = filter_azimuth_time(swathTiming['burstList']['burst'], start_utc_dt, start_utc_dt, 1)
                assert len(items) == 1, 'ERROR: unexpected bursts count, should be 1'
                # add TiFF file information
                items[0]['byteOffset'] = offset
                swathTiming['burstList'] = {'@count': len(items), 'burst': items}
                product = product   | {'swathTiming': swathTiming}

                geolocationGrid = annotation['geolocationGrid']
                items = filter_azimuth_time(geolocationGrid['geolocationGridPointList']['geolocationGridPoint'], start_utc_dt, stop_utc_dt, 1)
                # re-numerate line numbers for the burst
                for item in items: item['line'] = str(int(item['line']) - (length * burstIndex))
                geolocationGrid['geolocationGridPointList'] = {'@count': len(items), 'geolocationGridPoint': items}
                product = product   | {'geolocationGrid': geolocationGrid}

                product = product   | {'coordinateConversion': annotation['coordinateConversion']}
                product = product   | {'swathMerging': annotation['swathMerging']}

                with open(xml_file, 'w') as file:
                    file.write(xmltodict.unparse({'product': product}, pretty=True, indent='  '))

                # output noise xml
                content = xmltodict.parse(xml_content)['burst']['metadata']['noise'][subswathidx]
                assert polarization == content['polarisation'], 'ERROR: XML polarization differs from burst polarization'
                annotation = content['content']

                noise = {}

                adsHeader = annotation['adsHeader']
                adsHeader['startTime'] = start_utc
                adsHeader['stopTime'] = stop_utc
                adsHeader['imageNumber'] = '001'
                noise = noise   | {'adsHeader': adsHeader}

                if 'noiseVectorList' in annotation:
                    noiseRangeVector = annotation['noiseVectorList']
                    items = filter_azimuth_time(noiseRangeVector['noiseVector'], start_utc_dt, stop_utc_dt)
                    # re-numerate line numbers for the burst
                    for item in items: item['line'] = str(int(item['line']) - (length * burstIndex))
                    noiseRangeVector = {'@count': len(items), 'noiseVector': items}
                    noise = noise   | {'noiseVectorList': noiseRangeVector}

                if 'noiseRangeVectorList' in annotation:
                    noiseRangeVector = annotation['noiseRangeVectorList']
                    items = filter_azimuth_time(noiseRangeVector['noiseRangeVector'], start_utc_dt, stop_utc_dt)
                    # re-numerate line numbers for the burst
                    for item in items: item['line'] = str(int(item['line']) - (length * burstIndex))
                    noiseRangeVector = {'@count': len(items), 'noiseRangeVector': items}
                    noise = noise   | {'noiseRangeVectorList': noiseRangeVector}

                if 'noiseAzimuthVectorList' in annotation:
                    noiseAzimuthVector = annotation['noiseAzimuthVectorList']
                    items = noiseAzimuthVector['noiseAzimuthVector']['line']['#text'].split(' ')
                    items = [int(item) for item in items]
                    lowers = [item for item in items if item <= burstIndex * length] or items[0]
                    uppers = [item for item in items if item >= (burstIndex + 1) * length - 1] or items[-1]
                    mask = [True if item>=lowers[-1] and item<=uppers[0] else False for item in items]
                    items = [item - burstIndex * length for item, m in zip(items, mask) if m]
                    noiseAzimuthVector['noiseAzimuthVector']['firstAzimuthLine'] = lowers[-1] - burstIndex * length
                    noiseAzimuthVector['noiseAzimuthVector']['lastAzimuthLine'] = uppers[0] - burstIndex * length
                    noiseAzimuthVector['noiseAzimuthVector']['line'] = {'@count': len(items), '#text': ' '.join([str(item) for item in items])}
                    items = noiseAzimuthVector['noiseAzimuthVector']['noiseAzimuthLut']['#text'].split(' ')
                    items = [item for item, m in zip(items, mask) if m]
                    noiseAzimuthVector['noiseAzimuthVector']['noiseAzimuthLut'] = {'@count': len(items), '#text': ' '.join(items)}
                    noise = noise   | {'noiseAzimuthVectorList': noiseAzimuthVector}

                with open(xml_noise_file, 'w') as file:
                    file.write(xmltodict.unparse({'noise': noise}, pretty=True, indent='  '))

                # output calibration xml
                content = xmltodict.parse(xml_content)['burst']['metadata']['calibration'][subswathidx]
                assert polarization == content['polarisation'], 'ERROR: XML polarization differs from burst polarization'
                annotation = content['content']

                calibration = {}

                adsHeader = annotation['adsHeader']
                adsHeader['startTime'] = start_utc
                adsHeader['stopTime'] = stop_utc
                adsHeader['imageNumber'] = '001'
                calibration = calibration   | {'adsHeader': adsHeader}

                calibration = calibration   | {'calibrationInformation': annotation['calibrationInformation']}

                calibrationVector = annotation['calibrationVectorList']
                items = filter_azimuth_time(calibrationVector['calibrationVector'], start_utc_dt, stop_utc_dt)
                # re-numerate line numbers for the burst
                for item in items: item['line'] = str(int(item['line']) - (length * burstIndex))
                calibrationVector = {'@count': len(items), 'calibrationVector': items}
                calibration = calibration   | {'calibrationVectorList': calibrationVector}

                with open(xml_calib_file, 'w') as file:
                    file.write(xmltodict.unparse({'calibration': calibration}, pretty=True, indent='  '))

        # prepare authorized connection
        if session is None:
            session = self._get_asf_session()

        with tqdm(desc=f'ASF Downloading Bursts Catalog', total=1) as pbar:
            results = asf_search.granule_search(bursts_missed)
            pbar.update(1)

        if n_jobs is None or debug == True:
            print ('Note: sequential joblib processing is applied when "n_jobs" is None or "debug" is True.')
            joblib_backend = 'sequential'

        def download_burst_with_retry(result, basedir, session, retries, timeout_second):
            for retry in range(retries):
                try:
                    download_burst(result, basedir, session)
                    return True
                except Exception as e:
                    print(f'ERROR: download attempt {retry+1} failed for {result}: {e}')
                    if retry + 1 == retries:
                        return False
                time.sleep(timeout_second)

        # download bursts
        with self.tqdm_joblib(tqdm(desc='ASF Downloading Sentinel-1 SLC Bursts', total=len(bursts_missed))) as progress_bar:
            statuses = joblib.Parallel(n_jobs=n_jobs, backend=joblib_backend)(joblib.delayed(download_burst_with_retry)\
                                    (result, basedir, session, retries=retries, timeout_second=timeout_second) for result in results)

        failed_count = statuses.count(False)
        if failed_count > 0:
            raise Exception(f'Bursts downloading failed for {failed_count} items.')
        # parse processed bursts and convert to dataframe
        bursts_downloaded = pd.DataFrame(bursts_missed, columns=['burst'])
        # return the results in a user-friendly dataframe
        return bursts_downloaded

    @staticmethod
    def search(geometry, startTime=None, stopTime=None, flightDirection=None,
               platform='SENTINEL-1', processingLevel='auto', polarization='VV', beamMode='IW'):
        import geopandas as gpd
        import asf_search
        import shapely

        # cover defined time interval
        if len(startTime)==10:
            startTime=f'{startTime} 00:00:01'
        if len(stopTime)==10:
            stopTime=f'{stopTime} 23:59:59'

        if flightDirection == 'D':
            flightDirection = 'DESCENDING'
        elif flightDirection == 'A':
            flightDirection = 'ASCENDING'

        # convert to a single geometry
        if isinstance(geometry, (gpd.GeoDataFrame, gpd.GeoSeries)):
            geometry = geometry.geometry.union_all()
        # convert closed linestring to polygon
        if geometry.type == 'LineString' and geometry.coords[0] == geometry.coords[-1]:
            geometry = shapely.geometry.Polygon(geometry.coords)
        if geometry.type == 'Polygon':
            # force counterclockwise orientation.
            geometry = shapely.geometry.polygon.orient(geometry, sign=1.0)
        #print ('wkt', geometry.wkt)

        if isinstance(processingLevel, str) and processingLevel=='auto' and platform == 'SENTINEL-1':
            processingLevel = asf_search.PRODUCT_TYPE.BURST

        # search bursts
        results = asf_search.search(
            start=startTime,
            end=stopTime,
            flightDirection=flightDirection,
            intersectsWith=geometry.wkt,
            platform=platform,
            processingLevel=processingLevel,
            polarization=polarization,
            beamMode=beamMode,
        )
        return gpd.GeoDataFrame.from_features([product.geojson() for product in results], crs="EPSG:4326")

    @staticmethod
    def plot(bursts):
        import pandas as pd
        import matplotlib
        import matplotlib.pyplot as plt

        bursts['date'] = pd.to_datetime(bursts['startTime']).dt.strftime('%Y-%m-%d')
        bursts['label'] = bursts.apply(lambda rec: f"{rec['flightDirection'].replace('E','')[:3]} {rec['date']} [{rec['pathNumber']}]", axis=1)
        unique_labels = sorted(bursts['label'].unique())
        unique_paths = sorted(bursts['pathNumber'].astype(str).unique())
        colors = {label[-4:-1]: 'orange' if label[0] == 'A' else 'cyan' for i, label in enumerate(unique_labels)}
        fig, ax = plt.subplots(figsize=(10, 8))
        for label, group in bursts.groupby('label'):
            group.plot(ax=ax, edgecolor=colors[label[-4:-1]], facecolor='none', linewidth=1, alpha=1, label=label)
        burst_handles = [matplotlib.lines.Line2D([0], [0], color=colors[label[-4:-1]], lw=1, label=label) for label in unique_labels]
        aoi_handle = matplotlib.lines.Line2D([0], [0], color='red', lw=1, label='AOI')
        handles = burst_handles + [aoi_handle]
        ax.legend(handles=handles, loc='upper right')
        ax.set_title('Sentinel-1 Burst Footprints')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
