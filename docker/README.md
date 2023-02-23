[![GMTSAR tests](https://github.com/gmtsar/gmtsar/actions/workflows/gmtsar.yml/badge.svg)](https://github.com/gmtsar/gmtsar/actions/workflows/gmtsar.yml)

GMTSAR is an open source (GNU General Public License) InSAR processing system designed for users familiar with Generic Mapping Tools (GMT). The code is written in C and will compile on any computer where GMT and NETCDF are installed. The system has three main components:
 - a preprocessor for each satellite data type (ERS-1/2, Envisat, ALOS-1, TerraSAR-X, COSMOS-SkyMed, Radarsat-2, Sentinel-1A/B, and ALOS-2) to convert the native format and orbital information into a generic format;
 - an InSAR processor to focus and align stacks of images, map topography into phase, and form the complex interferogram;
 - a postprocessor, mostly based on GMT, to filter the interferogram and construct interferometric products of phase, coherence, phase gradient, and line-of sight displacement in both radar and geographic coordinates;
GMT is used to display all the products as pdf files and KML images for Google Earth. A set of shell scripts has been developed for standard 2-pass processing as well as geometric image alignment for stacking and time series. Users are welcome to contribute to this effort. 

![](https://topex.ucsd.edu/gmtsar/images/landing-rotate/001_intf.jpg)

## Docker Image

Download and launch the Docker image to open "ubuntu" user terminal using the commands below:
```
docker pull mobigroup/gmtsar
docker run -it --name gmtsar docker.io/mobigroup/gmtsar
```

or built the image yourself locally using Docker file inside the repository and open the terminal: 
```
docker build --no-cache -t gmtsar . -f gmtsar.Dockerfile
docker run -it --name gmtsar docker.io/library/gmtsar
```
To download orbit data for ERS and Envisat (6 GB) enter this command:
```
download_orbits.sh
```

Note: use no-password "sudo" for all operations requiring "root" user access.

Note: to build the multi-arch images and share them to DockerHub  use the command below with your own repository name:

```
docker buildx create --name gmtsar
docker buildx use gmtsar
docker buildx inspect --bootstrap
docker buildx build . -f gmtsar.Dockerfile \
    --platform linux/amd64,linux/arm64 \
    --tag mobigroup/gmtsar:2022-10-11 \
    --tag mobigroup/gmtsar:latest \
    --pull --push --no-cache
```

## Run Examples

Run example "Sentinel-1 TOPS Time Series" from https://topex.ucsd.edu/gmtsar/downloads/ using the commands:

```
wget -c http://topex.ucsd.edu/gmtsar/tar/S1A_Stack_CPGF_T173.tar.gz
tar xvzf S1A_Stack_CPGF_T173.tar.gz -C .
./README_prep.txt
./README_proc.txt
./README_sbas.txt
```

## View Results

View the results in the terminal window using terminal image viewer:

```
chafa intf_all/2015020_2015068/phasefilt_mask_ll.png
```
<img src="https://user-images.githubusercontent.com/7342379/195099112-6bccc8e5-a57c-4189-ad94-c3fa9a10b891.jpg" width="60%" />

View the results in Docket Desktop CLI using terminal image viewer:

```
chafa intf_all/2015020_2015068/phasefilt_mask.pdf
```
<img src="https://user-images.githubusercontent.com/7342379/195099102-738f404d-8a80-426b-8c53-9cbe0d4a7ec2.jpg" width="60%" />


```
chafa SBAS/vel_ll.png
```
<img src="https://user-images.githubusercontent.com/7342379/195388458-80112ce8-3c87-470e-ae77-99f8a45894c4.png" width="60%" />

## Citations

Sandwell, D. ., R. . Mellors, X. Tong, M. Wei, and P. Wessel (2011), Open radar interferometry software for mapping surface deformation, Eos Trans. AGU, 92(28), doi:10.1029/2011EO280002. 

Sandwell, D., Mellors, R., Tong, X., Xu, X., Wei, M., Wessel, P. (2016) GMTSAR: An InSAR Processing System Based on Generic Mapping Tools. UC San Diego: Scripps Institution of Oceanography. Retrieved from: http://topex.ucsd.edu/gmtsar/tar/GMTSAR_2ND_TEX.pdf 

## Learn more

- Official site: https://topex.ucsd.edu/gmtsar/

- Issue tracker: https://github.com/gmtsar/gmtsar/issues

- Source code: https://github.com/gmtsar/gmtsar

- Example datasets: https://topex.ucsd.edu/gmtsar/downloads/