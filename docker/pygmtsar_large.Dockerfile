# place the Dockerfile into your data files directory
FROM mobigroup/pygmtsar:latest

ENV ORIG=./
ENV DATA=data_desc

# download example notebooks
RUN mkdir tmp \
&&  cd tmp \
&&  svn export https://github.com/mobigroup/YamchiDam/trunk/notebooks \
&&  cd .. \
&&  mv tmp/notebooks/*.ipynb . \
&&  rm -rf tmp

RUN mkdir data_desc
COPY --chown="${NB_UID}:${NB_GID}" $ORIG/DEM_WGS84.nc .
COPY --chown="${NB_UID}:${NB_GID}" $ORIG/*.EOF "$DATA"
COPY --chown="${NB_UID}:${NB_GID}" $ORIG/S1A_*.SAFE/annotation/*.xml   "$DATA"
# split layers
COPY --chown="${NB_UID}:${NB_GID}" $ORIG/S1A_*202201*.SAFE/measurement/*.tiff "$DATA"
COPY --chown="${NB_UID}:${NB_GID}" $ORIG/S1A_*202202*.SAFE/measurement/*.tiff "$DATA"
COPY --chown="${NB_UID}:${NB_GID}" $ORIG/S1A_*202203*.SAFE/measurement/*.tiff "$DATA"
COPY --chown="${NB_UID}:${NB_GID}" $ORIG/S1A_*202204*.SAFE/measurement/*.tiff "$DATA"
