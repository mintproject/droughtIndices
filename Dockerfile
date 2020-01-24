FROM continuumio/miniconda3

ADD environment.yml /tmp/environment.yml
RUN conda env update -f /tmp/environment.yml
ADD WM_climate_indices /usr/bin/
