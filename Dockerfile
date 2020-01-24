FROM continuumio/miniconda3

ADD environment.yaml /tmp/environment.yaml
RUN conda env update -f /tmp/environment.yaml
ADD WM_climate_indices /usr/bin/
