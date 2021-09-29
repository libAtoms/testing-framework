FROM python:3.8

RUN apt-get -y update && apt-get upgrade -y && apt-get install -y gfortran libopenblas-dev

RUN pip install numpy scipy matplotlib seaborn 
RUN pip install quippy-ase
RUN pip install --global-option=build_ext --global-option="-L/opt/OpenBLAS/lib" git+https://github.com/Atomistica/atomistica

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}

# include files from repository in the container
COPY --chown=${NB_USER} . ${HOME}
