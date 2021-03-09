# start from QUIP base container, i.e. without GAP (to avoid licence check)
FROM libatomsquip/quip-base:latest

# add some additional Python packages
RUN pip install seaborn atomistica

# set dynamic library path to include OpenBLAS
ENV LD_LIBRARY_PATH=/opt/OpenBLAS/lib

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