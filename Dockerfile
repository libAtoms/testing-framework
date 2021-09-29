FROM python:3.8

RUN pip install numpy scipy matplotlib seaborn 
RUN pip install atomistica quippy-ase

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
