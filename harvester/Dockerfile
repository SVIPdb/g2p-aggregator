FROM python:2

WORKDIR /app

RUN apt-get update && apt-get install -y postgresql-client jq curl

# # install vep as well
# RUN mkdir -p /tools/vep && cd /tools/vep && \
#     curl -L -O https://github.com/Ensembl/ensembl-vep/archive/release/104.zip && \
#     unzip 104.zip && \
#     cd ensembl-vep-release-104/ && \
#     perl INSTALL.pl

# =======================================================================================
# === VEP SETUP 
# === based on https://github.com/Ensembl/ensembl-vep/blob/release/104/docker/Dockerfile

# Update aptitude and install some required packages
# a lot of them are required for Bio::DB::BigFile
RUN apt-get update && apt-get -y install \
    build-essential \
    cpanminus \
    curl \
    default-libmysqlclient-dev \
    libpng-dev \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    locales \
    openssl \
    perl \
    perl-base \
    unzip \
    samtools \
    vim && \
    apt-get -y purge manpages-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Setup VEP environment
ENV OPT /opt/vep
ENV OPT_SRC $OPT/src
ENV PERL5LIB_TMP $PERL5LIB:$OPT_SRC/ensembl-vep:$OPT_SRC/ensembl-vep/modules
ENV PERL5LIB $PERL5LIB_TMP:$OPT_SRC/bioperl-live
ENV KENT_SRC $OPT/src/kent-335_base/src
ENV HTSLIB_DIR $OPT_SRC/htslib
ENV MACHTYPE x86_64
ENV DEPS $OPT_SRC
ENV PATH $OPT_SRC/ensembl-vep:$OPT_SRC/var_c_code:$PATH
ENV LANG_VAR en_US.UTF-8

# install vep from ensembl-vep image
COPY --from=ensemblorg/ensembl-vep:release_104.3 /opt/vep $OPT

# # Install bioperl-ext, faster alignments for haplo (XS-based BioPerl extensions to C libraries)
# WORKDIR $OPT_SRC/bioperl-ext/Bio/Ext/Align/
# RUN perl Makefile.PL && make && make install && rm -f Makefile*

# # Install ensembl-xs, faster run using re-implementation in C of some of the Perl subroutines
# WORKDIR $OPT_SRC/ensembl-xs
# RUN perl Makefile.PL && make && make install && rm -f Makefile* cpanfile

WORKDIR $OPT_SRC
# Install/compile more libraries
RUN cpanm --installdeps --with-recommends --notest Archive::Zip DBD::mysql DBI && \
    # Configure "locale", see https://github.com/rocker-org/rocker/issues/19
    echo "$LANG_VAR UTF-8" >> /etc/locale.gen && locale-gen en_US.utf8 && \
    /usr/sbin/update-locale LANG=$LANG_VAR && \
    # Copy htslib executables. It also requires the packages 'zlib1g-dev', 'libbz2-dev' and 'liblzma-dev'
    cp $HTSLIB_DIR/bgzip $HTSLIB_DIR/tabix $HTSLIB_DIR/htsfile /usr/local/bin/

# RUN cp $HTSLIB_DIR/bgzip $HTSLIB_DIR/tabix $HTSLIB_DIR/htsfile /usr/local/bin/

ENV LC_ALL $LANG_VAR
ENV LANG $LANG_VAR

ENV PERL5LIB $PERL5LIB_TMP

# Final steps
WORKDIR $OPT_SRC/ensembl-vep
# Update bash profile
RUN echo >> $OPT/.profile && \
    echo PATH=$PATH:\$PATH >> $OPT/.profile && \
    echo export PATH >> $OPT/.profile && \
    # Run INSTALL.pl and remove the ensemb-vep tests and travis
    ./INSTALL.pl -a a -l -n && rm -rf t travisci .travis.yml

# === end VEP setup ===
# =======================================================================================

# set up a volume for VEP so it has someplace to download its data
VOLUME $OPT/data

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD [ "/bin/bash", "entrypoint.sh" ]
