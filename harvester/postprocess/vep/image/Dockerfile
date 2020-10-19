FROM ensemblorg/ensembl-vep:latest

USER root

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

COPY ensure_cache.sh /opt/ensure_cache.sh
RUN chmod +x /opt/ensure_cache.sh

# note: these should match the host's definition for this user
RUN usermod -u 9099 vep
RUN groupmod -g 1005 vep

USER vep

# ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
