#!/usr/bin/env python
import os
import sys
import importlib, pkgutil

from utils_ex.iterables import grouper_flat
from utils_ex.tqdm_logger import TqdmLoggingHandler

sys.path.append('silos')  # NOQA
sys.path.append('harvesters')  # NOQA

import json
import xxhash
import argparse
import logging
import logging.config
import yaml

from utils_ex.instrumentation import DelayedOpLogger, show_runtime_stats
from normalizers import gene_enricher, disease_normalizer, oncogenic_normalizer, reference_genome_normalizer, \
    biomarker_normalizer, location_normalizer, myvariant_enricher
from filters import has_hgvs

# these silos are added on line 4, but adding them again allows for code navigation
# FIXME: the sys.path.append on line 4 should be removed in favor of a real import
from silos.elastic_silo import ElasticSilo
import silos.elastic_silo as elastic_silo
from silos.kafka_silo import KafkaSilo
import silos.kafka_silo as kafka_silo
from silos.file_silo import FileSilo
import silos.file_silo as file_silo
from silos.postgres_silo import PostgresSilo
import silos.postgres_silo as postgres_silo

import requests_cache
import hashlib

# add handler to global logger
logger = logging.getLogger()
logger.addHandler(TqdmLoggingHandler())

# # drop into a shell if we raise an uncaught exception
# sys.excepthook = ultratb.FormattedTB(mode='Verbose', color_scheme='Linux', call_pdb=1)


# announce each item that's added
VERBOSE_ITEMS = False

# we just need to check for membership, so a set is faster than a list
DUPLICATES = set()

# cache responses
# expire the contents of the cache after 72 hours (which should be long enough for the harvester to run once)
requests_cache.install_cache('harvester', allowable_codes=(200, 400, 404), expire_after=60*60*72)

# a list of normalizers to annotate each feature association (and optionally defines what to report on the console if
# they take more time than expected to run)
normalizers = [
    # FIXME: drug_normalizer removes drug combination info and destroys the capitalization returned by each source.
    #  it's currently disabled, but it should eventually be fixed.
    # (drug_normalizer, lambda dol:
    #     feature_association['association'].get('environmentalContexts', None)),

    (disease_normalizer, lambda dol, feature_association:
        feature_association['association']['phenotypes'][0]['description']
        if 'phenotypes' in feature_association['association']
           and len(feature_association['association']['phenotypes']) > 0 else None),
    (oncogenic_normalizer, None),  # functionality for oncogenic_normalizer already mostly in harvesters
    (location_normalizer, None),
    (reference_genome_normalizer, None),
    (biomarker_normalizer, None),
    (gene_enricher, None),
    (myvariant_enricher, None)
]

# a list of filters applied to feature associations immediately before yielding to the silo; if any filter returns
# false for that feature association, it's discarded
filters = [
    has_hgvs
]

args = None
silos = []


class NormalizerException(Exception):
    """
    Raised when a normalizer fails to complete for any reason.
    """
    pass


def is_duplicate(feature_association):
    """ return true if already harvested """
    is_dup = False
    m = hashlib.md5()
    try:
        m.update(json.dumps(feature_association, sort_keys=True))
        hexdigest = m.hexdigest()
        if hexdigest in DUPLICATES:
            is_dup = True
            logging.info('is duplicate {}'.format(
                feature_association['association']['evidence']))
        else:
            DUPLICATES.add(hexdigest)
    except Exception as e:
        logging.warn('duplicate {}'.format(e))
    return is_dup


def _make_silos(in_args):
    """ construct silos """
    instances = []
    for s in in_args.silos:
        if s == 'elastic':
            instances.append(ElasticSilo(in_args))
        if s == 'kafka':
            instances.append(KafkaSilo(in_args))
        if s == 'file':
            instances.append(FileSilo(in_args))
        if s == 'postgres':
            instances.append(PostgresSilo(in_args))
    return instances


def memoized(src, harvester, genes, phases):
    """
    Records the contents of src and saves them to a file. On the next run, returns the results of the file if it exists.
    :param src: the generator from which to acquire samples if they don't already exist
    :param harvester: the name of the harvester, used to identify the cache
    :param genes: whether we were just harvesting or harvesting+converting, also used to identify the cache
    :param phases: whether we were just harvesting or harvesting+converting, also used to identify the cache
    :return: a generator over either the contents of src (if no cache found) or the contents of the cache
    """
    # the cache key basically encodes the script's arguments as part of the cache's filename, so we can re-acquire if we
    # change the arguments.
    cache_key = xxhash.xxh32(json.dumps([genes, phases])).hexdigest()
    cache_path = os.path.join(".harvest_cache", "%s_%s.jsonl" % (cache_key, harvester))

    # check if there's something in the cache to use
    if os.path.exists(cache_path):
        logging.info("Using cached values for %s (in file %s)" % (harvester, cache_path))
        with open(cache_path, "r") as cache_fp:
            for record in cache_fp:
                yield json.loads(record)
    else:
        # populate the cache while yielding samples from source, then eventually save it
        logging.info("Writing output from %s (to file %s)" % (harvester, cache_path))
        with open(cache_path, "w") as cache_fp:
            for record in src:
                cache_fp.write("%s\n" % json.dumps(record))
                yield record


# cgi, jax, civic, oncokb
def harvest(genes):
    """ get evidence from all sources """
    for h in args.harvesters:
        logging.info("=> Initializing harvester: %s" % h)
        harvester = importlib.import_module("harvesters.%s" % h)
        if args.delete_source:
            for silo in silos:
                if h == 'cgi_biomarkers':
                    h = 'cgi'
                silo.delete_source(h)

        # if we're testing, we can avoid hitting the remote sources repeatedly if we can be reasonably sure that they
        # haven't changed. note that this *does not* perform cache validation, so use with caution. you also need to
        # manually clear the cache (the files in .harvest_cache) if you want to recreate the memoization.
        if args.memoize_harvest:
            assoc_source = memoized(harvester.harvest_and_convert(genes),
                                    harvester=h, genes=args.genes, phases=args.phases)
        else:
            assoc_source = harvester.harvest_and_convert(genes)

        for feature_association in assoc_source:
            if VERBOSE_ITEMS:
                logging.info(
                    '{} yielded feat for gene {}, {} w/evidence level {}'.format(
                        harvester.__name__,
                        feature_association['genes'],
                        feature_association['feature_names'],
                        feature_association['association']['evidence_level']
                    )
                )
            yield feature_association


def harvest_only(genes):
    """ get evidence from all sources """
    for h in args.harvesters:
        harvester = importlib.import_module("harvesters.%s" % h)
        if args.delete_source:
            for silo in silos:
                if h == 'cgi_biomarkers':
                    h = 'cgi'
                silo.delete_source(h)

        # see the notes in harvest() above
        if args.memoize_harvest:
            assoc_source = memoized(harvester.harvest(genes),
                                    harvester=h, genes=args.genes, phases=args.phases)
        else:
            assoc_source = harvester.harvest(genes)

        for evidence in assoc_source:
            yield {'source': h, h: evidence}


def normalize(feature_association):
    """ standard representation of drugs,disease etc. """

    for normalizer, more in normalizers:
        with DelayedOpLogger(normalizer.__name__) as d:
            try:
                normalizer.normalize_feature_association(feature_association)
                if more:
                    d.logdelayed(more(d, feature_association))
            except Exception, ex:
                # FIXME: this is probably too broad of an exception clause, but almost anything
                #  can throw an exception in the normalizer...review if we can tighten it up, though.
                # fyi, downstream bits *do* depend on normalizers running, so we might need to bail
                #  on the entire entry if something goes wrong here.
                raise NormalizerException(), None, sys.exc_info()[2]


# FIXME: _check_dup() does a lot more than just check for duplicates; it injects essential info into each entry
#  it should probably have its name changed to indicate that it's a more essential part of the pipeline, and
#  the "normalizers" should be renamed to reflect that they also add important annotations/corrections
def _check_dup(harvested):
    total_fas = 0
    ignored_fas = 0
    duplicated_fas = 0
    normalizer_failed_fas = 0

    for feature_association in harvested:
        try:
            total_fas += 1

            feature_association['tags'] = []
            feature_association['dev_tags'] = []
            normalize(feature_association)

            # ensure it's not a duplicate
            if is_duplicate(feature_association):
                duplicated_fas += 1
                continue

            # ensure it passes all the filters
            if any(not f.filter_feature_association(feature_association) for f in filters):
                ignored_fas += 1
                continue

            yield feature_association

        except NormalizerException as ex:
            normalizer_failed_fas += 1
            # logging.warn("Normalizer failed, skipping; reason: %s" % str(ex))
            continue

    if ignored_fas > 0:
        logging.info("FAs that failed to pass filters: %d (out of %d)" % (ignored_fas, total_fas))
    if duplicated_fas > 0:
        logging.info("Duplicated FAs ignored: %d (out of %d)" % (duplicated_fas, total_fas))
    if duplicated_fas > 0:
        logging.info("Failed normalizer FAs: %d (out of %d)" % (normalizer_failed_fas, total_fas))


def main():
    global args
    global silos
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--harvesters',  nargs='+',
                           help='''harvest from these sources. default:
                                   [cgi_biomarkers,jax,civic,oncokb,cosmic,clinvar,
                                   pmkb]''',
                           default=['cgi_biomarkers', 'jax', 'civic', 'cosmic', 'clinvar',
                                    'oncokb', 'pmkb', 'brca', 'jax_trials',
                                    'molecularmatch_trials'])

    argparser.add_argument('--silos',  nargs='+',
                           help='''save to these silos. default:[elastic]''',
                           default=['elastic'],
                           choices=['elastic', 'kafka', 'file', 'postgres'])

    argparser.add_argument('--delete_index', '-d',
                           help='''delete all from index''',
                           default=False, action="store_true")

    argparser.add_argument('--delete_source', '-ds',
                           help='delete all content for any harvester',
                           default=False, action="store_true")

    argparser.add_argument('--skip_existing_genes', '-skg',
                           help='''skips reinserting genes that already exist in the silo''',
                           default=False, action="store_true")

    argparser.add_argument('--memoize_harvest', '-mh',
                           help='memoize the harvested data and use that, if available',
                           default=False, action="store_true")

    argparser.add_argument('--genes',   nargs='+',
                           help='array of hugo ids, no value will harvest all',
                           default=None)

    # if set, splits gene processing into chunks so we don't run out of memory (FIXME: fix the memory issue)
    argparser.add_argument('--gene_chunk_size', type=int,
                           help='size of each chunk of genes to process, no value will disable chunking',
                           default=None)

    argparser.add_argument('--phases',   nargs='+',
                           help='array of harvest phases to run '
                                '[harvest,convert,enrich,all]. default is all',
                           default=['all'],
                           choices=['all', 'harvest'])

    elastic_silo.populate_args(argparser)
    kafka_silo.populate_args(argparser)
    file_silo.populate_args(argparser)
    postgres_silo.populate_args(argparser)

    args = argparser.parse_args()

    # get list of harvesters in harvesters package
    # ignores anything that's a package and not a module for now
    harvesters_available = [x[1] for x in pkgutil.iter_modules(path=['./harvesters']) if not x[2]]

    for h in args.harvesters:
        assert h in harvesters_available, "harvester is not a module: %r" % h

    path = 'logging.yml'
    with open(path) as f:
        config = yaml.load(f)
    logging.config.dictConfig(config)

    logging.info("harvesters: %r" % args.harvesters)
    logging.info("silos: %r" % args.silos)

    # silo-specific argument reporting (consult each silo's populate_args() method for its parameters)
    if 'elastic' in args.silos:
        logging.info("elastic_search: %r" % args.elastic_search)
        logging.info("elastic_index: %r" % args.elastic_index)
    if 'postgres' in args.silos:
        logging.info("postgres_host: %r" % args.pg_host)
        logging.info("postgres_user: %r" % args.pg_user)
        logging.info("postgres_db: %r" % args.pg_db)
    if 'file' in args.silos:
        logging.info("file_output_dir: %r" % args.file_output_dir)

    logging.info("delete_index: %r" % args.delete_index)
    logging.info("phases: %r" % args.phases)

    silos = _make_silos(args)

    if not args.genes:
        logging.info("genes: all")
    else:
        logging.info("genes: %r" % args.genes)

    if args.delete_index:
        for silo in silos:
            silo.delete_all()

    for silo in silos:
        if args.skip_existing_genes:
            # if this is the postgres silo, attempt to get the list of genes already in the db
            # so we can skip harvesting them
            skip_genes = silo.get_genes()

            # remove these genes from args.genes
            remaining_genes = (x for x in args.genes if x not in skip_genes)
        else:
            remaining_genes = args.genes

        if 'all' in args.phases:
            if args.gene_chunk_size:
                for gene_chunk in grouper_flat(remaining_genes, args.gene_chunk_size):
                    logging.info(" -> Processing gene chunk: %s" % (gene_chunk,))
                    silo.save_bulk(_check_dup(harvest(gene_chunk)))
            else:
                silo.save_bulk(_check_dup(harvest(remaining_genes)))
        else:
            # FIXME: this is really just for debugging...is this necessary
            silo.save_bulk(harvest_only(remaining_genes))

    # afterward, print some statistics on how long the normalizers are taking
    show_runtime_stats()


if __name__ == '__main__':
    main()
