#!/usr/bin/python
import sys
import importlib, pkgutil

sys.path.append('silos')  # NOQA
sys.path.append('harvesters')  # NOQA

import json
import argparse
import logging
import logging.config
import yaml

from normalizers import drug_normalizer, gene_enricher, disease_normalizer, reference_genome_normalizer, \
    oncogenic_normalizer, biomarker_normalizer, location_normalizer

# these silos are added on line 4, but adding them again allows for code navigation
# FIXME: the sys.path.append on line 4 should be removed in favor of a real import
from silos.elastic_silo import ElasticSilo
import silos.elastic_silo as elastic_silo
from silos.kafka_silo import KafkaSilo
import silos.kafka_silo as kafka_silo
from silos.file_silo import FileSilo
import silos.file_silo as file_silo

import requests_cache
import timeit
import hashlib

# we just need to check for membership, so a set is faster than a list
DUPLICATES = set()

# cache responses
requests_cache.install_cache('harvester', allowable_codes=(200, 400, 404))

args = None
silos = None


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


def _make_silos(args):
    """ construct silos """
    silos = []
    for s in args.silos:
        if s == 'elastic':
            silos.append(ElasticSilo(args))
        if s == 'kafka':
            silos.append(KafkaSilo(args))
        if s == 'file':
            silos.append(FileSilo(args))
    return silos


# cgi, jax, civic, oncokb
def harvest(genes):
    """ get evidence from all sources """
    for h in args.harvesters:
        harvester = importlib.import_module("harvesters.%s" % h)
        if args.delete_source:
            for silo in silos:
                if h == 'cgi_biomarkers':
                    h = 'cgi'
                silo.delete_source(h)

        for feature_association in harvester.harvest_and_convert(genes):
            logging.info(
                '{} {} {}'.format(
                    harvester.__name__,
                    feature_association['genes'],
                    feature_association['association']['evidence_label']
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
        for evidence in harvester.harvest(genes):
            yield {'source': h, h: evidence}


def normalize(feature_association):
    """ standard representation of drugs,disease etc. """
    start_time = timeit.default_timer()
    drug_normalizer.normalize_feature_association(feature_association)
    elapsed = timeit.default_timer() - start_time
    if elapsed > 1:
        environmentalContexts = feature_association['association'].get(
            'environmentalContexts', None)
        logging.info('drug_normalizer {} {}'.format(elapsed,
                                                    environmentalContexts))

    start_time = timeit.default_timer()
    disease_normalizer.normalize_feature_association(feature_association)
    elapsed = timeit.default_timer() - start_time
    if elapsed > 1:
        disease = feature_association['association']['phenotypes'][0]['description']
        logging.info('disease_normalizer {} {}'.format(elapsed, disease))

    start_time = timeit.default_timer()
    # functionality for oncogenic_normalizer already mostly in harvesters
    oncogenic_normalizer.normalize_feature_association(feature_association)
    elapsed = timeit.default_timer() - start_time
    if elapsed > 1:
        logging.info('oncogenic_normalizer {}'.format(elapsed))

    start_time = timeit.default_timer()
    location_normalizer.normalize_feature_association(feature_association)
    elapsed = timeit.default_timer() - start_time
    if elapsed > 1:
        logging.info('location_normalizer {}'.format(elapsed))

    start_time = timeit.default_timer()
    reference_genome_normalizer \
        .normalize_feature_association(feature_association)
    elapsed = timeit.default_timer() - start_time
    if elapsed > 1:
        logging.info('reference_genome_normalizer {}'.format(elapsed))

    start_time = timeit.default_timer()
    biomarker_normalizer.normalize_feature_association(feature_association)
    if elapsed > 1:
        logging.info('biomarker_normalizer {}'.format(elapsed))

    start_time = timeit.default_timer()
    gene_enricher.normalize_feature_association(feature_association)
    if elapsed > 1:
        logging.info('gene_enricher {}'.format(elapsed))


def main():
    global args
    global silos
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--harvesters',  nargs='+',
                           help='''harvest from these sources. default:
                                   [cgi_biomarkers,jax,civic,oncokb,
                                   pmkb]''',
                           default=['cgi_biomarkers', 'jax', 'civic',
                                    'oncokb', 'pmkb', 'brca', 'jax_trials',
                                    'molecularmatch_trials'])

    argparser.add_argument('--silos',  nargs='+',
                           help='''save to these silos. default:[elastic]''',
                           default=['elastic'],
                           choices=['elastic', 'kafka', 'file'])

    argparser.add_argument('--delete_index', '-d',
                           help='''delete all from index''',
                           default=False, action="store_true")

    argparser.add_argument('--delete_source', '-ds',
                           help='delete all content for any harvester',
                           default=False, action="store_true")

    argparser.add_argument('--genes',   nargs='+',
                           help='array of hugo ids, no value will harvest all',
                           default=None)

    argparser.add_argument('--phases',   nargs='+',
                           help='array of harvest phases to run '
                                '[harvest,convert,enrich,all]. default is all',
                           default=['all'],
                           choices=['all', 'harvest'])

    elastic_silo.populate_args(argparser)
    kafka_silo.populate_args(argparser)
    file_silo.populate_args(argparser)

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
    logging.info("elastic_search: %r" % args.elastic_search)
    logging.info("elastic_index: %r" % args.elastic_index)
    logging.info("delete_index: %r" % args.delete_index)
    logging.info("file_output_dir: %r" % args.file_output_dir)
    logging.info("phases: %r" % args.phases)

    silos = _make_silos(args)

    if not args.genes:
        logging.info("genes: all")
    else:
        logging.info("genes: %r" % args.genes)

    if args.delete_index:
        for silo in silos:
            silo.delete_all()

    def _check_dup(harvest):
        for feature_association in harvest:
            feature_association['tags'] = []
            feature_association['dev_tags'] = []
            normalize(feature_association)
            if not is_duplicate(feature_association):
                yield feature_association

    if 'all' in args.phases:
        silos[0].save_bulk(_check_dup(harvest(args.genes)))
    else:
        silos[0].save_bulk(harvest_only(args.genes))


if __name__ == '__main__':
    main()
