import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON
from time import gmtime, strftime
import sys
import os
import requests
import pprint


__author__ = 'timputman'


def link_disease2org(doid):
    """
    used NCBO Annotator from BioPortal to parse mine the disease organism relationships from the descriptions of Disease ontology
    entries.
    return: Returns a dictionary with data to link wikidata disease to organism
    {'doid': 'DOID:11265',
     'domain': 'Bacteria',
     'genspec': 'Chlamydia trachomatis',
     'rank': 'species',
     'taxid': '813'}
    """
    orgdo = {'doid': doid}
    # split the prefix off of DOID:XXXXX to return XXXXX for url formatting
    doid = doid.split(':')[1]
    try:
        # construct url to request disease ontology record by DOID
        url = 'http://www.ebi.ac.uk/ols/api/ontologies/doid/terms?iri=http://purl.obolibrary.org/obo/DOID_{}'.format(doid)
        r = requests.get(url)
        do_entry = r.json()
        # Parse description from DO record
        description = do_entry['_embedded']['terms'][0]['description'][0]
        # construct the formmatter url to use NCBO Annotator to extract relationship from DO description
        bpurl = 'http://data.bioontology.org/annotator?apikey=0a99c359-d2a2-483a-8dca-148c3bb4e8c1&format=json&ontologies=NCBITAXON&text={}'.format(description)
        tm_results = requests.get(url=bpurl).json()
        # get purl for taxon record to make request to BioPortal taxon and make
        tid_purl = tm_results[0]['annotatedClass']['@id']
        bpurl = 'http://data.bioontology.org/search?apikey=0a99c359-d2a2-483a-8dca-148c3bb4e8c1&q={}&include=properties&format=json'.format(tid_purl)
        tid_results = requests.get(bpurl)
        taxon = tid_results.json()
        rank = taxon['collection'][0]['properties']['http://purl.bioontology.org/ontology/NCBITAXON/RANK'][0]
        if rank == 'species':
            orgdo['rank'] = rank
            orgdo['domain'] = taxon['collection'][0]['properties']['http://purl.bioontology.org/ontology/NCBITAXON/DIV'][0]
            orgdo['taxid'] = tid_purl.split('/')[-1]
            orgdo['genspec'] = taxon['collection'][0]['properties']['http://www.w3.org/2004/02/skos/core#prefLabel'][0]
        else:
            pass
        #
    except Exception as e:
        return None
    return orgdo