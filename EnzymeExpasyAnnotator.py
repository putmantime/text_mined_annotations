import Bio.ExPASy.Enzyme as bee
import re
import pprint
import requests
import urllib.request
import sys
import pandas as pd
import tempfile

__author__ = 'timputman'

enzyme_out = open('enz_out.txt', 'w')

api_key = 'xxxx'


def link_compound2chebi(compound):
    """
    used NCBO Annotator from BioPortal to return ChEBI IDS
    for substrates and products of reactions from Expasy enzyme

    """
    url = 'http://data.bioontology.org/annotator'
    params = dict(apikey=api_key, text=compound, ontologies='CHEBI', longest_only='true',
                  include='properties', exlude_numbers='false', exclude_synonyms='false', mappins='all')
    tm_results = requests.get(url=url, params=params).json()

    for i in tm_results:
        prefLabel = i['annotatedClass']['properties']['http://data.bioontology.org/metadata/def/prefLabel']
        text_input = i['annotations'][0]['text']
        if prefLabel[0].lower() == text_input.lower():
            final_chebi = i['annotatedClass']['@id'].split("/")[-1]
            return final_chebi.split("_")[-1]


def replace_strings(compound):

    """
    takes list of compound names and replaces with names formatted for NCBO annotator
    """

    names = {'H(2)O': 'water',
             'CO(2)': 'carbon dioxide',
             'An alcohol': 'alcohol',
             'A secondary alcohol': 'secondary alcohol',
             'O(2)': 'dioxygen',
             'H(+)': 'hydron',
             'a ketone': 'ketone',
            }

    if compound in names:
        return names[compound]

    else:
        return compound


def get_expasy_enzyme():
    """

    """
    url = "ftp://ftp.expasy.org/databases/enzyme/enzyme.dat"
    print("Retrieving enzyme records from Expasy Enzyme")
    enzyme = urllib.request.urlretrieve(url)
    enzyme_p = bee.parse(open(enzyme[0], 'r'))
    chebiout = open('chebi_list.txt', 'w')
    annotations = open('annotations_out.txt', 'w')

    enz_records = []
    chebi_list = []
    count = 0
    tester = []
    for record in enzyme_p:
        enz_rec = {}
        count += 1
        print(count)
        enz_rec['ECNumber'] = record['ID']
        enz_rec['Reaction(s)'] = []
        enz_rec['Substrates'] = {}
        enz_rec['Products'] = {}
        #enz_records.append(enz_rec)

        # split split to seperate multiple reactions
        reaction1 = record['CA'].split('.')

        for rxn in reaction1:
            try:
                if len(reaction1) > 2:
                    rxn = rxn[3:]
                enz_rec['Reaction(s)'].append(rxn)
                #split reactions into [substrates, products]
                constituents = rxn.split('=')
                # split each side of reaction on '+' not '(+)'
                r = re.compile(r'(?:[^\+(]|\([^)]*\))+')
                subr = r.findall(constituents[0])
                for sub in subr:
                    sub = sub.lstrip().rstrip()
                    sub = replace_strings(sub)
                    schebi = link_compound2chebi(sub)
                    enz_rec['Substrates'][sub] = schebi

                    if schebi:
                        chebi_list.append(schebi)
                prodr = r.findall(constituents[-1])
                for prod in prodr:
                    prod = prod.lstrip().rstrip()
                    prod = replace_strings(prod)
                    pchebi = link_compound2chebi(prod)
                    enz_rec['Products'][prod] = pchebi
                    if pchebi:
                        chebi_list.append(pchebi)
            except Exception as e:
                print(e)
                continue

        enz_records.append(enz_rec)
    print(chebi_list, file=chebiout)
    print(enz_records, file=annotations)
    return enz_records

expasy = get_expasy_enzyme()
OUT = open('enzyme_parsed.txt', 'w')
print(expasy, file=OUT)

