import Bio.ExPASy.Enzyme as bee
import re
import pprint
import requests

__author__ = 'timputman'

enzyme_out = open('enz_out.txt', 'w')


def link_compound2chebi(compound):
    """
    used NCBO Annotator from BioPortal to return ChEBI IDS
    for substrates and products of reactions from Expasy enzyme

    """
    url = 'http://data.bioontology.org/annotator'
    params = dict(apikey='8b5b7825-538d-40e0-9e9e-5ab9274a9aeb', text=compound, ontologies='CHEBI', longest_only='true',
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
             'CO(2)': 'carbon dioxide'
            }

    for k, v in names.items():
        if compound == k:
            return compound.replace(k,v)
        else:
            return compound


enzyme = open('enzyme.dat', 'r')
enzyme_p = bee.parse(enzyme)
enz_rec= {}
count = 0

# use Bio.Expasy.enzyme to parse enzyme.dat record from ISB expasy enzyme
for record in enzyme_p:
    count += 1

    #create record for each enzyme with EC number as primary key
    enz_rec['PreferedName'] = record['DE']
    enz_rec['ECNumber'] = record['ID']
    enz_rec['Reaction(s)'] = []
    enz_rec['Substrates'] = {}
    enz_rec['Products'] = {}
    enz_rec['UniProt'] = {}
    # split split to seperate multiple reactions
    reaction1 = record['CA'].split('.')

    for rxn in reaction1:
        if len(reaction1) > 2:
            rxn = rxn[3:]
        enz_rec['Reaction(s)'].append(rxn)
        #split reactions into [substrates, products]
        constituents = rxn.split('=')
        # split each side of reaction on '+' not '(+)'
        r = re.compile(r'(?:[^\+(]|\([^)]*\))+')
        for sub in r.findall(constituents[0]):
            sub = replace_strings(sub.lstrip().rstrip())
            schebi = link_compound2chebi(sub)
            enz_rec['Substrates'][sub] = schebi

        for prod in r.findall(constituents[-1]):
            prod = replace_strings(prod.lstrip().rstrip())
            pchebi = link_compound2chebi(prod)
            enz_rec['Products'][prod] = pchebi

            # populate enz_rec['UniProt'] with dictionary of uniprotid:name key, value pairs for protein
        for unpid in record['DR']:
            enz_rec['UniProt'][unpid[0]] = unpid[1]
    print(count)
    print(enz_rec, file=enzyme_out)
