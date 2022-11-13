import json, io
import pandas as pd
from Bio import Entrez


def get_mesh_term(query_result):
    '''Output MeSH terms as a list from given MeSH query results'''
    mesh_terms = []
    for k, v in query_result.items():
        mesh_terms.append(v['heading'])
    return mesh_terms


def search_pubmed(mesh_terms):
    '''Query PubMed by MeSH terms and retrieve articles with MeSH drug therapy tags'''

    retrieved_results = {}
    for term in mesh_terms:
        print(f'|| Retrieving drug therapy for {term} from PubMed...')
        retrieved_list = []
        query = f'("{term}/drug therapy"[MeSH Major Topic]) \
                  AND ((booksdocs[Filter] OR meta-analysis[Filter] OR review[Filter] OR systematicreview[Filter]) \
                  AND (humans[Filter]) \
                  AND (medline[Filter]))'
        Entrez.email = 'yun-tzu.yin.21@ucl.ac.uk'
        handle_search = Entrez.esearch(db='pubmed', term=query, retmax='100', sort='relevance')
        record_search = Entrez.read(handle_search)
        handle = Entrez.efetch(db='pubmed', id=record_search['IdList'], retmode='xml')
        record = Entrez.read(handle)
        for article in record['PubmedArticle']:
            retrieved_result = {'pmid': str(),
                                'year': str(),
                                'title': str(),
                                'summary': str(),
                                'substances': list()}
            try:
                if article['MedlineCitation']['ChemicalList']:
                    retrieved_result['pmid'] = str(article['MedlineCitation']['PMID'])
                    retrieved_result['year'] = article['MedlineCitation']['DateCompleted']['Year']
                    retrieved_result['title'] = article['MedlineCitation']['Article']['ArticleTitle']
                    for chemical in article['MedlineCitation']['ChemicalList']:
                        retrieved_result['substances'].append(str(chemical['NameOfSubstance']))
                    if len(article['MedlineCitation']['Article']['Abstract']['AbstractText']) == 1:
                        retrieved_result['summary'] = article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
                    else:
                        for text in article['MedlineCitation']['Article']['Abstract']['AbstractText']:
                            if text.attributes['NlmCategory'].lower() == 'conclusions' \
                               or text.attributes['Label'].lower() == 'conclusion' \
                               or text.attributes['Label'].lower() == 'conclusions' \
                               or text.attributes['Label'].lower() == 'short conclusion' \
                               or text.attributes['Label'].lower() == "authors' conclusions" \
                               or text.attributes['Label'].lower() == 'conclusions and relevance' \
                               or text.attributes['Label'].lower() == 'summary' \
                               or text.attributes['Label'].lower() == 'findings' \
                               or text.attributes['Label'].lower() == 'recent findings' \
                               or text.attributes['Label'].lower() == 'interpretation' \
                               or text.attributes['Label'].lower() == 'expert opinion':
                                retrieved_result['summary'] = str(text)
                    retrieved_list.append(retrieved_result)
            except KeyError:
                pass
        retrieved_results[term] = retrieved_list

    with open('../outputs/0_query_results/drug_therapy_from_pubmed.txt', 'a') as file:
        file.write(json.dumps(retrieved_results))
        file.close()

    return retrieved_results


def tabulate_pubmed_result(query_result, drug_therapy):
    '''Process and tabulate PubMed results, and output CSV files for each disease'''

    mesh_terms_dict = {v['heading']:k for k,v in query_result.items()}
    table_column = ['MeSH Term', 'Year of Publication', 'Medication', 'Details', 'PMID', 'PubMed Link']
    summary_table_column = ['Disease Name', 'MeSH Term', 'Medication', 'PMID of Reference']
    summary_table_list = []
    number = 0
    for k, v in drug_therapy.items():
        print(f'|| Processing and tabulating PubMed results for {k}...')
        number += 1
        substance_count = {}
        for article in v:
            for substance in article['substances']:
                if substance in substance_count.keys():
                    substance_count[substance] += 1
                else:
                    substance_count[substance] = 1

        if len(v) > 10:            
            for sub, count in substance_count.copy().items():
                if count <= 10:
                    del substance_count[sub]

        table_list = []
        pmid_list = []
        for article in v:
            for substance in article['substances']:
                if substance in substance_count.keys() and len(article['summary']) != 0:
                    year = article['year']
                    pmid = article['pmid']
                    pmid_list.append(pmid)
                    medication = ', '.join(str(x) for x in article['substances'])
                    details = article['summary']
                    pubmed_link = f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'
                    row = pd.DataFrame([[k, year, medication, details, pmid, pubmed_link]], columns=table_column)
                    table_list.append(row)
                    break

        if len(table_list) != 0:
            table_frame = pd.concat(table_list).sort_values(by=['Year of Publication'])  
            year_range = f'{table_frame["Year of Publication"].min()}-{table_frame["Year of Publication"].max()}'
            medication_sum = ', '.join(str(x) for x in substance_count.keys())
            summary_frame = pd.DataFrame([[k, year_range, medication_sum, '', '', '']], columns=table_column)
            table = pd.concat([summary_frame, table_frame])

            with io.open(f'../outputs/2_therapy_suggestion/{number}_{k}.csv', 'w', encoding='utf-8') as output:
                table.to_csv(output)

            row_sum = pd.DataFrame([[mesh_terms_dict[k], k, medication_sum,', '.join(str(x) for x in pmid_list)]], columns=summary_table_column)
            summary_table_list.append(row_sum)
        else:
            print(f'No eligible medication was found for {k}.')
    
    summary_table = pd.concat(summary_table_list)  
    with io.open(f'../outputs/2_therapy_suggestion/medication_suggestion_summary.csv', 'w', encoding='utf-8') as output:
        summary_table.to_csv(output)



if __name__ == "__main__":

    with open('../outputs/0_query_results/mesh_query_result.txt') as file:
        query_result = json.loads(file.read())
        file.close()
    
    mesh_terms = get_mesh_term(query_result)
    drug_therapy = search_pubmed(mesh_terms)

    tabulate_pubmed_result(query_result, drug_therapy)