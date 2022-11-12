import json
from Bio import Entrez


def get_mesh_term(query_result):
    mesh_terms = []
    for k, v in query_result.items():
        mesh_terms.append(v['heading'])
    return mesh_terms


def search_pubmed(mesh_terms):

    # year_count_dict = {}

    for term in mesh_terms:
        print(f"|| Retrieving drug therapy for {term} from PubMed...")
        query = f'"{term}/drug therapy"[MeSH Major Topic]'
        Entrez.email = 'yun-tzu.yin.21@ucl.ac.uk'
        handle = Entrez.esearch(db="pubmed", term=query)
        record = Entrez.read(handle)
        print(record)

    # with open('../outputs/0_query_results/publication_year_count.txt', 'a') as file:
    #     file.write(json.dumps(year_count_dict))
    #     file.close()

    # return year_count_dict





def search_medline(mesh_terms):

    year_count_all = {}

    for term in mesh_terms:
        print(term)
        year_start = 2000
        year_end = 2001
        year_count_dict = {}
        for year in range(year_start, year_end):
            # query = f'({term}[MeSH Terms] AND (medline[Filter])'
            query = f'({term}[MeSH Terms]) AND ((medline[Filter]) AND ({year}/1/1:{year}/12/31[pdat]))'
            Entrez.email = 'yun-tzu.yin.21@ucl.ac.uk'
            handle = Entrez.esearch(db="pubmed", term=query)
            record = Entrez.read(handle)
            year_count_dict[year] = record['Count']
        
        year_count_all[term] = year_count_dict












if __name__ == "__main__":

    with open('../outputs/0_query_results/mesh_query_result.txt') as file:
        query_result = json.loads(file.read())
        file.close()
    
    mesh_terms = get_mesh_term(query_result)
    search_pubmed(mesh_terms)