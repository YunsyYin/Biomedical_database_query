import json
import numpy as np
import matplotlib.pyplot as plt
from Bio import Entrez


def get_mesh_term(query_result):
    mesh_terms = []
    for k, v in query_result.items():
        mesh_terms.append(v['heading'])
    return mesh_terms


def get_year_count_from_medline(mesh_terms):

    year_count_dict = {}

    for term in mesh_terms:
        print(f'|| Retrieving year count of {term} from PubMed...')
        year_start = 1900
        year_end = 2023
        year_count= {}
        for year in range(year_start, year_end):
            query = f'({term}[MeSH Terms]) AND ((medline[Filter]) AND ({year}/1/1:{year}/12/31[pdat]))'
            Entrez.email = 'yun-tzu.yin.21@ucl.ac.uk'
            handle = Entrez.esearch(db="pubmed", term=query)
            record = Entrez.read(handle)
            year_count[str(year)] = record['Count']
        
        year_count_dict[term] = year_count

    with open('../outputs/0_query_results/publication_year_count.txt', 'a') as file:
        file.write(json.dumps(year_count_dict))
        file.close()

    return year_count_dict

    
def plot_publication_trend(year_count_dict):
    plt.figure(figsize=(20,10))

    for k, v in year_count_dict.items():
        arr_x = np.array(list(map(int, list(v.keys()))))
        arr_y = np.array(list(map(int, list(v.values()))))
        plt.plot(arr_x, arr_y, label = f'{k}')

    plt.title('Temporal trend of publication on Medline')
    plt.xlabel('Year') 
    plt.ylabel('Number of publication') 
    plt.xlim([1900, 2022]) 
    plt.ylim([0, 7000]) 
    plt.xticks(range(1900, 2022, 10)) 
    plt.legend()
    plt.grid()
    plt.savefig('../outputs/3_trend_of_publications/trend_plot.jpg')
    plt.show()



if __name__ == "__main__":

    # with open('../outputs/0_query_results/mesh_query_result.txt') as file:
    #     query_result = json.loads(file.read())
    #     file.close()
    
    # mesh_terms = get_mesh_term(query_result)
    # year_count_dict = get_year_count_from_medline(mesh_terms)

    with open('../outputs/0_query_results/publication_year_count.txt') as file:
        year_count_dict = json.loads(file.read())
        file.close()

    plot_publication_trend(year_count_dict)

