import re, json
import xml.etree.ElementTree as ET


def query_mesh(disease_list):
    '''Search MeSH database with given disease name and return a dict containing MeSH id, heading, synonyms and category'''

    query_result = {}

    for disease_name in disease_list:
        print(f'|| Searching MeSH database by {disease_name}...')
        mesh_dict = {'id': str(), 
                     'heading': str(),
                     'synonyms': list(),
                     'category_code': list(),
                     'category': list()}
        disease_found = False
        tree = ET.parse('../mesh_data/desc2022.xml')

        for record in tree.findall('DescriptorRecord'):
            for entry in record.findall("./ConceptList/Concept/TermList/Term//"):
                if entry.text.lower() == disease_name.lower():
                    disease_found = True
                    mesh_dict['id'] = record.find('./DescriptorUI').text
                    mesh_dict['heading'] = record.find('./DescriptorName/String').text
                    mesh_dict['synonyms'] = [entry.text for entry in record.findall("./ConceptList/Concept/TermList/Term/String")]
                    mesh_dict['category_code'] = set(tree.text.split('.')[0] for tree in record.findall("./TreeNumberList/TreeNumber"))
                    mesh_dict['category'] = []
                
                    meshTreeFile = '../mesh_data/mtrees2022.bin'
                    with open(meshTreeFile, mode='rb') as file:
                        meshTree = file.readlines() 
                    for line in meshTree:
                        line = str(line, encoding = "utf-8")
                        for category_code in mesh_dict['category_code']:
                            if re.search(f'{category_code}$', line):
                                mesh_dict['category'].append(line.split(';')[0])
                    break
                
        if disease_found:
            mesh_dict['category_code'] = list(mesh_dict['category_code'])
        else:
            print(f'Disease not found in MeSH queries: {disease_name}')
        query_result[disease_name] = mesh_dict

    with open('../outputs/0_query_results/mesh_query_result.txt', 'a') as file:
        file.write(json.dumps(query_result))
        file.close()

    return query_result


def output_thesaurus_format(query_result, output_file_name):
    '''Divide MeSH search results by organ system and output them in thesaurus format'''

    categories = set()
    for k, v in query_result.items():
        for category in v['category']:
            categories.add(category)
    categories.remove(output_file_name)
    categories = sorted(list(categories))

    thesaurus_dict = {'Dnode1': {'PT': 'Disease Category', 'NT': list()}}
    for i in range(len(categories)):
        thesaurus_dict['Dnode1']['NT'].append(f'Dnode{i+2}')
        thesaurus_dict[f'Dnode{i+2}'] = {'PT': categories[i], 'NT': list(), 'Child_Dnode': dict()}

        for k, v in query_result.items():
            if categories[i] in v['category']:
                char = chr(97 + len(thesaurus_dict[f'Dnode{i+2}']['NT']))
                thesaurus_dict[f'Dnode{i+2}']['NT'].append(f'Dnode{i+2}{char}')
                thesaurus_dict[f'Dnode{i+2}']['Child_Dnode'][f'Dnode{i+2}{char}'] = {'PT': k, 'SYN': list()}
                for syn in v['synonyms']:
                    thesaurus_dict[f'Dnode{i+2}']['Child_Dnode'][f'Dnode{i+2}{char}']['SYN'].append(syn)

    with open(f'../outputs/1_disease_classification/{output_file_name}.txt', 'a') as file:
        for k1, v1 in thesaurus_dict.items():
            file.write(f'{k1}\n')
            for k2, v2 in v1.items():
                if k2 == 'PT':
                    file.write(f'\t{k2} {v2}\n')
                elif k2 == 'NT' or k2 == 'SYN':
                    for i in v2:
                        file.write(f'\t{k2} {i}\n')
                elif k2 == 'Child_Dnode':
                    for k3, v3 in v2.items():
                        file.write(f'{k3}\n')
                        for k4, v4 in v3.items():
                            if k4 == 'PT':
                                file.write(f'\t{k4} {v4}\n')
                            elif k4 == 'NT' or k4 == 'SYN':
                                for i in v4:
                                    file.write(f'\t{k4} {i}\n')
        file.close()



if __name__ == "__main__":

    disease_list = ['Acute Lymphoblastic Leukemia',
                    'Cancer of the bladder',
                    'Duodenal cancer',
                    'Gastric cancer',
                    'Hemangioblastoma',
                    'Hepatocellular carcinoma',
                    'Hodgkin disease',
                    'Melanoma',
                    'NSCLC',
                    'Pituitary adenoma']

    query_result = query_mesh(disease_list)

    # with open('../outputs/0_query_results/mesh_query_result.txt') as file:
    #     query_result = json.loads(file.read())
    #     file.close()

    output_file_name = 'Neoplasms'
    output_thesaurus_format(query_result, output_file_name)
    