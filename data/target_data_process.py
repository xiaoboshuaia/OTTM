'''
Date: 2022-10-31 17:40:44
LastEditTime: 2022-11-06 17:59:42
FilePath: \Project\OTTM\data\dataTodatabase.py

将药物相关的信息储存到数据库中

[{'_index': 'drug',
    '_id': 'UniprotId',
    '_source': {
        'UniprotId': 'Q9Y6K9',
        'UniprotName': 'NR4A1_HUMAN',
        'geneName': 'NR4A1',
        'target': {'T00254':'research'},
        'drug': {
            {'KOS-1815': 'Preclinical',
            'RG7800': 'Phase 1/2'
            } 需要target
        }
    },
    {},{}]
'''
from elasticsearch import Elasticsearch
from elasticsearch import helpers
from bs4 import BeautifulSoup
import pickle
import requests
import re
import os

# uniPortId 转化为geneName
def uniPortIdToGeneName(uniPortId):
    url = 'https://rest.uniprot.org/uniprotkb/' + uniPortId +'.txt'
    #r = requests.get('https://rest.uniprot.org/uniprotkb/Q9Y2T7.txt')
    r = requests.get(url)
    data = r.text
    # 获取geneName
    if data != '':
        geneName = data.split('GN   Name=')[1].split(';')[0]
        # 去除{}中的内容
        geneName = re.sub(r'\{.*?\}', '', geneName)
    else:
        geneName = ''
    return geneName

# geneName 转化为uniprotId
# 输入的格式需要geneName和uniprotName之间使用\t分割(tab键)
def geneNameToUniprotId(file_name):
    with open(file_name, 'r', encoding='utf-8') as f:
        f = f.readlines()
    data = {}
    for i in f:
        i = i.replace('\t', ' ').replace('\n', '').split(' ')
        data[i[0]] = i[1]
    return data

# uniprotName 转化为uniprotId
def uniprotNameToId():
    with open(r'data\UniProt_ID_2_Name.hao', 'r', encoding='utf-8') as f:
        f = f.readlines()
    data = {}
    for i in f:
        i = i.replace('\t', '').replace('\n', '').split('  ')
        data[i[1]] = i[0]
    return data

# uniprotId 转化为uniprotName
def uniprotIdToName():
    with open(r'data\UniProt_ID_2_Name.hao', 'r', encoding='utf-8') as f:
        f = f.readlines()
    data = {}
    for i in f:
        i = i.replace('\t', '').replace('\n', '').split('  ')
        data[i[0]] = i[1]
    return data

# uniprotName 转化为靶标
def uniprotNameToTarget():
    with open(r'data\TTD_phase_List.hao', 'r', encoding='utf-8') as f:
        f = f.readlines()
    data = {}
    mRNA_target = mRNATarget()
    for i in f:
        i = i.replace('\t', '  ').replace('\n', '').split('   ')# 根据文件中的数据格式进行分割
        # 根据情况不同添加字典的键值对
        if i[1] not in data.keys() and i[0] not in mRNA_target:
            data[i[1]] = {}
            data[i[1]][i[0]] = i[2]
    return data

# 将靶标转化为药物
def targetToDrug():
    with open(r'data\Drug_Target_Mapping.hao', 'r', encoding='utf-8') as f:
        f = f.readlines()
    data = {}
    for i in f:
        i = i.replace('|||','  ').replace('\n', '').split('  ')# 根据文件中的数据格式进行分割
        # 根据情况不同添加字典的键值对
        if i[0] not in data.keys():
            data[i[0]] = {}
            data[i[0]][i[1]] = i[2]
        else:
            data[i[0]][i[1]] = i[2]
    return data

# 一个蛋白应该有且只对应一个target，多个target的蛋白需要进行筛选，将其中是对应的mRna的target删除
def mRNATarget():
    with open(r'data\TTD_uniprot_all.txt', 'r', encoding='utf-8') as f:
        f = f.readlines()
    target = []
    for line in f:
        if 'mRNA' in line:
            target.append(line.split('\t')[0])
    return target

# 获得txt结尾的文件
def get_txt():
    for file_name in os.listdir():
        if file_name.endswith('.txt'):
            return file_name
        
# 获得输入文件的uniPortName
def get_file_uniPortName():
    # 获得所有输入的gene的名字
    geneName = list(geneName_to_uniprotId.keys())
    # 获得所有输入的gene的uniprotId
    uniprotId = []
    for gene in geneName:
        uniprotId.append(geneName_to_uniprotId[gene])
    # 获得所有输入的gene的uniprotName
    uniprotName = []
    uniprotId_to_name = uniprotIdToName()
    for i in uniprotId:
        uniprotName.append(uniprotId_to_name[i])
    return uniprotName      

# 获得输出的数据
def get_output(uniprotName):
    target_have_drug = []
    target_FDA_approved = []
    target_research = []
    target_clinical_trial = []
    for name in uniprotName:
        if name in uniPortName_to_target.keys():
            if ('Successful target' in [*uniPortName_to_target[name].values()]
                ) or (
                'Research target' in [*uniPortName_to_target[name].values()]
                ) or (
                'Clinical Trial target' in [*uniPortName_to_target[name].values()]
            ):
                target_have_drug.append([*uniPortName_to_target[name].keys()][0])
            if 'Successful target' in [*uniPortName_to_target[name].values()]:
                target_FDA_approved.append([*uniPortName_to_target[name].keys()][0])
            if 'Research target' in [*uniPortName_to_target[name].values()]:
                target_research.append([*uniPortName_to_target[name].keys()][0])
            if 'Clinical Trial target' in [*uniPortName_to_target[name].values()]:
                target_clinical_trial.append([*uniPortName_to_target[name].keys()][0])
    target_have_drug = len(target_have_drug)
    target_FDA_approved = len(target_FDA_approved)
    target_research = len(target_research)
    target_clinical_trial = len(target_clinical_trial)   
    return target_have_drug, target_FDA_approved, target_research, target_clinical_trial

# 获得转化数据
def get_translate_data():
    uniPortName_to_uniPortId = uniprotNameToId() # uniprotId 转化为 uniprotName
    uniPortId_to_uniPorName = uniprotIdToName() # uniprotName 转化为 uniprotId
    uniPortName_to_target = uniprotNameToTarget() # uniprotName 转化为 靶标
    target_to_drug = targetToDrug() # 将 靶标 转化为 药物
    with open(r'data\uniPortIdToGeneName.json','rb') as fp:
        uniprotId_to_geneName = pickle.load(fp)# uniprotId 转化为 geneName
    return (uniPortName_to_uniPortId, uniPortId_to_uniPorName, 
            uniPortName_to_target, target_to_drug, 
            uniprotId_to_geneName)
    
# 将html中的数据进行更改后输出
def output_html():
    text_html = open(r'template.html', 'r', encoding='utf-8').read()
    text_html = text_html.replace(
                    'Compound data',str(target_have_drug)).replace(
                        'No-drug data',str(target_no_drug)).replace(
                            'FDA Approved data',str(target_FDA_approved)).replace(
                                'Research data',str(target_research)).replace(
                                    'Clinical data',str(target_clinical_trial))

    soup = BeautifulSoup(text_html, 'html.parser')
    with open(r'output.html', 'w', encoding='utf-8') as fp:
        fp.write(str(soup))

# 将所有蛋白和靶标相关的信息为elasticSearch的数据格式
def elasticSearchData():
    data = [] 
    for uniprotName in [*uniPortName_to_uniPortId.keys()]:
        protein = {'_index': 'uniprot',
                    '_id': '',
                    '_source': {
                        'uniPortId': '',
                        'uniPortName': '',
                        'geneName': '',
                        'target_name': '',
                        'target_phase': '',
                    }}
        protein['_source']['uniPortName'] = uniprotName
        protein['_source']['uniPortId'] = uniPortName_to_uniPortId[uniprotName]
        # 根据uniprotId作为数据库的索引
        protein['_id'] = uniPortName_to_uniPortId[uniprotName]
        # 添加geneName
        if protein['_source']['uniPortId'] in uniprotId_to_geneName.keys():
            protein['_source']['geneName'] = uniprotId_to_geneName[protein['_source']['uniPortId']]
        # 添加target
        if  uniprotName in uniPortName_to_target.keys():
            protein['_source']['target_name'] = [*uniPortName_to_target[uniprotName].keys()][0]
            protein['_source']['target_phase'] = [*uniPortName_to_target[uniprotName].values()][0]
        # 添加drug
        if  protein['_source']['target_name'] in target_to_drug.keys():
            drug_name_list = []
            drug_state_list = []
            for drug_name,drug_state in target_to_drug[protein['_source']['target_name']].items():
                drug_name_list.append(drug_name)
                drug_state_list.append(drug_state)
            protein['_source']['drug_name'] = drug_name_list
            protein['_source']['drug_state'] = drug_state_list
        data.append(protein)
        # 链接es数据库，并创建index为uniprot的数据库
        es = Elasticsearch(timeout=30, max_retries=10, retry_on_timeout=True)
        es.indices.create(index='uniprot', ignore=400)#
        # 将数据写入数据库
        helpers.bulk(es, data)   
''' 
处理数据需要从网上进行操作 2w个需要5个小时不推荐使用 但是可以使用本地保存的数据uniprotId_to_geneName.json
uniprotId_to_geneName = get_geneID()

'''
if __name__ == '__main__':
    (uniPortName_to_uniPortId, 
        uniPortId_to_uniPorName, 
        uniPortName_to_target,
        target_to_drug, 
        uniprotId_to_geneName ) = get_translate_data()  

    geneName_to_uniprotId = geneNameToUniprotId(get_txt()) # geneName 转化为 uniprotId
    file_uniprotName = get_file_uniPortName() # 获得输入文件的uniPortName
    
    target_have_drug, target_FDA_approved, target_research, target_clinical_trial = get_output(file_uniprotName)
    target_no_drug = len(file_uniprotName) - target_have_drug
    
    output_html()






