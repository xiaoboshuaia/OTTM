'''
Date: 2022-10-31 17:40:44
LastEditTime: 2022-11-12 15:41:05
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

import textwrap
from pyecharts.charts import Sunburst ,Tree ,Bar ,Grid ,Page
from pyecharts import options as opts
from operator import index
import traceback
from elasticsearch import Elasticsearch
from elasticsearch import helpers
from bs4 import BeautifulSoup
import pickle
import requests
import re
import os

# uniPortId 转化为geneName
def uniPortIdToGeneName(uniPortId):
    url = 'https://rest.uniprot.org/uniprotkb/' + uniPortId + '.txt'
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
    target_to_uniprotName = {}
    mRNA_target = mRNATarget()
    for i in f:
        i = i.replace('\t', '  ').replace(
            '\n', '').split('   ')  # 根据文件中的数据格式进行分割
        # 根据情况不同添加字典的键值对
        if i[1] not in data.keys() and i[0] not in mRNA_target:
            data[i[1]] = {}
            data[i[1]][i[0]] = i[2]
            target_to_uniprotName[i[0]] = i[1]
    return data, target_to_uniprotName

# 将靶标转化为药物
def targetToDrug():
    with open(r'data\Drug_Target_Mapping.hao', 'r', encoding='utf-8') as f:
        f = f.readlines()
    data = {}
    for i in f:
        i = i.replace('|||', '  ').replace(
            '\n', '').split('  ')  # 根据文件中的数据格式进行分割
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
    target_clinical_trial = []
    for name in uniprotName:
        if name in uniPortName_to_target.keys():
            target_have_drug.append([*uniPortName_to_target[name].keys()][0])
            if 'Successful target' in [*uniPortName_to_target[name].values()]:
                target_FDA_approved.append(
                    [*uniPortName_to_target[name].keys()][0])
            if 'Clinical Trial target' in [*uniPortName_to_target[name].values()]:
                target_clinical_trial.append(
                    [*uniPortName_to_target[name].keys()][0])
    return target_have_drug, target_FDA_approved,  target_clinical_trial

# 获得转化数据
def get_translate_data():
    uniPortName_to_uniPortId = uniprotNameToId()  # uniprotId 转化为 uniprotName
    uniPortId_to_uniPorName = uniprotIdToName()  # uniprotName 转化为 uniprotId
    uniPortName_to_target, target_to_uniprotName = uniprotNameToTarget()  # uniprotName 转化为 靶标
    target_to_drug = targetToDrug()  # 将 靶标 转化为 药物
    with open(r'data\uniPortIdToGeneName.json', 'rb') as fp:
        uniprotId_to_geneName = pickle.load(fp)  # uniprotId 转化为 geneName
    return (uniPortName_to_uniPortId, uniPortId_to_uniPorName,
            uniPortName_to_target, target_to_drug,
            uniprotId_to_geneName, target_to_uniprotName)

# 将html中的数据进行更改后输出
def output_html():
    text_html = open(r'Template/target_pie_template.html',
                    'r', encoding='utf-8').read()
    text_html = text_html.replace(
        'Compound data', str(len(target_have_drug))).replace(
        'No-drug data', str(target_no_drug)).replace(
        'FDA Approved data', str(len(target_FDA_approved))).replace(
        'Others data', str(target_others)).replace(
        'Clinical data', str(len(target_clinical_trial)))

    soup = BeautifulSoup(text_html, 'html.parser')
    with open(r'Target/targets_info.html', 'w', encoding='utf-8') as fp:
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
        if uniprotName in uniPortName_to_target.keys():
            protein['_source']['target_name'] = [
                *uniPortName_to_target[uniprotName].keys()][0]
            protein['_source']['target_phase'] = [
                *uniPortName_to_target[uniprotName].values()][0]
        # 添加drug
        if protein['_source']['target_name'] in target_to_drug.keys():
            drug_name_list = []
            drug_state_list = []
            for drug_name, drug_state in target_to_drug[protein['_source']['target_name']].items():
                drug_name_list.append(drug_name)
                drug_state_list.append(drug_state)
            protein['_source']['drug_name'] = drug_name_list
            protein['_source']['drug_state'] = drug_state_list
        data.append(protein)
        # 链接es数据库，并创建index为uniprot的数据库
        es = Elasticsearch(timeout=30, max_retries=10, retry_on_timeout=True)
        es.indices.create(index='uniprot', ignore=400)
        # 将数据写入数据库
        helpers.bulk(es, data)


''' 
处理数据需要从网上进行操作 2w个需要5个小时不推荐使用 但是可以使用本地保存的数据uniprotId_to_geneName.json
uniprotId_to_geneName = get_geneID()

'''
'''
通过输入一个target的名称,返回该target的所有相关的靶标数据和药物数据

'''

# 通过循环来进行查询 同时将uniprotId name 以及geneName 作为检索的关键字
def query(uni_ID, uni_Name):
    gene_name = uniprotId_to_geneName[uni_ID]
    body = {
        'query': {
            'query_string': {
                'query': uni_ID + ' OR ' + uni_Name + ' OR ' + gene_name
                + ' AND hepatocellular carcinoma ',
                'fields': ['abstract']
            },
        }
    }
    res = es.search(index='abstract', body=body, scroll='5m')
    hit_number = res['hits']['total']['value']
    '''
    在翻页查询完成后 释放scroll 避免内存溢出
    '''
    es.clear_scroll(scroll_id=res['_scroll_id'])
    return hit_number

''' 
在已知的abstractid的abstract中进行gene和肝癌相关性的查询
只要摘要中出现肝癌 就说明两者存在相关性
'''
def query2(uni_ID, uni_Name):
    gene_name = uniprotId_to_geneName[uni_ID].replace(
        ' ', '').replace(
            '\n', ''
    ).replace(
            '\r', ''
    )
    gene_name_re = re.sub(r'\{.*?}', '', gene_name)
    print(gene_name_re)
    pubMedId = gene_to_pubMedID[gene_name_re]
    sql = {
        'query': {
            'bool': {
                'must': [
                    {
                        'terms': {
                            'pubMedId': pubMedId
                        }
                    },
                    {
                        'query_string': {
                            'query': 'hepatocellular AND carcinoma',
                            'fields': ['abstract']
                        }
                    },
                ]
            }
        }
    }
    res = es.search(index='abstract', body=sql, scroll='5m')
    hit_number = res['hits']['total']['value']
    es.clear_scroll(scroll_id=res['_scroll_id'])

    abs_id = []
    for i in res['hits']['hits']:
        abs_id.append(i['_source']['pubMedId'])

    return hit_number

# 获得FDA批准的靶标并且没有研究过肝癌
def FDA_approved_drug():
    target_FDA_name = []
    for key, value in uniPortName_to_target.items():
        if [*value.keys()][0] in target_FDA_approved:
            target_FDA_name.append(key)

    re_target = []
    unre_target = []
    for i in target_FDA_name:
        print(i)
        number = query2(uniPortName_to_uniPortId[i], i)
        if number == 0:
            unre_target.append(i)
        else:
            re_target.append(i)
    dict1 = []
    dict2 = []
    for i in unre_target:
        dict1.append({'name': i})
    for i in re_target:
        dict2.append({'name': i})
    return len(unre_target + re_target), unre_target, dict1, dict2

# 获得Clinical trial阶段的靶标且没有研究过肝癌
def Clinical_trial_drug():
    target_clinical_name = []
    for key, value in uniPortName_to_target.items():
        if [*value.keys()][0] in target_clinical_trial:
            target_clinical_name.append(key)

    cl_re_target = []
    cl_unre_target = []
    for i in target_clinical_name:
        print(i)
        number = query2(uniPortName_to_uniPortId[i], i)
        if number == 0:
            cl_unre_target.append(i)
        else:
            cl_re_target.append(i)

    dict1 = []
    dict2 = []
    for i in cl_unre_target:
        dict1.append({'name': i})
    for i in cl_re_target:
        dict2.append({'name': i})
    return len(cl_unre_target + cl_re_target), cl_unre_target, dict1, dict2

# 通过phase将药物进行分类
def drug_classify(target_name):
    drug_phase = {'Approved': [], 'Clinical_trial': [], 'Others': []}
    # 'T00140'
    drug_ap_cl = []
    drug_ap = []
    drug_cl = []
    for key, value in target_to_drug[target_name].items():
        if value == 'Approved':
            drug_phase['Approved'].append({
                'name': key
            })
            drug_ap_cl.append(key)
            drug_ap.append(key)
        elif value.startswith('Phase'):
            drug_phase['Clinical_trial'].append({
                'name': key
            })
            drug_ap_cl.append(key)
            drug_cl.append(key)
        else:
            drug_phase['Others'].append({
                'name': key
            })
    return drug_phase, drug_ap_cl, drug_ap, drug_cl

# 通过循环来进行查询药物是否有研究过肝癌 若没有则返回数字0
def queryDrug(drug_name):
    drug_name = drug_name.replace(
        '+/-', ' ').replace(
            '/', ' ').replace(
                '[', '').replace(
                    ']', '').replace(
                        '-', ' ')
    body = {
        'query': {
            'query_string': {
                'query':  drug_name + ' AND hepatocellular AND carcinoma',
                'fields': ['abstract'],
                'default_operator': 'AND',  # 这里很关键将默认的OR改为AND可以明显的提升查询的质量
            }
        }
    }
    res = es.search(index='abstract', body=body, scroll='5m')
    hit_number = res['hits']['total']['value']
    es.clear_scroll(scroll_id=res['_scroll_id'])
    return hit_number

# 将药物通过approved和clinical和是否有研究过肝癌进行分类
def drug_classify_final():
    drug_classify_final = {'Approved': {'Not relevant': [], 'Relevant': []},
                            'Clinical_trial': {'Not relevant': [], 'Relevant': []}}
    ok_drug = []
    no_drug = []
    for i in drug_ap_cl:
        number = queryDrug(i)
        # print(number)
        if number == 0:
            ok_drug.append(i)
        else:
            no_drug.append(i)

    for i in drug_ap_cl:
        if i in ok_drug and i in drug_ap:
            drug_classify_final['Approved']['Not relevant'].append(i)
        elif i in ok_drug and i in drug_cl:
            drug_classify_final['Clinical_trial']['Not relevant'].append(i)
        elif i in no_drug and i in drug_ap:
            drug_classify_final['Approved']['Relevant'].append(i)
        elif i in no_drug and i in drug_cl:
            drug_classify_final['Clinical_trial']['Relevant'].append(i)
    return drug_classify_final

# 通过循环来进行查询药物频率(热度)
def query_Drug_frequency(drug_name):
    drug_name = drug_name.replace('+/-', ' ').replace('/', ' ')
    body = {
        'query': {
            'query_string': {
                'query':  drug_name,
                'fields': ['abstract']
            }
        }
    }
    res = es.search(index='abstract', body=body, scroll='5m')
    hit_number = res['hits']['total']['value']
    es.clear_scroll(scroll_id=res['_scroll_id'])
    return hit_number

# 获得所有未研究过肝癌的药物
def get_drug_frequency():
    drug_classify_final_data = drug_classify_final()
    drug_a_not_relevant = []
    drug_c_not_relevant = []
    drug_not_relevant = []
    for i in drug_classify_final_data['Approved']['Not relevant']:
        drug_a_not_relevant.append(i)
        drug_not_relevant.append(i)
    for i in drug_classify_final_data['Clinical_trial']['Not relevant']:
        drug_c_not_relevant.append(i)
        drug_not_relevant.append(i)
    drug_frequency = []
    # 有的靶标对应的药物全部都报告过跟肝癌有关 这里需要进行判断
    # 如果全都报道过 则从报道过的药物中 来进行排名和热度的查询
    if drug_not_relevant == []:
        for i in drug_ap_cl:
            frequency = query_Drug_frequency(i)
            drug_frequency.append(frequency)
    else:
        for i in drug_not_relevant:
            frequency = query_Drug_frequency(i)
            drug_frequency.append(frequency)

    return drug_a_not_relevant, drug_c_not_relevant, drug_not_relevant, drug_frequency

# 将数据处理成为echarts所需要的格式
def drug_frequency_list():
    data = []
    for i in range(len(drug_frequency)):
        data.append([drug_not_relevant[i], drug_frequency[i]])
    return data

# 通过genename 寻找pubmedid
def gene_to_pubMedID():
    with open(r'data/Symbol_To_PubMed_ID.txt') as file:
        file = file.read().splitlines()
    gene_to_pubMedID = {}
    for i in file:
        i = i.split('|||')
        gene_to_pubMedID[i[0]] = [j for j in i[1].split(',') if j != 'XXXXX']

    return gene_to_pubMedID

# 将药物转化为tree图需要的数据类型
def drug_treetype_data(drugs):
    drug = []
    for i in drugs:
        drug.append({'name': wrap_text(i)})
    return drug

# 将字符串处理，过长的字符串换行
def wrap_text(text, max_length=20):
    if len(text) > max_length:
        return textwrap.fill(text, max_length)
    return text

# 生成靶标推荐药物的Sunburst图
def all_targets_tree():
    data = [
    {
        "children": [
            
            {
                "children": [
                    {
                        "children": fda_targets_unre,
                        "name": "Not Relevant",
                        'value': len(fda_targets_unre)
                    },
                    {
                        "children":fda_targets_re,
                        "name": "Relevant",
                        'value': len(fda_targets_re)
                    }],
                "name": "FDA Approved",
                'value': fda_all_nmb,    
            },
            {
                "children": [
                    {
                        "children":clinical_targets_unre,
                        "name": "Not Relevant",
                        'value': len(clinical_targets_unre)
                    },
                    {
                        "children":clinical_targets_re,
                        "name": "Relevant",
                        'value': len(clinical_targets_re)
                    }
                ],
                "name": "Clinical",
                'value': cli_all_nmb
            },
        ],
        "name": "Target",
        'value': fda_all_nmb + cli_all_nmb
    }
]
    
    c = (
    Tree(init_opts=opts.InitOpts(width="1200px", height="800px"))
    .add("", data, collapse_interval=2,
        symbol_size=10,
        symbol="emptyCircle",
        leaves_label_opts=opts.LabelOpts(position="right"),
        itemstyle_opts=opts.ItemStyleOpts(border_width=1, border_color="#48466d"),
        )
    .set_global_opts(title_opts=opts.TitleOpts(title="All Targets"))
    .set_series_opts(label_opts=opts.LabelOpts(
                                            font_size=20,
                                            font_weight='bold',
                                            color="#48466d"
                                            ),
                    )
    
    .render( 'Target/' + "ALL_Targets.html")
)

# 生成每个靶标对应药物的树状图以及柱状图
def target_tree_bar():
    tree_data = [
    {
        "children": [
            
            {
                "children": [
                    {
                        "children": drug_a_not_relevant,
                        "name": "Not Relevant",
                        'value': len(drug_a_not_relevant)
                    },
                    {
                        "children":drug_a_relevant,
                        "name": "Relevant",
                        'value': len(drug_a_relevant)
                    }],
                "name": "FDA Approved",
                'value': len(drug_ap),    
            },
            {
                "children": [
                    {
                        "children":drug_c_not_relevant,
                        "name": "Not Relevant",
                        'value': len(drug_c_not_relevant)
                    },
                    {
                        "children":drug_c_relevant,
                        "name": "Relevant",
                        'value': len(drug_c_relevant)
                    }
                ],
                "name": "Clinical",
                'value': len(drug_cl)
            },
        ],
        "name": target,
        'value': len(drug_ap_cl)
    }
]

    tree = (
    Tree(init_opts=opts.InitOpts(width="1600px", height="800px"))
    .add("", tree_data, collapse_interval=2,
        symbol_size=10,
        symbol="emptyCircle",
        leaves_label_opts=opts.LabelOpts(position="right"),
        itemstyle_opts=opts.ItemStyleOpts(border_width=1, border_color="#48466d"),
        edge_fork_position="100%",
        )
    .set_global_opts(title_opts=opts.TitleOpts(title=target + ' corresponding drugs'))
    .set_series_opts(label_opts=opts.LabelOpts(
                                            font_size=15,
                                            font_weight='bold',
                                            color="#48466d"
                                            ),
                    )
)
    if drug_not_relevant == []:
        drug_data = drug_ap_cl
    else:
        drug_data = drug_not_relevant
        
    bar = (
    Bar(init_opts=opts.InitOpts(width="2000px", height="800px"))
    .add_xaxis(drug_data)
    .add_yaxis("Drug Frequency", drug_frequency, color='#617bdb' )
    .set_global_opts(
        xaxis_opts=opts.AxisOpts(
                                is_show=False,
                                axislabel_opts=opts.LabelOpts(font_size=15,
                                        font_weight='bold',
                                        color="#48466d",
                                        ),
),
        yaxis_opts=opts.AxisOpts(
            
            axislabel_opts=opts.LabelOpts(font_size=15,
                                        font_weight='bold',
                                        color="#48466d"),
            
        ),
        legend_opts=opts.LegendOpts(is_show=False),
        
    )
    .set_series_opts(
        label_opts=opts.LabelOpts(position="right",
                                color="#48466d",
                                font_size=15,
                                font_weight='bold',),
    )
    .reversal_axis()
    )
    
    (
        Page()
        .add(tree,bar)
        ).render(
            'Target/' + target + '/'+ target + '.html'
            )

# 制作sunburst图
def get_sunburst(un_relevant_targets_recommend_drug):
    
    data2 = [
        {
            'name': 'FDA approve',
            "itemStyle": {"color": '#fac858'},
            'children': []

        },
            {
            'name': 'Clinical trial',
            "itemStyle": {"color": '#73c0de'},
            'children': []

        },
    ]

    for key, value in un_relevant_targets_recommend_drug.items():
        if key in fda_targets_li:
            
            children = {
                "name": key,
                'value': 1,
                "itemStyle": {"color":'#fac858' },
                'children': [
                        {'name': value,
                        'value': 1,
                        "itemStyle": {"color":'#fac858' },
                        }
                ]
            }
            data2[0]['children'].append(children)    
        else:
            children = {
                "name": key,
                'value': 1,
                "itemStyle": {"color": '#73c0de'},
                'children': [
                        {'name': value,
                        'value': 1,
                        "itemStyle": {"color": '#73c0de'}
                        }
                ]
            }
            data2[1]['children'].append(children)


    c = (
        Sunburst(init_opts=opts.InitOpts(width="1200px", height="1200px"))
        .add(
            "",
            data_pair=data2,
            highlight_policy="ancestor",
            sort_="null",
            radius=[0, "95%"],
            # center=["55%", "55%"],
            # 居中
            
            levels=[
                {},
                {
                    "r0": "15%",
                    "r": "35%",
                    "itemStyle": {"borderWidth": 2},
                    "label": {"rotate": "tangential",},
                },
                {"r0": "35%", "r": "60%", "label": {"align": "right"}},
                {
                    "r0": "60%",
                    "r": "62%",
                    "label": {"position": "outside", "padding": 3, "silent": False},
                    "itemStyle": {"borderWidth": 3},
                },
            ],
        )
        .set_global_opts(title_opts=opts.TitleOpts(title="Suggestions for drug targets",
                                                    pos_left='center',
                                                    ))
        .set_series_opts(label_opts=opts.LabelOpts(formatter="{b}",
                                                color='black',
                                                font_weight='bold',
                                                font_size=15,
                                                font_family='Microsoft YaHei',                   
                                                ))
        
        .render("Target/drug_suggestion.html")
    )

if __name__ == '__main__':
    (   uniPortName_to_uniPortId,
        uniPortId_to_uniPorName,
        uniPortName_to_target,
        target_to_drug,
        uniprotId_to_geneName,
        target_to_uniprotName) = get_translate_data()

    gene_to_pubMedID = gene_to_pubMedID()
    geneName_to_uniprotId = geneNameToUniprotId(
        get_txt())  # geneName 转化为 uniprotId
    file_uniprotName = get_file_uniPortName()  # 获得输入文件的uniPortName

    target_have_drug, target_FDA_approved,  target_clinical_trial = get_output(
        file_uniprotName)
    target_no_drug = len(file_uniprotName) - len(target_have_drug)
    target_others = len(target_have_drug) - \
        len(target_FDA_approved) - len(target_clinical_trial)

# 获得靶标的数量信息
    output_html()


# hepatocellular carcinoma(HCC)
# 将有药物的target对于早期肺癌在进行一次筛选 找出没有研究过跟该疾病有关的靶标
    es = Elasticsearch(timeout=30, max_retries=10, retry_on_timeout=True)

# 获取没有报道过肝癌的靶标和数量 报道过肝癌的靶标
    fda_all_nmb, fda_targets_li, fda_targets_unre, fda_targets_re = FDA_approved_drug()
    cli_all_nmb, clinical_targets_li, clinical_targets_unre, clinical_targets_re = Clinical_trial_drug()


# 获取所有没有被报道过的靶标
    un_relevant_uniportname = fda_targets_li + clinical_targets_li
    un_relevant_targets = []

    for uniport_name in un_relevant_uniportname:
        # 获取基因名对应的靶标名字
        target_name = [*uniPortName_to_target[uniport_name].keys()][0]
        un_relevant_targets.append(target_name)
        

# 获得所有没有报道过跟肝癌有关的靶标 及其推荐的药物
    un_relevant_targets_recommend_drug = {}
    for target in un_relevant_targets:
        '''
        drug_phase 药物的研究阶段
        drug_ap_cl 药物处于FDA批准和临床试验的阶段
        drug_ap 药物处于FDA批准的阶段
        drug_cl 药物处于临床试验的阶段
        drug_a_not_relevant FDA批准的药物中没有报道过肝癌的药物
        drug_c_not_relevant 临床试验的药物中没有报道过肝癌的药物
        drug_not_relevant FDA批准和临床试验的药物中没有报道过肝癌的药物
        drug_frequency 没有报道过肝癌药物的药物频率
        '''
        drug_phase, drug_ap_cl, drug_ap, drug_cl = drug_classify(target)
        drug_a_not_relevant, drug_c_not_relevant, drug_not_relevant, drug_frequency = get_drug_frequency()
        
        number_index = drug_frequency.index(
            max(drug_frequency))  # 获取最大值的索引
        if drug_not_relevant == []:
            drug_suggest = drug_ap_cl[number_index]
        else:
            drug_suggest = drug_not_relevant[number_index]  # 获取热度最高对应的药物名字
        print(target, drug_suggest)

        uniprotName = target_to_uniprotName[target]

        un_relevant_targets_recommend_drug[uniprotName] = drug_suggest

    get_sunburst(un_relevant_targets_recommend_drug) # 所有靶标的药物推荐
    all_targets_tree() # 所有靶标具体信息

# 生成单个靶标的药物推荐
    for target in un_relevant_targets:
        drug_phase, drug_ap_cl, drug_ap, drug_cl = drug_classify(target)
        drug_a_not_relevant, drug_c_not_relevant,drug_not_relevant, drug_frequency = get_drug_frequency()
        
        drug_c_relevant = list(set(drug_cl) - set(drug_c_not_relevant))
        drug_a_relevant = list(set(drug_ap) - set(drug_a_not_relevant))
        
        # 数据类型转化
        drug_ap_cl = [wrap_text(i) for i in drug_ap_cl]
        if drug_not_relevant != []:
            drug_not_relevant = [wrap_text(i) for i in drug_not_relevant]
        
        drug_ap = drug_treetype_data(drug_ap)
        drug_cl = drug_treetype_data(drug_cl)
        drug_c_not_relevant = drug_treetype_data(drug_c_not_relevant)
        drug_a_not_relevant = drug_treetype_data(drug_a_not_relevant)
        drug_a_relevant = drug_treetype_data(drug_a_relevant)
        drug_c_relevant = drug_treetype_data(drug_c_relevant)
        
        target = target_to_uniprotName[target]
        os.makedirs('Target/'+target)
        target_tree_bar()
    print('The program is over')







































# %%
