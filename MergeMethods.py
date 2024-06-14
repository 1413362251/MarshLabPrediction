from lxml import etree
import pandas as pd
from methods import *
import re

def clinvar_big_orignal_sep(Path_clinvarfile,GeneName,OutPath_clinvar_uniprot):
    df1 = pd.read_csv(Path_clinvarfile, sep='\t', skiprows=19)
    df_gene = df1[df1['gene'] == GeneName]
    df_gene.to_excel(OutPath_clinvar_uniprot, index = False )

#处理df
def simple_to_full(text):

    text_new = text.split('.')[-1]
    if not any(char.islower() for char in text_new):
        parts = re.split(r'(\d+)', text_new)
        original = convert_sequence(parts[0], conversion_type='1to3')
        result = convert_sequence(parts[-1], conversion_type='1to3')
        final_str = f'p.{original}{parts[1]}{result}'
        return final_str
    else:
        return text

def getTagDetails(feature,ns,tagname):
    tags = feature.findall(f'.//ns:{tagname}', namespaces=ns)
    tagList = []#({},text)
    for tag in tags:
        # 提取所有属性
        attribute = tag.attrib
        text = tag.text.strip() if tag.text is not None else 'N/A'
        tagList.append((attribute,text))
        # 打印提取的数据
    return tagList

def getdbReference(feature,ns):
    list = getTagDetails(feature,ns,'dbReference')
    list_only_attribute = [i[0] for i in list]
    return  list_only_attribute


def print_element_details(element, indent=0):
    # 打印当前元素的标签和属性
    #print(' ' * indent + f"Tag: {element.tag.replace('{http://www.ebi.ac.uk/proteins/api/doc/xsd/feature}', 'ns:')}")
    #print(' ' * indent + f"Attributes: {element.attrib}")
    if element.text is not None and element.text.strip():
        print()
        #print(' ' * indent + f"Text: {element.text.strip()}")

    # 递归打印所有子元素
    for child in element:
        print()
        #print_element_details(child, indent + 4)



#initilise a df
def main_Uniprot(XML_filepath,outputpath):
    data = {
        'Position': [],
        'HGVS_P': [],
        'TranscriptId': [],
        'HGVS_D': [],
        'NCReseq': [],
        'CytogeneticBand': [],
        'ConsequenceType': [],
        'VariantClinicalSignificance': [],
        'Phenotype': [],
        'References': []
    }

    df = pd.DataFrame(data)


    # 读取XML文件

    tree = etree.parse(XML_filepath)
    root = tree.getroot()

    # 定义命名空间
    ns = {'ns': 'http://www.ebi.ac.uk/proteins/api/doc/xsd/feature'}
    # 获取前十个 feature 元素
    top_ten_features = root.findall('.//ns:feature[@type="Variant"]', namespaces=ns)
    #遍历是10个元素
    ##printDetail(1)


    x = 1
    for feature in top_ten_features:
        #print('Element:',x)
        x+=1
        if getTagDetails(feature,ns,'name'):
            phenotypes = getTagDetails(feature,ns,'name')
            phenotype_list = []
            for i in phenotypes:
                phenotype_list.append(i[1])
            phenotype = '|'.join(phenotype_list)
        else:
            phenotype = 'not Provided'
        #print( phenotype)




        references = getdbReference(feature,ns)
        if getTagDetails(feature,ns,'position'):
            position = getTagDetails(feature,ns,'position')[0][0]['position']
        else:
            begin= getTagDetails(feature,ns,'begin')[0][0]['position']
            end= getTagDetails(feature,ns,'end')[0][0]['position']
            if begin == end:
                position = begin
            else:
                position = f'{begin}-{end}'


        for location in getTagDetails(feature,ns,'variantLocation'):
            if location[0]['loc'].startswith('p'):
                HGVS_P =location[0]['loc']
                TranscriptId = location[0]['seqId']

        if getTagDetails(feature,ns,'genomicLocation'):
            NCReseq,HGVS_D= getTagDetails(feature,ns,'genomicLocation')[0][1].split(':')
        else:
            NCReseq, HGVS_D=None,None
        cytogeneticBand =getTagDetails(feature,ns,'cytogeneticBand')[0][1]
        consequenceType = getTagDetails(feature,ns,'consequenceType')[0][1]

        if getTagDetails(feature,ns,'variantClinicalSignificance'):
            variantClinicalSignificance = getTagDetails(feature,ns,'variantClinicalSignificance')[0][0]['type']
        else:
            variantClinicalSignificance=''

        prediction = ''

        #print(position)
        #print(HGVS_P)
        #print(TranscriptId)
        #print(HGVS_D)
        #print(NCReseq)
        #print(cytogeneticBand)
        #print(consequenceType)
        #print(variantClinicalSignificance)
        #print(prediction)
        #print(references)
        # 创建新行的DataFrame
        new_row = pd.DataFrame({
            'Position': [position],
            'HGVS_P': [HGVS_P],
            'TranscriptId': [TranscriptId],
            'HGVS_D': [HGVS_D],
            'NCReseq': [NCReseq],
            'CytogeneticBand': [cytogeneticBand],
            'ConsequenceType': [consequenceType],
            'VariantClinicalSignificance': [variantClinicalSignificance],
            'Phenotype': [phenotype],
            'References': [references]
        })

        # 使用concat将新行添加到现有DataFrame
        df = pd.concat([df, new_row], ignore_index=True)

    #printDetail(1363)



    df['HGVS_P'] = df['HGVS_P'].apply(simple_to_full)

    #精确访问多少列的多少行
    #print(df.iloc[176] )

    output_path = outputpath
    #df.to_csv(r"D:\file\毕业设计\Data\Uniprot_data.xlsx",sep="\t",header=True)
    df.to_excel(output_path, index=False)



###Uniprot_Clinvar
def nomalizeprotein(str):
    start = convert_sequence(str[0])
    end = convert_sequence(str[-1])
    position = str[1:-1]
    #print(start,end,position)
    new_HGVS = f'p.{start}{position}{end}'
    return new_HGVS

def getHGVS_p(str):
    lst = str.split()
    HGVS_p = lst[-1].strip()[1:-1]
    return HGVS_p

# 定义一个提取数字的函数
def extract_position(text):
    # 分割字符串并取最后一个元素
    last_part = text.split()[-1].strip()
    # 使用正则表达式查找数字
    match = re.search(r'\d+', last_part)
    if match:
        return int(match.group())  # 返回找到的数字
    return None  # 如果没有找到数字，返回None
def extract_transcriptId(text):
    first_part = text.split('(')[0].strip()
    ##print(text.split('('))
    return first_part
def extract_NCReseq(text):
    first_part = text.split(':')[0].strip()
    return first_part
def extract_HGVS_D(text):
    lst = text.split(':')
    HGVS_D = f'g.{lst[1]}{lst[2]}>{lst[3]}'
    return HGVS_D

def extract_Source(text):
    return text+'(Clinvar)'
def Uniprot_Clinvar(ClinvarPath,UniprotPath,OutputPath):
    # 使用pandas的read_excel函数读取数据
    path = ClinvarPath
    df = pd.read_excel(path)

    HGVS_p= df['Name'].apply(getHGVS_p)
    Position = df['Name'].apply(extract_position)
    TranscriptId =df['Name'].apply(extract_transcriptId)
    NCReseq =df['Canonical SPDI'].apply(extract_NCReseq)
    HGVS_D = df['Canonical SPDI'].apply(extract_HGVS_D)
    CytogeneticBand = '20q13.33'
    ConsequenceType = df['Molecular consequence']
    VariantClinicalSignificance = df['Germline classification']
    Phenotype = df['Condition(s)']
    References = df['Accession'].apply(extract_Source)

    new_df = pd.DataFrame({
        'Position': Position,
        'HGVS_P': HGVS_p,
        'TranscriptId': TranscriptId,
        'HGVS_D': HGVS_D,
        'NCReseq': NCReseq,
        'CytogeneticBand': CytogeneticBand,
        'ConsequenceType': ConsequenceType,
        'VariantClinicalSignificance': VariantClinicalSignificance,
        'Phenotype': Phenotype,
        'References': References
    })


    path_uniprot = UniprotPath
    df_uniprot = pd.read_excel(path_uniprot)



    # 检查new_df中的每个HGVS_P是否存在于df_uniprot的HGVS_P列中
    for index, row in new_df.iterrows():
        if row['HGVS_P'] not in df_uniprot['HGVS_P'].values:
            # 如果不存在，则将该行添加到rows_to_add中
            df_uniprot = pd.concat([df_uniprot, pd.DataFrame([row])], ignore_index=True)

    # 创建一个新列'SortablePosition'用于排序，处理范围中的第一个值
    df_uniprot['SortablePosition'] = df_uniprot['Position'].apply(lambda x: int(str(x).split('-')[0]))
    # 根据'SortablePosition'列对df_uniprot进行升序排序
    df_uniprot.sort_values('SortablePosition', inplace=True)
    # 删除'SortablePosition'列
    df_uniprot.drop(columns=['SortablePosition'], inplace=True)
    # 保存更新后的df_uniprot到Excel文件（可选）
    path_out = OutputPath
    df_uniprot.to_excel(path_out, index=False)















def Uniprot_Clinvar_Cleaned(GeneName,Uniprot_ClinvarPath,Uniprot_Clinvar_CleanedPath):
    path = Uniprot_ClinvarPath
    df = pd.read_excel(path)


    #异常处理Y374D
    if GeneName == 'KCNQ2':
        df['Position'] = df['Position'].astype(str)
        df.loc[df['Position'] == '172107', 'HGVS_P'] = 'p.Tyr374Asp'
        df.loc[df['Position'] == '172107', 'Position'] = '374'
        df.drop(df[df['Position'] == '374'].index, inplace=True)



    #missense variant>missense
    conditions  = (df['ConsequenceType'] == 'missense variant')
    df.loc[conditions, 'ConsequenceType'] = 'missense'
    #print(df['ConsequenceType'])

    #仅保留missense的数据
    #print(len(df))
    df = df[df['ConsequenceType'].str.contains('missense', na=False)]
    #print(len(df))

    #去除p
    def remove_p(str):
        return str.split('.')[-1]
    df['HGVS_P'] = df['HGVS_P'].apply(remove_p)
    #print(df['HGVS_P'])

    #取出插入缺失fs等
    conditions = df['HGVS_P'].str.contains('delins|del|ins|fs', na=False)
    df = df[~conditions]
    ##print(df['HGVS_P'])
    ##print(df.iloc[359:372]['HGVS_P'])
    #print(df.shape)

    #取出Ter(停止)
    conditions = df['HGVS_P'].str.contains('Ter', na=False)
    df = df[~conditions]
    #print(df.shape)

    #转换为1code
    def change_3_to_1(column_str):
        parts = re.split('(\d+)', column_str)
        one_code_parts = [convert_sequence(i, conversion_type='3to1') if not i.isdigit() else i for i in parts]
        return ''.join(one_code_parts)


    df['HGVS_P'] = df['HGVS_P'].apply(change_3_to_1)


    #delete Met1
    conditions = df['HGVS_P'].str.contains(r'M1(?!\d)', na=False)
    df = df[~conditions]
    #print(df.shape)



    df['Position'] = pd.to_numeric(df['Position'], errors='coerce')
    df = df.sort_values(by='Position', ascending=True)

    path_out = Uniprot_Clinvar_CleanedPath
    df.to_excel(path_out,index=False)





def Uniprot_Clinvar_MarshP(Pathclinvar_uniprot,GeneName,Uniprot_Clinvar_CleanedPath,Uniprot_Clinvar_MarshPPath):

    path1 = Pathclinvar_uniprot
    df1 = pd.read_excel(path1)

    Position = df1['uniprot_start']
    HGVS_p= df1['wild_type'] + df1['uniprot_start'].astype(str) + df1['mutant']

    TranscriptId = ''
    NCReseq = ''
    HGVS_D = df1['id']
    CytogeneticBand = '20q13.33'
    ConsequenceType = df1['consequence']
    VariantClinicalSignificance = df1['clinsig']
    Phenotype = df1['phenotype']
    References = df1['allele_id'].astype(str)


    new_df = pd.DataFrame({
        'Position': Position,
        'HGVS_P': HGVS_p,
        'TranscriptId': TranscriptId,
        'HGVS_D': HGVS_D,
        'NCReseq': NCReseq,
        'CytogeneticBand': CytogeneticBand,
        'ConsequenceType': ConsequenceType,
        'VariantClinicalSignificance': VariantClinicalSignificance,
        'Phenotype': Phenotype,
        'References': 'Clinvar(allele_id)'+References
    })

    new_df=new_df[new_df['ConsequenceType'].str.contains('missense', na=False)]


    path_uniprot_clinvar = Uniprot_Clinvar_CleanedPath
    df_uniprot_clinvar = pd.read_excel(path_uniprot_clinvar)



    # 检查new_df中的每个HGVS_P是否存在于df_uniprot的HGVS_P列中
    for index, row in new_df.iterrows():
        if row['HGVS_P'] not in df_uniprot_clinvar['HGVS_P'].values:
            # 如果不存在，则将该行添加到rows_to_add中
            df_uniprot_clinvar = pd.concat([df_uniprot_clinvar, pd.DataFrame([row])], ignore_index=True)

    # 创建一个新列'SortablePosition'用于排序，处理范围中的第一个值
    df_uniprot_clinvar['SortablePosition'] = df_uniprot_clinvar['Position'].apply(lambda x: int(str(x).split('-')[0]))

    # 根据'SortablePosition'列对df_uniprot进行升序排序
    df_uniprot_clinvar.sort_values('SortablePosition', inplace=True)

    # 删除'SortablePosition'列
    df_uniprot_clinvar.drop(columns=['SortablePosition'], inplace=True)

    # 保存更新后的df_uniprot到Excel文件（可选）
    df_uniprot_clinvar.to_excel(Uniprot_Clinvar_MarshPPath, index=False)

def RSA_split(GeneName,UniprotID,Path_PDB_RSA,Path_AF_RSA,OutPath_PDB_RSA,OutPath_AF_RSA):
    df1 = pd.read_csv(Path_AF_RSA)
    df2 = pd.read_csv(Path_PDB_RSA)
    df_filter1 = df1['Uniprot'] == f'{UniprotID}'
    df_af = df1[df_filter1]
    df_af.to_excel(OutPath_AF_RSA, index=False)
    df_split2 = df2['ID'].str.split('|', expand=True)
    df_pdb = df2[df_split2[1] == f'{UniprotID}']
    df_pdb.to_excel(OutPath_PDB_RSA, index=False)





def GnomAD_split(GeneName,Path_MarshLaab_gnomAD,OutPath_GnomAD):
    df1 = pd.read_csv(Path_MarshLaab_gnomAD, sep='\t', skiprows=0)
    df_Prot = df1[df1['gene'] == f'{GeneName}']
    df_Prot.to_excel(OutPath_GnomAD, index=False)



def Uniprot_Clinvar_MarshP_MarshG(GenomADPath,Uniprot_Clinvar_MarshPPath,Uniprot_Clinvar_MarshP_MarshGPath):
    path1 = GenomADPath
    df1 = pd.read_excel(path1)

    Position = df1['variant'].str.extract('(\d+)').iloc[:, 0]
    HGVS_p= df1['variant']
    TranscriptId = ''
    NCReseq = ''
    HGVS_D = df1['id']
    CytogeneticBand = '20q13.33'
    ConsequenceType = df1['consequence']
    VariantClinicalSignificance = 'Benign(gnomAD)'
    Phenotype = 'Benign(gnomAD)'
    References = df1['id'].astype(str)


    new_df = pd.DataFrame({
        'Position': Position,
        'HGVS_P': HGVS_p,
        'TranscriptId': TranscriptId,
        'HGVS_D': HGVS_D,
        'NCReseq': NCReseq,
        'CytogeneticBand': CytogeneticBand,
        'ConsequenceType': ConsequenceType,
        'VariantClinicalSignificance': VariantClinicalSignificance,
        'Phenotype': Phenotype,
        'References': 'gnomAD(id)'+References
    })

    new_df=new_df[new_df['ConsequenceType'].str.contains('missense', na=False)]
    new_df['Position'] = pd.to_numeric(new_df['Position'], errors='coerce')

    path_uniprot_clinvar_marshP = Uniprot_Clinvar_MarshPPath
    df_uniprot_clinvar_marshP = pd.read_excel(path_uniprot_clinvar_marshP)


    # 检查new_df中的每个HGVS_P是否存在于df_uniprot的HGVS_P列中
    for index, row in new_df.iterrows():
        if row['HGVS_P'] not in df_uniprot_clinvar_marshP['HGVS_P'].values:
            # 如果不存在，则将该行添加到rows_to_add中
            df_uniprot_clinvar_marshP = pd.concat([df_uniprot_clinvar_marshP, pd.DataFrame([row])], ignore_index=True)

    # 创建一个新列'SortablePosition'用于排序，处理范围中的第一个值
    df_uniprot_clinvar_marshP['SortablePosition'] = df_uniprot_clinvar_marshP['Position'].apply(lambda x: int(str(x).split('-')[0]))
    # 根据'SortablePosition'列对df_uniprot进行升序排序
    df_uniprot_clinvar_marshP.sort_values('SortablePosition', inplace=True)
    # 删除'SortablePosition'列
    df_uniprot_clinvar_marshP.drop(columns=['SortablePosition'], inplace=True)
    # 保存更新后的df_uniprot到Excel文件（可选）


    path_out =Uniprot_Clinvar_MarshP_MarshGPath
    df_uniprot_clinvar_marshP.to_excel(path_out, index=False)




def Uniprot_Clinvar_MarshP_MarshG_Cleaned(Uniprot_Clinvar_MarshP_MarshGPath,Uniprot_Clinvar_MarshP_MarshG_CleanedPath):
    path = Uniprot_Clinvar_MarshP_MarshGPath
    df = pd.read_excel(path)

    df.drop(columns=['TranscriptId','HGVS_D','NCReseq','CytogeneticBand'],inplace=True)

    #处理所有的missense列
    conditions  = (df['ConsequenceType'] == 'missense variant')
    df.loc[conditions, 'ConsequenceType'] = 'missense'
    conditions  = (df['ConsequenceType'] == 'missense_variant')
    df.loc[conditions, 'ConsequenceType'] = 'missense'

    #处理VariantClinicalSignificance列
    df.loc[df['VariantClinicalSignificance'] == 'Variant of uncertain significance', 'VariantClinicalSignificance'] = 'Uncertain_significance'
    df.loc[df['VariantClinicalSignificance'] == 'Likely benign', 'VariantClinicalSignificance'] = 'Likely_benign'
    df.loc[df['VariantClinicalSignificance'] == 'not_provided', 'VariantClinicalSignificance'] = ''
    df.loc[df['VariantClinicalSignificance'] == 'Likely pathogenic', 'VariantClinicalSignificance'] = 'Likely_pathogenic'
    df.loc[df['VariantClinicalSignificance'] == 'Pathogenic/Likely pathogenic', 'VariantClinicalSignificance'] = 'pathogenic/Likely_pathogenic'
    df.loc[df['VariantClinicalSignificance'] == 'Pathogenic', 'VariantClinicalSignificance'] = 'pathogenic'




    out_path = Uniprot_Clinvar_MarshP_MarshG_CleanedPath
    df.to_excel(out_path,index=False)










def merge_data(Uniprot_Clinvar_MarshP_MarshG_CleanedPath,
               af_foldxPath,
               pdb_foldxPath,
               af_RSAPath,
               pdb_RSAPath,
               VEPSPath,
               ConsurfPath,
               DESCRIBEPROPath,
               FULLPath):
    #originaldata
    path1 = Uniprot_Clinvar_MarshP_MarshG_CleanedPath
    df1 = pd.read_excel(path1)


    #foldx_af
    path2 = af_foldxPath
    df2 = pd.read_csv(path2)
    #foldx_pdb
    path3 = pdb_foldxPath
    df3 = pd.read_csv(path3)



    #merge_foldx_af
    df_merge_af = pd.merge(df1, df2[['variant', 'ddG_fold', 'fold rank', 'abs(fold) rank']], left_on='HGVS_P', right_on='variant', how='left')
    df_merge_af.rename(columns={'ddG_fold': 'ddG_fold_af', 'fold rank': 'fold_rank_af', 'abs(fold) rank': 'abs(fold)_rank_af'}, inplace=True)

    #merge_foldx_af with pdb
    df_merge_af_pdb = pd.merge(df_merge_af, df3[['variant', 'ddG_fold','ddG_full', 'fold rank','full rank', 'abs(fold) rank','abs(full) rank']], left_on='HGVS_P', right_on='variant', how='left')
    df_merge_af_pdb.rename(columns={'ddG_fold': 'ddG_fold_pdb','ddG_full':'ddG_full_pdb', 'fold rank': 'fold_rank_pdb','full rank':'full_rank_pdb', 'abs(fold) rank': 'abs(fold)_rank_pdb','abs(full) rank':'abs(full)_rank_pdb'}, inplace=True)
    df_merge_af_pdb.drop(columns=['variant_x','variant_y'], inplace=True)

    #print(df_merge_af)
    #print(df_merge_af_pdb)

    #df_merge_af_pdb.to_excel(r"D:\file\毕业设计\Data\KCNQ2_withFoldX.xlsx", index=False)



    df_merge1 = df_merge_af_pdb
    #RSA_af
    path4 = af_RSAPath
    df4 = pd.read_excel(path4)
    #RSA_pdb
    path5 = pdb_RSAPath
    df5 = pd.read_excel(path5)

    #merge_RSA_af
    df_merge1_af = pd.merge(df_merge1, df4, left_on='Position', right_on='residue', how='left')
    #df_merge_af.rename(columns={'ddG_fold': 'ddG_fold_af', 'fold rank': 'fold_rank_af', 'abs(fold) rank': 'abs(fold)_rank_af'}, inplace=True)
    df_merge1_af.drop(columns=['Gene','Uniprot','residue','amino acid','number of structure'], inplace=True)
    df_merge1_af.rename(columns={'secondary structure': 'secondary_structure_af', 'SASA': 'SASA_af', 'RSA': 'RSA_af','score':'score_af'}, inplace=True)

    #merge_RSA_pdb
    df_merge1_af_pdb=pd.merge(df_merge1_af, df5, left_on='Position', right_on='residue', how='left')
    df_merge1_af_pdb.drop(columns=['ID','residue'], inplace=True)
    df_merge1_af_pdb.rename(columns={'PDB residue': 'PDB_structure', 'location': 'location_pdb', 'secondary structure': 'secondary_structure_pdb','SASA full':'SASA_full_pdb','SASA chain':'SASA_chain_pdb','RSA full':'RSA_full_pdb','RSA chain':'RSA_chain_pdb'}, inplace=True)
    #print(df_merge1_af_pdb)

    #split the structure
    df_split = df_merge1_af_pdb['PDB_structure'].str.split(' ', expand=True)
    df_merge1_af_pdb['PDB_structure'] = df_split[0]
    #df_merge1_af_pdb['PDB_structure'] = df_merge1_af_pdb['PDB_structure'].fillna('alphafold')

    #print(df_merge1_af_pdb['PDB_structure'])
    #df_merge1_af_pdb.to_excel(r"D:\file\毕业设计\Data\KCNQ2_withFoldX_withRSA.xlsx", index=False)



    #start to merge veps
    df_merge2 = df_merge1_af_pdb

    #Veps
    path6 = VEPSPath
    df6 = pd.read_csv(path6)
    #print(df6.shape)


    #Conservartoin
    #consurf
    path7 = ConsurfPath
    df7 = pd.read_excel(path7)

    #DESCRIBEPRO
    path8 = DESCRIBEPROPath
    df8 = pd.read_csv(path8,skiprows=47)


    df_merge2_DESCRIBEPRO = pd.merge(df_merge2, df8[['Index','MMseqs2_conservation_propensity','MMseqs2_conservation_decile']], left_on='Position', right_on='Index', how='left')
    df_merge2_DESCRIBEPRO.drop(columns=['Index'], inplace=True)

    df_merge2_DESCRIBEPRO_consurf = pd.merge(df_merge2_DESCRIBEPRO, df7[['POS',' SCORE','COLOR','CONFIDENCE INTERVAL','RESIDUE VARIETY']], left_on='Position', right_on='POS', how='left')
    df_merge2_DESCRIBEPRO_consurf.drop(columns=['POS'], inplace=True)
    df_merge2_DESCRIBEPRO_consurf.rename(columns={'CONFIDENCE INTERVAL': 'CONFIDENCE_INTERVAL', 'RESIDUE VARIETY': 'RESIDUE_VARIETY'}, inplace=True)

    #print(df_merge2_DESCRIBEPRO_consurf.columns)

    #VEPS
    df_merge3 = df_merge2_DESCRIBEPRO_consurf
    df_merge3_veps = pd.merge(df_merge3, df6, left_on='HGVS_P', right_on='variant', how='left')
    df_merge3_veps.drop(columns=['variant'], inplace=True)
    #print(df_merge3_veps.columns)


    df_merge3_veps.to_excel(FULLPath, index=False)


















import shutil
def split_full(FullPath,GeneName,OutPath_P,OutPath_B):
    #print(FullPath)
    # 原始文件路径
    original_path = FullPath
    #print(original_path)
    # 目标文件路径，原路径基础上添加后缀 '_copy'
    target_path = original_path[:-5]+'_copy.xlsx'
    # 使用shutil.copyfile复制文件
    shutil.copyfile(original_path, target_path)

    #根据致病性拆分
    df = pd.read_excel(target_path)

    df_benign = df[df['VariantClinicalSignificance'].str.contains('Benign',case=False,na=False)]
    #print(df_benign)
    df_pathogenic = df[df['VariantClinicalSignificance'].str.contains('pathogenic',case=False,na=False)]
    #print(df_pathogenic)


    out_path_benign =rf"D:\file\毕业设计\DataHandle\OutputData\{GeneName}_benign.xlsx"
    out_path_pathogenic =rf"D:\file\毕业设计\DataHandle\OutputData\{GeneName}_pathogenic.xlsx"
    df_benign.to_excel(OutPath_B,index=False)
    df_pathogenic.to_excel(OutPath_P,index=False)




    def split_sheets(filepath):
        file = pd.ExcelFile(filepath)
        # 从sheet1读取数据
        df = file.parse('Sheet1')
        # Bacis Info
        # 选择A和B列
        sheet2 = df[['HGVS_P', 'ConsequenceType','VariantClinicalSignificance','Phenotype','References']]
        # 将选定的列复制到新的DataFrame中（如果sheet2不存在，将创建一个新的sheet2）
        with pd.ExcelWriter(file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            if 'Basic_info' not in writer.book.sheetnames:
                # 创建一个新的sheet2，如果它不存在
                writer.book.create_sheet('Basic_info')
            sheet2.to_excel(writer, sheet_name='Basic_info', index=False)
        # ddG
        sheet3 = df.iloc[:, [1] + list(range(6, 15))]  # Python中的索引是从0开始，所以K列索引是10，S列索引是18
        with pd.ExcelWriter(file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            if 'ddG' not in writer.book.sheetnames:
                # 创建一个新的sheet2，如果它不存在
                writer.book.create_sheet('ddG')
            sheet3.to_excel(writer, sheet_name='ddG', index=False)
        # RSA
        sheet4 = df.iloc[:, [1] + list(range(15, 26))]
        with pd.ExcelWriter(file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            if 'RSA' not in writer.book.sheetnames:
                # 创建一个新的sheet2，如果它不存在
                writer.book.create_sheet('RSA')
            sheet4.to_excel(writer, sheet_name='RSA', index=False)
        #conservation
        sheet5  = df.iloc[:, [1] + list(range(26, 32))]
        with pd.ExcelWriter(file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            if 'Conservation' not in writer.book.sheetnames:
                # 创建一个新的sheet2，如果它不存在
                writer.book.create_sheet('Conservation')
            sheet5.to_excel(writer, sheet_name='Conservation', index=False)
        # VEPS
        sheet6  = df.iloc[:, [1] + list(range(32, len(df.columns)))]
        with pd.ExcelWriter(file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            if 'VEPS' not in writer.book.sheetnames:
                # 创建一个新的sheet2，如果它不存在
                writer.book.create_sheet('VEPS')
            sheet6.to_excel(writer, sheet_name='VEPS', index=False)

    # 加载Excel文件

    path_full = target_path
    path_benign = out_path_benign
    path_pathogenic = out_path_pathogenic

    split_sheets(path_full)
    split_sheets(path_benign )
    split_sheets(path_pathogenic)