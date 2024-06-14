import os

def check_folders():
    #get parent_dir of current pwd
    parent_dir = os.path.dirname(os.getcwd())
    # Set DataHandle dir
    data_handle_path = os.path.join(parent_dir, 'DataHandle')
    # check wether DataHandle exist
    if not os.path.exists(data_handle_path):
        os.makedirs(data_handle_path)
        for sub_dir in ['RawData', 'OutPutData', 'IntermediateData']:
            os.makedirs(os.path.join(data_handle_path, sub_dir))


def convert_sequence(sequence, conversion_type='1to3'):
    """转换蛋白质序列的缩写类型。

    参数:
        sequence (str): 输入的蛋白质序列。
        conversion_type (str): 转换类型，'1to3' 或 '3to1'。

    返回:
        str: 转换后的蛋白质序列。
    """
    # 定义1字母缩写到3字母缩写的映射
    one_to_three = {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'E': 'Glu',
        'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
        'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp',
        'Y': 'Tyr', 'V': 'Val'
    }
    # 定义3字母缩写到1字母缩写的映射
    three_to_one = {v: k for k, v in one_to_three.items()}
    if conversion_type == '1to3':
        return ''.join(one_to_three.get(res, 'X') for res in sequence)
    elif conversion_type == '3to1':
        return ''.join(three_to_one.get(sequence[i:i + 3], 'X') for i in range(0, len(sequence), 3))