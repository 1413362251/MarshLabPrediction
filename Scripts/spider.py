from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from webdriver_manager.microsoft import EdgeChromiumDriverManager
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.edge.options import Options
import time
import pandas as pd
import os
import requests
import sys



def get_clinvar_data(GeneName):
    print('Spider is acquiring Clinvar Data...')
    #get parent_dir of current pwd
    parent_dir = os.path.dirname(os.getcwd())
    # Set DataHandle dir
    download_dir = os.path.join(parent_dir, 'DataHandle','RawData')
    Path_file = os.path.join(parent_dir, 'DataHandle', 'RawData','clinvar_result.txt')
    # 创建一个Options实例
    options = Options()
    options.add_experimental_option("prefs", {
        "download.default_directory": download_dir,  # 设置默认下载路径
        "download.prompt_for_download": False,  # 不提示下载窗口
        "download.directory_upgrade": True,  # 使用系统下载目录
        "safebrowsing.enabled": True  # 启用安全浏览
    })
    options.add_argument("--headless")  # 启用无头模式

    # 初始化WebDriver
    service = Service(EdgeChromiumDriverManager().install())
    driver = webdriver.Edge(service=service, options=options)


    clinvar_url = fr'https://www.ncbi.nlm.nih.gov/clinvar/?term={GeneName}%5Bgene%5D'
    # 打开网站
    driver.get(clinvar_url)  # 替换为实际的网址

    # 定位按钮并点击
    button = driver.find_element(By.XPATH, '//*[@id="_Properties"]/li[1]/ul[1]/li[5]/input[1]')
    button.click()

    # 定位按钮并点击
    button = driver.find_element(By.XPATH, '//*[@id="_Properties"]/li[1]/ul[1]/li[6]/input[1]')
    button.click()
    # 定位按钮并点击
    button = driver.find_element(By.XPATH, '//*[@id="_molConseq"]/li[1]/ul[1]/li[2]/input[1]')
    button.click()

    element = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.XPATH, '//*[@id="submenu_File"]/button[1]'))
    )
    # 使用 JavaScript 执行点击
    driver.execute_script("arguments[0].click();", element)
    time.sleep(5)

    # 完成后关闭浏览器
    driver.quit()


    # 读取TXT文件
    df = pd.read_csv(Path_file, delimiter='\t')  # 假设数据是用制表符分隔的

    # 将DataFrame导出到Excel文件
    df.to_excel(f'{download_dir}/{GeneName}_clinvar_result.xlsx', index=False, engine='openpyxl')
    # 删除原始的TXT文件
    os.remove(Path_file)

def get_uniprot_data(GeneName,UniportID):
    print('Spider is acquiring Uniprot Data...')
    #get parent_dir of current pwd
    parent_dir = os.path.dirname(os.getcwd())
    # Set DataHandle dir
    download_dir = os.path.join(parent_dir, 'DataHandle','RawData')
    url = f'https://www.ebi.ac.uk/proteins/api/variation/{UniportID}?format=xml'
    response = requests.get(url)
    if response.status_code == 200:
        # 指定要保存的文件路径
        file_path = fr"{download_dir}/{UniportID}.xml"
        # 打开文件进行写入
        with open(file_path, 'wb') as file:
            file.write(response.content)

#P56696 · KCNQ4_HUMAN P51787 · KCNQ1_HUMAN
KCNQ_family = [['KCNQ1','P51787'],['KCNQ2','O43526'],['KCNQ3','O43525'],['KCNQ4','P56696'],['KCNQ5','Q9NR82']]
GeneName = sys.argv[1]
UniportID = sys.argv[2]

get_clinvar_data(GeneName)
get_uniprot_data(GeneName,UniportID)
print("Spider Success")