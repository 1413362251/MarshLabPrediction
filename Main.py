from methods import *
import subprocess


def menu():
    subprocess.run(["python", "Menu.py",GeneName,UniprotID])
def spider():
    subprocess.run(["python", "spider.py",GeneName,UniprotID])
def merge():
    subprocess.run(["python", "MergeProcess.py",GeneName,UniprotID])
def check_missing():
    parent_dir = os.path.dirname(os.getcwd())
    download_dir = os.path.join(parent_dir, 'DataHandle', 'RawData')
    Web_data_lst = [f"{download_dir}/{GeneName}_clinvar_result.xlsx",
                    f"{download_dir}/{UniprotID}.xml",
                    f"{download_dir}/{UniprotID}_consurf.xlsx",
                    f"{download_dir}/{UniprotID}_testresult.csv",
                    ]
    MarshLab_data_lst = [f"{download_dir}/MarshLab_clinvar_uniprot.tsv",
                         f"{download_dir}/MarshLab_gnomad.tsv",
                         f"{download_dir}/MarshLab_PDB_RSA.csv",
                         f"{download_dir}/MarshLab_AF_RSA.csv",
                         f"{download_dir}/{UniprotID}_pdbstruct_foldx.csv",
                         f"{download_dir}/{UniprotID}_alphafold_foldx.csv",
                         f"{download_dir}/{UniprotID}_VEPs.csv",
                         ]
    # a space list to store their names
    missing_files = []
    for file in Web_data_lst:
        if not os.path.exists(file):
            missing_files.append(os.path.basename(file))
    for file in MarshLab_data_lst:
        if not os.path.exists(file):
            missing_files.append(os.path.basename(file))
    return missing_files

GeneName = 'KCNQ2'
UniprotID = 'O43526'
check_folders()

if not GeneName:
    GeneName = input('Please set the GeneName (eg:KCNQ2) :')
    UniprotID = input('Please set the UniprotID (eg:O43526) :')




while True:
    missing_files = check_missing()
    menu()
    command = input('Command :')
    match command:
        case '1':
            GeneName = input('Please set the GeneName (eg:KCNQ2) :')
            UniprotID = input('Please set the UniprotID (eg:O43526) :')
            next
        case '2':
            print()
            print()
            print("Spider Start".center(100,'='))
            spider()
            print()

        case '3':
            if missing_files:
                print("Lack these files :")
                print(missing_files)
                print("Pls put them into the folder RawData")
                print()
            else:
                print()
                print("All files are prepared, can be merged")
                print()

        case '4':
            if missing_files:
                print()
                print("Some files are missing ,choose 3 to print the list and prepare")
                print()
            else:
                print()
                print("Merge Start".center(100,'='))
                merge()
                print()
        case _:
            print("Wrong command, input again")
            next