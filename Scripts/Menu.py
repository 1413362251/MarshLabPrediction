from colorama import init, Fore, Style
import sys
import os


init(autoreset=True)

GeneName = sys.argv[1]
UniprotID = sys.argv[2]

parent_dir = os.path.dirname(os.getcwd())
download_dir = os.path.join(parent_dir, 'DataHandle', 'RawData')
Flag_WebData = False
Flag_MarshLabData = False
Web_data_lst = [f"{download_dir}/{GeneName}_clinvar_result.xlsx",
                f"{download_dir}/{UniprotID}.xml",
                f"{download_dir}/{UniprotID}_A_consurf_grades.txt",
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

Flag_WebData = all(os.path.exists(file) for file in Web_data_lst)
Flag_MarshLabData = all(os.path.exists(file) for file in MarshLab_data_lst)

print('Meun'.center(100,'='))
print()
print(f'GeneName : {GeneName}\tUniprotID : {UniprotID}'.center(100)  )
print(f'Web Data : {Flag_WebData}'.center(100)  )
print(f'MarshLab Data : {Flag_MarshLabData}'.center(100)  )
print()
print('1.Change the GeneName and UniprotID'.ljust(50).center(100))
print('2.Acquire web data'.ljust(50).center(100))
print('3.Print the list of missing files'.ljust(50).center(100))
print('4.Merge data'.ljust(50).center(100))
print()
print(''.center(100,'='))
