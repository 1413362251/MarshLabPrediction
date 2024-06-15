from MergeMethods import *
import sys

#开始之前,需要准备:
# 1.Uniprot下载的变体数据,{UniprotID}.xml在RawData文件夹里
# 2.Gene的Clinvar数据 \RawData\{GeneName}_clinvar_result.xlsx

def data_collectting_process(GeneName,UniprotID):
    parent_dir = os.path.dirname(os.getcwd())
    Rawdata_dir = os.path.join(parent_dir, 'DataHandle', 'RawData')
    IntermediateData_dir = os.path.join(parent_dir, 'DataHandle', 'IntermediateData')
    OutPutData_dir = os.path.join(parent_dir, 'DataHandle', 'OutPutData')

    #从Uniprot数据处理开始
    print("Handling webdata(Uniprot)...")
    XML_filePath = rf"{Rawdata_dir}\{UniprotID}.xml"
    UniprotPath = fr"{IntermediateData_dir}\{GeneName}_Uniprot.xlsx"
    main_Uniprot(XML_filePath,UniprotPath )

    #Uniprot_clinvar
    print("Handling webdata(Clinvar)...")
    UniprotPath = UniprotPath
    Path_clinvarfile = fr"{Rawdata_dir}\MarshLab_clinvar_uniprot.tsv"
    OutPath_clinvar_uniprot =fr"{IntermediateData_dir}\{GeneName}_clinvar_uniprot.xlsx"
    clinvar_big_orignal_sep(Path_clinvarfile, GeneName, OutPath_clinvar_uniprot)
    ClinvarPath = rf"{Rawdata_dir}\{GeneName}_clinvar_result.xlsx"
    Uniprot_ClinvarPath =fr"{IntermediateData_dir}\{GeneName}_Uniprot_Clinvar.xlsx"
    Uniprot_Clinvar(ClinvarPath,UniprotPath,Uniprot_ClinvarPath)

    #Uniprot_Clinvar_cleaned
    print("Normalising webdata...")
    Uniprot_ClinvarPath = Uniprot_ClinvarPath
    Uniprot_Clinvar_CleanedPath = fr"{IntermediateData_dir}\{GeneName}_Uniprot_Clinvar_cleaned.xlsx"
    GeneName = GeneName
    Uniprot_Clinvar_Cleaned(GeneName,Uniprot_ClinvarPath,Uniprot_Clinvar_CleanedPath,)


    #Start to handle MarshLab Data
    print("Start to process MarshLabData...")
    #Seperate GnomAD
    print("Spliting MarshLabData(GnomAD)...")
    Path_MarshLaab_gnomAD = rf"{Rawdata_dir}\MarshLab_gnomad.tsv"
    OutPath_GnomAD = rf"{IntermediateData_dir}\{GeneName}_GenomAD.xlsx"
    GnomAD_split(GeneName, Path_MarshLaab_gnomAD, OutPath_GnomAD)

    #Seperate RSA
    print("Spliting MarshLabData(RSA)...")
    Path_PDB_RSA = rf"{Rawdata_dir}\MarshLab_PDB_RSA.csv"
    Path_AF_RSA = rf"{Rawdata_dir}\MarshLab_AF_RSA.csv"
    OutPath_PDB_RSA = rf"{Rawdata_dir}\{GeneName}_PDB_RSA.xlsx"
    OutPath_AF_RSA =rf"{Rawdata_dir}\{GeneName}_AF_RSA.xlsx"
    RSA_split(GeneName, UniprotID, Path_PDB_RSA, Path_AF_RSA, OutPath_PDB_RSA, OutPath_AF_RSA)


    #Uniprot_Clinvar_MarshP
    print("Updating Webdata with MarshLabData(Pathogenic)...")
    GeneName = GeneName
    Uniprot_Clinvar_CleanedPath = Uniprot_Clinvar_CleanedPath
    Uniprot_Clinvar_MarshPPath = fr"{IntermediateData_dir}\{GeneName}_Uniprot_Clinvar_MarshP.xlsx"
    Pathclinvar_uniprot = rf"{IntermediateData_dir}\{GeneName}_clinvar_uniprot.xlsx"
    Uniprot_Clinvar_MarshP(Pathclinvar_uniprot, GeneName, Uniprot_Clinvar_CleanedPath, Uniprot_Clinvar_MarshPPath)

    #Uniprot_Clinvar_MarshP_MarshG
    print("Updating Webdata with MarshLabData(Benign)...")
    GenomADPath = OutPath_GnomAD
    Uniprot_Clinvar_MarshPPath = Uniprot_Clinvar_MarshPPath
    Uniprot_Clinvar_MarshP_MarshGPath = rf"{IntermediateData_dir}\{GeneName}_Uniprot_Clinvar_MarshP_MarshG.xlsx"
    Uniprot_Clinvar_MarshP_MarshG(GenomADPath,Uniprot_Clinvar_MarshPPath,Uniprot_Clinvar_MarshP_MarshGPath)


    #Uniprot_Clinvar_MarshP_MarshG_Cleaned
    print("Normaling ...")
    Uniprot_Clinvar_MarshP_MarshGPath = Uniprot_Clinvar_MarshP_MarshGPath
    Uniprot_Clinvar_MarshP_MarshG_CleanedPath = rf"{IntermediateData_dir}\{GeneName}_Uniprot_Clinvar_MarshP_MarshG_Cleaned.xlsx"
    Uniprot_Clinvar_MarshP_MarshG_Cleaned(Uniprot_Clinvar_MarshP_MarshGPath, Uniprot_Clinvar_MarshP_MarshG_CleanedPath)


    #Merge with other data RSA,Conversation,VEPS....
    print("Merging [RSA,Conversation,VEPs] data...")
    Uniprot_Clinvar_MarshP_MarshG_CleanedPath = Uniprot_Clinvar_MarshP_MarshG_CleanedPath
    af_foldxPath = rf"{Rawdata_dir}\{UniprotID}_alphafold_foldx.csv"
    pdb_foldxPath = rf"{Rawdata_dir}\{UniprotID}_pdbstruct_foldx.csv"
    af_RSAPath =  rf"{Rawdata_dir}\{GeneName}_AF_RSA.xlsx"
    pdb_RSAPath = rf"{Rawdata_dir}\{GeneName}_PDB_RSA.xlsx"
    VEPSPath = rf"{Rawdata_dir}\{UniprotID}_VEPs.csv"
    ConsurfPath = fr"{Rawdata_dir}\{UniprotID}_A_consurf_grades.txt"
    DESCRIBEPROPath = rf"{Rawdata_dir}\{UniprotID}_testresult.csv"
    FULLPath = rf"{OutPutData_dir}\{GeneName}_full.xlsx"
    merge_data(Uniprot_Clinvar_MarshP_MarshG_CleanedPath,
                   af_foldxPath,
                   pdb_foldxPath,
                   af_RSAPath,
                   pdb_RSAPath,
                   VEPSPath,
                   ConsurfPath,
                   DESCRIBEPROPath,
                   FULLPath)


    #Split_Sheet
    print("Spliting the final sheet...")
    FULLPath = FULLPath
    GeneName = GeneName
    OutPath_P = fr"{OutPutData_dir}\{GeneName}_benign.xlsx"
    OutPath_B = fr"{OutPutData_dir}\{GeneName}_pathogenic.xlsx"
    split_full(FULLPath, GeneName, OutPath_P, OutPath_B)
    print("Merging success.")
#P56696 · KCNQ4_HUMAN P51787 · KCNQ1_HUMAN
GeneName = sys.argv[1]
UniprotID = sys.argv[2]
data_collectting_process(GeneName,UniprotID)