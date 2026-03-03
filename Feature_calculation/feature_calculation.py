# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 10:07:55 2026
Program to calculate features of Protein sequences using 'protr' package
@author: anju
"""

import pandas as pd
from Bio import SeqIO

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# ==========================
# USER INPUT
# ==========================
fasta_file = "sample_input.fasta"
output_file = "Features.csv"

# ==========================
# STEP 1: Read & Clean FASTA
# ==========================
print("Reading FASTA file...")

valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

sequence_ids = []
clean_sequences = []

for record in SeqIO.parse(fasta_file, "fasta"):
    seq = str(record.seq).strip().upper()
    if set(seq).issubset(valid_amino_acids) and len(seq) >= 3:
        sequence_ids.append(record.id.strip())
        clean_sequences.append(seq)
    else:
        print(f"Removed invalid sequence: {record.id}")

print("Valid sequences:", len(clean_sequences))
if len(clean_sequences) == 0:
    raise ValueError("No valid sequences remaining after cleaning.")
print("Minimum sequence length:", min(len(s) for s in clean_sequences))

# ==========================
# STEP 2: Load protr package
# ==========================
print("Loading protr package in R...")
protr = importr("protr")

# ==========================
# Functions to check quality of data and data conversions
# ==========================
def _safe_list(x):
    if x is None or x == ro.NULL:
        return None
    try:
        return list(x)
    except Exception:
        return None

def r_obj_to_df_and_dimnames(r_obj):
    
    # detect dim()
    try:
        dim = ro.r["dim"](r_obj)
        has_dim = not (dim is None or dim == ro.NULL)
    except Exception:
        has_dim = False

    if not has_dim:
        r_names = _safe_list(ro.r["names"](r_obj))
        values = list(r_obj)
        df = pd.DataFrame([values])
        if r_names is not None and len(r_names) == df.shape[1]:
            df.columns = r_names
        return df, None, None, r_names

    r_row = _safe_list(ro.r["rownames"](r_obj))
    r_col = _safe_list(ro.r["colnames"](r_obj))

    with localconverter(ro.default_converter + pandas2ri.converter):
        py_obj = ro.conversion.rpy2py(r_obj)
    df = pd.DataFrame(py_obj)

    if r_row is not None and len(r_row) == df.shape[0]:
        df.index = r_row
    if r_col is not None and len(r_col) == df.shape[1]:
        df.columns = r_col

    return df, r_row, r_col, None

def _looks_numeric_columns(cols):
    return all(isinstance(c, int) or (isinstance(c, str) and c.isdigit()) for c in cols)

def finalize_feature_columns(df, prefix, feature_names=None):
    out = df.copy()
    if feature_names is not None and len(feature_names) == out.shape[1]:
        out.columns = feature_names
    return out

def extract_named(func, sequences, prefix, force_per_sequence=False):

    n_seq = len(sequences)

    def one_run(r_obj):
        df, r_row, r_col, r_names = r_obj_to_df_and_dimnames(r_obj)

      
        if r_names is not None:
            df = df.reset_index(drop=True)
            return df, r_names

        
        if df.shape[0] != n_seq and df.shape[1] == n_seq:
            # features x seq -> transpose
            df = df.T

           
            if _looks_numeric_columns(df.columns) and r_row is not None and len(r_row) == df.shape[1]:
                df.columns = r_row

        
        if df.shape[0] == n_seq and _looks_numeric_columns(df.columns):
            
            if r_row is not None and len(r_row) == df.shape[1]:
                df.columns = r_row

        df = df.reset_index(drop=True)
        return df, list(df.columns)

    if not force_per_sequence:
        try:
            r_obj = func(ro.StrVector(sequences))
            df, feat_names = one_run(r_obj)
            if df.shape[0] == n_seq:
                return finalize_feature_columns(df, prefix, feat_names)
            else:
                print(f"[{prefix}] Batch gave {df.shape[0]} rows (expected {n_seq}). Falling back to per-sequence...")
        except Exception as e:
            print(f"[{prefix}] Batch failed. Falling back to per-sequence. Reason: {e}")

    rows = []
    ref_feat_names = None

    for seq in sequences:
        r_obj_i = func(ro.StrVector([seq]))
        df_i, feat_names_i = one_run(r_obj_i)

        # ensure single row
        if df_i.shape[0] != 1 and df_i.shape[1] == 1:
            df_i = df_i.T
        if df_i.shape[0] != 1:
            df_i = pd.DataFrame([df_i.to_numpy().ravel()])

        if ref_feat_names is None:
            ref_feat_names = feat_names_i
        else:
            # enforce same feature naming
            if len(ref_feat_names) == df_i.shape[1]:
                df_i.columns = ref_feat_names

        rows.append(df_i)

    df = pd.concat(rows, axis=0).reset_index(drop=True)
    return finalize_feature_columns(df, prefix, list(df.columns))

# ==========================
# STEP 3: SET 1 (Composition)
# ==========================
print("Generating Set1 (AAC + DPC + TPC)...")

aac_df = extract_named(protr.extractAAC, clean_sequences, "AAC")
dpc_df = extract_named(protr.extractDC,  clean_sequences, "DPC")
tpc_df = extract_named(protr.extractTC,  clean_sequences, "TPC")

set1 = pd.concat([aac_df, dpc_df, tpc_df], axis=1)
print("Set1 shape:", set1.shape)

# ==========================
# STEP 4: SET 2 (CTriad)
# ==========================
print("Generating Set2 (CTriad)...")

set2 = extract_named(protr.extractCTriad, clean_sequences, "CTriad")
print("Set2 shape:", set2.shape)

# ==========================
# STEP 5: SET 3 (CTD)
# ==========================
print("Generating Set3 (CTDC + CTDT + CTDD)...")

ctdc_df = extract_named(protr.extractCTDC, clean_sequences, "CTDC")
ctdt_df = extract_named(protr.extractCTDT, clean_sequences, "CTDT", force_per_sequence=True)
ctdd_df = extract_named(protr.extractCTDD, clean_sequences, "CTDD", force_per_sequence=True)

set3 = pd.concat([ctdc_df, ctdt_df, ctdd_df], axis=1)
print("Set3 shape:", set3.shape)

# ==========================
# STEP 6: Combine All Features
# ==========================
print("Combining all features...")

all_features = pd.concat([set1, set2, set3], axis=1)

if len(sequence_ids) != all_features.shape[0]:
    raise ValueError(f"Row mismatch: IDs={len(sequence_ids)} vs features={all_features.shape[0]}")

all_features.insert(0, "Sequence_ID", sequence_ids)

# ==========================
# STEP 7: Keep ONLY selected features
# ==========================

wanted_str = """GC,C,CP,CC,YVC,TC,LCY,TCP,VK,NT,CNT,VS537,prop4.G3.residue0,IN,VS236,DM,EMR,prop7.G2.residue0,RG,LHK,VS735,VS436,prop3.G1.residue0,PY,VS372,prop1.G2.residue100,NL,prop6.G1.residue0,KCN,VS273,prop7.G3.residue0,CPK,prop2.G2.residue50,KR,CN,prop5.G2.residue0,LVK,VS517,prop7.G1.residue25,CS,CI,KAS,prop2.G2.residue0,K,VS725,VKY,VC,prop2.G2.residue25,prop7.G2.residue25,VS427,VS574,KRG,VS244,prop7.G3.residue100,LC,prop7.G1.residue0,normwaalsvolume.Group1,VS433,prop2.G2.residue75,CCN,normwaalsvolume.Group2,prop7.G2.residue100,secondarystruct.Group2,VS743,REI,KD,ASI,NQK,LH,prop5.G3.residue0,EQV,GP,PYE,prop6.G2.residue0,TTK,QVE,NTL,FPY,YK,prop3.Tr1331,VS745,LTT,DLN,GW,prop7.Tr2332,KY,YEP,NM,LSG,PA,LY,prop5.G3.residue50,LKC,IG,CDS,KV,KM,TLS,prop5.G3.residue75,prop2.G2.residue100,NSS,prop5.G2.residue25,VS364,YGG,VS611,NCK,G,prop4.G3.residue75,VS442,CCK,prop6.Tr2332,KC,PV,prop1.G2.residue75,PTM,CNP,VS635,prop2.Tr2332,solventaccess.Group1,VS762,LML,RQ,VS361,GCP,HN,SSK,prop5.G1.residue0,PSG,W,VS142,SG,GT,prop7.G1.residue50,prop3.Tr2332,KN,AV,prop3.G2.residue0,EA,CCE,VS362,GPS,NLY,A,SMN,VS276,QET,NKT,Y,prop4.G2.residue100,VS552,VS563,prop4.G3.residue25,DGP,DS,ES,VS123,prop4.G2.residue0,KTC,prop5.Tr2332,prop3.G1.residue25,TI,VS776,I,VS775,RIY,FN,prop7.G2.residue50,VS774,VS623,YDD,TL,VS723,VDS,GN,LFS,KVK,VS632,YV,prop2.G1.residue75,VS625,CLC,GPT,prop3.G1.residue100,secondarystruct.Group1,VS626,polarity.Group2,SGT,KA,DLG,V,N,NR,IGW,DN,prop3.G2.residue100,INL,IEN,VS326,VVV,prop5.Tr1221,AC,VS522,prop4.G1.residue100,DA,prop6.G1.residue50,HAR,CAG,VS261,prop2.G1.residue100,TA,VS172,YN,ELD,VS213,S,DSR,KL,prop5.G3.residue100,VS355,SDG,IDV,VS721,VS144,PE,KYV,KT,prop4.G3.residue50,VS226,FGA,prop5.G1.residue100,prop3.G2.residue50,AM,NYY,prop4.G1.residue50,GY,KVI,R,VS757,HID,SP,FG,polarity.Group1,LV,NV,QR,prop3.G1.residue50,prop3.G2.residue25,NN,prop6.G3.residue100,VS314,NI,PKS,DNK,polarizability.Group3,prop3.G1.residue75,FW,SQ,CSY,VS532,EL,DW,VS627,PQ,prop2.Tr1331,MKK,VS111,VLM,QSK,DT,NKV,KLY,CK,VS554,VS166,DC,LD,prop6.G2.residue100,FVL,ECM,prop6.Tr1331,prop5.G3.residue25,VS131,VS315,VS662,NQ,QP,VS443,VS511,prop2.Tr1221,NYT,DY,CPS,prop4.G1.residue75,AA,VS633,ESK,prop6.G3.residue0,GWC,GCI,IY,FV,CLD,FK,VS211,QS,VCL,VS444,VS413,SR,VS726,polarizability.Group1,prop3.G2.residue75,VS434,VS773,EC,VS132,EKI,prop6.G1.residue100,VS312,VS432,VS222,NC,VS641,VSV,MKT,KP,prop4.G1.residue0,PT,VR,WCC,secondarystruct.Group3,prop4.G2.residue25,FC,EP,LLK,prop4.G2.residue50,CQ,VFV,RQL,LR,VS266,RT,WY,RP,VS317,DV,VS646,VS242,TM,VS252,VS336,PP,VS126,IIN,VS262,prop4.G2.residue75,prop1.G2.residue25,CNK,LA,P,MID,LYT,KIA,EI,prop2.G1.residue25,AD,solventaccess.Group3,FYK,VS565,NE,DE,NS,MA,VS663,prop6.G2.residue50,prop7.G3.residue25,VS143,VT,prop4.Tr2332,prop6.G1.residue25,LQ,VS112,prop5.G2.residue75,prop2.G1.residue0,VS127,IL,LIN,WL,M,YCC,prop5.G1.residue50,VS125,VS256,VS356,VS544,prop6.G2.residue75,VS573,LL,SS,prop7.Tr1221,VS731,D,prop6.G2.residue25,VS621,VS331,VS343,VS323,VS713,CE,PF,VS421,YD,EH,VS115,prop6.G3.residue50,EF,VS513,LSV,GI,VS137,VS363,LST,prop6.G3.residue25,AY,CY,VI,QF,VS245,CV,VS271,AF,prop5.G1.residue25,VS415,II,AQ,VS157,VS263,CT,VS411,IQ,KE,QL,VS255,VS763,VS546,AG,GQ,VS122,PL,VS171,VS311,KQ,PAC,CL,prop5.G2.residue50,AN,CD,prop5.G1.residue75,prop6.G1.residue75,PG,VS212,CH,SFG,GS,LS,AP,prop7.G1.residue100,VS441,VPC,HG,HI,F,VS512,VDV,ESQ,VS366,VLS,ER,VS551,YT,WT,YND,DTN,END,EW,DVN,LW,VW,TPN,NY,EFN,polarizability.Group2,ACD,YKD,TT,MSD,LEC,VCC,LCC,ET,DRC,QT,TVD,LYD,prop3.Tr1221,NYD,IT,LT,ST,DSD,QPD,MT,QMD,FT,PKN,TKD,EID,CID,LHD,LCD,LMN,PIA,KKN,TV,TDR,MDR,MV,DNR,RNR,MRR,MVA,SV,VTA,ETA,VSA,NSA,CPA,QKN,VV,NPA,GMA,HKA,CKA,IAA,ACA,WLA,KLA,YIA,DHA,NIA,IV,GV,LKR,QV,CLN,EY,TIN,PIN,EIN,SIA,NGN,MY,FY,RGN,SEN,SY,RCN,HDN,GDN,KNN,YY,VY,INN,solventaccess.Group2,SVR,RV,WYR,TTR,EV,NIN,VS446,AT,ID,IE,GE,QE,EE,RE,AE,YC,WC,SC,PC,prop6.G3.residue75,IC,RC,VD,TD,SD,PD,FD,MD,LE,ME,FE,NG,prop7.Tr1331,MG,KG,LG,GG,EG,CG,DG,VQ,SE,FQ,EQ,DQ,prop1.G2.residue50,VE,YE,WE,TE,prop7.G1.residue75,HD,VS,GD,CR,DR,RR,AR,VA,YA,SA,FA,IA,GA,QA,CA,NA,RA,T,L,H,Q,prop7.G3.residue50,GR,IR,prop7.G2.residue75,MN,ED,DD,ND,RD,VN,TN,SN,PN,LN,MR,prop2.G1.residue50,QN,EN,RN,YR,TR,PR,FR,TG,YG,VG,AH,FF,MF,KF,LF,IF,HF,GF,CF,DF,NF,RF,YM,SM,FM,prop5.Tr1331,LM,GM,QM,EM,SF,TF,YF,prop4.Tr1331,YS,TS,PS,FS,KS,IS,HS,prop4.Tr1221,RS,VF,AS,VP,TP,FP,LP,IP,DP,NP,CM,prop6.Tr1221,WK,KI,DL,RL,AL,YI,SI,PI,FI,MI,LI,GL,QI,DI,RI,AI,YH,WH,IH,QH,RGC,HL,TK,QK,SK,PK,MK,KK,LK,IK,HK,GK,EK,ML,DK,NK,RK,AK,VL,YL,SL,FL,AGC,EQI,FGC,PSV,VS231,VS521,VS321,VS221,VS121,VS711,VS217,LVV,RVV,IYV,LTV,ASV,VS531,VPV,EPV,YMV,FMV,YKV,PKV,KKV,TLV,KLV,LLV,ALV,VS431,VS631,VS732,VS561,VS232,VS722,VS622,VS422,VS322,VS712,VS612,VS412,VS771,VS371,VS661,VS461,VS141,VS117,VS161,VS751,VS651,VS451,VS351,VS251,VS151,VS541,VS341,VS241,TIV,KIV,CIV,KLT,VST,IPT,KFT,LFT,KMT,VKT,YKT,SKT,LKT,VS327,VLT,LLT,FGV,ALT,VIT,KIT,IIT,NQT,GDT,NDT,FNT,INT,GNT,ENT,ITT,VS227,GYT,EVT,GWY,WPY,IPY,NPY,LFY,CMY,VS617,YIY,QIY,LGY,AGY,GCY,SDY,NDY,SNY,LNY,MAY,EPW,PFW,SGW,DEW,VVT,IVT,VS332,VS342,ELC,VS664,VS435,VS335,VS235,VS135,VS525,VS425,VS325,VS225,VS515,VS215,VS274,VS564,VS156,VS464,VS264,VS164,VS654,VS454,VS354,VS254,VS154,VS644,VS344,VS734,VS535,VS345,VS666,VS216,VS246,VS146,VS636,VS536,VS136,VS526,VS426,VS616,VS516,VS416,VS316,VS116,VS445,VS475,VS665,VS465,VS365,VS265,VS165,VS655,VS555,VS455,VS155,VS645,VS634,VS534,VS334,VS272,VS533,VS333,VS233,VS133,VS523,VS423,VS223,VS613,VS313,VS113,prop4.G3.residue100,VS562,VS234,VS462,VS466,VS162,VS652,VS566,VS452,VS352,VS152,VS742,VS642,VS542,VS733,VS243,VS543,VS643,VS134,VS624,VS424,VS324,VS224,VS124,VS614,VS514,VS414,VS214,VS114,VS373,VS456,VS173,VS556,VS463,VS163,VS653,VS553,VS453,VS353,VS253,VS153,NNT,VAT,IAT,MKI,YYI,PYI,GYI,VTI,PTI,NTI,KSI,prop5.G2.residue100,TPI,LFI,EMI,KKI,SVI,GLI,DLI,KII,III,SGI,VS346,RQI,KEI,DEI,VS567,YCI,DVI,IAL,YVS,ELL,QPL,PFL,QFL,VS267,KKL,LKL,NKL,VLL,TLL,SLL,LLL,AHL,SAL,GGL,VEL,TEL,QEL,CEL,DEL,PCL,TDL,LDL,HDL,ADL,KRI,TLH,KCH,GQE,PSE,KSE,NSE,CPE,VKE,YKE,SKE,TLE,ILE,QLE,EIE,PDE,KAH,KDE,LDE,YNE,ENE,TRE,IVC,KYC,DWC,SKC,DKC,NKC,KTE,YTE,EWE,LYE,MYG,IYG,VS177,VS277,DPG,CFG,TIG,AHG,VS477,LGG,PEG,VRG,PAG,LAG,TVQ,CKQ,SGQ,KEQ,IEQ,INQ,NRQ,RRQ,VS677,VPL,SSL,DTL,prop4.G1.residue25,GVP,GTP,ATP,LPP,YFP,PLP,SIP,LIP,HGP,VS237,VS337,KNP,LIF,VS147,IVF,EYF,AYF,SSF,LSF,PPF,KMF,VLF,PLF,DLF,LVP,KVP,TVP,IAS,MVS,EVS,KTS,KSS,VS527,SFS,DMS,KKS,VS727,ILS,FIS,NGS,TQS,TES,AES,KCS,LCS,ADS,FNS,KNS,INS,NNS,TAS,FIF,DGF,ITL,IDK,VIK,YIK,EGK,AGK,prop7.G3.residue75,SEK,KEK,IEK,VS457,VS657,TDK,NDK,YEF,TNK,QNK,NNK,IAK,TVL,SVL,LWL,YTL,FTL,KTL,LTL,ALK,TLK,NKK,DKK,ACF,VS447,ASM,TLM,PGM,KEM,LEM,YDM,IDM,VVK,SYK,KYK,GYK,VS257,TSK,VS357,PSK,DSK,TMK,NMK,SKK,EKK,CKK,DKL"""

wanted = [x.strip() for x in wanted_str.split(",") if x.strip()]

final_cols = ["Sequence_ID"] + wanted

missing = [c for c in final_cols if c not in all_features.columns]
if missing:
    print(f"WARNING: {len(missing)} requested columns not found. Example:", missing[:20])

final_cols_existing = [c for c in final_cols if c in all_features.columns]
all_features_filtered = all_features.loc[:, final_cols_existing].copy()

print("Filtered shape:", all_features_filtered.shape)

all_features_filtered.to_csv("Features.csv", index=False)
print("Saved file: Features.csv")
print("Final feature dimension:", all_features.shape)

# ==========================
# STEP 7: Save CSV
# ==========================
all_features.to_csv(output_file, index=False)
print("Saved file:", output_file)