import pandas as pd


# define dirs
dataset_dir = r"Q:\coding_dir\Rep_dir\VScode_Lib\testersession\porj\SPC707P\labs\Week12\Protein_Solubility_ML_Model\Dataset"
log_dir = r"Q:\coding_dir\Rep_dir\VScode_Lib\testersession\porj\SPC707P\labs\Week12\Protein_Solubility_ML_Model\logs"



d1 = pd.read_csv(os.path.join(dataset_dir, "esol.csv"))


cols = ["Solubility (%)", "Calculated MW (kDa)", "Calculated pI", "Yield (uM)", "Yield (ug/ml)", "Type of gene product", "Cell location"]
crit_col = ["Solubility (%)", "Calculated MW (kDa)", "Calculated pI"]
optional_cols = ["Yield (uM)", "Yield (ug/ml)", "Type of gene product", "Cell location"]


cdf = d1.dropna(subset=["Solubility (%)"]).copy()
cdf_woptional = cdf[optional_cols].copy()