import requests
import json
import pandas as pd
import os 
import time


dataset_dir = r"Q:\coding_dir\Rep_dir\VScode_Lib\testersession\porj\SPC707P\labs\Week12\Protein_Solubility_ML_Model\Dataset"
log_dir = r"Q:\coding_dir\Rep_dir\VScode_Lib\testersession\porj\SPC707P\labs\Week12\Protein_Solubility_ML_Model\logs"

# single file json test for identification of key name 
Targ_ULR = "https://rest.uniprot.org/uniprotkb/search"
params = {
    "query": "JW2231",
    "format": "json",   # full JSON
}

response = requests.get(Targ_ULR, params=params)
print("Status:", response.status_code)

data = response.json()
print(json.dumps(data, indent=4))

save_path = "JW2231_uniprot.json"

# with open(save_path, "w", encoding="utf-8") as f:
#     json.dump(data, f, indent=4)




def fetch_uniprot_seq_by_jw(jw_id):
    if pd.isna(jw_id):
        return None, "JW_ID is NaN"

    query = f"{jw_id} AND organism_id:83333"
    
    params = {
        "query": query,
        "fields": "accession,id,gene_names,organism_id,sequence",
        "format": "tsv",
        "size": 5
    }

    try:
        r = requests.get(Targ_URL, params=params)

        if r.status_code != 200:
            return None, f"HTTP {r.status_code}: {r.text[:200]}"

        text = r.text.strip()
        if not text:
            return None, "Empty response"

        lines = text.splitlines()
        if len(lines) <= 1:
            return None, "No hits"

        header = lines[0].split("\t")
        if "Sequence" not in header:
            return None, f"'Sequence' column missing. Header={header}"

        seq_idx = header.index("Sequence")
        first = lines[1].split("\t")
        seq = first[seq_idx]

        if seq == "":
            return None, "Empty sequence"

        return seq, None

    except Exception as e:
        return None, f"Exception: {str(e)}"
    

def add_sequences_with_progress(cdf, id_col="JW_ID", sleep_sec=0.2):
    df = cdf.copy()
    df["UniProt_seq"] = None

    unique_ids = df[id_col].dropna().unique()
    total = len(unique_ids)

    id_to_seq = {}
    error_log = []

    for i, jw in enumerate(unique_ids, start=1):

        seq, err = fetch_uniprot_seq_by_jw(jw)
        id_to_seq[jw] = seq        
        if err is None:
            print(f"[{i}/{total}] {jw}: OK")                # keeping  track of the progression
        else:
            print(f"[{i}/{total}] {jw}: FAIL → {err}")
            error_log.append({"JW_ID": jw, "Error": err})

        time.sleep(sleep_sec)

    df["UniProt_seq"] = df[id_col].map(id_to_seq)
    err_df = pd.DataFrame(error_log)
    return df, id_to_seq, err_df


if __name__ == "__main__":
    cdf_seq, seq_map, err_df = add_sequences_with_progress(cdf)
    cdf_seq.to_csv(os.path.join(dataset_dir, "esol_with_sequences.csv"), index=False)
    err_df.to_csv(os.path.join(log_dir, "uniprot_sequence_errors.csv"), index=False)

    print("DONE!")