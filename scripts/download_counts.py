import os
import urllib.request
import tarfile
import pandas as pd

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--gse_id', default='GSE116239')
args = parser.parse_args()
GSE_ID = args.gse_id
RAW_TAR = f"{GSE_ID}_RAW.tar"
OUT_DIR = f"{GSE_ID}_RAW"

os.makedirs(OUT_DIR, exist_ok=True)

print(f"[INFO] Working with series: {GSE_ID}")
print(f"[INFO] Archive: {RAW_TAR}")
print(f"[INFO] Output directory for extraction: {OUT_DIR}")
print(f"[INFO] Current working directory: {os.getcwd()}")

series_prefix = GSE_ID[:-3] + "nnn"
url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_prefix}/{GSE_ID}/suppl/{GSE_ID}_RAW.tar"

print(f"[INFO] URL for download: {url}")

if not os.path.exists(RAW_TAR):
    print(f"[INFO] Downloading {url} -> {RAW_TAR}")
    urllib.request.urlretrieve(url, RAW_TAR)
else:
    print(f"[INFO] {RAW_TAR} already exists")

print(f"[INFO] Extracting archive {RAW_TAR} to {OUT_DIR}")

with tarfile.open(RAW_TAR, "r") as tar:
    tar.extractall(path=OUT_DIR)

files = sorted(os.listdir(OUT_DIR))

gz_files = [f for f in files if f.endswith(".gz")]
print(f"[INFO] Found {len(gz_files)} .gz files")

expression_data = {}
meta_rows = []

for filename in gz_files:
    file_path = os.path.join(OUT_DIR, filename)
    print(f"[INFO] Reading {file_path}")

    df = pd.read_csv(file_path, sep="\t", compression="gzip")

    if "Feature" in df.columns:
        df = df.set_index("Feature")
    else:
        df = df.set_index(df.columns[0])

    sample_name = filename.replace(".txt.gz", "").replace(".gz", "")
    condition = "non_foam" if "non_foam" in filename else "foam"
    parts = filename.split("_")
    replicate = parts[-3] if len(parts) >= 3 else "1"
    expression_data[sample_name] = df
    meta_rows.append({
        "sample": sample_name,
        "filename": filename,
        "condition": condition,
        "replicate": replicate
    })
pd.to_pickle(expression_data, "expression_data.pkl")
meta_df = pd.DataFrame(meta_rows)
meta_df.to_csv("samples.tsv", sep="\t", index=False)

print("[INFO] Output files generated: expression_data.pkl and samples.tsv")
