import scanpy as sc
import pandas as pd
import requests
from pathlib import Path

def download_and_convert_to_h5ad(tsv_path, output_dir):

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(
        tsv_path,
        sep="\t"
    )

    print("Downloading files...")
    for _, row in df.iterrows():
        url = row["url"]
        filename = row["name"]
        out_file = output_dir / filename

        if out_file.exists():
            print(f"{filename} already exists, skipping download.")
            continue

        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(out_file, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

        print(f"Downloaded {filename}")

        rename_map = {
        "barcodes": "barcodes.tsv.gz",
        "features": "features.tsv.gz",
        "matrix": "matrix.mtx.gz",
    }

    for f in output_dir.iterdir():
        for key, new_name in rename_map.items():
            if key in f.name and f.name != new_name:
                target = output_dir / new_name
                if not target.exists():
                    f.rename(target)

    print(list(output_dir.iterdir()))
    
    print("Reading 10x matrix...")
    adata = sc.read_10x_mtx(
        output_dir,
        var_names="gene_symbols",
        cache=True
    )

    adata.var_names_make_unique()
    adata.write(output_dir / "filtered_feature_bc_matrix.h5ad")

    print("Done!")

if __name__ == "__main__":
    download_and_convert_to_h5ad(
        "data/references/GEO_GSM8766920/file.tsv", # Path(s) to TSV file
        "data/references/GEO_GSM8766920/filtered_feature_bc_matrix" # Output dir
    )
