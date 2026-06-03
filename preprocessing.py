import pandas as pd
import numpy as np
import os

def main():
        # Read data
    df = pd.read_csv("data/raw/extracted_data.csv")

    # Bioclim covariate columns
    bioclim_cols = [c for c in df.columns if c.startswith("Bioclim_")]

    # Create a unique site ID from lon/lat
    df["site_id"] = df["lon"].astype(str) + "_" + df["lat"].astype(str)

    # Store results per class
    species_matrices = {}
    covariates = {}

    # create processed dir
    os.makedirs("data/processed", exist_ok=True)

    ### THIS PART IS IMPORTANT as it defines what species would be used in terms of intersection
    # union, minimum counts, it should be clearly stated in the paper
    # species list (This could be later be changed so that the union of species are used)
    # for now we will use only those present on PAC
    species_list = df[df["class"] == "PAC"]["scientificName"].unique()
    # also filter at least 10 occurences for PAC
    species_counts = df[df["class"] == "PAC"]["scientificName"].value_counts()
    species_list = species_counts[species_counts >= 10].index.tolist()


    # # filter so that only species from the list are used # this would either keep the sites with species from the target or not
    # df = df[df["scientificName"].isin(species_list)]

    for cls in ["PAC", "PAU", "POV", "POU"]:
        subset = df[df["class"] == cls].copy()
        
        if subset.empty:
            print(f"No data for class {cls}")
            continue

        # --- Presence-Absence matrix ---
        # One row per site, one column per species, 1 if species was observed at site
        species_matrix = (
            subset
            .groupby(["site_id", "scientificName"])
            .size()
            .unstack(fill_value=0)
            .clip(upper=1)
            .astype(np.int8)
            .reindex(columns=species_list, fill_value=0)  # ← add this
        )
        species_matrices[cls] = species_matrix

        # --- Covariates (one row per site, mean across records at that site) ---
        cov = (
            subset
            .groupby("site_id")[bioclim_cols]
            .mean()
        )
        covariates[cls] = cov

        print(f"\n[{cls}]  {species_matrix.shape[0]} sites × {species_matrix.shape[1]} species")
        print(f"       {cov.shape[1]} Bioclim covariates")

        # save results
        species_matrix.to_csv(f"data/processed/species_matrix_{cls}.csv")
        cov.to_csv(f"data/processed/covariates_{cls}.csv")

    


if __name__ == "__main__":
    main()







