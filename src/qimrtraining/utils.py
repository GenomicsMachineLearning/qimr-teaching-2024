import pandas as pd
import anndata as ad
from PIL import Image
import numpy as np
from anndata import AnnData

def preprocess_fluo(adata: AnnData) -> pd.DataFrame:
    """
    Preprocesses fluorescence data in the given AnnData object.

    Parameters:
        adata (AnnData): Annotated data object containing fluorescence data.

    Returns:
        pd.DataFrame: Preprocessed fluorescence data.

    """
    df = adata.obsm["fluo"]

    divider = 5 * np.quantile(df, 0.2, axis=0)
    divider[divider == 0] = df.max(axis=0)[divider == 0]

    scaled = np.arcsinh(df / divider)
    return (scaled - scaled.mean(0)) / scaled.std(0)

def higher_z_score(adata: AnnData, marker_cell_dict: dict, cell_type_key: str = "cell_type"):
    """
    Assigns cell types to cells based on the highest z-score of fluorescence intensity for given markers.

    Parameters:
        adata (AnnData): Annotated data object containing fluorescence intensity information.
        marker_cell_dict (dict): Dictionary mapping marker names to cell types.
        cell_type_key (str, optional): Key to store the assigned cell types in the `adata.obs` attribute. 
            Defaults to "cell_type".

    Returns:
        None
    """
    
    adata.obsm["fluo_scale"] = preprocess_fluo(adata)
    # adata.obsm["fluo_scale"] = adata.to_df()

    markers, cell_types = list(marker_cell_dict.keys()), np.array(list(marker_cell_dict.values()))
    ct_indices = adata.obsm["fluo_scale"][markers].values.argmax(1)

    adata.obs[cell_type_key] = cell_types[ct_indices]



def Read_CODEX_QuPath(path, parent=None, library_id=None):
    '''
    Read and preprocess CODEX data from a QuPath output file.

    Parameters:
    path (str): The path to the QuPath output file.
    parent (str, optional): The parent cell ID to filter the data. Default is None.
    library_id (str, optional): The library ID for the data. Default is None.

    Returns:
    tuple: A tuple containing the library ID and the processed AnnData object.

    '''

    df = pd.read_csv(path, sep="\t")
    if parent:
        df = df.loc[df["Parent"] == parent,:]
    df.index = df.index.map(lambda x:f"cell_{x}")
    df_mtx = df.loc[:,df.columns.str.contains("Cell: Mean")]
    df_mtx.columns = df_mtx.columns.map(lambda x:x.split(":")[0])
    del df_mtx["DAPI"]

    df_obs = df.iloc[:,0:22]
    df_obs["imagecol"] = df_obs["Centroid X µm"]
    df_obs["imagerow"] = df_obs["Centroid Y µm"]

    adata = ad.AnnData(X=df_mtx, obs=df_obs)
    adata.uns["spatial"] = dict()
    if library_id:
        library_id=library_id
    elif parent:
        library_id=parent
    else:
        library_id=path.stem
        
    adata.obsm["spatial"] = df_obs[["imagecol", "imagerow"]].values
    adata.obsm["fluo"] = adata.to_df()
    # Create image
    max_size = np.max([adata.obs["imagecol"].max(), adata.obs["imagerow"].max()])
    max_size = int(max_size + 0.1 * max_size)
    image = Image.new("RGBA", (max_size, max_size))
    imgarr = np.array(image)
    # Create spatial dictionary
    adata.uns["spatial"] = {}
    adata.uns["spatial"][library_id] = {}
    adata.uns["spatial"][library_id]["images"] = {}
    adata.uns["spatial"][library_id]["images"]["hires"] = imgarr
    adata.uns["spatial"][library_id]["use_quality"] = "hires"
    adata.uns["spatial"][library_id]["scalefactors"] = {}
    adata.uns["spatial"][library_id]["scalefactors"][
        "tissue_" + "hires" + "_scalef"
    ] = 1.0
    adata.uns["spatial"][library_id]["scalefactors"][
        "spot_diameter_fullres"
    ] = 15
    return library_id, adata