import pandas as pd
import anndata as ad
from PIL import Image
import numpy as np
from anndata import AnnData
from matplotlib import pyplot as plt
from typing import Optional, Union

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


def QC_plot_Protein(
    adata: AnnData,
    library_id: str = None,
    name: str = None,
    data_alpha: float = 0.8,
    tissue_alpha: float = 1.0,
    cmap: str = "Spectral_r",
    spot_size: tuple = (5, 40),
    show_color_bar: bool = True,
    show_size_legend: bool = True,
    show_axis: bool = False,
    cropped: bool = True,
    margin: int = 100,
    dpi: int = 150,
    output: str = None,
) -> Optional[AnnData]:
    """\
        QC plot for sptial transcriptomics data.

        Parameters
        ----------
        adata
            Annotated data matrix.
        library_id
            Library id stored in AnnData.
        data_alpha
            Opacity of the spot.
        tissue_alpha
            Opacity of the tissue.
        cmap
            Color map to use.
        spot_size
            Size of the spot (min, max).
        show_color_bar
            Show color bar or not.
        show_axis
            Show axis or not.
        show_size_legend
            Show size legend or not.
        name
            Name of the output figure file.
        output
            Save the figure as file or not.
        copy
            Return a copy instead of writing to adata.
        Returns
        -------
        Nothing
        """

    imagecol = adata.obs["imagecol"]
    imagerow = adata.obs["imagerow"]
    from sklearn.preprocessing import MinMaxScaler

    reads_per_spot = adata.to_df().sum(axis=1)
    scaler = MinMaxScaler(feature_range=spot_size)
    reads_per_spot_size = scaler.fit_transform(reads_per_spot.to_numpy().reshape(-1, 1))
    genes_per_spot = adata.to_df().astype(bool).sum(axis=1)

    # plt.rcParams['figure.dpi'] = dpi

    # Option for turning off showing figure
    plt.ioff()

    # Initialize matplotlib
    fig, a = plt.subplots()

    vmin = min(genes_per_spot)
    vmax = max(genes_per_spot)
    # Plot scatter plot based on pixel of spots
    plot = a.scatter(
        adata.obs["imagecol"],
        adata.obs["imagerow"],
        edgecolor="none",
        alpha=data_alpha,
        s=reads_per_spot_size,
        marker="o",
        vmin=vmin,
        vmax=vmax,
        cmap=plt.get_cmap(cmap),
        c=genes_per_spot,
    )

    if show_color_bar:
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        axins = inset_axes(
            a,
            width="100%",
            height="100%",
            loc="upper left",
            bbox_to_anchor=(1.0, 0.73, 0.05, 0.35),
            bbox_transform=a.transAxes,
            borderpad=4.3,
        )
        cb = plt.colorbar(plot, cax=axins)
        cb.ax.set_xlabel("Number of Proetein", fontsize=10)
        cb.ax.xaxis.set_label_coords(0.98, 1.20)
        cb.outline.set_visible(False)

    if show_size_legend:
        size_min, size_max = spot_size
        markers = [
            size_min,
            size_min + 1 / 3 * (size_max - size_min),
            size_min + 2 / 3 * (size_max - size_min),
            size_max,
        ]
        legend_markers = [plt.scatter([], [], s=i, c="grey") for i in markers]
        labels = [
            str(int(scaler.inverse_transform(np.array(i).reshape(1, 1))))
            for i in markers
        ]
        a.legend(
            handles=legend_markers,
            labels=labels,
            loc="lower left",
            bbox_to_anchor=(1, 0.05),
            scatterpoints=1,
            frameon=False,
            title="Intensity",
        )

    if not show_axis:
        a.axis("off")
    if library_id is None:
        library_id = list(adata.uns["spatial"].keys())[0]

    image = adata.uns["spatial"][library_id]["images"][
        adata.uns["spatial"][library_id]["use_quality"]
    ]
    # Overlay the tissue image
    a.imshow(
        image,
        alpha=tissue_alpha,
        zorder=-1,
    )

    if cropped:
        a.set_xlim(imagecol.min() - margin, imagecol.max() + margin)

        a.set_ylim(imagerow.min() - margin, imagerow.max() + margin)

        a.set_ylim(a.get_ylim()[::-1])
        # plt.gca().invert_yaxis()

    # fig.tight_layout()
    if output is not None:
        fig.savefig(output + "/" + name, dpi=dpi, bbox_inches="tight", pad_inches=0)

    plt.show()
