import anndata as ad

def sampleFromClusters(adata):
    combined_data.obs["original_index"] = combined_data.obs.index
    xs= []
    ys = []
    chosen_idx = []
    colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive']
    clusters = []
    for cluster in list(set(combined_data.obs["leiden"])):
        clusters.append(cluster)
        subset = combined_data.obs[combined_data.obs["leiden"] == cluster]
        
        random_element = subset.sample(n=1)
        random_index = combined_data[random_element.index].obs["original_index"]
        chosen_idx.append(random_index)
        x, y = combined_data[random_index].obsm["X_umap"][0]
    #     plt.imshow(np.array(im_dataset[int(random_index)].squeeze(dim=0)))
    #     plt.title(f"{cluster} - {combined_data[random_index].obs['pattern'][0]} x:{x:.2f} y:{y:.2f}")
        xs.append(x)
        ys.append(y)
    return xs, ys

