import scanpy as sc
import pandas as pd
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt

soupx=sc.read_mtx("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/RNA_soupX_data.mtx")
soupx=soupx.transpose()
soupx = remove_sparsity(soupx)
anno=pd.read_csv("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/metadata_clustering_w_header_upd.csv")
anno=anno.iloc[1:, :]
anno.set_index("NAME")
features=pd.read_csv("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/features_RNA_soupX_data.csv")
cells=pd.read_csv("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/cells_RNA_soupX_data.csv")
#print("Features DF")
#print(features.head(n=5))
#print("Cells DF")
#print(cells.head(n=5))
soupx.obs_names=cells["cells"].values.tolist()
soupx.var_names=features["Genes"].values.tolist()
#print("obs_names")
#print(soupx.obs_names)
#print("var_names")
#print(print(soupx.var_names))
#soupx.var_names=cells
#soupx.obs_names=features
anno.reindex(soupx.obs_names)
print(anno.head(n=10))
print(soupx.obs_names)
soupx.obs["inflammation_group"]=pd.Categorical(anno["inflammation_group"])
soupx.obs["Broad_cell_identity"]=pd.Categorical(anno["Broad_cell_identity"])

tmp=anno[anno['malignant'].isin(['microenvironment','Control'])]
subset_name=tmp['NAME'].values.tolist()
#print(subset_name)
subset_soupx=soupx[soupx.obs_names.isin(subset_name)].copy()
tmp2=anno[anno['malignant']=='malignant']
ref_name=tmp2['NAME'].values.tolist()
#print(ref_name)
ref_soupx=soupx[soupx.obs_names.isin(ref_name)].copy()

#print("query shape:")
#print(subset_soupx.shape)
#print("reference shape:")
#print(ref_soupx.shape)
condition_key='inflammation_group'
cell_type_key='Broad_cell_identity'
#sc.pp.normalize_total(ref_soupx)
sca.models.SCVI.setup_anndata(ref_soupx, batch_key=condition_key)
vae = sca.models.SCVI(
    ref_soupx,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)
vae.train()


reference_latent = sc.AnnData(vae.get_latent_representation())
reference_latent.obs["Broad_cell_identity"] = ref_soupx.obs[cell_type_key].tolist()
reference_latent.obs['inflammation_group'] = ref_soupx.obs[condition_key].tolist()
sc.pp.neighbors(reference_latent, n_neighbors=8)
sc.tl.leiden(reference_latent)
sc.tl.umap(reference_latent)
sc.pl.umap(reference_latent,
           color=['inflammation_group', 'Broad_cell_identity'],
           frameon=False,
           wspace=0.6,
           save="soupx_reference_umap.pdf")

ref_path = 'ref_model_soupx/'
vae.save(ref_path, overwrite=True)

model = sca.models.SCVI.load_query_data(
    subset_soupx,
    ref_path,
    freeze_dropout = True,
)
model.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['Broad_cell_identity'] = subset_soupx.obs[cell_type_key].tolist()
query_latent.obs['inflammation_group'] = subset_soupx.obs[condition_key].tolist()

sc.pp.neighbors(query_latent)
sc.tl.leiden(query_latent)
sc.tl.umap(query_latent)
plt.figure()
sc.pl.umap(
    query_latent,
    color=["inflammation_group", "Broad_cell_identity"],
    frameon=False,
    wspace=0.6,
    save = "soupx_query_umap.pdf"
)

query_path = 'query_model_soupx'
model.save(query_path, overwrite=True)

soupx_full = ref_soupx.concatenate(subset_soupx)
full_latent = sc.AnnData(model.get_latent_representation(adata=soupx_full))
full_latent.obs['Broad_cell_identity'] = soupx_full.obs[cell_type_key].tolist()
full_latent.obs['inflammation_group'] = soupx_full.obs[condition_key].tolist()

sc.pp.neighbors(full_latent)
sc.tl.leiden(full_latent)
sc.tl.umap(full_latent)
plt.figure()
sc.pl.umap(
    full_latent,
    color=["inflammation_group", "Broad_cell_identity"],
    frameon=False,
    wspace=0.6,
    save="soupx_full_umap.pdf"
)

