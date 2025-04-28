import scanpy as sc
import pandas as pd
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import scvi
from scvi.data.fields import LayerField
from scvi.data import AnnDataManager

soupx=sc.read_mtx("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/RNA_soupX1.mtx")
soupx=soupx.transpose()
soupx = remove_sparsity(soupx)
anno=pd.read_csv("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/metadata_clustering_w_header_upd.csv")
anno=anno.iloc[1:, :]
anno.set_index("NAME")
features=pd.read_csv("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/features_RNA_soupX1.csv")
cells=pd.read_csv("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/cells_RNA_soupX1.csv")
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
sca.models.SCVI.setup_anndata(ref_soupx, batch_key=condition_key,labels_key=cell_type_key)
vae = sca.models.SCVI(
    ref_soupx,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)
vae.train()

scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")
print("Labelled Indices: ", len(scanvae._labeled_indices))
print("Unlabelled Indices: ", len(scanvae._unlabeled_indices))
scanvae.train(max_epochs=20)

reference_latent = sc.AnnData(scanvae.get_latent_representation())
reference_latent.obs["Broad_cell_identity"] = ref_soupx.obs[cell_type_key].tolist()
reference_latent.obs['inflammation_group'] = ref_soupx.obs[condition_key].tolist()
sc.pp.neighbors(reference_latent, n_neighbors=8)
sc.tl.leiden(reference_latent)
sc.tl.umap(reference_latent)
sc.pl.umap(reference_latent,
           color=['inflammation_group', 'Broad_cell_identity'],
           frameon=False,
           wspace=0.6,
           save="reference_soupx1_scanvi_umap.pdf")

reference_latent.obs['predictions'] = scanvae.predict()
print("Acc: {}".format(np.mean(reference_latent.obs.predictions == reference_latent.obs.Broad_cell_identity)))


df = reference_latent.obs.groupby(["Broad_cell_identity", "predictions"]).size().unstack(fill_value=0)
norm_df = df / df.sum(axis=0)

plt.figure(figsize=(20, 20))
_ = plt.pcolor(norm_df)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("Predicted")
plt.ylabel("Observed")
plt.savefig("accuracy_scanvi_soupx1_reference.pdf")

ref_path = 'ref_model_soupx1_scanvi/'
scanvae.save(ref_path, overwrite=True)

subset_soupx.obs['orig_cell_types'] = subset_soupx.obs[cell_type_key].copy() 
subset_soupx.obs[cell_type_key] = scanvae.unlabeled_category_

model = sca.models.SCANVI.load_query_data(
    subset_soupx,
    ref_path,
    freeze_dropout = True
)

model._unlabeled_indices = np.arange(subset_soupx.n_obs)
model._labeled_indices = []

print("Labelled Indices: ", len(model._labeled_indices))
print("Unlabelled Indices: ", len(model._unlabeled_indices))

model.train(
    max_epochs=100,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['Broad_cell_identity'] = subset_soupx.obs['orig_cell_types'].tolist()
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
    save = "query_soupx1_scanvi_umap.pdf"
)

query_path = 'query_model_soupx1_scanvi'
model.save(query_path, overwrite=True)

query_latent.obs['predictions'] = model.predict()
print("Acc: {}".format(np.mean(query_latent.obs.predictions == query_latent.obs.Broad_cell_identity)))

df = query_latent.obs.groupby(["Broad_cell_identity", "predictions"]).size().unstack(fill_value=0)
norm_df = df / df.sum(axis=0)

plt.figure(figsize=(20, 20))
_ = plt.pcolor(norm_df)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("Predicted")
plt.ylabel("Observed")
plt.savefig("accuracy_scanvi_soupx1_query.pdf")


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
    save="full_soupx1_scanvi_umap.pdf"
)

full_latent.obs['predictions'] = pd.concat([reference_latent.obs['predictions'],query_latent.obs['predictions']]).to_list()
df = full_latent.obs.groupby(["Broad_cell_identity", "predictions"]).size().unstack(fill_value=0)
norm_df = df / df.sum(axis=0)

plt.figure(figsize=(20, 20))
_ = plt.pcolor(norm_df)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("Predicted")
plt.ylabel("Observed")
plt.savefig("accuracy_scanvi_soupx1_full.pdf")

print("Acc: {}".format(np.mean(full_latent.obs.predictions == full_latent.obs.Broad_cell_identity)))

