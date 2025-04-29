library(Seurat)
library(ggplot2)
library(data.table)
library(pheatmap)
library(reshape2)

#purpose: build integrated reference, then use it to annotate new query dataset
#load("/home/exacloud/gscratch/mcweeney_lab/resources/sc_references/GSE185381_scrna_ref.RData")
#read in files for query matrix and reference matrix
query_matrix<-ReadMtx("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/GSE185381/query_matrix.mtx",cells="/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/GSE185381/query_barcodes.tsv",features="/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/GSE185381/query_features.tsv")
ref_matrix<-ReadMtx("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/GSE185381/ref_matrix.mtx",cells="/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/GSE185381/ref_barcodes.tsv",features="/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/GSE185381/ref_features.tsv")

query_metadata<-read.csv("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/GSE185381/query_metadata.csv")
query<-CreateSeuratObject(query_matrix)
query<-AddMetaData(query,query_metadata)

ref_metadata<-read.csv("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/GSE185381/ref_metadata.csv")
ref<-CreateSeuratObject(ref_matrix)
ref<-AddMetaData(ref,ref_metadata)



#query.s<- merged.s[,merged.s$samples == "Control4"]
#ref.s<-merged.s[,merged.s$samples != "Control4"]
#query<-query.s
#ref<-ref.s
#print(dim(query))
#print(dim(ref))
#preprocess without integration
ref <- NormalizeData(ref)
ref <- FindVariableFeatures(ref)
ref <- ScaleData(ref)
ref <- RunPCA(ref)
ref <- FindNeighbors(ref, dims = 1:50)
ref <- FindClusters(ref) 

#run umap and save plot
ref <- RunUMAP(ref, dims = 1:50)
p<-DimPlot(ref, group.by = c("Broad_cell_identity"))
ggsave(p,file="ref_dimplot.pdf")

#integrate datasets into shared reference and run umap
#ref <- IntegrateLayers(object = ref, method = CCAIntegration, orig.reduction = "pca",
#    new.reduction = "integrated.cca", verbose = FALSE)
#ref <- FindNeighbors(ref, reduction = "integrated.cca", dims = 1:30)
#ref <- FindClusters(ref)
#ref <- RunUMAP(ref, reduction = "integrated.cca", dims = 1:30)
#p2<-DimPlot(ref, group.by = c("celltype"))
#ggsave(p2,file="ref_integrated_dimplot.pdf")

#normalize query
#identify anchors between reference and query
#make predictions for query
query <- NormalizeData(query)
anchors <- FindTransferAnchors(reference = ref, query = query, dims = 1:30,
    reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = ref$Broad_cell_identity, dims = 1:30)
query <- AddMetaData(query, metadata = predictions)
#get the accuracy of predicted cell type annotations
query$prediction.match <- query$predicted.id == query$Broad_cell_identity
#see how well the prediction was for cell type
table(query$prediction.match)

#can create a violin plot by gene for each cell type
#VlnPlot(query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")

#unimodal umap projection
#projection of query onto reference umap 
#mapquery performs three functions: transferdata (transfers cell type labels and impute antibody derived tags values), 
#integrateembeddings (integrate reference query by correcting query's project low dimensional embeddings)
# and projectumap (project query data onto umap of reference)
ref <- RunUMAP(ref, dims = 1:50, return.model = TRUE)
query <- MapQuery(anchorset = anchors, reference = ref, query = query,
    refdata = list(celltype = "Broad_cell_identity"), reference.reduction = "pca", reduction.model = "umap")

#create visualization of query next to reference
p1 <- DimPlot(ref, reduction = "umap", group.by = "Broad_cell_identity", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
ggsave(p1 + p2,file="transferred_labels.pdf")

save(query,file="query_s.rda")
print("saved query")

#query_pred<-query$predicted.id
#query_true<-query$Broad_cell_identity

#create heatmap to show predicted versus observed cell type annotation
query_df<-cbind(query$predicted.id,query$Broad_cell_identity,query$prediction.score.B)
print("created query df")
query_df<-as.data.frame(query_df)
colnames(query_df)<-c("predicted.id","Broad_cell_identity","prediction.score.B")
query_df[,"prediction.score.B"]<-as.numeric(query_df[,"prediction.score.B"])
query_dt<-data.table(query_df)
query_hm<-query_dt[,.(mean(prediction.score.B)),by=.(Broad_cell_identity,predicted.id)]
print("created query_hm")
query_hm<-dcast(query_hm,Broad_cell_identity ~ predicted.id)
print("dcast query hm")
rownames(query_hm)<-query_hm$Broad_cell_identity
print("assigned rownames")
query_hm<-query_hm[,-1]

png("query_seurat_heatmap.png")
pheatmap(query_hm,cluster_rows=F,cluster_cols=F)
dev.off()


#split up query data frame to different broad cell identity types for plotting
cell_ids<-unique(query_df[,"Broad_cell_identity"])
query_df1<-query_df[which(query_df[,"Broad_cell_identity"] %in% cell_ids[1:9]),]
query_df2<-query_df[which(query_df[,"Broad_cell_identity"] %in% cell_ids[10:18]),]
query_df3<-query_df[which(query_df[,"Broad_cell_identity"] %in% cell_ids[19:27]),]

box1<-ggplot(query_df1,aes(x=predicted.id,y=prediction.score.B)) + geom_boxplot() + facet_wrap(vars(Broad_cell_identity)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(box1,file="seurat_boxplot1.pdf")


box2<-ggplot(query_df2,aes(x=predicted.id,y=prediction.score.B)) + geom_boxplot() + facet_wrap(vars(Broad_cell_identity)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(box2,file="seurat_boxplot2.pdf")

box3<-ggplot(query_df3,aes(x=predicted.id,y=prediction.score.B)) + geom_boxplot() + facet_wrap(vars(Broad_cell_identity)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(box3,file="seurat_boxplot3.pdf")

