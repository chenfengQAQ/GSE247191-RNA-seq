library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(pheatmap)
library(org.Hs.eg.db)
library(clusterProfiler)
library(readxl)
library(magrittr)
library(enrichplot)
library(msigdbr)
library(GSVA)
library(GSEABase)
library(limma)

setwd('D:/GEO_test/GSE247191/')

options(scipen = 200)
windowsFonts(A = windowsFont('Times New Roman'),
             B = windowsFont('Arial'))

data1<-read.table('GSE247191_raw_counts_GRCh38.p13_NCBI.tsv',sep = '\t',
                  header = TRUE)
data2<-data1[,-1]
row.names(data2)<-data1$GeneID
data3<-data.frame(read_xlsx('sample_info.xlsx',sheet = 'Sheet1'))
name_map<-setNames(data3[,1],data3[,2])
colnames(data2)<-name_map[colnames(data2)]
data4<-data2[,c(1,2,4,5,7,8)]
head(data4)
names(data4)<-c('HIV1','Control1','HIV2','Control2','HIV3','Control3')
data5<-data4[c('Control1','Control2','Control3','HIV1','HIV2','HIV3')]

coldata<-data.frame(sample = colnames(data5),
                    group = as.factor(rep(c('Control','HIV'),each=3)))
dds<-DESeqDataSetFromMatrix(countData = data5,colData = coldata,
                            design = ~group)
dds<-dds[rowSums(counts(dds))>10,]
dds<-DESeq(dds)
result<-results(dds,contrast = c('group','HIV','Control'))
norm_count<-counts(dds,normalized=TRUE)
result_df<-cbind.data.frame(as.data.frame(result),
                            as.data.frame(norm_count))
fwrite(result_df,file = 'DESeq2_all_gene.csv',row.names = TRUE)

#volcano
for(i in which(is.na(result_df$pvalue))){
  result_df$pvalue[i]<-1
}
volcano_df1<-result_df[c(row.names(subset(result_df,result_df$log2FoldChange > 0 & 
                                          result_df$pvalue < 0.05)),
                       row.names(subset(result_df,result_df$log2FoldChange < 0 & 
                                          result_df$pvalue < 0.05)),
                       row.names(subset(result_df,result_df$pvalue > 0.05 & 
                                          result_df$pvalue <= 1))),]
volcano_df1$group<-rep(c('up-regulated','down-regulated','none_diff'),
                       time = c(length(row.names(subset(result_df,result_df$log2FoldChange > 0 &
                                                          result_df$pvalue < 0.05))),
                                length(row.names(subset(result_df,result_df$log2FoldChange < 0 &
                                                          result_df$pvalue < 0.05))),
                                length(row.names(subset(result_df,result_df$pvalue > 0.05 &
                                                          result_df$pvalue <= 1)))))
volcano_df1$group<-factor(volcano_df1$group,levels = unique(volcano_df1$group),
                          ordered = TRUE)

ggplot(volcano_df1,aes(x=log2FoldChange,y=-1*log10(pvalue)))+
  geom_point(aes(color=group),shape=16,size=2)+
  scale_color_manual(values = c('#DB0000','#723097','black'),
                     labels = c(paste0('up-regulated',' ',
                                       length(row.names(subset(volcano_df1,volcano_df1$group == 'up-regulated')))),
                                paste0('down-regulated',' ',
                                       length(row.names(subset(volcano_df1,volcano_df1$group == 'down-regulated')))),
                                paste0('none_diff',' ',
                                       length(row.names(subset(volcano_df1,volcano_df1$group == 'none_diff'))))))+
  geom_hline(yintercept = -log10(0.05),linetype=2)+
  labs(y = '-log10(pvalue)',title = 'Control vs HIV Volcano')+
  # geom_vline(xintercept = c(-1,1),linetype=3)+
  theme_bw()+
  theme(axis.line = element_line('black'),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        text = element_text(family = 'A'),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = 'top',
        plot.title = element_text(size = 20,hjust = 0.5))

save(volcano_df1,file = 'volcano.Rdata')
load('volcano.Rdata')


#heatmap
result_diff<-subset(result_df,abs(result_df$log2FoldChange) > 2 &
                      result_df$padj < 0.05)
head(result_diff)
result_diff$expression_state<-factor(ifelse(result_diff$log2FoldChange>0,'up','down'),
                                     levels = c('up','down'))
result_diff_arrange<-result_diff%>%arrange(expression_state)
heatmap_df1<-result_diff_arrange[,c(7:12)]
heatmap_df2<-data.frame(t(scale(t(heatmap_df1),center = TRUE,scale = TRUE)))

pheatmap(heatmap_df2,show_rownames = TRUE,show_colnames = TRUE,
         color = colorRampPalette(c("#54B9C5", "#FFE5DA", "#E61B1B"))(100),
         cluster_cols = FALSE,cluster_rows = FALSE,
         border=FALSE,fontsize = 10,cellwidth = 15,cellheight = 10,
         annotation_names_col = FALSE,annotation_legend = TRUE,fontfamily = 'A',
         legend_breaks = c(-1,0.5,2))

gene_symbol_1<-bitr(geneID = as.numeric(row.names(data5)),
                    OrgDb = org.Hs.eg.db,
                    fromType = 'ENTREZID',
                    toType = 'SYMBOL')
gene_symbol_2<-na.omit(gene_symbol_1)
row.names(gene_symbol_2)<-gene_symbol_2$ENTREZID
gene_symbol_3<-gene_symbol_2[row.names(heatmap_df2),]

heatmap_df3<-data.frame(heatmap_df2,gene_symbol = gene_symbol_3$SYMBOL)
heatmap_df4<-heatmap_df3

for(i in 1:length(row.names(heatmap_df4))){
  heatmap_df4$gene_symbol[i]<-ifelse(is.na(heatmap_df4$gene_symbol[i]),
                                     paste0('Unknown',row.names(heatmap_df4)[i]),
                                     heatmap_df4$gene_symbol[i])
}
heatmap_df5<-heatmap_df4[,-length(colnames(heatmap_df4))]
row.names(heatmap_df5)<-heatmap_df4$gene_symbol

pheatmap(heatmap_df5[1:50,],show_rownames = TRUE,show_colnames = TRUE,
         color = colorRampPalette(c("#54B9C5", "#FFE5DA", "#E61B1B"))(100),
         cluster_cols = FALSE,cluster_rows = FALSE,
         border=FALSE,fontsize = 10,cellwidth = 15,cellheight = 10,
         annotation_names_col = FALSE,annotation_legend = TRUE,fontfamily = 'A',
         legend_breaks = c(-1,0.5,2))

save(data5,result_df,result_diff,result_diff_arrange,heatmap_df2,heatmap_df5,
     gene_symbol_2,gene_symbol_3,file = 'pheatmap.Rdata')
load('pheatmap.Rdata')

#KEGG,GO
kegg_enrich<-enrichKEGG(as.numeric(row.names(result_diff)),
                     organism = 'hsa',
                     keyType = 'kegg')
kegg_df1<-data.frame(kegg_enrich)

gobp_enrich<-enrichGO(as.numeric(row.names(result_diff)),
                   keyType = 'ENTREZID',
                   ont = 'BP',
                   OrgDb = org.Hs.eg.db)
gocc_enrich<-enrichGO(as.numeric(row.names(result_diff)),
                      keyType = 'ENTREZID',
                      ont = 'CC',
                      OrgDb = org.Hs.eg.db)
gomf_enrich<-enrichGO(as.numeric(row.names(result_diff)),
                      keyType = 'ENTREZID',
                      ont = 'MF',
                      OrgDb = org.Hs.eg.db)
gobp_df1<-data.frame(gobp_enrich)
gocc_df1<-data.frame(gocc_enrich)
gomf_df1<-data.frame(gomf_enrich)

kegg_df2<-kegg_df1[order(kegg_df1$pvalue,decreasing = TRUE),]
kegg_df2$Description<-factor(kegg_df2$Description,
                             levels = unique(kegg_df2$Description),
                             ordered = TRUE)
fwrite(kegg_df2,file = 'kegg_enrich.csv',row.names = TRUE)

ggplot(kegg_df2,aes(x=-1*log10(pvalue),y=Description))+
  geom_point(aes(size=Count,color=-1*log10(pvalue)))+
  scale_color_gradient(low="#BEB8DC",high = "red")+
  labs(x="EnrichmentScore(-log10(pvalue))",
       size = 'Numbers of genes',
       title = 'KEGG Enrichment')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust = 0.5,size = 20),
        axis.title.y = element_blank(),
        text = element_text(family = 'A'),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text.y.left = element_text(size = 18),
        title = element_text(size = 15),
        panel.grid = element_blank())

gobp_df2<-gobp_df1[order(gobp_df1$pvalue,decreasing = TRUE),]
gocc_df2<-gocc_df1[order(gocc_df1$pvalue,decreasing = TRUE),]
gomf_df2<-gomf_df1[order(gomf_df1$pvalue,decreasing = TRUE),]
fwrite(gobp_df2,file = 'gobp.csv',row.names = TRUE)
fwrite(gocc_df2,file = 'gocc.csv',row.names = TRUE)
fwrite(gomf_df2,file = 'gomf.csv',row.names = TRUE)

gobp_df3<-gobp_df2[c((length(row.names(gobp_df2))-4):length(row.names(gobp_df2))),]
gocc_df3<-gocc_df2[c((length(row.names(gocc_df2))-4):length(row.names(gocc_df2))),]
gomf_df3<-gomf_df2[c((length(row.names(gomf_df2))-4):length(row.names(gomf_df2))),]

go_merge1<-rbind.data.frame(gomf_df3,gocc_df3,gobp_df3)
go_merge1$go_type<-rep(c('MF','CC','BP'),
                       time=c(length(row.names(gomf_df3)),
                              length(row.names(gocc_df3)),
                              length(row.names(gobp_df3))))
go_merge1$Description<-factor(go_merge1$Description,levels = unique(go_merge1$Description),
                  ordered = TRUE)

ggplot(go_merge1,aes(x=Description,y=-log10(pvalue),fill=go_type))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  scale_fill_brewer(palette = 'Set2')+
  labs(y="EnrichmentScore(-log10(pvalue))",
       title = "GO Enrichment",
       fill = 'GO_TYPE')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        text = element_text(family = 'A'),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.title.x  = element_text(size = 20),
        axis.text.y.left = element_text(size = 15),
        title = element_text(size = 15),
        panel.grid = element_blank())

save(kegg_df2,gobp_df2,gocc_df2,gomf_df2,go_merge1,
     file = 'kegg_go_enrichment.Rdata')
load('kegg_go_enrichment.Rdata')

#GSEA
gsea_df1<-data.frame(log2FoldChange = result_df$log2FoldChange)
row.names(gsea_df1)<-row.names(result_df)
gsea_list1<-gsea_df1[order(gsea_df1$log2FoldChange,decreasing = TRUE),]
names(gsea_list1)<-as.character(row.names(gsea_df1))

gsea_kegg<-gseKEGG(gsea_list1)
gsea_kegg_df1<-as.data.frame(gsea_kegg)

fwrite(gsea_kegg_df1,file = 'gsea_kegg_all.csv')

gsea_kegg_df2<-subset(gsea_kegg_df1,abs(gsea_kegg_df1$NES)>1 & 
                        gsea_kegg_df1$pvalue < 0.05 & 
                        gsea_kegg_df1$qvalue < 0.25)

for(i in 1:length(row.names(gsea_kegg_df2))){
  p1<-gseaplot2(gsea_kegg,geneSetID = gsea_kegg_df2$ID[i],pvalue_table = TRUE,
                title = gsea_kegg_df2$Description[i])
  ggsave(filename = paste0('D:/GEO_test/GSE247191/GSEA/',
                           gsea_kegg_df2$ID[i],'.png'),plot = p1)
}

hallmark_dataset1<-msigdbr(species = 'Homo sapiens',category = 'H')
colnames(hallmark_dataset1)
hallmark_dataset2<-data.frame(gs_name = hallmark_dataset1$gs_name,
                              entrez_gene = hallmark_dataset1$entrez_gene)

gsea_hallmark<-GSEA(gsea_list1,TERM2GENE = hallmark_dataset2,
                    pvalueCutoff = 1)
gsea_hallmark_df1<-data.frame(gsea_hallmark)
gsea_hallmark_df2<-subset(gsea_hallmark_df1,abs(gsea_hallmark_df1$NES)>1 & 
                            gsea_hallmark_df1$pvalue<0.05 &
                            gsea_hallmark_df1$qvalue<0.25)

for(i in 1:length(row.names(gsea_hallmark_df2))){
  p1<-gseaplot2(gsea_hallmark,geneSetID = gsea_hallmark_df2$ID[i],
                pvalue_table = TRUE,
                title = gsea_hallmark_df2$Description[i])
  ggsave(filename = paste0('D:/GEO_test/GSE247191/GSEA/',
                           gsea_hallmark_df2$ID[i],'.png'),plot = p1)
}

save(gsea_df1,gsea_list1,gsea_kegg,gsea_kegg_df1,gsea_kegg_df2,
     hallmark_dataset1,hallmark_dataset2,gsea_hallmark,gsea_hallmark_df1,
     gsea_hallmark_df2,
     file = 'gsea.Rdata')
load('gsea.Rdata')

#GSVA
hallmark_dataset3<-getGmt('h.all.v2024.1.Hs.entrez.gmt')
gsva_hallmark<-gsva(expr = as.matrix(data5),gset.idx.list = hallmark_dataset3,
                    kcdf = 'Poisson',method = 'gsva',parallel.sz = 5)
gsva_hallmark_df1<-data.frame(gsva_hallmark)
gsva_hallmark_df2<-gsva_hallmark_df1[c('HIV1','HIV2','HIV3','Control1',
                                       'Control2','Control3')]

group_list<-c(rep('HIV',3),rep('Control',3))%>%factor(.,levels = c('Control','HIV'),
                                                    ordered = FALSE)
group_list<-model.matrix(~factor(group_list)+0)
colnames(group_list)<-c('Control','HIV')           
gsva_hallmark_df2.fit<-lmFit(gsva_hallmark_df2,group_list)
gsva_hallmark_df2.matrix<-makeContrasts(HIV - Control,levels = group_list)
fit<-contrasts.fit(gsva_hallmark_df2.fit,gsva_hallmark_df2.matrix)
fit2<-eBayes(fit)
gsva_hallmark_output<-topTable(fit2,number = length(row.names(gsva_hallmark_df2)))
fwrite(gsva_hallmark_output,'gsva_limma_all.csv')

# save(hallmark_dataset3,gsva_hallmark,gsva_hallmark_df1,gsva_hallmark_df2,
#      gsva_hallmark_output,
#      group_list,file = 'gsva.Rdata')
# load('gsva.Rdata')


gsva_hallmark_df3<-gsva_hallmark_output[c(row.names(subset(gsva_hallmark_output,
                                                           gsva_hallmark_output$logFC > 0 & gsva_hallmark_output$P.Value < 0.05)),
                                          row.names(subset(gsva_hallmark_output,
                                                           gsva_hallmark_output$P.Value > 0.05)),
                                          row.names(subset(gsva_hallmark_output,
                                                           gsva_hallmark_output$logFC < 0 & gsva_hallmark_output$P.Value < 0.05))),]
gsva_hallmark_df3$group<-factor(rep(c('up-regulated','none_diff','down-regulated'),
                                    time = c(length(row.names(subset(gsva_hallmark_df3,
                                                                     gsva_hallmark_df3$logFC > 0 & gsva_hallmark_df3$P.Value < 0.05))),
                                             length(row.names(subset(gsva_hallmark_df3,
                                                                     gsva_hallmark_df3$P.Value > 0.05))),
                                             length(row.names(subset(gsva_hallmark_df3,
                                                                     gsva_hallmark_df3$logFC < 0 & gsva_hallmark_df3$P.Value < 0.05))))),
                                levels = c('up-regulated','none_diff','down-regulated'))

gsva_hallmark_df4<-gsva_hallmark_df3 %>% arrange(group,desc(logFC))

gsva_hallmark_df4$hallmark_pathway<-factor(row.names(gsva_hallmark_df4),
                                           levels = unique(row.names(gsva_hallmark_df4)),
                                           ordered = TRUE)
gsva_hallmark_df5<-gsva_hallmark_df4[c(length(row.names(gsva_hallmark_df4)):1),]

gsva_hallmark_df5$hallmark_pathway<-factor(row.names(gsva_hallmark_df5),
                                           levels = unique(row.names(gsva_hallmark_df5)),
                                           ordered = FALSE)

ggplot(gsva_hallmark_df5,aes(x = logFC,y = hallmark_pathway))+
  geom_bar(aes(fill = group),stat = 'identity')+
  scale_fill_manual(values = c('red','grey','grey','blue'))+
  labs(x = 'Log2FoldChange',title = 'Control vs HIV GSVA')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(family = 'A'),
        title = element_text(size = 15),
        panel.grid = element_blank(),
        axis.ticks.y.left = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_blank())

save(hallmark_dataset3,gsva_hallmark,gsva_hallmark_df1,gsva_hallmark_df2,
     group_list,gsva_hallmark_output,gsva_hallmark_df3,gsva_hallmark_df4,
     gsva_hallmark_df5,file = 'gsva.Rdata')
load('gsva.Rdata')  