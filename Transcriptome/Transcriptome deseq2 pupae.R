# Library Loading
packs = c("tximport", "rtracklayer", "DESeq2", "GenomicFeatures", "ggbio", "pheatmap", "dplyr", "tidyr",
         "clusterProfiler", "org.Dm.eg.db", "pathview", "apeglm", "ggplot2", "DEGreport", "ggrepel", "vsn")
lapply(packs, require, character.only = TRUE)
rm(packs)

setwd("C:/R/R drosophila") #press 'set as working directory' in right 'files' tap
getwd()

# Files Loading
load("./0.txdb.RData")
load("1.txisalmon.RData") #txi.salmon
load("3.u.gene.DF.RData") #u.gene.DF
load("./3_1.res.RData")

# txdb processing
txdb <- makeTxDbFromGFF(file="Drosophila_melanogaster.BDGP6.28.100.chr.gtf.gz")
transcripts(txdb)
keytypes(txdb)

k = keys(txdb, keytype = "TXNAME")
head(k)
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME") 
#normally, function 'select' might be crashed with dplyr::select
head(tx2gene)
dim(tx2gene)

dir <- getwd()
list.files(getwd())
folder.names  <- list.files(getwd())[-c(1:4)]
folder.names
files <- file.path(dir, folder.names, "quant.sf")
files
names(files) <- folder.names
names(files)

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon$counts)
head(txi.salmon$abundance)
save(txi.salmon, file="1.txisalmon.RData")

#Annotation gene name of all count
totalcount = read.csv("count.csv")
totalcount$gene = mapIds(org.Dm.eg.db,
                         keys = totalcount$X, 
                         column = "SYMBOL", 
                         keytype = "FLYBASE",
                         multiVals = "first")
write.csv(totalcount, "count_genename.csv")

# txi processing
load("1.txisalmon.RData") #txi.salmon
name.temp <- colnames(txi.salmon$counts)
name.temp

aa <- import.gff3("Drosophila_melanogaster.BDGP6.28.100.chr.gff3.gz")
head(aa)
head(aa$gene_id)
head(aa$Name) #just checking step
gene.DF <- data.frame(aa$Name, aa$gene_id)
u.gene.DF<- unique(gene.DF)
head(u.gene.DF)
save(u.gene.DF, file="3.u.gene.DF.RData")
load("3.u.gene.DF.RData")

u.gene.DF2 = na.omit(u.gene.DF)
row.names(u.gene.DF2) = 1:nrow(u.gene.DF2)
u.gene.DF2$aa.Name = as.character(u.gene.DF2$aa.Name)
u.gene.DF2$aa.gene_id = as.character(u.gene.DF2$aa.gene_id)

# letter preprocessing
tmp_table = data.frame(r.num.input = c(1:4), letter_1 = NA, letter_2 = NA)
for (i in 1:nrow(tmp_table)){
  if(i == 1){ tmp_table[i,c(2,3)] =  c("pupae-AP-F", "pupae-AP-M") }
  if(i == 2){ tmp_table[i,c(2,3)] =  c("pupae-APX-F", "pupae-DF-F") }
  if(i == 3){ tmp_table[i,c(2,3)] =  c("pupae-GF-F", "pupae-GF-M") }
  if(i == 4){ tmp_table[i,c(2,3)] =  c("pupae-AP-F", "pupae-GF-F") }
}

# define deseq_run
txi.salmon.input <- txi.salmon

#preprocessing
    res05.list = list()
    res.list = list()
    
    letter = as.character(tmp_table[r.num.input, c(2,3)])
    target.num1 = base::grep(letter[1], name.temp)
    target.num2 = base::grep(letter[2], name.temp)
    
    tt = c(rep(letter[1], length(target.num1)), rep(letter[2], length(target.num2)))
    print(paste0(r.num.input,"_tt complete") )
    
    txi.salmon.input$abundance = txi.salmon$abundance[,c(target.num1, target.num2)]
    txi.salmon.input$counts = txi.salmon$counts[,c(target.num1, target.num2)]
    txi.salmon.input$length = txi.salmon$length[,c(target.num1, target.num2)]
    
    sampleTable = data.frame(Treatment = factor(tt))
    rownames(sampleTable) = colnames(txi.salmon.input$counts)
    
    #deseq2_function
    dds = DESeqDataSetFromTximport(txi.salmon.input, sampleTable, ~Treatment)  
    dds.results = DESeq(dds)
    print(paste0(r.num.input,"_DEseq complete"))
    
    res = results(dds.results)
    summary(res)
    normalizSizeFact <- normalizationFactors(dds.results)
    filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2], "_normalized sidefactors.csv")
    write.csv(normalizSizeFact, filename)
      
    res05 = results(dds.results, alpha=0.05)
    res05.list[[r.num.input]] = res05
    print(paste0(r.num.input,"_res05.list complete"))
    res.list[[r.num.input]] = results(dds.results)
    print(paste0(r.num.input,"_res.list complete"))
    
    temp = rownames(res05)[res05$padj < 0.05]
    target.gene = temp[!is.na(temp)]
    target.gene.frame = data.frame(txi = as.character(target.gene), gene.name = NA)
    
    for (i in 1:nrow(target.gene.frame)){
      target.gene.frame[i,2] = u.gene.DF2[which(u.gene.DF2$aa.gene_id == target.gene.frame[i,1]),1]
    }
    
    #Volcano plot
    cut_lfc = 1
    cut_pvalue = 0.05
    topT = as.data.frame(na.omit(res05))
    topT$gene = mapIds(org.Dm.eg.db,
                                       keys = row.names(topT), 
                                       column = "SYMBOL", 
                                       keytype = "FLYBASE",
                                       multiVals = "first")
    topT05 = topT[topT$padj<0.05,]
    topT05$gene.name = target.gene.frame$gene.name
    tiff(filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_Volcano.tiff"), width = 10, height = 8, units = "in", res = 300)
    with(topT, plot(log2FoldChange, -log10(padj), pch=20, main=paste0(letter[1]," vs ",letter[2]), col='grey', cex=1.0, 
                    xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)), ylim=c(-25,10))
    with(subset(topT, padj<cut_pvalue & log2FoldChange>cut_lfc), points(log2FoldChange, -log10(padj), pch=20, col='red', cex=1.5))
    with(subset(topT, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), pch=20, col='blue', cex=1.5))
    abline(v=0, col='black', lty=3, lwd=1.0)
    abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)
    abline(v=cut_lfc, col='black', lty=4, lwd=2.0)
    abline(h=-log10(max(topT$padj[topT$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)## Add lines for FC and P-value cut-off for volcano plot
    dev.off()
    
    #Top20 multiple genes plot
    top20_topT05 <- topT05 %>% 
      arrange(padj) %>% 	#Arrange rows by padj values
      pull(gene.name) %>% 		#Extract character vector of ordered genes
      head(n=20) 		#Extract the first 20 genes
    top20log_topT05 <- topT05 %>% 
      arrange(padj) %>% 	#Arrange rows by padj values
      pull(log2FoldChange) %>% 		#Extract character vector of ordered genes
      head(n=20) 		#Extract the first 20 genes
    top20padj_topT05 <- topT05 %>% 
      arrange(padj) %>% 	#Arrange rows by padj values
      pull(padj) %>% 		#Extract character vector of ordered genes
      head(n=20) 		#Extract the first 20 genes
    
    top20DEG_table <- cbind(top20_topT05, top20log_topT05)
    top20DEG_table <- cbind(top20DEG_table, top20padj_topT05)
    top20DEG = as.data.frame(top20DEG_table)
    
    filename2 = paste0(r.num.input,"_",letter[1],"Vs",letter[2], "_Top20DEG.csv")
    write.csv(top20DEG, filename2)
    
    topT05_top20DEG = subset(totalcount, gene %in% top20_topT05) ##subset(totalcount, gene %in% c("Muc96D", "CG7567", "CG34282", "CG5084", "CG14300"))
    filename3 = paste0(r.num.input,"_",letter[1],"Vs",letter[2], "_Top20DEG count.csv")
    write.csv(topT05_top20DEG, filename3)
    ?pull
    topT05_20DEG_Plot = DEGreport::degPlot(dds = dds.results, genes = NULL, res = res, n = 20, xs = "Treatment", group = "Treatment") # dds object is output from DESeq2
    dev.off()
   
    #pHeatmap
    rld = rlog(dds.results)
    select = order(results(dds.results)$pvalue, decreasing=TRUE)[1:20]
    df = as.data.frame(colData(dds.results))
    if(ncol(df) != 1){if(colnames(df)[4] == "replaceable"){df = df[,-4]}}
    tiff(filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_pHeatmap rlognorm.tiff"),
         width = 10, height = 6, units = "in", res = 300)
    pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
             cluster_cols=TRUE, annotation_col=df)
    dev.off()
    

    #MAplot-evaluate the magnitude of fold changes and 
    #how they are distributed relative to mean expression
    
    resultsNames(dds.results)
    tiff(filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_MAplot.tiff"),
         width = 6, height = 6, units = "in", res = 300)
    plotMA(res, alpha = 0.05, ylim = c(-4, 4), ylab = expression(MLE ~ log[2] ~ fold ~
                                                                   change), cex = 0.6, colNonSig = rgb(0, 0, 0, 0.3), colSig = rgb(1, 0, 0, 1), colLine = NULL)
    abline(h = 0, col = "dodgerblue", lwd = 2)
    
    dev.off()
    
    resLFC <- lfcShrink(dds.results, coef=2) #log fold change shrinkage
    tiff(filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_MAplot_lfcshrink.tiff"),
         width = 6, height = 6, units = "in", res = 300)
    plotMA(resLFC, alpha = 0.05, ylim = c(-4, 4), ylab = expression(MLE ~ log[2] ~ fold ~
                                                                      change), cex = 0.6, colNonSig = rgb(0, 0, 0, 0.3), colSig = rgb(1, 0, 0, 1), colLine = NULL)
    abline(h = 0, col = "dodgerblue", lwd = 2)
    legend("bottomright", "adj. p < .05", pch = 16, col = "red", cex = 0.9, bg = "white")
    
    dev.off()
    
    #pathway analysis
    if(r.num.input %in% c(1:4)){
      try({
        #over_representation analysis Kegg pathway enrichment
        target.gene.frame$uniprot = mapIds(org.Dm.eg.db,
                                           keys = as.character(target.gene.frame$txi), 
                                           column = "UNIPROT", 
                                           keytype = "FLYBASE",
                                           multiVals = "first")
        
        eKG = enrichKEGG(target.gene.frame$uniprot, organism="dme", keyType = "uniprot",
                         pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
        head(eKG)
        print(paste0(r.num.input,"_eKG complete"))
        
        if(is.null(eKG) == FALSE){
          ggsave(plot = clusterProfiler::dotplot(eKG), filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_enrichKegg.tiff"),
                 width = 9, height = 9, units = "in", dpi = 300)
        }
        if(is.null(eKG) == FALSE){
          ggsave(plot = barplot(eKG), filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_enrichKegg_bar.tiff"),
                 width = 9, height = 9, units = "in", dpi = 300)
        }
        filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2], "_enrichKegg.csv")
        write.csv(eKG, filename)
        print(paste0(r.num.input,"_enrichKegg dot&bar plot complete"))
      })
      
      try({
        #Over_representation analysis GO pathway enrichment
        target.gene.frame$ENTREZID = mapIds(org.Dm.eg.db,
                                            keys = as.character(target.gene.frame$txi), 
                                            column = "ENTREZID", 
                                            keytype = "FLYBASE",
                                            multiVals = "first")
        
        groupgo <- groupGO(gene = target.gene.frame$ENTREZID, OrgDb = org.Dm.eg.db, ont = "BP", level = 3
                           , readable = T)
        
        ego <- enrichGO(gene = target.gene.frame$ENTREZID, OrgDb = org.Dm.eg.db, ont = "ALL",
                         readable = T, pAdjustMethod = "BH",
                         pvalueCutoff = 0.05, qvalueCutoff = 0.05, pool = TRUE)
        filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2], "_enrichgo_ALL.csv")
        write.csv(ego, filename)
        print(paste0(r.num.input,"_enrichgo complete"))
        
        if(is.null(ego) == FALSE){
          ggsave(plot = clusterProfiler::dotplot(ego, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free"), filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_enrichgo_ALL.tiff"),
                 width = 12, height = 15, units = "in", dpi = 300)
        }
        print(paste0(r.num.input,"_enrichGO dotplot complete"))
        
        #GO pathway enrichment bar plot
        if(is.null(ego) == FALSE){
          ggsave(plot = barplot(ego, drop = T, title = "GO Biological Pathways"),
                 filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_enrichgo_bar.tiff"),
                 width = 12, height = 9, units = "in", dpi = 300)
        }
        print(paste0(r.num.input,"_enrichGO barplot complete"))
      })
      try({
        #GO pathway enrichment enrichmap plot
        if(is.null(ego) == FALSE){
          ggsave(plot = emapplot(ego),
                 filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_enrichgo_enrichmap.tiff"),
                 width = 12, height = 9, units = "in", dpi = 300)
        }
        print(paste0(r.num.input,"_enrichGO enrichmap complete"))
      })
  
      try({
        #Gene set enrichment preprocessing
        topT = as.data.frame(na.omit(res05))
        topT05 <- topT[topT$padj<0.05,]
        topT05_table <- cbind(target.gene.frame$txi, topT05)
        gene_list <- topT05_table$log2FoldChange
        names(gene_list) <- topT05_table$`target.gene.frame$txi`
        is.character(topT05_table$`target.gene.frame$txi`)
        gene_list = sort(gene_list, decreasing = TRUE)
        #Gene set enrichment GO
        gse <- gseGO(geneList = gene_list, ont = "ALL", keyType = "FLYBASE",
                     minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 0.05, verbose = TRUE,
                     OrgDb = org.Dm.eg.db, pAdjustMethod = "BH")
        filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2], "_gse_ALL.csv")
        write.csv(gse, filename)
        print(paste0(r.num.input,"_gse complete"))
      })
      try({
        #gse GO Dot plot
        if(is.null(gse) == FALSE){
          ggsave(plot = dotplot(gse, showCategory=15, split=".sign") + facet_grid(ONTOLOGY~.sign), filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_gseGO_dot_All_category.tiff"),
                 width = 13, height = 15, units = "in", dpi = 300)
        }
        if(is.null(gse) == FALSE){
          ggsave(plot = dotplot(gse, showCategory=15, split=".sign") + facet_grid(~.sign), filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_gseGO_dot_All.tiff"),
                 width = 13, height = 9, units = "in", dpi = 300)
        }
        print(paste0(r.num.input,"_gene set enrichment GO dotplot complete"))
      })
       
      try({
        #gse GO enrichmap plot
        if(is.null(gse) == FALSE){
          ggsave(plot = emapplot(gse, showCategory = 10), filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_gseGO_emap_all.tiff"),
                 width = 13, height = 9, units = "in", dpi = 300)
        }
        print(paste0(r.num.input,"_gene set enrichment GO enrichmap complete"))
      })
      
      try({
        #gse GO category netplot
      if(is.null(gse) == FALSE){
        ggsave(plot = cnetplot(gse, categorySize="pvalue", foldChange=gene_list), filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_gseGO_cnet_all.tiff"),
               width = 13, height = 9, units = "in", dpi = 300)
      }
      print(paste0(r.num.input,"_gene set enrichment GO category netplot complete"))
      })
      
      try({
        #gse GO ridgeplot
      if(is.null(gse) == FALSE){
        ggsave(plot = ridgeplot(gse) + labs(x = "enrichment distribution"), filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_gseGO_ridge_all.tiff"),
               width = 13, height = 9, units = "in", dpi = 300)
      }
      print(paste0(r.num.input,"_gene set enrichment GO ridgeplot complete"))
      })
    }

    
    #PCA type 2
    if (r.num.input %in% c(1:4)){
      par(mar=c(4,4,3,3))
      plotPCA(rld, intgroup = c("Treatment"))#regularised log transformation
      pcaData = plotPCA(rld, intgroup = c("Treatment"), returnData=TRUE)
      percentVar = round(100 * attr(pcaData, "percentVar"))
      pca_plot_rld = ggplot(pcaData, aes(x = PC1, y = PC2, color=Treatment)) + 
        scale_alpha_manual(values=c(0.4,1)) +
        scale_color_manual(values=c("red","blue")) +
        geom_point(size=3) + 
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        coord_fixed(ratio=1) #ratio=y/x, alpha means transparency
      ggsave(filename = paste0(r.num.input,"_",letter[1],"Vs",letter[2],"_PCA_rld.tiff"), 
             plot = pca_plot_rld, dpi = 300, height = 6, width = 6, units = "in")
      
      dev.off()
    }
    print(paste0(r.num.input,"_PCA plot complete"))

    #save
    write.csv(target.gene.frame, paste0(r.num.input,"_",letter[1],"Vs",letter[2],".csv"), quote = F, row.names = F)
    print(paste0(r.num.input,"_save.csv complete"))

#RUN
for (i in c(1:4)){
deseq_run_with_letter(i,txi.salmon.input)
}
dev.off()

# Extract deseq2 geneDB
getwd()
file.n <- list.files(pattern="*.csv")
data.frame(file.n)
for(i in 1:length(file.n)){
  temp <- read.csv(file.n[i], stringsAsFactors = F)
  temp$order <- i
  if(i == 1){
    temp2 <- temp
  } else {
    temp2 = rbind(temp2, temp)
  }
}

total.gene <- temp2  

load("1.txisalmon.RData") #txi.salmon
str(txi.salmon)

for(k in 1:length(total.gene[,1])){
  t1 <- (1:length(rownames(txi.salmon$counts)))[total.gene[k,1] == rownames(txi.salmon$counts)]
  if(k == 1){
    t2 <- t1
  } else {
    t2 <- c(t2, t1)
  }
}
gene.order <- t2

table(total.gene[,1] == rownames(txi.salmon$counts)[gene.order]) # =All true

F.DB <- cbind(total.gene, txi.salmon$counts[gene.order,])
head(F.DB)
tail(F.DB)
unique(F.DB$order)

## Check sub DB
file.n
targe.num = 6
file.n[targe.num]

Gene.DB.case <- subset(F.DB, F.DB$order == targe.num)
Gene.DB.case

## Extract all sub geneDB
for(i in 1:length(file.n)){
  Gene.DB.case <- subset(F.DB, F.DB$order == i)
  write.csv(Gene.DB.case, paste0("GeneDB_case_new",file.n[i]))
} # must delete un-involved groups in the csv file # must delete column 'order' in the csv file

## heatmap version 2
setwd("C:/R/R drosophila/GeneDB_case and heatmap")
getwd()
dir()
file.GeneDB <- list.files(pattern="*.csv")
data.frame(file.GeneDB)

library(RColorBrewer)
coul <- colorRampPalette(c("blue", "black", "yellow"))(100)
tiff(filename = paste0("heatmapcolor.tiff"), width = 10, height = 10, units = "in", res = 300)
plot(1:length(coul), rep(1,length(coul)), col=coul, pch=19, cex=10,)
dev.off()
#coul <- colorRampPalette(c("forestgreen", "black", "red"))(100)
#plot(1:length(coul), rep(1,length(coul)), col=coul, pch=19, cex=10,)

for(i in 1:length(file.GeneDB)){
  af <- read.csv(file.GeneDB[i], header = TRUE)
  row.names(af) <- af$gene.name
  af2 <- af[,-c(1:5)]
  
  tiff(filename = paste0(file.GeneDB[i], "_heatmap_hierarchical.tiff"), width = 10, height = 12, units = "in", res = 300)
  #heatmap(as.matrix(af2), Colv=NA, Rowv=NA, col=coul) #labRow = ""
  #heatmap(as.matrix(af2), Colv=NA, col=coul)
  heatmap(as.matrix(af2), col=coul)
  dev.off()
  
}