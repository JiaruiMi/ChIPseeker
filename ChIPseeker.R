#==========================================================================================================
#
#                                       ChIPseeker包的使用
#
#==========================================================================================================

# ChIP-seq的分析流程
# quality control using FastQC
# short reads mapping using Bowtie2
# peak calling using MACS
# peak annotation using R package ChIPseeker


#-------------------------------# 使用ChIPseeker包对peak（from MACS2）进行注释 #--------------------------#

require(ChIPseeker)

library('dplyr')
library('clusterProfiler')
## ChIPseeker自带了5个BED文件，用getSampleFiles()可以拿到文件的全路径，它返回的是个named list
getSampleFiles()
f = getSampleFiles()[[4]]    # 这里取第4个文件来演示，获得bed.gz文件的路径
f

## 载入基因注释文件
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
### 这里需要注意的是，启动子区域是没有明确的定义的，所以你可能需要指定tssRegion，把基因起始转录位点的上下游区域来
### 做为启动子区域。有了这两个输入（BED文件和TxDb对象），你就可以跑注释了，然后就可以出结果了。
### 注意在这里我们载入的是hg19的注释文件，这是因为在之前运行BWA和Bowtie2的时候我们使用的是hg19的基因组版本
### 你最初用了哪个版本，后面也应该用相应的版本。

## 为什么说ChIPseeker支持所有物种？
### ChIPseeker支持所有已经有注释的物种信息，这是基于已有的TxDb包的，目前Bioconductor提供了30个TxDb包，可以供我们使用。
### 这当然只能覆盖到一小部分物种，我们的物种基因组信息，多半要从UCSC或者Ensembl获得，我敢说支持所有物种，就是因为UCSC
### 和ensembl上所有的基因组都可以被ChIPseeker支持。

## 利用GenomicFeatures包函数来制作TxDb对象：
### makeTxDbFromUCSC： 通过UCSC在线制作TxDb；
### makeTxDbFromBiomart: 通过ensembl在线制作TxDb；
### makeTxDbFromGRanges：通过GRanges对象制作TxDb；
### makeTxDbFromGFF：通过解析GFF文件制作TxDb；
### 具体参见：https://mp.weixin.qq.com/s/_OPzvaEAbiMolCA2mqJXLw

## 对peak进行基本的注释
x = annotatePeak(f, tssRegion=c(-1000, 1000), TxDb=txdb) 
x   # 告诉我们ChIPseq的位点落在基因组上什么样的区域，分布情况如何

as.GRanges(x) %>% head(3)   # 如果我想看具体的信息呢？你可以用as.GRanges方法，这里我只打印前三行
x1 <- as.data.frame(x)   # 只想输出个Excel表格，用as.data.frame就可以转成data.frame，然后你就可以用write.table输出表格了
x1  # 建议使用View进行查看，理解表格中每一列的输出的含义
### 两种不同的注释：
### 1, genomic annotation (也就是annotation这一列）：注释的是peak的位置（可以是UTR，可以是内含子或者外显子）。这种peak
### 所在的位置可能是转录调控的根本，例如可变剪切。
### 2, nearest gene annotation (也就是多出来的其它列): peak相对于转录起始位点的距离，在做基因表达调控的，promoter区域
### 是重点，离结合位点最近的gene更有可能被调控。


## 查看peak的上下游某个范围内（在这里是-5k到5k的距离）都有什么基因。
x = annotatePeak(f, tssRegion=c(-1000, 1000), TxDb=txdb, addFlankGeneInfo=TRUE, flankDistance=5000)
as.GRanges(x) %>% head(3) 
### 输出中多三列: flank_txIds, flank_geneIds和flank_gene_distances，在指定范围内所有的基因都被列出。

## 让得到的基因ID更加友好：只需要给annotatePeak传入annoDb参数就行了
## 如果你的TxDb的基因ID是Entrez，它会转成ENSEMBL，反之亦然，当然不管是那一种，都会给出SYMBOL，还有描述性的GENENAME
x = annotatePeak(f, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb = 'org.Hs.eg.db')
as.GRanges(x) %>% head(3) 


## 其他常用参数
### 参见： https://mp.weixin.qq.com/s/_OPzvaEAbiMolCA2mqJXLw


#-------------------------------# 使用ChIPseeker包对peak进行可视化 #--------------------------#
files <- getSampleFiles()
peak <- readPeakFile(files[[4]])
### heatmap：我们可以限定在一个固定的窗口里，然后把这些peak全部align起来，看看在某个窗口上的结合谱图。
### 比如说启动子区域，使用转录起始位点，然后指定上下游，使用getPromoters函数，就可以为我们准备好窗口。
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
### 再使用getTagMatrix函数，把peak比对到这个窗口，并生成矩阵，供我们分析、可视化。
tagMatrix <- getTagMatrix(peak, windows=promoter)


### 这时候，我们可以用tagHeatmap来画热图，颜色随意指定：
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
### 同时绘制多个图：
### 你可以直接从文件出发，ChIPseeker提供了peakHeatmap函数，你给文件名，就可以直接出图
### 两个函数都支持多种数据同时展示，方便比较，这里可以看出前三个样品结合位点不在转录起始区域，而后两个则是。
peakHeatmap(files, weightCol=NULL, TxDb=txdb,   # 关于weightCol，请参考CS2:BED文件
            upstream=3000, downstream=3000, 
            color=rainbow(length(files)))


### 另外我们还可以以谱图的形式来展示结合的强度(只是针对第4个测试数据的)：
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", 
            ylab = "Read Count Frequency") # 上面这个图可以非常清晰地看到，这个结合位点是集中在转录起始位点的
### ChIPseeker也提供了plotAvgProf2函数，支持以文件为输入一步做图；下面这个代码，会产生上面一样的图。
plotAvgProf2(files[[4]], TxDb=txdb, 
             upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", 
             ylab = "Read Count Frequency")
### 另外这个函数还支持使用bootstrap估计置信区间：
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), 
            conf = 0.95, resample = 1000)   # 现实95%可信区间
### 这个函数也支持多个数据的比较：
tagMatrixList <- lapply(files, getTagMatrix, 
                        windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
### 我们可以对其进行分面展示，当然也支持估计置信区间：
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), 
            conf=0.95,resample=500, facet="row")  # 这和前面的热图是一致的，前三个样品显然不是调控转录的。
### ChIPseeker还提供了getBioRegion函数，你可以指定’exon’或’intron’，它会像getPromoters一样为你准备好数据，
### 即使这不是你感兴趣的，你也依然可以用它来准备窗口，然后可视化看一下，或许你能意外地发现你研究的蛋白竟
### 然在外显子/内含子起始位置上有很强的结合谱，做为前期data exploration

plotAnnoPie(x)
plotAnnoPie.csAnno(x)
plotAnnoBar(x)
### peak落在了什么地方，启动子？外显子？内含子？还是基因间区？饼图给出了比例，当然一个peak所在的位置可能是一个基因的
### 外显子而同时又是另一个基因的内含子，所以annotatePeak有个参数genomicAnnotationPriority让你指定优先顺序，默认顺序如下:
### Promoter
### 5' UTR
### 3' UTR
### Exon
### Intron
### Downstream
### Intergenic

### 一个结合位点可以是一个基因的Exon和另一个基因的Intron，除了按照优先级别给出唯一的注释之外，ChIPseeker其实对所有的注释都有记录
### 这个图我称之为vennpie取维恩图venn plot画交集和饼图pie plot画分布的结合。
vennpie(x)

### 不过上面的图依然无法展示所有信息，所以可以综合使用UpSetR包
upsetplot(x)

### 这个可以很好的把信息全部展示出来，然而第一不够vennpie直观，第二右上角空白太多，所以我就想把vennpie强插上去，
### 最终通过参数实现这个效果
upsetplot(x, vennpie=TRUE)

### ChIPseeker还注释了最近的基因，peak离最近基因的距离分布是什么样子的？ChIPseeker提供了plotDistToTSS函数来画这个分布：
plotDistToTSS(x, title="Distribution of transcription factor-binding loci relative to TSS")

### plotAnnoBar和plotDistToTSS这两个柱状图都支持多个数据同时展示，方便比较
peakAnnoList <- lapply(files, annotatePeak, 
                       TxDb=txdb,tssRegion=c(-3000, 3000))
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

### ChIPseeker还提供了一个vennplot函数，比如我想看注释的最近基因在不同样本中的overlap
peakAnnoList
genes <- lapply(peakAnnoList, function(i) 
  as.data.frame(i)$geneId)
genes
library('Vennerable')
vennplot(genes[2:4], by='Vennerable')


#----------------------------# 使用ChIPseeker包对Genomic Coordination做富集性分析 #-----------------------#
bedfile=getSampleFiles()
seq=lapply(bedfile, readPeakFile)
seq
### seq2gene的函数，套用seq2pathway的思路，把一个基因位置上所有关联的基因全部返回来，我们可以使用它去把基因位置
### 信息转换成基因列表，然后用于富集分析
genes=lapply(seq, function(i) 
  seq2gene(i, c(-1000, 3000), 3000, TxDb=txdb))  # seq, tssRegion, flankDistance, TxDb, sameStrand

cc = compareCluster(geneClusters = genes,   # 来自clusterProfiler包
                    fun="enrichKEGG", organism="hsa")

dotplot(cc, showCategory=10)
