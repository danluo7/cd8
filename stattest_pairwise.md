# comparing CD8 to non-CD8, using stattest function, and output a heatmap

#Need to make a new header file and recode the "type" header,since this will determine what gets compared in a stattest. the stattest function will use type as a covariate and use fpkm as a meansurement. since this function can't compare multiple things, need to make another file called GBM049_all_stattest.csv and make the "type" tissue vs. non-tissue. Then heatmaps can be done comparing invitro to everything else.


	mkdir -p $cd8/de/ballgown/ref_only/stattest
	cd $cd8/de/ballgown/ref_only/stattest
	
printf "\"ids\",\"type\",\"path
\"\n\"1_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/1_WTCD8
\"\n\"2_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/2_WTCD8
\"\n\"3_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/3_WTCD8

\"\n\"4_SK2KOCD8\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/4_SK2KOCD8
\"\n\"5_SK2KOCD8\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/5_SK2KOCD8
\"\n\"6_SK2KOCD8\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/6_SK2KOCD8

\"\n\"7_MDSCWT\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/7_MDSCWT
\"\n\"8_MDSCWT\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/8_MDSCWT
\"\n\"9_MDSCWT\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/9_MDSCWT

\"\n\"10_MDSCSK2KO\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/10_MDSCSK2KO
\"\n\"11_MDSCSK2KO\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/11_MDSCSK2KO
\"\n\"12_MDSCSK2KO\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/12_MDSCSK2KO


	printf "\"ids\",\"type\",\"path\"\n\"1_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/1_WTCD8\"\n\"2_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/2_WTCD8\"\n\"3_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/3_WTCD8\"\n\"4_SK2KOCD8\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/4_SK2KOCD8\"\n\"5_SK2KOCD8\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/5_SK2KOCD8\"\n\"6_SK2KOCD8\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/6_SK2KOCD8\"\n\"7_MDSCWT\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/7_MDSCWT\"\n\"8_MDSCWT\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/8_MDSCWT\"\n\"9_MDSCWT\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/9_MDSCWT\"\n\"10_MDSCSK2KO\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/10_MDSCSK2KO\"\n\"11_MDSCSK2KO\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/11_MDSCSK2KO\"\n\"12_MDSCSK2KO\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/12_MDSCSK2KO\"\n" > SK2_all_stattest.csv
	
	
#now rerun all the R scripts to do stattest and heat map.



	R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	library(ggplot2)
	library(gplots)
	library(GenomicRanges)

	pheno_data = read.csv("SK2_all_stattest.csv")  

	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg

	bg_table = texpr(bg, 'all')

	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)

	save(bg, file='bg.rda')

	bg

	pdf(file="SK2_R_output_stattest.pdf")


	dir()

	load('bg.rda')

	bg_table = texpr(bg, 'all')
	bg_gene_names = unique(bg_table[, 9:10])

	gene_expression = as.data.frame(gexpr(bg))
	head(gene_expression)

	colnames(gene_expression)
	dim(gene_expression)

	i = row.names(gene_expression) == "BRD4"
	gene_expression[i,]

	transcript_gene_table = indexes(bg)$t2g
	head(transcript_gene_table)

	length(row.names(transcript_gene_table)) #Transcript count
	length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count

	counts=table(transcript_gene_table[,"g_id"])
	c_one = length(which(counts == 1))
	c_more_than_one = length(which(counts > 1))
	c_max = max(counts)
	hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
	legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
	legend("topright", legend_text, lty=NULL)

	full_table <- texpr(bg , 'all')
	hist(full_table$length, breaks=500, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

	min(gene_expression[,"FPKM.1_WTCD8"])
	max(gene_expression[,"FPKM.2_WTCD8"])

	min_nonzero=1
	
	
#make sure the columns AND short names are correct, this is for the names of columns in the heatmap 


	data_columns=c(1:12) 		
	short_names=c("CD8_WT","CD8_WT","CD8_WT","CD8_SK2KO","CD8_SK2KO","CD8_SK2KO","MDSC_WT","MDSC_WT","MDSC_WT","MDSC_SK2KO","MDSC_SK2KO","MDSC_SK2KO")


#colors()

data_colors=c("black","black","black","green","green","green","red","red","red","blue","blue","blue")


boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for libraries of all 12 samples")





gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

i = which(gene_expression[,"sum"] > 10)

r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
r



## Calculate the differential expression results including significance

	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))


## View the distribution of differential expression values as a histogram

#Display only those that are significant according to Ballgown

	sig=which(results_genes$pval<0.05)
	results_genes[,"de"] = log2(results_genes[,"fc"])
	hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) CD8_WT vs non-CD8WT", main="Distribution of differential expression values")
	abline(v=-2, col="black", lwd=2, lty=2)
	abline(v=2, col="black", lwd=2, lty=2)
	legend("topleft", "Fold-change > 4", lwd=2, lty=2)
	
	#Display the grand expression values from WTCD8 vs non-WTCD8 and mark those that are significantly differentially expressed. Make sure all the columns are correctly indicated. 

	gene_expression[,"CD8_WT"]=apply(gene_expression[,c(1:3)], 1, mean)
	gene_expression[,"non-CD8_WT"]=apply(gene_expression[,c(4:12)], 1, mean)

	x=log2(gene_expression[,"CD8_WT"]+min_nonzero)
	y=log2(gene_expression[,"non-CD8_WT"]+min_nonzero)
	plot(x=x, y=y, pch=16, cex=0.25, xlab="FPKM (log2)", ylab="FPKM (log2)", main="tissue vs invitro FPKMs")
	abline(a=0, b=1)
	xsig=x[sig]
	ysig=y[sig]
	points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
	legend("topleft", "Significant", col="magenta", pch=16)



Get the gene symbols for the top 25 genes (according to corrected p-value) and display them on the plot
	
	topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
	topn = order(results_genes[sig,"qval"])[1:25]
	text(x[topn], y[topn], results_genes[topn,"gene_name"], col="black", cex=0.75, srt=45)


Write a simple table of differentially expressed transcripts to an output file

Each should be significant with a log2 fold-change >= 2

	sigpi = which(results_genes[,"pval"]<0.05)
	sigp = results_genes[sigpi,]
	sigde = which(abs(sigp[,"de"]) >= 2)
	sig_tn_de = sigp[sigde,]


	o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE) #Order the output by or p-value and then break ties using fold-change

	output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
	write.table(output, file="SigDE_R_ballgown_CD8_WT_vs_non-CD8_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)

#View selected columns of the first 25 lines of output

	output[1:25,c(1,4,5)]
	
	
	
	
## Create a heatmap to vizualize expression differences between the 16 samples
	

#Define custom dist and hclust functions for use with heatmaps

	mydist=function(c) {dist(c,method="euclidian")}
	myclust=function(c) {hclust(c,method="average")}

	main_title="sig DE Transcripts for CD8_WT vs. non-CD8_WT"
	par(cex.main=0.8)
	sig_genes_de=sig_tn_de[,"id"]
	sig_gene_names_de=sig_tn_de[,"gene_name"]

	data=log2(as.matrix(gene_expression[as.vector(sig_genes_de),data_columns])+1)
	heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(10,4), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names_de,labCol=short_names,col=rev(heat.colors(75)))


	dev.off()
	
	
	
	
	
	
	
	
	
	


# Then do stattest for alternative analysis: 1) CD8 wildtype vs. SK2 knockout, 2) MDSC wildtype vs. SK2 knockout

	mkdir -p $cd8/de/ballgown/ref_only/stattest_cd8
	cd $cd8/de/ballgown/ref_only/stattest_cd8
	
printf "\"ids\",\"type\",\"path
\"\n\"1_WTCD8\",\"CD8_WT\",\"$cd8/expression/stringtie/ref_only/1_WTCD8
\"\n\"2_WTCD8\",\"CD8_WT\",\"$cd8/expression/stringtie/ref_only/2_WTCD8
\"\n\"3_WTCD8\",\"CD8_WT\",\"$cd8/expression/stringtie/ref_only/3_WTCD8

\"\n\"4_SK2KOCD8\",\"CD8_SK2KO\",\"$cd8/expression/stringtie/ref_only/4_SK2KOCD8
\"\n\"5_SK2KOCD8\",\"CD8_SK2KO\",\"$cd8/expression/stringtie/ref_only/5_SK2KOCD8
\"\n\"6_SK2KOCD8\",\"CD8_SK2KO\",\"$cd8/expression/stringtie/ref_only/6_SK2KOCD8

\"\n" > cd8wt_vs_sk2ko.csv


	printf "\"ids\",\"type\",\"path\"\n\"1_WTCD8\",\"CD8_WT\",\"$cd8/expression/stringtie/ref_only/1_WTCD8\"\n\"2_WTCD8\",\"CD8_WT\",\"$cd8/expression/stringtie/ref_only/2_WTCD8\"\n\"3_WTCD8\",\"CD8_WT\",\"$cd8/expression/stringtie/ref_only/3_WTCD8\"\n\"4_SK2KOCD8\",\"CD8_SK2KO\",\"$cd8/expression/stringtie/ref_only/4_SK2KOCD8\"\n\"5_SK2KOCD8\",\"CD8_SK2KO\",\"$cd8/expression/stringtie/ref_only/5_SK2KOCD8\"\n\"6_SK2KOCD8\",\"CD8_SK2KO\",\"$cd8/expression/stringtie/ref_only/6_SK2KOCD8\"\n" > cd8wt_vs_sk2ko.csv
	
	#stattest and heat map.



	R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	library(ggplot2)
	library(gplots)
	library(GenomicRanges)

	pheno_data = read.csv("cd8wt_vs_sk2ko.csv")  

	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg

	bg_table = texpr(bg, 'all')

	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)

	save(bg, file='bg.rda')

	bg

	pdf(file="cd8wt_vs_sk2ko.pdf")


	dir()

	load('bg.rda')

	bg_table = texpr(bg, 'all')
	bg_gene_names = unique(bg_table[, 9:10])

	gene_expression = as.data.frame(gexpr(bg))
	head(gene_expression)

	colnames(gene_expression)
	dim(gene_expression)

	i = row.names(gene_expression) == "BRD4"
	gene_expression[i,]

	transcript_gene_table = indexes(bg)$t2g
	head(transcript_gene_table)

	length(row.names(transcript_gene_table)) #Transcript count
	length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count

	counts=table(transcript_gene_table[,"g_id"])
	c_one = length(which(counts == 1))
	c_more_than_one = length(which(counts > 1))
	c_max = max(counts)
	hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
	legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
	legend("topright", legend_text, lty=NULL)

	full_table <- texpr(bg , 'all')
	hist(full_table$length, breaks=500, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

	min(gene_expression[,"FPKM.1_WTCD8"])
	max(gene_expression[,"FPKM.2_WTCD8"])

	min_nonzero=1
	
	
#make sure the columns AND short names are correct, this is for the names of columns in the heatmap 


	data_columns=c(1:6) 		
	short_names=c("CD8_WT","CD8_WT","CD8_WT","CD8_SK2KO","CD8_SK2KO","CD8_SK2KO")


#colors()

data_colors=c("green","green","green","red","red","red")


boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for libraries of all 6 samples")





gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

i = which(gene_expression[,"sum"] > 10)

r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
r



## Calculate the differential expression results including significance

	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))


## View the distribution of differential expression values as a histogram

#Display only those that are significant according to Ballgown

	sig=which(results_genes$pval<0.05)
	results_genes[,"de"] = log2(results_genes[,"fc"])
	hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) CD8_WT vs non-CD8WT", main="Distribution of differential expression values")
	abline(v=-2, col="black", lwd=2, lty=2)
	abline(v=2, col="black", lwd=2, lty=2)
	legend("topleft", "Fold-change > 4", lwd=2, lty=2)
	
	#Display the grand expression values from CD8_WT vs CD8_SK2KO and mark those that are significantly differentially expressed. Make sure all the columns are correctly indicated. 

	gene_expression[,"CD8_WT"]=apply(gene_expression[,c(1:3)], 1, mean)
	gene_expression[,"non-CD8_WT"]=apply(gene_expression[,c(4:6)], 1, mean)

	x=log2(gene_expression[,"CD8_WT"]+min_nonzero)
	y=log2(gene_expression[,"non-CD8_WT"]+min_nonzero)
	plot(x=x, y=y, pch=16, cex=0.25, xlab="FPKM (log2)", ylab="FPKM (log2)", main="CD8_WT vs CD8_SK2KO FPKMs")
	abline(a=0, b=1)
	xsig=x[sig]
	ysig=y[sig]
	points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
	legend("topleft", "Significant", col="magenta", pch=16)



Get the gene symbols for the top 25 genes (according to corrected p-value) and display them on the plot
	
	topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
	topn = order(results_genes[sig,"qval"])[1:25]
	text(x[topn], y[topn], results_genes[topn,"gene_name"], col="black", cex=0.75, srt=45)


Write a simple table of differentially expressed transcripts to an output file

Each should be significant with a log2 fold-change >= 2

	sigpi = which(results_genes[,"pval"]<0.05)
	sigp = results_genes[sigpi,]
	sigde = which(abs(sigp[,"de"]) >= 2)
	sig_tn_de = sigp[sigde,]


	o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE) #Order the output by or p-value and then break ties using fold-change

	output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
	write.table(output, file="SigDE_R_ballgown_CD8_WT_vs_SK2KO.txt", sep="\t", row.names=FALSE, quote=FALSE)

#View selected columns of the first 25 lines of output

	output[1:25,c(1,4,5)]
	
	
	
	
## Create a heatmap to vizualize expression differences between the 6 samples
	

#Define custom dist and hclust functions for use with heatmaps

	mydist=function(c) {dist(c,method="euclidian")}
	myclust=function(c) {hclust(c,method="average")}

	main_title="sig DE Transcripts for CD8_WT vs. CD8_SK2KO"
	par(cex.main=0.8)
	sig_genes_de=sig_tn_de[,"id"]
	sig_gene_names_de=sig_tn_de[,"gene_name"]

	data=log2(as.matrix(gene_expression[as.vector(sig_genes_de),data_columns])+1)
	heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(10,4), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names_de,labCol=short_names,col=rev(heat.colors(75)))


	dev.off()
	
	







	
# 2) MDSC wildtype vs. SK2 knockout


	
	mkdir -p $cd8/de/ballgown/ref_only/stattest_mdsc
	cd $cd8/de/ballgown/ref_only/stattest_mdsc


printf "\"ids\",\"type\",\"path

\"\n\"7_MDSCWT\",\"MDSC_WT\",\"$cd8/expression/stringtie/ref_only/7_MDSCWT
\"\n\"8_MDSCWT\",\"MDSC_WT\",\"$cd8/expression/stringtie/ref_only/8_MDSCWT
\"\n\"9_MDSCWT\",\"MDSC_WT\",\"$cd8/expression/stringtie/ref_only/9_MDSCWT

\"\n\"10_MDSCSK2KO\",\"MDSC_SK2KO\",\"$cd8/expression/stringtie/ref_only/10_MDSCSK2KO
\"\n\"11_MDSCSK2KO\",\"MDSC_SK2KO\",\"$cd8/expression/stringtie/ref_only/11_MDSCSK2KO
\"\n\"12_MDSCSK2KO\",\"MDSC_SK2KO\",\"$cd8/expression/stringtie/ref_only/12_MDSCSK2KO

\"\n" > mdsc_wt_vs_sk2ko.csv


	printf "\"ids\",\"type\",\"path\"\n\"7_MDSCWT\",\"MDSC_WT\",\"$cd8/expression/stringtie/ref_only/7_MDSCWT\"\n\"8_MDSCWT\",\"MDSC_WT\",\"$cd8/expression/stringtie/ref_only/8_MDSCWT\"\n\"9_MDSCWT\",\"MDSC_WT\",\"$cd8/expression/stringtie/ref_only/9_MDSCWT\"\n\"10_MDSCSK2KO\",\"MDSC_SK2KO\",\"$cd8/expression/stringtie/ref_only/10_MDSCSK2KO\"\n\"11_MDSCSK2KO\",\"MDSC_SK2KO\",\"$cd8/expression/stringtie/ref_only/11_MDSCSK2KO\"\n\"12_MDSCSK2KO\",\"MDSC_SK2KO\",\"$cd8/expression/stringtie/ref_only/12_MDSCSK2KO\"\n" > mdsc_wt_vs_sk2ko.csv


	#stattest and heat map.



	R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	library(ggplot2)
	library(gplots)
	library(GenomicRanges)

	pheno_data = read.csv("mdsc_wt_vs_sk2ko.csv")  

	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg

	bg_table = texpr(bg, 'all')

	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)

	save(bg, file='bg.rda')

	bg

	pdf(file="mdsc_wt_vs_sk2ko.pdf")


	dir()

	load('bg.rda')

	bg_table = texpr(bg, 'all')
	bg_gene_names = unique(bg_table[, 9:10])

	gene_expression = as.data.frame(gexpr(bg))
	head(gene_expression)

	colnames(gene_expression)
	dim(gene_expression)

	i = row.names(gene_expression) == "BRD4"
	gene_expression[i,]

	transcript_gene_table = indexes(bg)$t2g
	head(transcript_gene_table)

	length(row.names(transcript_gene_table)) #Transcript count
	length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count

	counts=table(transcript_gene_table[,"g_id"])
	c_one = length(which(counts == 1))
	c_more_than_one = length(which(counts > 1))
	c_max = max(counts)
	hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
	legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
	legend("topright", legend_text, lty=NULL)

	full_table <- texpr(bg , 'all')
	hist(full_table$length, breaks=500, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

	min(gene_expression[,"FPKM.7_MDSCWT"])
	max(gene_expression[,"FPKM.8_MDSCWT"])

	min_nonzero=1
	
	
#make sure the columns AND short names are correct, this is for the names of columns in the heatmap 


	data_columns=c(1:6) 		
	short_names=c("MDSC_WT","MDSC_WT","MDSC_WT","MDSC_SK2KO","MDSC_SK2KO","MDSC_SK2KO")


#colors()

data_colors=c("blue","blue","blue","magenta","magenta","magenta")


boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for libraries of all 6 samples")





gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

i = which(gene_expression[,"sum"] > 10)

r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
r



## Calculate the differential expression results including significance

	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))


## View the distribution of differential expression values as a histogram

#Display only those that are significant according to Ballgown

	sig=which(results_genes$pval<0.05)
	results_genes[,"de"] = log2(results_genes[,"fc"])
	hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) CD8_WT vs non-CD8WT", main="Distribution of differential expression values")
	abline(v=-2, col="black", lwd=2, lty=2)
	abline(v=2, col="black", lwd=2, lty=2)
	legend("topleft", "Fold-change > 4", lwd=2, lty=2)
	
	#Display the grand expression values from MDSC_WT vs MDSC_Sk2KO and mark those that are significantly differentially expressed. Make sure all the columns are correctly indicated. 

	gene_expression[,"CD8_WT"]=apply(gene_expression[,c(1:3)], 1, mean)
	gene_expression[,"non-CD8_WT"]=apply(gene_expression[,c(4:6)], 1, mean)

	x=log2(gene_expression[,"CD8_WT"]+min_nonzero)
	y=log2(gene_expression[,"non-CD8_WT"]+min_nonzero)
	plot(x=x, y=y, pch=16, cex=0.25, xlab="FPKM (log2)", ylab="FPKM (log2)", main="MDSC_WT vs MDSC_Sk2KO FPKMs")
	abline(a=0, b=1)
	xsig=x[sig]
	ysig=y[sig]
	points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
	legend("topleft", "Significant", col="magenta", pch=16)



Get the gene symbols for the top 25 genes (according to corrected p-value) and display them on the plot
	
	topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
	topn = order(results_genes[sig,"qval"])[1:25]
	text(x[topn], y[topn], results_genes[topn,"gene_name"], col="black", cex=0.75, srt=45)


Write a simple table of differentially expressed transcripts to an output file

Each should be significant with a log2 fold-change >= 2

	sigpi = which(results_genes[,"pval"]<0.05)
	sigp = results_genes[sigpi,]
	sigde = which(abs(sigp[,"de"]) >= 2)
	sig_tn_de = sigp[sigde,]


	o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE) #Order the output by or p-value and then break ties using fold-change

	output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
	write.table(output, file="SigDE_R_ballgown_MDSC_WT_vs_SK2KO.txt", sep="\t", row.names=FALSE, quote=FALSE)

#View selected columns of the first 25 lines of output

	output[1:25,c(1,4,5)]
	
	
	
	
## Create a heatmap to vizualize expression differences between the 6 samples
	

#Define custom dist and hclust functions for use with heatmaps

	mydist=function(c) {dist(c,method="euclidian")}
	myclust=function(c) {hclust(c,method="average")}

	main_title="sig DE Transcripts for MDSC_WT vs. MDSC_Sk2KO"
	par(cex.main=0.8)
	sig_genes_de=sig_tn_de[,"id"]
	sig_gene_names_de=sig_tn_de[,"gene_name"]

	data=log2(as.matrix(gene_expression[as.vector(sig_genes_de),data_columns])+1)
	heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(10,4), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names_de,labCol=short_names,col=rev(heat.colors(75)))


	dev.off()
	
