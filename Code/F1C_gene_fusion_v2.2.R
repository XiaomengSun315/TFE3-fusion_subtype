rm(list=ls())

args <- commandArgs(T)
setwd(paste0(args[1], "Data/"))
result_dir <- paste0(args[1], "Results/")

#### import packages ####
# BiocManager::install("gggenes")
# BiocManager::install("ggplotify")
# BiocManager::install("eoffice")
# devtools::install_github('YuLab-SMU/ggtree')
# devtools::install_github('GuangchuangYu/yyplot')
# install.packages('randomcoloR')
library("openxlsx")
library("ggplot2")
library("gggenes")
library("ggplotify")
library("eoffice")
library("randomcoloR")

#### read data ####
gene_fusion <- read.xlsx("gene_fusion_1C.xlsx", colNames = TRUE, na.strings = c("","NA","Nan","#N/A"), rowNames = FALSE)
gene_exon <- read.xlsx("gene_exon.xlsx", colNames = TRUE, na.strings = c("","NA","Nan","#N/A"), rowNames = FALSE)
gene_domain <- read.xlsx("gene_domain_1C.xlsx", colNames = TRUE, na.strings = c("","NA","Nan","#N/A"), rowNames = FALSE)

# # set domains' color
# colors <- data.frame(matrix(nrow=length(unique(gene_domain$domain)), ncol=2))
# names(colors) <- c("domain", "color")
# colors$domain <- unique(gene_domain$domain)
# colors$color <- randomColor(count = length(unique(gene_domain$domain)))
# gene_domain$color <- colors[match(gene_domain$domain, colors$domain),"color"]

# gene_domain[gene_domain$domain=="AD","color"] <- "red"
# gene_domain[gene_domain$domain=="bHLH-LZ","color"] <- c("green", "orange", "purple")
# # gene_domain[gene_domain$domain=="AD","color"] <- "#E81123"
# # gene_domain[gene_domain$domain=="bHLH-LZ","color"] <- c("#6AA329", "#EF7B28", "#7C5A95")

#### output format ####
genes <- data.frame(matrix(nrow=0, ncol=8))
names(genes) <- c("id", "molecule", "gene", "start", "end", "strand", "orientation", "color")
subgenes <- data.frame(matrix(nrow=0, ncol=11))
names(subgenes) <- c("id", "molecule","gene","start","end","strand","subgene","from","to","orientation", "color")
current_line <- 1

###############################################################################
############################# data transformation

# get position of genes & exons (GRCh37.p13)
for(i in 1:nrow(gene_fusion)){
	# basic info
	fusion <- gene_fusion[i,"fusion"]
	gene1 <- gsub("-.*", "", fusion)
	gene1_color <- unlist(strsplit(gene_fusion[i, "first_gene_color"], split=","))
	gene1_exon_start <- 1
	gene1_exon_end <- as.numeric(gsub("^.*-", "", gene_fusion[i,"first_gene_exon"]))
	gene1_exon_count <- gene1_exon_end - gene1_exon_start + 1
	gene2 <- "TFE3"
	gene2_color <- unlist(strsplit(gene_fusion[i, "second_gene_color"], split=","))
	gene2_exon_start <- as.numeric(gsub("-.*$", "", gene_fusion[i,"second_gene_exon"]))
	gene2_exon_end <- 10
	gene2_exon_count <- gene2_exon_end - gene2_exon_start + 1
	
	# plot_subgene basic info
	subgenes[current_line:(current_line+gene1_exon_count+gene2_exon_count),"id"] <- paste0("fusion", i)
	subgenes[current_line:(current_line+gene1_exon_count+gene2_exon_count),"molecule"] <- fusion
	
	############### gene 1 ###############
	# basic info
	subgenes[current_line:(current_line+gene1_exon_count), "gene"] <- gene1
	subgenes[current_line:(current_line+gene1_exon_count), "subgene"] <- c("5'UTR", seq(gene1_exon_start, gene1_exon_end, 1))
	
	# gene data and transcripts
	gene_data <- gene_exon[gene_exon$Gene.name==gene1,]
	transcripts <- unique(gene_data$Transcript.stable.ID)
	if(length(transcripts) > 1){
		if(gene1=="MED15"){
			gene_data <- gene_data[gene_data$Transcript.stable.ID=="ENST00000263205",]
		}else if(gene1=="NONO"){
			gene_data <- gene_data[gene_data$Transcript.stable.ID=="ENST00000373856",]
		}else{
			gene_data <- gene_data[gene_data$Transcript.stable.ID==transcripts[1],]
		}
	}
	
	# fill in exon position info
	# exons
	temp <- merge(subgenes[(current_line+1):(current_line+gene1_exon_count),c("gene", "subgene")], gene_data[,c("Strand", "Exon.rank.in.transcript", "Exon.region.start.(bp)", "Exon.region.end.(bp)", "Genomic.coding.start", "Genomic.coding.end")], by.x="subgene", by.y="Exon.rank.in.transcript", sort = FALSE)
	subgenes[(current_line+1):(current_line+gene1_exon_count), "from"] <- temp$Genomic.coding.start
	subgenes[(current_line+1):(current_line+gene1_exon_count), "to"] <- temp$Genomic.coding.end
	subgenes[(current_line+1):(current_line+gene1_exon_count), "color"] <- rep(gene1_color, gene1_exon_count)[1:gene1_exon_count]
	# 5'UTR
	temp_5UTR <- na.omit(gene_data[,c("5'.UTR.start", "5'.UTR.end")])
	temp_5UTR$length <- temp_5UTR[,2] - temp_5UTR[,1]
	length_5UTR <- sum(temp_5UTR$length)

	temp_5UTR_whole_exon <- temp[is.na(temp$Genomic.coding.start),]
	if(nrow(temp_5UTR_whole_exon) > 0){
		# whole exon in 5'UTR
		hit_line <- intersect(grep("TRUE", is.na(subgenes$from)), current_line:(current_line+gene1_exon_count-1))
		subgenes[hit_line[1:nrow(temp_5UTR_whole_exon)], "from"] <- temp_5UTR_whole_exon[,"Exon.region.start.(bp)"]
		subgenes[hit_line[1:nrow(temp_5UTR_whole_exon)], "to"] <- temp_5UTR_whole_exon[,"Exon.region.end.(bp)"]
		subgenes[hit_line[1:nrow(temp_5UTR_whole_exon)], "subgene"] <- paste0(temp_5UTR_whole_exon$subgene, "-5'UTR")

		# partial exon in 5'UTR
		temp_5UTR_partial_exon <- temp[temp$subgene==max(as.numeric(temp_5UTR_whole_exon$subgene)) + 1,]
		if(unique(temp_5UTR_partial_exon$Strand)>0){
			subgenes[hit_line[length(hit_line)], "from"] <- temp_5UTR_partial_exon[,"Exon.region.start.(bp)"]
			subgenes[hit_line[length(hit_line)], "to"] <- temp_5UTR_partial_exon[,"Genomic.coding.start"] - 1
			subgenes[hit_line[length(hit_line)], "subgene"] <- paste0(temp_5UTR_partial_exon$subgene, "-5'UTR")
		}else{
			# ZC3H4: exons in UTR & negative strand
			subgenes[hit_line[length(hit_line)], "from"] <- temp_5UTR_partial_exon[,"Genomic.coding.end"] + 1
			subgenes[hit_line[length(hit_line)], "to"] <- temp_5UTR_partial_exon[,"Exon.region.end.(bp)"]			
			subgenes[hit_line[length(hit_line)], "subgene"] <- paste0(temp_5UTR_partial_exon$subgene, "-5'UTR")
		}
	}else{
		hit_line <- intersect(grep("5'UTR", subgenes$subgene), current_line:(current_line+gene1_exon_count))
		subgenes[hit_line, "to"] <- max(as.numeric(na.omit(gene_data[,"5'.UTR.end"])))
		subgenes[hit_line, "from"] <- subgenes[hit_line,"to"] - length_5UTR
	}
	subgenes[hit_line, "color"] <- "#d9d9d9"

	# fill in gene position info
	subgenes[current_line:(current_line+gene1_exon_count),"strand"] <- unique(gene_data$Strand)
	subgenes[current_line:(current_line+gene1_exon_count),"orientation"] <- unique(gene_data$Strand)
	subgenes[current_line:(current_line+gene1_exon_count),"start"] <- min(na.omit(subgenes[current_line:(current_line+gene1_exon_count),"from"]))
	subgenes[current_line:(current_line+gene1_exon_count),"end"] <- max(na.omit(subgenes[current_line:(current_line+gene1_exon_count),"to"]))
	
	############### gene 2: TFE3 ###############
	# basic info
	subgenes[(current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count), "gene"] <- gene2
	subgenes[(current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count), "subgene"] <- c(seq(gene2_exon_start, gene2_exon_end, 1))
	
	# gene data
	gene_data <- gene_exon[gene_exon$Gene.name==gene2,]
	transcripts <- unique(gene_data$Transcript.stable.ID)
	if(length(transcripts) > 1){
		gene_data <- gene_data[gene_data$Transcript.stable.ID==transcripts[1],]
	}
	
	# fill in exon position info
	# hit_line <- intersect(grep("3'UTR", subgenes$subgene), (current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count+1))
	# subgenes[hit_line, "from"] <- as.numeric(na.omit(gene_data[,"3'.UTR.start"]))
	# subgenes[hit_line, "to"] <- as.numeric(na.omit(gene_data[,"3'.UTR.end"]))
	# subgenes[hit_line, "color"] <- "#d9d9d9"
	temp <- merge(subgenes[(current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count),c("gene", "subgene")], gene_data[,c("Exon.rank.in.transcript", "Genomic.coding.start", "Genomic.coding.end")], by.x="subgene", by.y="Exon.rank.in.transcript", sort = FALSE, all.y = FALSE)
	subgenes[(current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count), "from"] <- temp$Genomic.coding.start
	subgenes[(current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count), "to"] <- temp$Genomic.coding.end
	subgenes[(current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count), "color"] <- rep(gene2_color, gene2_exon_count)[(gene2_exon_count+1):(gene2_exon_count*2)]
	
	# fill in gene position info
	subgenes[(current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count),"strand"] <- unique(gene_data$Strand)
	subgenes[(current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count),"orientation"] <- unique(gene_data$Strand)
	subgenes[(current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count),"start"] <- min(subgenes[(current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count),"from"])
	subgenes[(current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count),"end"] <- max(subgenes[(current_line+gene1_exon_count+1):(current_line+gene1_exon_count+gene2_exon_count),"to"])

	# update current_line
	current_line <- nrow(subgenes) + 1
	
	############################################################
	############### gene domains
	gene1_domain_exons <- intersect(subgenes[subgenes$id==paste0("fusion", i) & subgenes$gene==gene1, "subgene"], unique(unlist(strsplit(gene_domain[gene_domain$gene==gene1,"exons"], split=","))))
	gene2_domain_exons <- intersect(subgenes[subgenes$id==paste0("fusion", i) & subgenes$gene==gene2, "subgene"], unique(unlist(strsplit(gene_domain[gene_domain$gene==gene2,"exons"], split=","))))
	
	if(length(gene1_domain_exons) > 0){
		# domains included in the fusion gene
		current_line_temp <- current_line
		temp_domain <- gene_domain[gene_domain$gene==gene1,]
		j = 1
		for(j in 1:nrow(temp_domain)){
			overlap_exon <- unlist(strsplit(temp_domain[j,"exons"], split=","))
			if(length(grep("FALSE", overlap_exon %in% gene1_domain_exons))==0){
				subgenes[current_line_temp, "subgene"] <- temp_domain[j,"domain"]
				subgenes[current_line_temp, "from"] <- temp_domain[j,"start_gDNA"]
				subgenes[current_line_temp, "to"] <- temp_domain[j,"end_gDNA"]
				subgenes[current_line_temp, "color"] <- temp_domain[j,"color"]
				current_line_temp <- current_line_temp + 1
			}
		}

		# gene domains basic info
		subgenes[current_line:nrow(subgenes),"id"] <- paste0("fusion", i, "_domain")
		subgenes[current_line:nrow(subgenes),"molecule"] <- paste0(fusion, "_domain")
		subgenes[current_line:nrow(subgenes),"gene"] <- paste0(gene1, "_domain")
		subgenes[current_line:nrow(subgenes), "start"] <- unique(na.omit(subgenes[subgenes$id==paste0("fusion",i)&subgenes$gene==gene1, "start"]))
		subgenes[current_line:nrow(subgenes), "end"] <- unique(na.omit(subgenes[subgenes$id==paste0("fusion",i)&subgenes$gene==gene1, "end"]))
		subgenes[current_line:nrow(subgenes), "strand"] <- unique(na.omit(subgenes[subgenes$id==paste0("fusion",i)&subgenes$gene==gene1, "strand"]))
		subgenes[current_line:nrow(subgenes), "orientation"] <- unique(na.omit(subgenes[subgenes$molecule==fusion&subgenes$gene==gene1, "orientation"]))
	}else{
		# ASPSCR1 (without domain)
		subgenes[current_line, "id"] <- paste0("fusion", i, "_domain")
		subgenes[current_line, "molecule"] <- "ASPSCR1-TFE3_domain"
		subgenes[current_line, "gene"] <- "ASPSCR1_domain"
		subgenes[current_line, "subgene"] <- " "
		subgenes[current_line, "start"] <- unique(na.omit(subgenes[subgenes$id==paste0("fusion",i)&subgenes$gene=="ASPSCR1", "start"]))
		subgenes[current_line, "end"] <- unique(na.omit(subgenes[subgenes$id==paste0("fusion",i)&subgenes$gene=="ASPSCR1", "end"]))
		subgenes[current_line, "from"] <- subgenes[current_line, "start"]
		subgenes[current_line, "to"] <- subgenes[current_line, "end"]
		subgenes[current_line, "strand"] <- unique(na.omit(subgenes[subgenes$id==paste0("fusion",i)&subgenes$gene=="ASPSCR1", "strand"]))
		subgenes[current_line, "orientation"] <- unique(na.omit(subgenes[subgenes$id==paste0("fusion",i)&subgenes$gene=="ASPSCR1", "orientation"]))
		subgenes[current_line, "color"] <- "white"
	}
	# update current_line
	current_line <- nrow(subgenes) + 1

	if(length(gene2_domain_exons) > 0){
		# domains included in the fusion gene
		current_line_temp <- current_line
		temp_domain <- gene_domain[gene_domain$gene==gene2,]
		j = 1
		for(j in 1:nrow(temp_domain)){
			overlap_exon <- unlist(strsplit(temp_domain[j,"exons"], split=","))
			if(length(grep("FALSE", overlap_exon %in% gene2_domain_exons))==0){
				subgenes[current_line_temp, "subgene"] <- temp_domain[j,"domain"]
				subgenes[current_line_temp, "from"] <- temp_domain[j,"start_gDNA"]
				subgenes[current_line_temp, "to"] <- temp_domain[j,"end_gDNA"]
				subgenes[current_line_temp, "color"] <- temp_domain[j,"color"]
				current_line_temp <- current_line_temp + 1
			}
		}

		# gene domains basic info
		subgenes[current_line:nrow(subgenes),"id"] <- paste0("fusion", i, "_domain")
		subgenes[current_line:nrow(subgenes),"molecule"] <- paste0(fusion, "_domain")
		subgenes[current_line:nrow(subgenes),"gene"] <- paste0(gene2, "_domain")
		subgenes[current_line:nrow(subgenes), "start"] <- unique(na.omit(subgenes[subgenes$id==paste0("fusion",i)&subgenes$gene==gene2, "start"]))
		subgenes[current_line:nrow(subgenes), "end"] <- unique(na.omit(subgenes[subgenes$id==paste0("fusion",i)&subgenes$gene==gene2, "end"]))
		subgenes[current_line:nrow(subgenes), "strand"] <- unique(na.omit(subgenes[subgenes$id==paste0("fusion",i)&subgenes$gene==gene2, "strand"]))
		subgenes[current_line:nrow(subgenes), "orientation"] <- unique(na.omit(subgenes[subgenes$molecule==fusion&subgenes$gene==gene2, "orientation"]))
	}
	# update current_line
	current_line <- nrow(subgenes) + 1
}

############# convert positions into plot friendly version
#### output format ####
trans_dic <- data.frame(matrix(nrow=0, ncol=4))
names(trans_dic) <- c("fusion", "gene", "before", "after")
current_line <- nrow(trans_dic) + 1

new_subgenes <- subgenes
new_subgenes$gene_length <- new_subgenes$end - new_subgenes$start + 1
new_subgenes$exon_length <- new_subgenes$to - new_subgenes$from + 1
new_subgenes$from_trans <- new_subgenes$from
new_subgenes$to_trans <- new_subgenes$to

fusion_count <- unique(new_subgenes$id)[!grepl("domain", unique(new_subgenes$id))]
for(fusion_number in fusion_count){
	count <- as.numeric(gsub("fusion", "", fusion_number))
	fusion <- unique(new_subgenes[new_subgenes$id==fusion_number, "molecule"])
	line_range_fusion <- grep(paste0("^", fusion_number, "$"), new_subgenes$id, perl=TRUE)
	line_range_domain <- grep(paste0("^", fusion_number, "_domain$"), new_subgenes$id, perl=TRUE)

	############ fusion: transform exon from/to position to continious coordinates
	for(line in line_range_fusion){
		# transform to continous coordinates
		if(line != line_range_fusion[1]){
			new_subgenes[line, "from_trans"] <- max(na.omit(c(new_subgenes[line-1, "to_trans"], new_subgenes[line_range_fusion[1], "to_trans"]))) + 1
		}
		new_subgenes[line, "to_trans"] <- new_subgenes[line, "from_trans"] + new_subgenes[line, "exon_length"] - 1
		
		# set corresponding position dictionary
		exon_length <- new_subgenes[line, "exon_length"]
		trans_dic[current_line:(current_line+exon_length-1), "fusion"] <- new_subgenes[line, "id"]
		trans_dic[current_line:(current_line+exon_length-1), "gene"] <- new_subgenes[line, "gene"]
		if(new_subgenes[line, "strand"] == 1){
			# positive strand
			trans_dic[current_line:(current_line+exon_length-1), "before"] <- seq(new_subgenes[line, "from"], new_subgenes[line, "to"], 1)
		}else{
			# negative srand
			trans_dic[current_line:(current_line+exon_length-1), "before"] <- seq(new_subgenes[line, "to"], new_subgenes[line, "from"], -1)
		}
		trans_dic[current_line:(current_line+exon_length-1), "after"] <- seq(new_subgenes[line, "from_trans"], new_subgenes[line, "to_trans"], 1)
		current_line <- nrow(trans_dic) + 1
	}

	# gene info
	gene_names <- unique(new_subgenes[line_range_fusion, "gene"])
	
	genes[(4*count-3):(4*count), "id"] <- c(rep(fusion_number,2), rep(paste0(fusion_number, "_domain"), 2))
	genes[(4*count-3):(4*count), "molecule"] <- c(rep(fusion,2), rep(paste0(fusion,"_domain"),2))
	genes[(4*count-3):(4*count), "gene"] <- c(gene_names[1], gene_names[2], paste0(gene_names[1],"_domain"), paste0(gene_names[2],"_domain"))

	############ domain: transform domain position to continous coordiantes
	############ fusion: sum gene length
	for(gene in gene_names){
		line_range_gene <- intersect(line_range_fusion, grep(gene, new_subgenes$gene))
		line_range_gene_domain <- intersect(line_range_domain, grep(gene, new_subgenes$gene))
		# gene start & end
		gene_cds_start_trans <- min(na.omit(new_subgenes[line_range_gene, "from_trans"]))
		gene_cds_end_trans <- max(na.omit(new_subgenes[line_range_gene, "to_trans"]))
		new_subgenes[line_range_gene, "start_trans"] <- gene_cds_start_trans
		new_subgenes[line_range_gene, "end_trans"] <- gene_cds_end_trans
		
		# domain from & to & start & end
		trans_dic_specific <- trans_dic[trans_dic$fusion==fusion_number & trans_dic$gene==gene,]
		from_trans_matrix <- trans_dic_specific[match(new_subgenes[line_range_gene_domain, "from"], trans_dic_specific$before),]
		to_trans_matrix <- trans_dic_specific[match(new_subgenes[line_range_gene_domain, "to"], trans_dic_specific$before),]
		
		new_subgenes[line_range_gene_domain, "from_trans"] <- from_trans_matrix[,"after"]
		new_subgenes[line_range_gene_domain, "to_trans"] <- to_trans_matrix[,"after"]
		new_subgenes[line_range_gene_domain, "start_trans"] <- gene_cds_start_trans
		new_subgenes[line_range_gene_domain, "end_trans"] <- gene_cds_end_trans

		# for(line in 1:length(line_range_gene)){
		# 	if(line%%2 == 0){
		# 		new_subgenes[line_range_gene[line], "keep"] <- "yes"
		# 	}else{
		# 		new_subgenes[line_range_gene[line], "keep"] <- "no"
		# 	}
		# }

		# gene info
		if(gene != "TFE3"){
			# gene
			genes[4*count-3, "start_trans"] <- gene_cds_start_trans
			genes[4*count-3, "end_trans"] <- gene_cds_end_trans
			genes[4*count-3, "strand"] <- unique(new_subgenes[line_range_gene,"strand"])
			genes[4*count-3, "orientation"] <- unique(new_subgenes[line_range_gene,"orientation"])
			genes[4*count-3, "color"] <- paste0(gene1_color[1], ",", gene1_color[2])
			# domain
			genes[4*count-1, "start_trans"] <- gene_cds_start_trans
			genes[4*count-1, "end_trans"] <- gene_cds_end_trans
			genes[4*count-1, "strand"] <- unique(new_subgenes[line_range_gene,"strand"])
			genes[4*count-1, "orientation"] <- unique(new_subgenes[line_range_gene,"orientation"])
			genes[4*count-1, "color"] <- paste0(gene1_color[1], ",", gene1_color[2])
		}else{
			# gene
			genes[4*count-2, "start_trans"] <- gene_cds_start_trans
			genes[4*count-2, "end_trans"] <- gene_cds_end_trans
			genes[4*count-2, "strand"] <- unique(new_subgenes[line_range_gene,"strand"])
			genes[4*count-2, "orientation"] <- unique(new_subgenes[line_range_gene,"orientation"])
			genes[4*count-2, "color"] <- paste0(gene2_color[1], ",", gene2_color[2])
			# domain
			genes[4*count, "start_trans"] <- gene_cds_start_trans
			genes[4*count, "end_trans"] <- gene_cds_end_trans
			genes[4*count, "strand"] <- unique(new_subgenes[line_range_gene,"strand"])
			genes[4*count, "orientation"] <- unique(new_subgenes[line_range_gene,"orientation"])
			genes[4*count, "color"] <- paste0(gene2_color[1], ",", gene2_color[2])
		}
	}
}

write.xlsx(trans_dic, paste0(result_dir, "trans_dic.xlsx"))
write.xlsx(genes, paste0(result_dir, "genes_1C.xlsx"))
write.xlsx(new_subgenes, paste0(result_dir, "subgenes_1C.xlsx"))


# subgenes_plot <- new_subgenes[new_subgenes$keep=="yes",]

# Figure 1B
# ggplot(subset(gene_fusion, molecule == "U2AF2-TFE3" & gene == "TFE3"),
#        aes(xmin = start, xmax = end, y = strand)
#   ) +
#   geom_gene_arrow(arrowhead_width = unit(0, "mm")) +
#   geom_gene_label(aes(label = gene)) +
#   geom_subgene_arrow(arrowhead_width = unit(0, "mm"), 
#     data = subset(new_subgenes, molecule == "U2AF2-TFE3" & gene == "TFE3"),
#     aes(xsubmin = from, xsubmax = to, fill = subgene)
#   ) +
#   geom_subgene_label(
#     data = subset(new_subgenes, molecule == "U2AF2-TFE3" & gene == "TFE3"),
#     aes(xsubmin = from, xsubmax = to, label = subgene),
#     min.size = 0
#   )

###############################################################################
############################# Figure 1C (sperate plots)
fusions <- unique(genes$molecule)
dir.create(paste0(result_dir, "sperate/"))
for(fusion in fusions){
	plot_gene <- genes[genes$molecule==fusion,]
	plot_subgene <- new_subgenes[new_subgenes$molecule==fusion,]

	# gene color
	gene1_color <- unique(plot_subgene[plot_subgene$gene!="TFE3"&plot_subgene$subgene!="5'UTR","color"])
	gene2_color <- unique(plot_subgene[plot_subgene$gene=="TFE3"&plot_subgene$subgene!="3'UTR","color"])

	# gggene
	plot_gene$id <- factor(plot_gene$id, level=unique(plot_gene$id))
	plot_subgene$id <- factor(plot_subgene$id, level=unique(plot_subgene$id))
	dummies <- make_alignment_dummies(plot_gene, aes(xmin = start_trans, xmax = end_trans, y = id, id=gene), on = "TFE3", side = "right")
	names(dummies)[1] <- "id"

	gggene <- ggplot(plot_gene, aes(xmin = start_trans, xmax = end_trans, y = id)) +
		geom_blank(data = dummies) +
		facet_wrap(~ id, scales = "free", ncol = 1) +
		geom_gene_arrow(arrowhead_width = grid::unit(0, "mm"), arrow_body_height = grid::unit(8, "mm"), color="grey") +
		geom_subgene_arrow(data = plot_subgene, arrowhead_width = grid::unit(0.00001, "mm"), arrowhead_height = grid::unit(0.00001, "mm"), arrow_body_height = grid::unit(8, "mm"), aes(xmin = start_trans, xmax = end_trans, y = id, fill = color, xsubmin = from_trans, xsubmax = to_trans), color="grey", alpha=.7) +
		geom_subgene_label(data = plot_subgene, aes(xsubmin = from_trans, xsubmax = to_trans, label = subgene), min.size=4) +
		scale_fill_identity() +
		theme_genes() +
		theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(), axis.text.y=element_blank(), title=element_blank())

	# PDF version
	pdf(paste0(result_dir, "sperate/", fusion, ".pdf"))
	plot(gggene)
	dev.off()

	# TIFF version
	tiff(paste0(result_dir, "sperate/", fusion, ".tif"))
	plot(gggene)
	dev.off()

	# PDF version
	pptx <- paste0(result_dir, "sperate/", fusion, ".pptx")
	topptx(gggene, pptx)
}

############################# Fig 1C (all in one)
# fusions <- unique(genes$molecule)

# # shorten 5'UTR region length of gene NONO, RBM10, ZC3H4, U2AF2(deleted)
# for(number in c(3, 6, 9)){
# 	fusion <- fusions[number]
# 	gene1 <- gsub("-.*", "", fusion)
# 	modify_nono <- unique(plot_subgene[plot_subgene$gene==gene1&plot_subgene$subgene=="5'UTR","to"]) - 500
# 	plot_subgene[plot_subgene$gene==gene1&plot_subgene$subgene=="5'UTR","from_trans"] <- modify_nono
# 	plot_subgene[plot_subgene$gene==gene1,"start_trans"] <- modify_nono
# 	plot_gene[plot_gene$gene==gene1, "start_trans"] <- modify_nono
# }

fusion_plot <- function(plot_gene, plot_subgene, plot_type, aligan_target){
	# gggene plot
	plot_gene$id <- factor(plot_gene$id, level=unique(plot_gene$id))
	plot_subgene$id <- factor(plot_subgene$id, level=unique(plot_subgene$id))
	dummies <- make_alignment_dummies(plot_gene, aes(xmin = start_trans, xmax = end_trans, y = id, id=gene, forward = orientation), on = aligan_target, side = "right")
	names(dummies)[1] <- "id"
	gggene <- ggplot(plot_gene, aes(xmin = start_trans, xmax = end_trans, y = id)) +
		facet_wrap(~ id, scales = "free", ncol = 1) +
		geom_gene_arrow(arrowhead_width=unit(0,"mm"), arrow_body_height=unit(10,"mm"), color="#bdbdbd") +
		geom_subgene_arrow(data = plot_subgene, color="#bdbdbd", arrowhead_width=unit(0,"mm"), arrow_body_height=unit(10,"mm"), aes(xmin = start_trans, xmax = end_trans, y = id, fill = color, xsubmin = from_trans, xsubmax = to_trans), alpha=.7) +
		geom_subgene_label(data = plot_subgene, aes(xsubmin = from_trans, xsubmax = to_trans, label = subgene), min.size=4) +
		geom_blank(data = dummies) +
		scale_fill_identity() +
		theme_genes() +
		theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(), axis.text.y=element_blank(), title=element_blank())
	
	# PDF version
	pdf(paste0(result_dir, plot_type, "_ALL_pre.pdf"))
	plot(gggene)
	dev.off()
	
	# TIFF versio
	tiff(paste0(result_dir, plot_type, "_ALL_pre.tif"))
	plot(gggene)
	dev.off()
	
	# PPT version
	pptx <- paste0(result_dir, plot_type, "_ALL_pre.pptx")
	topptx(gggene, pptx)
}

# plot 
plot_gene <- genes[!grepl("domain", genes$id),]
plot_subgene <- new_subgenes[!grepl("domain", new_subgenes$id),]
fusion_plot(plot_gene, plot_subgene, "fusion", "TFE3")

plot_gene <- genes[grepl("domain", genes$id),]
plot_subgene <- new_subgenes[grepl("domain", new_subgenes$id),]
fusion_plot(plot_gene, plot_subgene, "domain", "TFE3_domain")
