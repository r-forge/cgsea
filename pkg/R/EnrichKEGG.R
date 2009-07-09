EnrichKEGG <-
function(supplyID,univerID=names(gene2map),p.adjust.methods=c("holm", 
	"hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))
{
	if(missing(p.adjust.methods)) p.adjust.methods="bonferroni"
	univerID <- intersect(univerID,names(gene2map))
	ids <- intersect(supplyID,univerID)
	totAnnNum <- length(univerID)
	annoDiffNum <- length(ids)
	if(! all(univerID==names(gene2map))){
		map2id <- lapply(map2gene,function(x) intersect(x,univerID))
		map2id <- map2id[sapply(map2id,length)>0]
	}else{
		map2id <- map2gene
	}
	mapNum <- length(map2id)
	diffKegg.tbl <-data.frame("MapID"=names(map2id),
		"MapTitle"=as.character(kegg.pthid2ti[names(map2id)]),
		"Pvalue"=rep(0,mapNum),"AdjPv"=rep(0,mapNum),
		"DiffExprGeneNumThisMap"=rep(0,mapNum),
		"TotalAnnotatedGeneNumThisMap"=rep(0,mapNum),
		"AnnoDiffExprGeneNum"=rep(0,mapNum),
		"AnnotatedGeneNum"=rep(0,mapNum),		
		"GeneIDs"=rep("a",mapNum),stringsAsFactors=FALSE)
	rownames(diffKegg.tbl) <- names(map2id)
	for(j in rownames(diffKegg.tbl)){
		knum <- length(intersect(ids,map2id[[j]]))
		if(knum <1) next
		dif.data <- c(knum,length(map2id[[j]]),annoDiffNum,totAnnNum)
		diffKegg.tbl[j,5:8] <- dif.data
		knum.rm <- length(map2id[[j]])-knum
		tot.rm <- totAnnNum-annoDiffNum-knum.rm
		dif.data <- c(knum,knum.rm,annoDiffNum-knum,tot.rm)
		if(any(dif.data<=5)) test.res <- fisher.test(matrix(dif.data,nr=2))
		if(all(dif.data>5)) test.res <- chisq.test(matrix(dif.data,nr=2))
		diffKegg.tbl[j,3] <- test.res$p.v	
		diffKegg.tbl[j,9] <- ifelse((dif.data[1]/dif.data[2])>(dif.data[3]/dif.data[4]),
						"Over","Under")
		diffKegg.tbl[j,10] <- paste(intersect(ids,map2id[[j]]),collapse=" ")
	}
	diffKegg.tbl <- diffKegg.tbl[diffKegg.tbl[,5]>0,]
	diffKegg.tbl[,4] <- p.adjust(diffKegg.tbl[,3],method=p.adjust.methods)
	diffKegg.tbl <- sort.data.frame(diffKegg.tbl,key="Pvalue");a <- diffKegg.tbl	
	return(diffKegg.tbl)
}

