DiffGOs <-
function(totNum,geneNum,goOnto,freq.supply,freq.ref,OnClass){
	dif.tbl <- matrix(0,nr=length(freq.supply),nc=9)
	dif.tbl <- as.data.frame.matrix(dif.tbl)
	ini <- 0
	for(i in names(freq.supply)){
		ini <- ini+1
		dif.data <- c(as.numeric(freq.supply[i]),as.numeric(freq.ref[i] - freq.supply[i]),
			geneNum-as.numeric(freq.supply[i]),
			totNum - geneNum-(as.numeric(freq.ref[i])-as.numeric(freq.supply[i])))
		dif.tbl[ini,c(3:8)] <- c(i,freq.supply[i],freq.ref[i],geneNum,totNum,goOnto[[i]][2])
		if(any(dif.data<=5)) test.res <- fisher.test(matrix(dif.data,nr=2))
		if(all(dif.data>5)) test.res <- chisq.test(matrix(dif.data,nr=2))
		dif.tbl[ini,1] <- test.res$p.v
	}
	dif.tbl[,9] <- OnClass
	return(dif.tbl)
}

