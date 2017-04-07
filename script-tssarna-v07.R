# rvencio, fela, pikachu, ema
library(quantmod)
library(Biostrings)

minL = 5 # tamanho mínimo aceitavel de read
maxL = 151 # tamanho do read na tecnologia vigente
minCount = 10 # número mínimo de reads que o pico tem que ter para ir em frente 
maxPeaks = 3 # número de picos válidos e ordenados por tamanho que permitimos ir em frente
maxUpstream = 50 # distância a partir do início da CDS upstream máxima que pode conter um TSS experimental.
maxDownstream = 10 # distância a partir do início da CDS downstream máxima que pode conter um TSS experimental.
filenameOut = "numeros-do-histograma.txt"
filenameTSSpositions = "tss-1e-10_2_10.txt"
filenameGenepositions = "gene_annotation.txt"
filenameRNAseqdataSAM = "minus_mapped_chr_R1.sam"
filenameGenomeFASTA = "genome_Halobacterium_salinarum_NRC-1.fasta"

D = read.delim(filenameRNAseqdataSAM, comment.char="@", header=FALSE)
tss = read.delim(filenameTSSpositions)
genes = read.delim(filenameGenepositions)
genome = readDNAStringSet(filenameGenomeFASTA)


# filtrar só TSSs que de fato são inícios de genes codificantes 
k = genes[,"start"] != genes[,"end"] 
genes = genes[k, ]

chromossomes = unique(genes[,1]) # coluna 3 é onde está a informação de cromossomo
nc = length(chromossomes) # quantos cromossomos tem?
write.table("go...", filenameOut, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
for(i in 1:nc){
	# i vai indexar o cromossomo da vez.

	# atividade: filtrando TSS válidos (perto do início de um CDS)
	if(nrow(genes) > 0){
		staytss = NULL
		k = genes[,"sequence"] == as.vector(chromossomes[i])
		d = genes[k,]
		q = as.vector(tss[,"chrom"]) == as.vector(chromossomes[i])
		f = tss[q,]
		if(nrow(f) == 0){ next } # isso acontece se não tiver nenhum tss no cromossomo da vez.
		ng = nrow(d)
		for(j in 1:ng){
			if(genes[j, "strand"] == "+"){
				k = d[j,"start"] - maxUpstream  <  f[, "start"] & f[, "start"]  <  d[j,"start"] + maxDownstream
				staytss = c(staytss, which(k))
			}
			if(genes[j, "strand"] == "-"){
				k = d[j,"start"] + maxUpstream  <  f[, "start"] & f[, "start"]  <  d[j,"start"] - maxDownstream
				staytss = c(staytss, which(k))
			}
		}
	}

	# atividade: achar os picos
	k = D[,3] == as.vector(chromossomes[i])
	d = D[k,] # só o cromossomo da vez.
	f = f[staytss, ]
	k = f[,"chrom"] == as.vector(chromossomes[i])
	t = f[k,] # só o cromossomo da vez.
	nt = nrow(t)
	write.table(chromossomes[i], filenameOut, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
	for(j in 1:nt){
		k = d[ ,4] == t[j,2] # acha todos os reads que começam na posição
		if(sum(k) != 0){
			w = d[k, 6]
			x = as.numeric( gsub("M", "", w)) ### aqui assumimos que ia ser só no formato numernumeroM mas as vezes tem codigos nada a ver.
			h = hist(x, plot=FALSE, br=seq(minL-5,maxL+5,by=5)) # trocar isso por achar picos
			p = findPeaks(h$density, thresh=0)
			p = p - 1 # a função findPeaks aponta o índice seguinte ao pico que queremos.
			if(length(p) == 0){ next }	
			p = p[h$counts[p] > minCount]
			if(length(p) == 0){ next }	
			p = p[order(h$counts[p], decreasing=TRUE)]
			p = p[1:maxPeaks]
			p = p[!is.na(p)]
#			x11(); plot(h, col="pink"); abline(h=minCount, col="red"); abline(v=h$mids[p], col="blue")
			k = t[j,5] == genes[,"strand"]
			g = genes[k,]
			k = which.min(abs(t[j,2] - g[, "start"]))
			genename = as.vector(g[k, "name"])
			seqs = NULL
			k = which(names(genome) == chromossomes[i])
			for(s in 1:length(p)){ seqs = c(seqs, as.character(DNAStringSet(genome[k], start=t[j,2], end=t[j,2]+floor(h$mids[p])[s]))) }
			results = rbind( t[j,2] + floor(h$mids[p]) , seqs ,  h$counts[p] )
			results = as.vector(results)
			output = t(c(j, genename, t[j,2], results))
#			output = t(c(j, t[j,2], t[j,2] + floor(h$mids[p])))
			write.table(output, filenameOut, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
		}
	}
}



