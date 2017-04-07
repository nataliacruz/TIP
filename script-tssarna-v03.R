# rvencio, fela, pikachu, ema
library(quantmod)

minL = 5 # tamanho mínimo aceitavel de read
maxL = 151 # tamanho do read na tecnologia vigente
minCount = 5 # número mínimo de reads que o pico tem que ter para ir em frente 
maxPeaks = 3 # número de picos válidos e ordenados por tamanho que permitimos ir em frente


D = read.delim("minus_mapped_chr_R1.sam", skip=5, header=FALSE)
tss = read.delim("tss-1e-10_2_10.txt")
genes = read.delim("gene_annotation.txt")


# filtrar só TSSs que de fato são inícios de genes codificantes 
k = genes[,"start"] != genes[,"end"] 
genes = genes[k, ]

for(i in 



chromossomes = unique(D[,3]) # coluna 3 é onde está a informação de cromossomo
nc = length(chromossomes) # quantos cromossomos tem?
write.table("go...", "numeros-do-histograma.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
for(i in 1:nc){
	k = D[,3] == chromossomes[i]
	d = D[k,] # só o cromossomo da vez.
	k = tss[,1] == chromossomes[i]
	t = tss[k,] # só o cromossomo da vez.
	nt = nrow(t)
	write.table(chromossomes[i], "numeros-do-histograma.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
	for(j in 1:nt){
		k = d[ ,4] == t[j,2] # acha todos os reads que começam na posição
		if(sum(k) != 0){
			w = d[k, 6]
			x = as.numeric( gsub("M", "", w))
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
			teste = t(c(j, t[j,2], floor(h$mids[p])))
			write.table(teste, "numeros-do-histograma.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
		}
	}
}



