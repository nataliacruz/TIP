# rvencio, fela, pikachu, ema

D = read.delim("minus_mapped_chr_R1.sam", skip=5, header=FALSE)
tss = read.delim("tss-1e-10_2_10.txt")

minL = 5 # tamanho mínimo aceitavel de read
maxL = 151 # tamanho do read na tecnologia vigente


chromossomes = unique(D[,3]) # coluna 3 é onde está a informação de cromossomo
nc = length(chromossomes) # quantos cromossomos tem?
for(i in 1:nc){
	k = D[,3] == chromossomes[i]
	d = D[k,] # só o cromossomo da vez.
	k = tss[,1] == chromossomes[i]
	t = tss[k,] # só o cromossomo da vez.
	nt = nrow(t)
	write.table(t(c(0,seq(minL-5,maxL+5,by=5))), "numeros-do-histograma.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
	for(j in 1:10){   ###nt){
		k = d[ ,4] == t[j,2] # acha todos os reads que começam na posição
		if(sum(k) != 0){
			w = d[k, 6]
			x = as.numeric( gsub("M", "", w))
			h = hist(x, plot=FALSE, br=seq(minL-5,maxL+5,by=5)) # trocar isso por achar picos
			p = findPeaks(h$density, thresh=0)
			x11(); plot(h, col="pink"); abline(v=h$mids[p - 1], col="red")
			#teste = t(c(j, t[j,2], h$counts))
			#write.table(teste, "numeros-do-histograma.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
		}
	}
}



