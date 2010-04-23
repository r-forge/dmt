### Tommi Suvitaival
### 30.9.2008

# Tallentaa klusteroinnin tulokset tiedostoon siten, että tallennetaan
# järjestyksessä kuhunkin klusteriin yleisimmän kuuluneen metaboliitin
# nimi tai järjestysnumero sekä sen esiintymisfrekvenssi tässä klusterissa.

# 'met_names' on vektori, joka sisältää metaboliittien nimet
#				  samassa järjestyksessä, kuin missä ne ilmenevät datassa.
# 'fname' on tiedoston nimi, johon funktio kirjoittaa.
#			 Voi sisältää myös polun.

write_clust = function(data,met_names,V_mode,V_posterior,fname) {

	K = ncol(V_mode)
	#tiedosto = paste(path,"klusterit2.txt",sep="")
	write("Metabolites in their most probable cluster.\nFirst column: posterior probability of metabolite m in cluster k, second column: number of the metabolite m, third column: name of the metabolite (if available).",file=fname)
	
	if (!is.null(met_names)) {
		for (k in 1:K) {
			write(paste("\nCluster ",k,"\n",sep=""),file=fname,append=T)
			for (m in which(V_mode[,k]==1)) {
				write(paste(round(V_posterior[m,k],digits=2),m,met_names[m],sep="\t\t"),
						file=fname,append=T)
			}	
		}
	} else {
		for (k in 1:K) {
			write(paste("\nCluster ",k,"\n",sep=""),file=fname,append=T)
			for (m in which(V_mode[,k]==1)) {
				write(paste(round(V_posterior[m,k],digits=2),m,sep="\t\t"),
						file=fname,append=T)
			}
		}
	}

}