# Copied 22.4.09

# Funktio 'show_clusters'
# Tehtävä:
# Piirtää kutakin klusteria vastaavan datan korrelaatiomatriisin.
# Punainen väri vastaa korkeaa positiivista korrelaatiota,
# sininen väri korkeaa negatiivista korrelaatiota ja
# vihreä väri nollakorrelaatiota.
# Parametrit:
# xxx - data (nrow='näytteiden lkm', ncol='näytteiden dimensio')
# VVV - klusterointi eli binäärinen matriisi, joka kertoo, 
# 		  mihin klusteriin xxx:n osoittama vastaava datapiste kuuluu.
#		  Datapiste voi kuulua ainoastaan yhteen klusteriin.

show_clusters <- function(xxx,VVV) {

	# Näytteiden lukumäärä
	m <- nrow(VVV)
	# Klusterien lukumäärä
	K <- ncol(VVV)
	corre <- array()
	klusteri <- array()

	# Etsitään kunkin rivin 'j' suurimman arvon indeksi.
	# Tätä ei kuitenkaan tarvita myöhemmin!?
	for (j in 1:m) {
		prob <- which.max(VVV[j,])
	}

	# Käydään läpi kaikki 'K' klusteria.
	for (k in 1:K) {
		# Etsitään datamatriisista klusteriin 'k' kuuluvat näytteet.
		klusteri <- xxx[VVV[,k]==1,]
		#print("klusteri ja corre")
		#print(dim(klusteri)) # Tulostaa klusteriin kuuluvien näytteiden lukumäärän.
		
		# Lasketaan klusterin näytteiden korrelaatio.
		if (length(klusteri)>0) {
			corre <- cor(t(klusteri))
		} else {
			#print("do this else")
			rm(corre)
			corre <- array(0,dim=c(2,2))
		}
		
		#corre[1,1] <- 0

		# Piirretään klusterin 'k' korrelaatiomatriisi.
		# Tätä funktiota kutsuva taho (funktio) voi määritellä, millaiseen kehykseen kuva tulee.
		screen(k)
		image(1:nrow(corre),1:nrow(corre),-corre, 
				zlim=c(-1,1), col=rainbow(100,start=0,end=4/6), xlab="",ylab="",main=sprintf('K=%d',k))
	}

}