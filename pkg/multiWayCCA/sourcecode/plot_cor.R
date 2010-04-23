# Copied 22.4.09

plot_cor = function(x,case,gender,VV,Vs,Vsum,path) {

	K = ncol(Vs) # Number of clusters -16.3.09
	#VV=t(as.matrix(VV))

	x_h = x[,case==0] # kaikki terveet
	x_s = x[,case==1] # kaikki sairaat
	x_1 = x[,gender==0] # sukupuoli 1
	x_2 = x[,gender==1] # sukupuoli 2
	x_h1 = x[,case==0 & gender==0] # terveet, sukupuoli 1
	x_h2 = x[,case==0 & gender==1] # terveet, sukupuoli 2
	x_s1 = x[,case==1 & gender==0] # sairaat, sukupuoli 1
	x_s2 = x[,case==1 & gender==1] # sairaat, sukupuoli 2

	mod = 5
	#print(dim(VV))
	corres_length = floor(dim(VV)[3]/mod)

	# Time series printing is performed only if there are enough time points. If 'corres_length=0', this function does not print any figure at all. -Tommi 30.1.09
	if (corres_length > 0) {
		#print(corres_length)
		# Matriisit eri ryhmien sisäisten korrelaatioiden tallentamiseen.
		corres_p = matrix(0,nrow=corres_length,ncol=K)
		corres_h = matrix(0,nrow=corres_length,ncol=K)
		corres_s = matrix(0,nrow=corres_length,ncol=K)
		corres_1 = matrix(0,nrow=corres_length,ncol=K)
		corres_2 = matrix(0,nrow=corres_length,ncol=K)
		corres_h1 = matrix(0,nrow=corres_length,ncol=K)
		corres_h2 = matrix(0,nrow=corres_length,ncol=K)
		corres_s1 = matrix(0,nrow=corres_length,ncol=K)
		corres_s2 = matrix(0,nrow=corres_length,ncol=K)
	
		Vsums = matrix(nrow=corres_length,ncol=K)
	
		#print(dim(VV[,,50]))
		for (i in 1:corres_length) {
	
			j = i*mod
			#print(dim(as.matrix(corres_p[i,])))
			#print(dim(VV[,,j]))
			#print(dim(corre_means(x,as.matrix(VV[,,j]))))
			corres_p[i,] = corre_means(x,VV[,,j])
			corres_h[i,] = corre_means(x_h,VV[,,j])
			if (sum(case)>0) {
				corres_s[i,] = corre_means(x_s,VV[,,j])
				corres_s1[i,] = corre_means(x_s1,VV[,,j])
			}
			corres_1[i,] = corre_means(x_1,VV[,,j])
			if (sum(gender)>0) {
				corres_2[i,] = corre_means(x_2,VV[,,j])
				corres_h2[i,] = corre_means(x_h2,VV[,,j])
			}
			corres_h1[i,] = corre_means(x_h1,VV[,,j])
			if (sum(case&gender)>0) {
				corres_s2[i,] = corre_means(x_s2,VV[,,j])
			}
	
			#Vsums[i,] =  colSums(samples$VV[,,j])
			if (K==1) { # If only one cluster, only sum needs to be calculated.
				Vsums[i,] =  sum(VV[,,j])
			} else { # If more than one cluster, colSums need to be calculated.
				Vsums[i,] =  colSums(VV[,,j]) # 'samples' not passed to the function -16.3.09
			}
		}
	
	# 	print(colMeans(corres_p))
	# 	print(colMeans(corres_h))
	# 	print(colMeans(corres_s))
	# 	print(colMeans(corres_1))
	# 	print(colMeans(corres_2))
	# 	print(colMeans(corres_h1))
	# 	print(colMeans(corres_h2))
	# 	print(colMeans(corres_s1))
	# 	print(colMeans(corres_s2))
	
# 		cor_labels = c("pooled","healthy","sick","gender1","gender2",
# 							"healthy_gender1","healthy_gender2",
# 							"sick_gender1","sick_gender2")
		cor_labels = c("pooled","a0","a1","b0","b1",
							"a0b0","a0b1",
							"a1b0","a1b1")
	
		V_posterior = Vs/dim(VV)[3]
		#print(V_posterior)
		#print(t(VV[,,50]))
		#print(V_posterior)
		V_mode = (V_posterior==apply(V_posterior,1,max))*1
		#print("V_mode")
		#print(dim(V_mode))
		#print(V_mode)
		#print(dim(x))
	
		# Tallennetaan kunkin ryhmän 'j' kuvaajat.
		for (j in 1:length(cor_labels)) {
	
			# Avataan tiedosto, johon kuvaaja tallennetaan.
			# 23.1.09: Kuvan tallentaminen ei toimi versiossa R-2.7.2 ilman
			# "type='cairo1'"-määrittelyä.
			png(paste(path,"cor-",cor_labels[j],".png",sep=""),
				width=840,height=840,type="cairo1")
			# Luodaan 'subplot' kuvaajia varten.
			layout(t(matrix(1:30,nrow=6)))
	
			# Piirretään klusterien kovarianssimatriisit.
			if (j==1)
				show_clusters(x,V_mode)
			if (j==2)
				show_clusters(x_h,V_mode)
			if (j==3 & sum(case)>0)
				show_clusters(x_s,V_mode)
			if (j==4)
				show_clusters(x_1,V_mode)
			if (j==5 & sum(gender>0))
				show_clusters(x_2,V_mode)
			if (j==6)
				show_clusters(x_h1,V_mode)
			if (j==7 & sum(gender)>0)
				show_clusters(x_h2,V_mode)
			if (j==8 & sum(case)>0)
				show_clusters(x_s1,V_mode)
			if (j==9 & sum(case&gender)>0)
				show_clusters(x_s2,V_mode)
		
			# Piirretään 'subplotiin' em. aikasarjat kullekin klusterille 'k'.
			for (k in 1:K) {
				screen(K+k)
	
				plot(Vsums[,k], type="l",col='blue',ylim=c(0,max(Vsums)),main=sprintf('N, K=%d',k))
				screen(2*K+k)
	
				if (j==1) {
					plot(corres_p[,k],type="l",col='magenta',ylim=c(0,1),
						main=sprintf('Cor, K=%d',k))
				}
				if (j==2) {
					plot(corres_h[,k],type="l",col='magenta',ylim=c(0,1),
						main=sprintf('Cor, K=%d',k))
				}
				if (j==3) {
					plot(corres_s[,k],type="l",col='magenta',ylim=c(0,1),
						main=sprintf('Cor, K=%d',k))
				}
				if (j==4) {
					plot(corres_1[,k],type="l",col='magenta',ylim=c(0,1),
						main=sprintf('Cor, K=%d',k))
				}
				if (j==5) {
					plot(corres_2[,k],type="l",col='magenta',ylim=c(0,1),
						main=sprintf('Cor, K=%d',k))
				}
				if (j==6) {
					plot(corres_h1[,k],type="l",col='magenta',ylim=c(0,1),
						main=sprintf('Cor, K=%d',k))
				}
				if (j==7) {
					plot(corres_h2[,k],type="l",col='magenta',ylim=c(0,1),
						main=sprintf('Cor, K=%d',k))
				}
				if (j==8) {
					plot(corres_s1[,k],type="l",col='magenta',ylim=c(0,1),
						main=sprintf('Cor, K=%d',k))
				}
				if (j==9) {
					plot(corres_s2[,k],type="l",col='magenta',ylim=c(0,1),
						main=sprintf('Cor, K=%d',k))
				}
			}
	
		# Suljetaan tiedosto.
		dev.off()
	
		}
	}

}

##########
##########

# Funktio 'corre_means'
# Tommi Suvitaival, 11.06.08
# Tehtävä: 
# Laskee klusteroidun datan klusterikohtaisten korrelaatioiden itseisarvon odotusarvon.
# Voidaan käyttää apuna konvergenssidiagnostiikassa.
# Parametrit:
# xxx - data (nrow='näytteiden lkm', ncol='näytteiden dimensio')
# VVV - klusterointi eli binäärinen matriisi, joka kertoo,
# 		  mihin klusteriin xxx:n osoittama vastaava datapiste kuuluu.
#		  Datapiste voi kuulua ainoastaan yhteen klusteriin.
#
# History:
# 27.3.09 - Handling for only one cluster added. In that case, the clusterin matrix is in fact a vector, and thus, only one correlation value is obtained.

corre_means = function(x,V) {
	
	#print(dim(V))
	#print("K");print(K)
	
	#print("m");print(m)
	klusteri <- array()
	#print(dim(V))
	# Käydään läpi kaikki 'K' klusteria.
	if (is.vector(V)) { # only one cluster
		M = length(V)
		# The data needs to be submitted to function 'cor' in matrix form, when auto-correlations are to be computed.
		corre_means = (sum(abs(cor(matrix(V))))-M)/(M^2-M)
	} else {
		K <- ncol(V) # Klusterien lukumäärä
		corre_means <- matrix(nrow=K,ncol=1)
		for (k in 1:K) {
			# Etsitään datamatriisista klusteriin 'k' kuuluvat näytteet.
			# Jotta voidaan laskea korrelaatio, on klusterissa oltava
			# enemmän kuin yksi näyte.
			if (sum(V[,k]==1)>1) {
				if (k==6) {
					#print(sum(V[,k]==1))
				}
				klusteri <- t(x[(V[,k]==1),])
				#print(dim(klusteri))
				# Lasketaan klusterin näytteiden
				# korrelaation itseisarvon odotusarvo.
				# 1. Muodostetaan klusterin näytteiden korrelaatiomatriisi.
				# 2. Otetaan korrelaatioista itseisarvot,
				#    koska ei haluta tehdä eroa 
				#	  positiivisen ja negatiivisen korrelaation välillä.
				# 3. Lasketaan korrelaatioiden klusterikohtainen keskiarvo.
				# 4. Poistetaan diagonaalilla oleva 1-korrelaatio keskiarvosta.
				#	  Otetaan siis huomioon vain diagonaalin ulkopuoliset näytteet.
				# 5. Kaava on sievennetty.

				# Kaava ilman diagonaalialkioiden poistoa keskiarvosta.
				# corre_means[k] = mean(abs(cor(klusteri)))

				M = ncol(klusteri)
				#print(dim(klusteri))
				corre_means[k] = (sum(abs(cor(klusteri)))-M)/(M^2-M)
			} else {
				# Jos klusterissa ei ole yhtään näytettä, asetetaan korrelaatioksi 0.
				#print("nollaa")
				corre_means[k] = 0
			}
		}
	}
	return(corre_means)
	
}
