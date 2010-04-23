### Leave-one-out cross-validation
# 30.7.08
# Tommi Suvitaival

# Generoidaan satunnaisdata.
# Datan ominaisuutena on, että yhden klusterin latenttimuuttujan odotusarvo
# on eri terveille ja sairaille.

# Tämä tiedosto on peräisin Ilkalta 22.1.09.

nod_cv = function() {

	for (j in 1:J) {
	
		gen = gener(26,9,67,j) #gen = gener(26,5,67)
		x_gen = gen$x
		case_gen = gen$case
	
		K = 10 # klusterien lkm.
		N = dim(x_gen)[2] # havaintojen lkm.
	
		# Leave-one-out-silmukka, joka opettaa mallin 25 näytteellä ja testaa
		# oppimista yhdellä pois jätetyllä.
		N_iterations = 50
		N_burnin = 50

		# luokittelijan laskemat luokkatodennäköisyydet
		p_z_te = matrix(nrow=J,ncol=2) # 'nrow' on luokittelijoiden lkm.
	
		p_sick = matrix(nrow=N,ncol=1)
		success = matrix(nrow=N,ncol=1)
	
		d_mu_c = matrix(nrow=N,ncol=J)
		weights = matrix(nrow=N,ncol=J)
	
		for (i in 1:N) { # Kukin näyte 'i' on vuorollaan 'left-out'.
			if (i==1) { # jos ensimmäinen näyte
				x_tr = x_gen[,2:N]
				case_tr = case_gen[2:N]
				x_te = x_gen[,1]
				case_te = case_gen[1]
			} else {
				if (i==N) { # jos viimeinen näyte
					x_tr = x_gen[,1:(N-1)]
					case_tr = case_gen[1:(N-1)]
					x_te = x_gen[,N]
					case_te = case_gen[N]
				} else { # jos välillä oleva näyte
					x_tr = cbind(x_gen[,1:(i-1)],x_gen[,(i+1):N])
					case_tr = rbind(as.matrix(case_gen[1:(i-1)]),
								as.matrix(case_gen[(i+1):N]))
					x_te = x_gen[,i]
					case_te = case_gen[i]
				}
			}
			#aika = proc.time()
			samples = aja_annettudata(x_tr,t(case_tr),K,N_iterations,N_burnin)
			#aika = (proc.time()-aika)[1]
			#print("Naytteistysaika");print(aika)
	
			#aika = proc.time()
			# Kunkin metaboliitin 'm' posterioritodennäköisyys klusterissa 'k'
			# on 'pV[m,k]'.
			pV = samples$VS / (matrix(1,nrow=dim(samples$VS)[1],ncol=1) %*%
									colSums(samples$VS))
			# Järjestetään klusterit aiempia kierroksia vastaavaan järjestykseen.
			# Jos kyseessä on ensimmäinen CV-kierros, tätä ei tarvitse tehdä.
			if (i==1) {
				pV_avg = pV
				idx = 1:K
			} else {
				# Verrataan uusinta klusterointimatriisia aiempien kierrosten
				# keskimääräiseen klusterointimatriisiin.
				idx = compare_clusters(pV,pV_avg)
				#print("idx");print(idx);print("dim(idx)");print(dim(idx))
				pV = pV[,idx]
				# Päivitetään keskimääräinen klusterointi.
				#print(dim(pV));print(dim(pV_avg))
				pV_avg = (pV+(i-1)*pV_avg)/i
				# Skaalataan klusterointimatriisi jakaumaksi.
				pV_avg = pV_avg / 
							(matrix(1,nrow=dim(pV_avg)[1],ncol=1) %*% colSums(pV_avg))
			}
	
			# Etsitään suurimman käsittelykohtaisen eroavuuden omaavat klusterit
			# Lasketaan keskimääräiset hyperparametrin arvot
			# kussakin klusterissa molemmille käsittelyille.
			mu_c_s = apply(samples$Mu_cs[,idx],2,mean)
			mu_c_h = apply(samples$Mu_cs[,(idx+K)],2,mean)
	
			sorted = sort(abs(mu_c_s-mu_c_h),decreasing=TRUE,
									index.return=TRUE)
			classifier.idx = sorted$ix[1:J]
			classifier.weights = matrix(nrow=J,ncol=1)
			classifier.weights = sorted$x[1:J]/sum(sorted$x[1:J])
	
			d_mu_c[i,] = sorted$x[1:J]
			weights[i,] = classifier.weights
	
			x_te = t(as.matrix(x_te))
	
			z_te = x_te %*% pV[,classifier.idx]
	
			# Lasketaan uuden havainnon todennäköisyys terveiden ja sairaiden
			# havaintojen joukossa.
			p_z_te[1,] = dnorm(z_te,mean=mu_c_s[classifier.idx],
											sd=sd(samples$Mu_cs[,classifier.idx]))
			p_z_te[2,] = dnorm(z_te,mean=mu_c_h[classifier.idx],
											sd=sd(samples$Mu_cs[,(classifier.idx+K)]))
	
			p_sick[i] = sum(classifier.weights*p_z_te[1,]) / 
							sum(classifier.weights*(p_z_te[1,]+p_z_te[2,]))
			success[i] = (round(p_sick[i])==case_te)
	
			#aika = (proc.time()-aika)[1]
			#print("Analysointiaika");print(aika)
	
			print("sample");print(i)
			#print("classifier clusters");print(classifier.idx)
			print("p_sick");print(p_sick[i])
			print("case");print(case_te)
	
		}
		E_success[j] = mean(success)
		E_p_sick[j] = mean(p_sick)
	}

	# Laske onnistumisprosentti

	list(case=case_gen,p_sick=p_sick,success=success,weights=weights,d_mu_c=d_mu_c)
}

# 05.08.08
# Klusterien identifiointifunktio.
# Etsii toisiaan vastaavat klusterit matriiseista 'V1' ja V2'.
# Funktiossa oli bugi mutta se on nyt korjattu. Siitä, toimiiko funktio
# sillä tapaa kuin toivotaan, ei kuitenkaan ole vielä varmuutta. Bugeja
# ei pitäisi kuitenkaan enää olla.
compare_clusters = function(V_new,V,scale=TRUE) {

	if (ncol(V)>1) {
		# Skaalataan pystyvektorit ykköseen.
		# If the clustering contains a zero-vector, scaling is not a good idea, and it should be performed before adding any zero-vectors. -Tommi 12.3.09
		if (scale) {
			V1 = V_new / 
					t(as.matrix(colSums(V_new))%*%matrix(1,ncol=dim(V_new)[1],nrow=1))
			V2 = V / t(as.matrix(colSums(V))%*%matrix(1,ncol=dim(V)[1],nrow=1))
		} else {
			V1 = V_new
			V2 = V
		}

		K = dim(V1)[2] # klusterien lkm.
		# 'found' on jo löydetyt klusterit ilmaiseva vektori.
		V1_ind = c(1:K) #matrix(1:K,nrow=1,ncol=K)
		V2_ind = c(1:K) #matrix(1:K,nrow=1,ncol=K)
		assigns = matrix(nrow=1,ncol=K)

		for (k in 1:(K-1)) {
			C = dim(V1)[2] # jäljellä olevien klusterien lkm.

			# Etsitään maksimiarvo neliömuotoisesta 'CxC' 
			# "sisätulotulomatriisista"
			# ja valitaan maksimeista ensimmäinen.
			prod = t(V1)%*%V2
			d = which(max(prod)==prod)[1]
			r = ((d-1)%%C) + 1 # vaakarivi, kertoo 'V1':n klusterin
			c = floor((d-1)/C)+1 # pystyrivi, kertoo 'V2':n klusterin

			# Tallennetaan 
			assigns[V1_ind[r]] = V2_ind[c]
			if (r==1) { # 'V1':n ensimmäinen pystyrivi
				V1 = V1[,2:C]
				V1_ind = V1_ind[2:C]
			} else {
				if (r==C) { # 'V1':n viimeinen pystyrivi
					V1 = V1[,1:(C-1)]
					V1_ind = V1_ind[1:(C-1)]
				} else { # 'V1':n keskellä oleva pystyrivi
					V1 = cbind(V1[,1:(r-1)],V1[,(r+1):C])
					V1_ind = c(V1_ind[1:(r-1)],V1_ind[(r+1):C])
				}
			}
			if (c==1) { # 'V2':n ensimmäinen pystyrivi
				V2 = V2[,2:C]
				V2_ind = V2_ind[2:C]
			} else {
				if (c==C) { # 'V2':n viimeinen pystyrivi
					V2 = V2[,1:(C-1)]
					V2_ind = V2_ind[1:(C-1)]
				} else { # 'V2':n keskellä oleva pystyrivi
					V2 = cbind(V2[,1:(c-1)],V2[,(c+1):C])
					V2_ind = c(V2_ind[1:(c-1)],V2_ind[(c+1):C])
				}
			}
		}
		# Viimeinen klusteri
		assigns[V1_ind] = V2_ind
	} else { # If only one cluster, no comparison needs to be performed.
		assigns = 1
	}
	#print(assigns)
	assigns
}

log_p_values = function(data1,data2) {

	M = dim(data1)[1]
	pvals = matrix(nrow=M,ncol=1)

	# Lasketaan p-arvo kullekin metaboliitille 'm'.
	for (m in 1:M) {
		pvals[m] = log(t.test(x=data1[m,],y=data2[m,],
							alternative="greater",paired=F)$p.value) -
						log(t.test(x=data1[m,],y=data2[m,],
								alternative="less",paired=F)$p.value)
	}

	pvals[1:M]

}