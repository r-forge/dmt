# 4.8.09 - for plotting 2-level 3-way designs

plotMuIntervals = function(posteriorN,Ngen,file,stds=FALSE,mu.gen=NULL,plotQuantiles=TRUE) {
	K = ncol(posteriorN$mu_c)
	N = length(Ngen)
	if (!is.null(posteriorN$mu_s)) {
		png(paste(file,".png",sep=""),width=5000,height=1600,type="cairo1")
		layout(matrix(1:(7*K),nrow=K,byrow=T))
	} else {
		png(paste(file,".png",sep=""),width=2400,height=1600,type="cairo1")
		layout(matrix(1:(3*K),nrow=K,byrow=T))
	}
	for (k in 1:K) {
# 		if (sum(Vmode$mat[,k])>0) {
# 			mu_c.gen[n,k] = mu_c.gen[n,k]+mean(X[Vmode$mat[,k]==1,(case==1&gender==0)])/R
# 			mu_g.gen[n,k] = mu_g.gen[n,k]+mean(X[Vmode$mat[,k]==1,(case==0&gender==1)])/R
# 			mu_cg.gen[n,k] = mu_cg.gen[n,k]+mean(X[Vmode$mat[,k]==1,(case==1&gender==1)])/R
# 		}
		screen(3*k-2)
		#mu_intervals(Mu_c,k,Ngen,paste("mu_c",k,sep=""),path,mu_c.gen[,k])
		#mu_intervals(Mu_c,k,stdGen,effN=1,TRUE,mu_c.gen[,k])
		mu_intervals(posteriorN$mu_c[,,1:N,drop=F],k,Ngen,effN=1,stds=stds,mu.gen=NULL,plotQuantiles=plotQuantiles)
		screen(3*k-1)
		#mu_intervals(Mu_g,k,Ngen,paste("mu_g",k,sep=""),path,mu_g.gen[,k])
		#mu_intervals(Mu_g,k,stdGen,effN=2,TRUE,mu_g.gen[,k])
		mu_intervals(posteriorN$mu_g[,,1:N,drop=F],k,Ngen,effN=2,stds=stds,mu.gen=NULL,plotQuantiles=plotQuantiles)
		screen(3*k)
		#mu_intervals(Mu_cg,k,Ngen,paste("mu_cg",k,sep=""),path,mu_cg.gen[,k])
		#mu_intervals(Mu_cg,k,stdGen,effN=12,TRUE,mu_cg.gen[,k])
		mu_intervals(posteriorN$mu_cg[,,1:N,drop=F],k,Ngen,effN=12,stds=stds,mu.gen=NULL,plotQuantiles=plotQuantiles)
		
		# For plotting 2-level 3-way analysis -4.8.09
		if (!is.null(posteriorN$mu_s)) {
			print("Plotting 3-way")
			mu_intervals(posteriorN$mu_s[,,1:N,drop=F],k,Ngen,effN=3,stds=stds,mu.gen=NULL,plotQuantiles=plotQuantiles)	
			mu_intervals(posteriorN$mu_sc[,,1:N,drop=F],k,Ngen,effN=13,stds=stds,mu.gen=NULL,plotQuantiles=plotQuantiles)
			mu_intervals(posteriorN$mu_sg[,,1:N,drop=F],k,Ngen,effN=23,stds=stds,mu.gen=NULL,plotQuantiles=plotQuantiles)
			mu_intervals(posteriorN$mu_scg[,,1:N,drop=F],k,Ngen,effN=123,stds=stds,mu.gen=NULL,plotQuantiles=plotQuantiles)
			
		}
	}
	dev.off()
}
