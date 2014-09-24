## Load Legendre-Gauss quadrature function with Fortran wrapper for numerical integration: 
##    all provided in GitHub repo
dyn.load('legauss.so')
source('legauss.R')  
#######################################################
##      Function for initiations of P cells			###
#######################################################
initiations = function(tau,mu1,screen_age){
	n2 = rpois(length(tau),1*mu1*(screen_age-tau))
	if (length(which(n2!=0))>0){
		ind = which(n2!=0)
		n2 = n2[ind]
		taunew = tau[ind]
		tau2 = rep(0,sum(n2))
		ptimes = function(start,n){
			t = runif(n,start,screen_age)
			return(t)
		}
		count = 1
		for ( i in 1:length(n2)){
			tau2[count:(count+n2[i]-1)] = ptimes(taunew[i],n2[i])		
			count = count + n2[i]
		}
	}
	else {return(0)}
	return(tau2)	
}

###################################################################################################
##      Function for stochastic growth of P clone using SSA/tau-leap hybrid: b-d-m process		###
###################################################################################################
pclone = function(params,n0,M,tmax){
	## tau leap variables, decide accuracy of tau leap algorithm : see Supplementary info
	epsilon = .005
	tau = epsilon/(params[5]-params[3])
	leaprates= params[3:5]
	v= c(-1,0,1) 
	Mold=M
	rates=params[3:5]
	probs = rates/sum(rates)
	dn = sample(-1:1,M,prob=probs,replace=T)
	n = c(n0,cumsum(dn)+n0)
	totalmaligsize = 0
	## sojourn time in years to EAC STARTING with first P cell that creates the EAC
	p_malig_soj = NULL
	EAdetection =0
	if(!all(n>0)) {
	  M = min(which(n==0))-1   ## >=1 since n0 >=1
	  dn  = dn[1:M]
	  n   =  n[1:(M+1)]  ## n has length M+1 !!!
	}
	## event time
	Dt=cumsum(rexp(M,n[1:M]*sum(rates)))
	## This is for a realization that goes past tmax, whether it eventually goes extinct or not, and needs to be truncated
	if (Dt[M]>tmax && Dt[1]<=tmax){
		last_time = max(which(Dt<=tmax))
		Dt=Dt[1:last_time]
		n = n[1:(last_time+1)]
		dn=dn[1:last_time]
		if (length(which(dn==0))>0){
			mu2ind = which(dn==0) 		
			tmaxmalig = tmax - Dt[mu2ind]
			for (j in 1: length(tmaxmalig)){
				malout = mclone_times(params,tmaxmalig[j])
				#  print(malout$finalsize)
			    totalmaligsize = totalmaligsize + malout$finalsize
			    if ( malout$EA==1){ 
			      	## keep track of EA detections to condition
			      	EAdetection = 1
			      	p_malig_soj = min(p_malig_soj,(tmax-tmaxmalig[j]+malout$soj_time))
			    }
			}
		}	
	}
  	## This is for a clone that dies before tmax
  	else if (Dt[M]<=tmax && length(dn)<Mold){
  		if (length(which(dn==0))>0){
			mu2ind = which(dn==0) 
			tmaxmalig = tmax - Dt[mu2ind]
			for (j in 1: length(tmaxmalig)){
				malout = mclone_times(params,tmaxmalig[j])
			    totalmaligsize = totalmaligsize + malout$finalsize
			    if ( malout$EA==1){ 
			    	# keep track of EA detections to condition and time of mu2 trasformation that caused the EAC
			      	EAdetection = 1
			      	p_malig_soj = min(p_malig_soj,(tmax-tmaxmalig[j]+malout$soj_time))
			    }
			}
		}
  	}
  	## This is for a first event that happens later than tmax ie. no event occurs 
    else if (Dt[1]>tmax){
		n = n[1]
		dn=1
		Dt=0
 	}
 	## This is if M events in Gillespie does not reach tmax. eg. a very large clone piles up the exponentials
 	## Now we will implement a tau leap algorithm to speed up the growth of the larger clones in the simulation.	
 	else if (length(dn)==Mold){
 		if (length(which(dn==0))>0){
			mu2ind = which(dn==0) 	
			tmaxmalig = tmax - Dt[mu2ind]
			for (j in 1: length(tmaxmalig)){
				malout = mclone_times(params,tmaxmalig[j])
			    totalmaligsize = totalmaligsize + malout$finalsize
			    if ( malout$EA==1){ 
			      	## keep track of EA detections to condition 
			      	p_malig_soj = min(p_malig_soj,(tmax-tmaxmalig[j]+malout$soj_time))
			      	EAdetection = 1
			    }
		    }
		}
		## only continue with tau leap if no EA is detected in first M Gillespie steps
		if (EAdetection==0 && n[length(n)]!=0){
			## need to create new time left for tau leap subtracting Gillespie times from full P clone possible alive time, tmax
			tfinal = tmax-Dt[Mold]
			time0 = 0   #running time since start of tau leap algorithm
			## create new variable x to track P clone size
			x = n[length(n)]     
			while (time0 < tfinal && EAdetection==0){
				lambda = x*leaprates
				N= rep(0,3)
				N = rpois(3,lambda*tau)
				if (N[2]>0){
					for (j in 1: N[2]){
						tmaxmalig = tfinal-time0
					  	malout = mclone_times(params,tmaxmalig)
				 	  	totalmaligsize = totalmaligsize + malout$finalsize
				 	  	### Check if individual develops EA from malignant clone		
			     		if ( malout$EA==1){ 
			     			p_malig_soj = min(p_malig_soj,(Dt[Mold]+time0+malout$soj_time))
			    			EAdetection = 1
			     		}
			    	}
				}
				add=0
				add = N*v
				xtemp = x + sum(add)
				## If count goes negative, discard this change and try smaller tau
				if (xtemp <0){
					tau=tau/2
				}
				else{	
					x =xtemp
					time0 = time0 + tau
				}
				if (x==0){ break}
			}
			finalsize = x
			n[length(n)]=finalsize
		}		
 	}  
	# keept track of mu2 transformation times that cause EAC for EAC death distribution times
	return(list(events=dn,etime=Dt,n=n,tmax=tmax,finalsize=n[length(n)], maligclone=totalmaligsize, EAdetect=EAdetection, p_malig_soj=p_malig_soj))
}

################################################################################################################
##      Functions for stochastic size of M clone using analytic distribution: b-d-m process					 ###
## First computes if M clone becomes a clinically incident EAC case then, if not, size at time of screening  ###
################################################################################################################

mclone = function(params,tmax){	
	alpham= params[8]
	pm = params[9]
	qm = params[10]
	finalsize=0
	u = tmax
	zetat = exp(-pm*u)-exp(-qm*u)
	zetab  = (qm+alpham)*exp(-pm*u)-(pm+ alpham)*exp(-qm*u)
	zetacond = zetat/zetab
	p0t = zetacond*(alpham+pm)*(alpham+qm)*(qm*exp(-pm*u)-pm*exp(-qm*u))
	p0b = qm*(alpham+pm)*exp(-pm*u)-pm*(alpham+qm)*exp(-qm*u)
	p0 = p0t/p0b
	## survival function to check if detection rho events occurs by in t-s: Pr[Z(s,t)=0 | M(s,s)=1]
	survival = 1+(1/alpham)*((pm*qm*exp(-pm*tmax)-qm*pm*exp(-qm*tmax))/(qm*exp(-pm*tmax)-pm*exp(-qm*tmax)))
	EAcheck = sample(c(0,1),1,prob=c(survival,1-survival))
	if (EAcheck==0){
		x=runif(1,0,1)
	 	if (p0 >= x){
	 		finalsize = 0
		}
	 	if (p0 < x){
	 		temp = log((1-x)/(1-p0))/log(alpham*zetacond)
	 		finalsize = round(temp)
	    }
	}	
	return(list(tmax=tmax,finalsize=finalsize, EA = EAcheck))
}


## Samples sojourn times of EAC clinical detection beginning with 1 malignant cell at time tau, conditional on non-extinction
malig_sojourn = function(params, tau=0){
 	alphaM = params[8]
 	betaM = params[6]
 	pM = params[9]
 	qM = params[10]
 	# conditional on non-extinction
 	r = runif(1,0,1)
 	if ( r >= (1-betaM/alphaM)){
 		u =Inf
 	}
 	else{
 	u = -(1/(pM-qM))*log((pM*(alphaM*r+qM))/(qM*(alphaM*r+pM)))
 	}
 	return(list(sojourn = u, age_EAC = u+tau))
}

## Samples sojourn times starting with m number of M cells: important for surveillance screening
malig_sojourn_m = function(params, tau=0,mtotal){
 	alphaM = params[8]
 	betaM = params[6]
 	pM = params[9]
 	qM = params[10]
 	# conditional on non-extinction
 	print((1-(betaM/alphaM)^mtotal))
 	r = runif(1,0,1)
 	if( r >= (1-(betaM/alphaM)^mtotal)){
 		u =Inf
 	}
 	else{
 	u = -(1/(pM-qM))*log((pM*(alphaM*(-r+1)^(1/mtotal)-alphaM-qM))/(qM*(alphaM*(-r+1)^(1/mtotal)-alphaM-pM)))
 	}
 	return(list(sojourn = u, age_EAC = u+tau))
}

mclone_times = function(params,tmax){
	alpham= params[8]
	pm = params[9]
	qm = params[10]
	finalsize=0
	u = tmax
	zetat = exp(-pm*u)-exp(-qm*u)
	zetab  = (qm+alpham)*exp(-pm*u)-(pm+ alpham)*exp(-qm*u)
	zetacond = zetat/zetab
	EAcheck =0
	p0t = zetacond*(alpham+pm)*(alpham+qm)*(qm*exp(-pm*u)-pm*exp(-qm*u))
	p0b = qm*(alpham+pm)*exp(-pm*u)-pm*(alpham+qm)*exp(-qm*u)
	p0 = p0t/p0b
	## times of rho detection drawn and saved as EAC case if sojourn time < tmax or discarded and conditional size drawn if sojourn > tmax
	soj_time = malig_sojourn(params,tau=0)
	if (soj_time$sojourn >0 && soj_time$sojourn <= u){
		EAcheck = 1
		# NOTE: don't need malignant size at time tmax if rho event occurs, clinical detection takes place
	}
	if (EAcheck==0){
		x=runif(1,0,1)
	 	if (p0 >= x){
	 		finalsize = 0
		}
	 	if (p0 < x){
	 		temp = log((1-x)/(1-p0))/log(alpham*zetacond)
	 		finalsize = round(temp)
	    }
	}	
	return(list(tmax=tmax,finalsize=finalsize, EA = EAcheck, soj_time=soj_time$sojourn))
}

################################################################################################################
##      Functions for intersection of circular clone assumption with rectangular biopsies					 ###
## 			uses Legendre-Gauss quadrature to integrate intersection of circle with rectangle  				 ###
################################################################################################################

## This radius computed for assumption of 1000 crypts/ 15 mm^2 biopsy.. this assumption can be changed and will also change kstem_biop if you want rho to stay the same
cloneradii = function(kclonesize,kstem){
	radii = rep(0,length(kclonesize))
	for (i in 1:length(kclonesize)){
		crypts = kclonesize[i]/kstem
		one_crypt_r= sqrt(15/(1000*sqrt(12)))
		area = crypts*(one_crypt_r)^2*sqrt(12)
		radii[i]= sqrt(area/pi)
	}
	return(radii)
}


fplus = function(xc,yc,radius,point){
 # upper half of circle
 return(yc + sqrt(radius^2-(point-xc)^2))	
}
fminus = function(xc,yc,radius,point){
 	# lower half of circle
	return(yc - sqrt(radius^2-(point-xc)^2))
}
upperg = function(xc,yc,radius,point){
	## top function is either circle or top of rectangle (or discontinuous mix of both )
	Lbottom = -5/2
	Ltop = 5/2
	Kleft = -3/2
	Kright = 3/2
	uppcircle = fplus(xc,yc,radius,point)
	temp = min(Ltop,uppcircle)
	return(max(Lbottom,temp))
}

lowerg = function(xc,yc,radius,point){
	## lower function is either circle or bottom of rectangle (or discontinuous mix of both )
	Lbottom = -5/2
	Ltop = 5/2
	
	lowcircle = fminus(xc,yc,radius,point)
	temp = max(Lbottom,lowcircle)
	return(min(Ltop,temp))
}

leftbound = function(xc,r){
	## left circle or left side of rectangle
	Kleft = -3/2
	Kright = 3/2
	temp = max(Kleft,(xc-r))
	return(min(temp,Kright))	
}
rightbound = function(xc,r){
	## right circle or right side of rectangle
	Kleft = -3/2
	Kright = 3/2
	temp = min(Kright,(xc+r))
	return(max(temp,Kleft))
}

### Intersection function assumes SEATTLE PROTOCOL for biopsy sampling but can be adjusted to any protocol  ###
intersect.area = function(dims,xc,yc,radius){
	## biopsy quadrant dimensions for 4 quadrants af 75 cm circumference esophagus
	a = -75/8
	b = 75/8
	## heights of quadrant depend on patient-specific BE segment length
	c = dims[1]
	d = dims[2]
	## biopsy rectangle dimensions
	Lbottom = -5/2
	Ltop = 5/2
	Kleft = -3/2
	Kright = 3/2
	xc =xc
	yc=yc
	r=radius
	xlower = leftbound(xc,r)
	xupper = rightbound(xc,r)
	# return 0 if there is no intersection with biopsy rectangle
	if (xupper==xlower) { return(0)}
	## Legendre- Gauss quadrature for numerical integration of area intersection: ng = number of Gauss points
	ng=100
	GLquad=legauss(xlower,xupper,ng)
	xg = GLquad$mesh
	wg = GLquad$weights
	integrand = rep(0,ng)
	for (i in 1:ng){
		gplus=upperg(xc,yc,radius,xg[i])
		gminus=lowerg(xc,yc,radius,xg[i])
		integrand[i] = gplus-gminus
	}	 
	area = sum(wg*integrand)
	return(area)
}

### BIOPSY SAMPLING PROCEDURE FOR CHECKING HGD AND MALIGNANT CLONES: Assumes Seattle protocol, may be changed to any protocol
biop.sample = function(numbers, Psizes, Msizes, totalpop, BEdist, BEdistmm, extramalig,kstem){
	one_crypt_r= sqrt(15/(1000*sqrt(12)))
	## outputs results for biopsy sensitivity range of 10%-95%
	coverage= c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9)
	mcoverage = c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9)
	## used double sensitivity for malignant clones since biopsy is under closer inspection because of necessary detected neoplasia: can change this
	mcoverage=mcoverage/2
	positive.biop= matrix(0,10*totalpop,nrow=10,ncol=totalpop)
	mpositive.biop= matrix(0,10*totalpop,nrow=10,ncol=totalpop)
	missedm = matrix(0,10*totalpop,nrow=10,ncol=totalpop)
	extramaligpop = extramalig[2,]
	numbers_new= rep(0,totalpop)
	Psizes_new = rep(0,length(Psizes))
	Msizes_new = rep(0,length(Psizes))
	## CHECK BIOPSIES FOR neoplastic (P + M) cells then decide if any of biopsy was M cells
	sizes = Psizes + Msizes
	a = -75/8
	b = 75/8
	## for 5mm X 3mm biopsy	
	bioparea = 15
   	if (max(Msizes[1:numbers[1]])<(BEdistmm[1]*kstem*5000) && BEdistmm[1]>0){		
		if (numbers[1]>0){
			clonesize1 = sizes[1:numbers[1]]
			n = ceiling(BEdist[1]/2)
			kradii = cloneradii(clonesize1,kstem)
			quadrant = rep(0,length(kradii))
			### dimensions of quadrant for appropriate BE length for individual
			dims=rep(0,2)
			dims[1] = -BEdistmm[1]/(2*n)
			dims[2]= BEdistmm[1]/(2*n)
			area = rep(0,length(kradii))
			mclonesize1 = Msizes[1:numbers[1]]
			## usually size 0 for malignant compartment in each P clone
			mkradii = cloneradii(mclonesize1,kstem)
			marea = rep(0,length(mkradii))
			for (j in 1:length(kradii)){
			    quadrant[j] = sample(1:(n*4),1,prob=rep(1/(n*4),n*4)) 
			    ## x and y are coordinates of center for premalignant clones in BE segment
				u1 = runif(1)
				x= (b-a)*u1+a
				u2 = runif(1)
				y = (dims[2]-dims[1])*u2+dims[1]	
				area[j]= intersect.area(dims,x,y,kradii[j])	
				if (mkradii[j]!=0){
					## place center of malignant clone s.t. its entire area is contained in P+M clone
					r_crit = kradii[j]-mkradii[j]
					# pick M center within square of sides r_criterion
					aM = x-r_crit
					bM = x+r_crit
					cM= y-r_crit
					dM = y+r_crit		 
					u1m = runif(1)
					xm= (bM-aM)*u1m+aM
					u2m = runif(1)
					ym = (dM-cM)*u2m+cM	
					while (sqrt((xm-x)^2+(ym-y)^2)>(kradii[j]-mkradii[j])){		
						u1m = runif(1)
						xm= (bM-aM)*u1m+aM
						u2m = runif(1)
						ym = (dM-cM)*u2m+cM	
					}
					marea[j]= intersect.area(dims,xm,ym,mkradii[j])
				}
			}
			moverlap_cells = round((kstem*marea)/((one_crypt_r)^2*sqrt(12)))
			neo_overlap_cells = round((kstem*area)/((one_crypt_r)^2*sqrt(12)))
			#deplete neoplastic cells in biopsies from main counts
			Msizes[1:numbers[1]]= mclonesize1-moverlap_cells
			Psizes[1:numbers[1]] = clonesize1 - neo_overlap_cells - mclonesize1+moverlap_cells
			psize_temp = Psizes[1:numbers[1]]
			p_nonextinct_ind=0
			if (sum(psize_temp>0)){
				p_nonextinct_ind = which(psize_temp>0)
			}
			if (any(p_nonextinct_ind>0)){
				numbers_new[1]= length(psize_temp[psize_temp>0])
				Psizes_new[1:numbers_new[1]]= psize_temp[psize_temp>0]
				Msizes_new[1:numbers_new[1]] = Msizes[p_nonextinct_ind] 
			}
			total.area = rep(0,(n*4))
			mtotal.area = rep(0,(n*4))
			for (t in 1:(n*4)){
				for (l in 1:length(kradii)){
					if (quadrant[l]==t){
						total.area[t] = total.area[t] + area[l]
						mtotal.area[t] = mtotal.area[t] + marea[l] 
						
					}
				}
			}
			for (s in 1:10){
				for (m in 1:length(total.area)){
					if (total.area[m] >= (coverage[s]*bioparea)){
						hit=1
						positive.biop[s,1] = positive.biop[s,1] + hit
					}
					else { hit = 0}
					## check if malignant clone is detected in biopsy
					if (mtotal.area[m] >= (mcoverage[s]*bioparea) && hit==1){
						mhit = 1				
						mpositive.biop[s,1] = mpositive.biop[s,1]+mhit
					}
				}	
				if (sum(mtotal.area)>0){
					if(positive.biop[s,1]>0 && mpositive.biop[s,1]==0){
						missedm[s,1] = missedm[s,1]+1
					}
				}
			}
		}
	}
	else if(max(Msizes[1:numbers[1]])>=(BEdistmm[1]*kstem*5000) && BEdistmm[1]>0) {
		positive.biop[,1]= rep(1,length(coverage))
		mpositive.biop[,1]= rep(1,length(mcoverage))
	}
 	for (k in 2:totalpop){
 		if (numbers[k]>0){
 			kclonesize= sizes[(sum(numbers[1:(k-1)])+1):(sum(numbers[1:(k-1)])+numbers[k])]
 			kmclonesize =Msizes[(sum(numbers[1:(k-1)])+1):(sum(numbers[1:(k-1)])+numbers[k])]
 			if (max(kmclonesize)<(BEdistmm[k]*kstem*5000) && BEdistmm[k]>0){
				## usually size 0 for malignant compartment in each P clone
				kradii = cloneradii(kclonesize,kstem)
				mkradii = cloneradii(kmclonesize,kstem)
				n = ceiling(BEdist[k]/2)
				quadrant = rep(0,length(kradii))
				### dimensions of quadrant for appropriate BE length for individual
				dims=rep(0,2)
				dims[1] = -BEdistmm[k]/(2*n)
				dims[2]= BEdistmm[k]/(2*n)
				## area of individual clone intersect with biopsy
				area = rep(0,length(kradii))
				marea = rep(0,length(mkradii))
				for (j in 1:length(kradii)){
					quadrant[j] = sample(1:(n*4),1,prob=rep(1/(n*4),n*4)) 
					u1 = runif(1)
					x= (b-a)*u1+a
					u2 = runif(1)
					y = (dims[2]-dims[1])*u2+dims[1]
					area[j]= intersect.area(dims,x,y,kradii[j])
					if (mkradii[j]!=0){
						## place center of malignant clone s.t. its entire area is contained in P+M clone 
						r_crit = kradii[j]-mkradii[j]
						## pick M center within square of sides r_criterion
						aM = x-r_crit
						bM = x+r_crit
						cM= y-r_crit
						dM = y+r_crit	
						u1m = runif(1)
						xm= (bM-aM)*u1m+aM
						u2m = runif(1)
						ym = (dM-cM)*u2m+cM	
						while (sqrt((xm-x)^2+(ym-y)^2)>(kradii[j]-mkradii[j])){
							u1m = runif(1)
							xm= (bM-aM)*u1m+aM
							u2m = runif(1)
							ym = (dM-cM)*u2m+cM
						}
						marea[j]= intersect.area(dims,xm,ym,mkradii[j])
					}
				}
				moverlap_cells = round((kstem*marea)/((one_crypt_r)^2*sqrt(12)))
				neo_overlap_cells = round((kstem*area)/((one_crypt_r)^2*sqrt(12)))
				#deplete neoplastic cells in biopsies from main counts
				Msizes[(sum(numbers[1:(k-1)])+1):(sum(numbers[1:(k-1)])+numbers[k])]= kmclonesize-moverlap_cells
				Psizes[(sum(numbers[1:(k-1)])+1):(sum(numbers[1:(k-1)])+numbers[k])] = kclonesize - neo_overlap_cells - kmclonesize + moverlap_cells
				psize_temp = Psizes[(sum(numbers[1:(k-1)])+1):(sum(numbers[1:(k-1)])+numbers[k])]
				msize_temp = Msizes[(sum(numbers[1:(k-1)])+1):(sum(numbers[1:(k-1)])+numbers[k])]
				p_nonextinct_ind=0
				if (sum(psize_temp>0)){
					p_nonextinct_ind = which(psize_temp>0)
				}
				if (any(p_nonextinct_ind>0)){
					numbers_new[k]= length(psize_temp[psize_temp>0])
					Psizes_new[(sum(numbers_new[1:(k-1)])+1):(sum(numbers_new[1:(k-1)])+numbers_new[k])]= psize_temp[psize_temp>0]
					Msizes_new[(sum(numbers_new[1:(k-1)])+1):(sum(numbers_new[1:(k-1)])+numbers_new[k])] = msize_temp[p_nonextinct_ind] 
				}				
				total.area = rep(0,(n*4))
				mtotal.area = rep(0,(n*4))
				for (t in 1:(n*4)){
					for (l in 1:length(kradii)){
						if (quadrant[l]==t){
							total.area[t] = total.area[t] + area[l]
							mtotal.area[t] = mtotal.area[t] + marea[l]
						}
					}
				}
				for (s in 1:10){
					for (m in 1:length(total.area)){
						if (total.area[m] >= (coverage[s]*bioparea)){
							hitcheck=1
							hit=1
							positive.biop[s,k] = positive.biop[s,k] + hit
						}
						else {hitcheck = 0}
						## check if malignant clone is detected in biopsy
						if (mtotal.area[m] >= (mcoverage[s]*bioparea) && hitcheck==1){
							hitcheckm=1
							mhit = 1
							mpositive.biop[s,k] = mpositive.biop[s,k]+mhit
						}
						else {hitcheckm=0}	
					}	
					if (sum(mtotal.area)>0){
						if(positive.biop[s,k]>0 && mpositive.biop[s,k]==0){
							missedm[s,k] = missedm[s,k]+1
						}
					}
				}
			}
			else if(max(kmclonesize)>=(BEdistmm[k]*kstem*5000) && BEdistmm[k]>0){
				positive.biop[,k]= rep(1,length(coverage))
				mpositive.biop[,k]= rep(1,length(mcoverage))
			}
		}
	}
	for (f in 1:totalpop){
		if (length(extramaligpop[extramaligpop==f])>0 && BEdistmm[f]>0){
			maligind = which(extramaligpop==f)
			maligs = extramalig[1,maligind]
			mkradiiextra = cloneradii(maligs,kstem)
			n = ceiling(BEdist[f]/2)
			quadrant = rep(0,length(mkradiiextra))
			## dimensions of quadrant for appropriate BE length for individual
			dims= rep(0,2)
			dims[1] = -BEdistmm[f]/(2*n)
			dims[2]= BEdistmm[f]/(2*n)
			## area of individual clone intersect with biopsy
			marea = rep(0,length(mkradiiextra))
			for (j in 1:length(mkradiiextra)){
				quadrant[j] = sample(1:(n*4),1,prob=rep(1/(n*4),n*4)) 
				u1 = runif(1)
				x= (b-a)*u1+a
				u2 = runif(1)
				y = (dims[2]-dims[1])*u2+dims[1]
				marea[j]= intersect.area(dims,x,y,mkradiiextra[j])			
			}
			moverlap_cells = round((kstem*marea)/((one_crypt_r)^2*sqrt(12)))
			## deplete neoplastic cells in biopsies from main counts
			extramalig[1,maligind]= maligs-moverlap_cells
			for (check_malig in 1:length(maligind)){
				if (extramalig[1,maligind[check_malig]]<=0){
					extramalig[,maligind[check_malig]]=0
				}
			}
			mtotal.area = rep(0,(n*4))
			for (t in 1:(n*4)){
				for (l in 1:length(mkradiiextra)){
					if (quadrant[l]==t){
						mtotal.area[t] = mtotal.area[t] + marea[l]
					}
				}
			}
			for (s in 1:10){
				for (m in 1:length(mtotal.area)){
					## greater than HGD type coverage threshold 
					if (mtotal.area[m] >= (coverage[s]*bioparea)){
						#cat(mtotal.area[m], 'M biop cov', mcoverage[s]*bioparea, 'threshold \n')
						mhit = 1
						mpositive.biop[s,f] = mpositive.biop[s,f]+mhit
						## Also count as HGD
						positive.biop[s,f] = positive.biop[s,f]+mhit
					}					
				}
				if (sum(mtotal.area)>0){
					if(positive.biop[s,f]>0 && mpositive.biop[s,f]==0){
						missedm[s,f] = missedm[s,f]+1
					}
				}	
			}	
		}	
	}
	Psizes_new=Psizes_new[1:sum(numbers_new)]
	Msizes_new= Msizes_new[1:sum(numbers_new)]
 	return(list(pbiop=positive.biop,mbiop=mpositive.biop, missedmalig=missedm, Pnumbers=numbers_new,Psizes=Psizes_new, Msizes=Msizes_new, extramalig=extramalig))
}

## FOR AN INDIVIDUAL PERSON'S BIOPSY.  Not the originial method for prevalences of a totalpop in Screening paper results 
## but used for more intensive surveillance screening when each patient followed individually		
biop.sample_ind = function(numbers, Psizes, Msizes, BEdist, BEdistmm, extramalig, kstem){
	coverage= c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9)
	mcoverage = c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9)
	mcoverage=mcoverage/2
	numbers_new= 0
	Psizes_new = rep(0,length(Psizes))
	Msizes_new = rep(0,length(Psizes))
	extramalig_new =NULL
	positive.biop= matrix(0,10,nrow=10,ncol=1)
	mpositive.biop= matrix(0,10,nrow=10,ncol=1)
	missedm = matrix(0,10,nrow=10,ncol=1)
	one_crypt_r= sqrt(15/(1000*sqrt(12)))
	## CHECK BIOPSIES FOR P + M cells then decide if any of biopsy was M cells
	sizes = Psizes + Msizes
	a = -75/8
	b = 75/8
	### for 5mm X 3mm biopsy	
	bioparea = 15
   	if (max(Msizes[1:numbers])<(BEdistmm*kstem*5000) && BEdistmm>0){		
		if (numbers[1]>0){
			clonesize1 = sizes[1:numbers[1]]
			n = ceiling(BEdist[1]/2)
			kradii = cloneradii(clonesize1,kstem)
			quadrant = rep(0,length(kradii))
			### dimensions of quadrant for appropriate BE length for individual
			dims=rep(0,2)
			dims[1] = -BEdistmm/(2*n)
			dims[2]= BEdistmm/(2*n)
			area = rep(0,length(kradii))
			mclonesize1 = Msizes[1:numbers[1]]
			## usually size 0 for malignant compartment in each P clone
			mkradii = cloneradii(mclonesize1,kstem)
			marea = rep(0,length(mkradii))
			for (j in 1:length(kradii)){
			    quadrant[j] = sample(1:(n*4),1,prob=rep(1/(n*4),n*4)) 
			    ## x and y are coordinates of center for premalignant clones in BE segment
				u1 = runif(1)
				x= (b-a)*u1+a
				u2 = runif(1)
				y = (dims[2]-dims[1])*u2+dims[1]	
				area[j]= intersect.area(dims,x,y,kradii[j])	
				if (mkradii[j]!=0){
					## place center of malignant clone s.t. its entire area is contained in P+M clone
					r_crit = kradii[j]-mkradii[j]
					# pick M center within square of sides r_criterion
					aM = x-r_crit
					bM = x+r_crit
					cM= y-r_crit
					dM = y+r_crit		 
					u1m = runif(1)
					xm= (bM-aM)*u1m+aM
					u2m = runif(1)
					ym = (dM-cM)*u2m+cM	
					while (sqrt((xm-x)^2+(ym-y)^2)>(kradii[j]-mkradii[j])){		
						u1m = runif(1)
						xm= (bM-aM)*u1m+aM
						u2m = runif(1)
						ym = (dM-cM)*u2m+cM	
					}
					marea[j]= intersect.area(dims,xm,ym,mkradii[j])
				}
			}
			moverlap_cells = round((kstem*marea)/((one_crypt_r)^2*sqrt(12)))
			neo_overlap_cells = round((kstem*area)/((one_crypt_r)^2*sqrt(12)))
			#deplete neoplastic cells in biopsies from main counts
			Msizes[1:numbers[1]]= mclonesize1-moverlap_cells
			Psizes[1:numbers[1]] = clonesize1 - neo_overlap_cells - mclonesize1+moverlap_cells
			psize_temp = Psizes[1:numbers[1]]
			p_nonextinct_ind=0
			p_extinct_ind=0
			if (sum(psize_temp)>0){
				p_nonextinct_ind = which(psize_temp>0)
			}
			else if (sum(psize_temp)==0){
				p_extinct_ind= which(psize_temp==0)
			}
			if (any(p_nonextinct_ind>0)){
				numbers_new[1]= length(psize_temp[psize_temp>0])
				Psizes_new[1:numbers_new[1]]= psize_temp[psize_temp>0]
				Msizes_new[1:numbers_new[1]] = Msizes[p_nonextinct_ind] 
			}
			else if(all(p_nonextinct_ind==0)){
				numbers_new[1]= 0
				Psizes_new =0
				Msizes_new=0
			}
			if (any(p_extinct_ind>0)){
				extramalig_new = Msizes[p_extinct_ind]
			}	
			total.area = rep(0,(n*4))
			mtotal.area = rep(0,(n*4))
			for (t in 1:(n*4)){
				for (l in 1:length(kradii)){
					if (quadrant[l]==t){
						total.area[t] = total.area[t] + area[l]
						mtotal.area[t] = mtotal.area[t] + marea[l] 
						
					}
				}
			}
			for (s in 1:10){
				for (m in 1:length(total.area)){
					if (total.area[m] >= (coverage[s]*bioparea)){
						hit=1
						positive.biop[s,1] = positive.biop[s,1] + hit
					}
					else { hit = 0}
					## check if malignant clone is detected in biopsy
					if (mtotal.area[m] >= (mcoverage[s]*bioparea) && hit==1){
						mhit = 1				
						mpositive.biop[s,1] = mpositive.biop[s,1]+mhit
					}
				}	
				if (sum(mtotal.area)>0){
					if(positive.biop[s,1]>0 && mpositive.biop[s,1]==0){
						missedm[s,1] = missedm[s,1]+1
					}
				}
			}
		}
	}
	else if(max(Msizes[1:numbers])>=(BEdistmm*kstem*5000) && BEdistmm>0) {
		positive.biop[,1]= rep(1,length(coverage))
		mpositive.biop[,1]= rep(1,length(mcoverage))
	}
	if (length(extramalig)>0 && BEdistmm>0){
		maligs = extramalig
		mkradiiextra = cloneradii(maligs,kstem)
		n = ceiling(BEdist/2)
		quadrant = rep(0,length(mkradiiextra))
		### dimensions of quadrant for appropriate BE length for individual
		dims= rep(0,2)
		dims[1] = -BEdistmm/(2*n)
		dims[2]= BEdistmm/(2*n)
		## area of individual clone intersect with biopsy
		marea = rep(0,length(mkradiiextra))
		for (j in 1:length(mkradiiextra)){
			quadrant[j] = sample(1:(n*4),1,prob=rep(1/(n*4),n*4)) 
			u1 = runif(1)
			x= (b-a)*u1+a
			u2 = runif(1)
			y = (dims[2]-dims[1])*u2+dims[1]
			marea[j]= intersect.area(dims,x,y,mkradiiextra[j])			
		}
		moverlap_cells = round((kstem*marea)/((one_crypt_r)^2*sqrt(12)))
		## deplete neoplastic cells in biopsies from main counts
		extramalig= extramalig-moverlap_cells
		for (check_malig in 1:length(maligs)){
			if (extramalig[check_malig]<0){
				extramalig[check_malig]=0
			}
		}
		if (all(extramalig==0)){
			extramalig=NULL
		}
		mtotal.area = rep(0,(n*4))
		for (t in 1:(n*4)){
			for (l in 1:length(mkradiiextra)){
				if (quadrant[l]==t){
					mtotal.area[t] = mtotal.area[t] + marea[l]
				}
			}
		}
		for (s in 1:10){
			for (m in 1:length(mtotal.area)){
				## greater than HGD type coverage threshold 
				if (mtotal.area[m] >= (coverage[s]*bioparea)){
					#cat(mtotal.area[m], 'M biop cov', mcoverage[s]*bioparea, 'threshold \n')
					mhit = 1
					mpositive.biop[s,1] = mpositive.biop[s,1]+mhit
					## Also count as HGD
					positive.biop[s,1] = positive.biop[s,1]+mhit
				}					
			}
			if (sum(mtotal.area)>0){
				if(positive.biop[s,1]>0 && mpositive.biop[s,1]==0){
					missedm[s,1] = missedm[s,1]+1
				}
			}	
		}	
	}
	Psizes_new=Psizes_new[1:sum(numbers_new)]
	Msizes_new= Msizes_new[1:sum(numbers_new)]
	extramalig=c(extramalig_new,extramalig)
	return(list(pbiop=positive.biop,mbiop=mpositive.biop, missedmalig=missedm, Pnumbers=numbers_new,Psizes=Psizes_new, Msizes=Msizes_new, extramalig=extramalig))
}


#######################################################################################################################
##      Functions for growth of diffusive clones on a BE segment and Seattle protocol procedure				 		###
## 	counts crypts of certain size on hexagonal grid for intersection with biopsies for different sensitivities		###
#######################################################################################################################


source("MSCE_EAC_HexGrid.R")
source("MSCE_EAC_neighborList.R")
source("MSCE_EAC_genShape.R")


biop.sample_abm = function(numbers, Psizes, Msizes, totalpop, BEdist, BEdistmm, extramalig,kstem,gamma){
	coverage= c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9)
	mcoverage = c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9)
	mcoverage=mcoverage/2
	positive.biop= matrix(0,10*totalpop,nrow=10,ncol=totalpop)
	mpositive.biop= matrix(0,10*totalpop,nrow=10,ncol=totalpop)
	missedm = matrix(0,10*totalpop,nrow=10,ncol=totalpop)
	extramaligpop = extramalig[2,]
	## diffusivity of agent based grown clones 
	## CHECK BIOPSIES FOR P + M cells then decide if any of biopsy was M cells
	sizes = Psizes + Msizes
	a = -75/8
	b = 75/8
	### for 5mm X 3mm biopsy (can change this)	
	biop_crypts = 1012
	n = ceiling(BEdist[1]/2)
	## If malignancy is larger than number of stem cells in biopsy quadrant, 
	## consider it either covering entire area or large enough to be seen by naked endoscpoist eye
	totalcrypts_quad1 = BEdistmm[1]*1250/n
	totalcells_quad1 = totalcrypts_quad1*kstem
   	if (max(Msizes[1:numbers[1]])<(totalcells_quad1) && BEdistmm[1]>0){		
		if (numbers[1]>0){
			clonesize1 = sizes[1:numbers[1]]
			n = ceiling(BEdist[1]/2)
			clone_crypts1 = round(clonesize1/kstem)
			crypt_intersect1 =rep(0,length(clone_crypts1))
			quadrant = rep(0,length(clonesize1))
			### dimensions of quadrant for appropriate BE length for individual
			dims=rep(0,2)
			c=dims[1] = -BEdistmm[1]/(2*n)
			d=dims[2]= BEdistmm[1]/(2*n)
			xy = genHexGrid(c(a,c),c(b,d),totalcrypts_quad1)
			nbr = neighborList(xy)
			mclonesize1 = Msizes[1:numbers[1]]
			mcrypts1 = round(mclonesize1/kstem)
			## usually size 0 for malignant compartment in each P clone
			mcrypt_intersect1 = rep(0,length(mcrypts1))
			for (j in 1:length(clonesize1)){
			    quadrant[j] = sample(1:(n*4),1,prob=rep(1/(n*4),n*4)) 
				if (clonesize1[j]<totalcells_quad1){
					grow_clone1 =genShape(xy,nbr,size=clone_crypts1[j], gamma=gamma, col="lightpink2")
					x_premalig1 = grow_clone1$xpremalig
					y_premalig1 = grow_clone1$ypremalig
					for (p in 1:length(x_premalig1)){
						if (abs(x_premalig1[p])<3/2 && abs(y_premalig1[p])<5/2){
							crypt_intersect1[j]=crypt_intersect1[j]+1
						}
					}
					if (mcrypts1[j]!=0){
						## place center of malignant clone s.t. its entire area is contained in P+M clone
						N=length(x_premalig1)
						## Makig malig clone
						## pick random point in premalig clone
						m_center = sample(1:N,1,prob=rep(1/N,N))
						dist = rep(0,length(x_premalig1))
						for (k in 1:length(x_premalig1)){
							xd = x_premalig1[m_center]-x_premalig1[k]
						    yd = y_premalig1[m_center]-y_premalig1[k]
						    dist[k] = sqrt(xd*xd + yd*yd)
						}
						cutoff_size=sort(dist)[mcrypts1[j]]
						x_malig1= x_premalig1[dist<=cutoff_size]
						y_malig1 =y_premalig1[dist<=cutoff_size]
						for (p in 1:length(x_malig1)){
							if (abs(x_malig1[p])<3/2 && abs(y_malig1[p])<5/2){
								mcrypt_intersect1[j]=mcrypt_intersect1[j]+1
							}
						}
					}
				}
				else {
					positive.biop[,1]= rep(1,length(coverage))
					if (mcrypts1[j]!=0){
						mgrow_clone1 =genShape(xy,nbr,size=mcrypts1[j], gamma=gamma, col="lightpink2")
						x_malig1 = mgrow_clone1$xpremalig
						y_malig1 = mgrow_clone1$ypremalig
						for (p in 1:length(x_malig1)){
							if (abs(x_malig1[p])<3/2 && abs(y_malig1[p])<5/2){
								mcrypt_intersect1[j]=mcrypt_intersect1[j]+1
							}
						}
					}
				}	
			}
			total.area = rep(0,(n*4))
			mtotal.area = rep(0,(n*4))
			for (t in 1:(n*4)){
				for (l in 1:length(clonesize1)){
					if (quadrant[l]==t){
						total.area[t] = total.area[t] + crypt_intersect1[l]
						mtotal.area[t] = mtotal.area[t] + mcrypt_intersect1[l] 
						
					}
				}
			}
			for (s in 1:10){
				for (m in 1:length(total.area)){
					if (total.area[m] >= (coverage[s]*biop_crypts)){
						hit=1
						positive.biop[s,1] = positive.biop[s,1] + hit
					}
					else { hit = 0}
					## check if malignant clone is detected in biopsy
					if (mtotal.area[m] >= (mcoverage[s]*biop_crypts) && hit==1){
						mhit = 1				
						mpositive.biop[s,1] = mpositive.biop[s,1]+mhit
					}
				}	
				if (sum(mtotal.area)>0){
					if(positive.biop[s,1]>0 && mpositive.biop[s,1]==0){
						missedm[s,1] = missedm[s,1]+1
					}
				}
			}
		}
	}
	else if(max(Msizes[1:numbers[1]])>=(totalcells_quad1) && BEdistmm[1]>0) {
		positive.biop[,1]= rep(1,length(coverage))
		mpositive.biop[,1]= rep(1,length(mcoverage))
	}
 	for (k in 2:totalpop){
 		n = ceiling(BEdist[k]/2)
		## If malignancy is larger than number of stem cells in biopsy quadrant, 
		## consider it either covering entire area or large enough to be seen by naked endoscpoist eye
		totalcrypts_quadk = BEdistmm[k]*1250/n
		totalcells_quadk = totalcrypts_quadk*kstem
 		if (numbers[k]>0){
 			kclonesize= sizes[(sum(numbers[1:(k-1)])+1):(sum(numbers[1:(k-1)])+numbers[k])]
 			kmclonesize =Msizes[(sum(numbers[1:(k-1)])+1):(sum(numbers[1:(k-1)])+numbers[k])]
			
 			if (max(kmclonesize)<(totalcells_quadk) && BEdistmm[k]>0){
 				clone_cryptsk = round(kclonesize/kstem)
				crypt_intersectk =rep(0,length(clone_cryptsk))
				### dimensions of quadrant for appropriate BE length for individual
				dims=rep(0,2)
				c=dims[1] = -BEdistmm[k]/(2*n)
				d=dims[2]= BEdistmm[k]/(2*n)
				mcryptsk = round(kmclonesize/kstem)
				## usually size 0 for malignant compartment in each P clone
				mcrypt_intersectk = rep(0,length(mcryptsk))
				quadrant = rep(0,length(kclonesize))
 				xy = genHexGrid(c(a,c),c(b,d),totalcrypts_quadk)
				nbr = neighborList(xy)
 				for (j in 1:length(kclonesize)){
	 				quadrant[j] = sample(1:(n*4),1,prob=rep(1/(n*4),n*4)) 
	 				if (kclonesize[j]<totalcells_quadk){
						grow_clonek =genShape(xy,nbr,size=clone_cryptsk[j], gamma=gamma, col="lightpink2")
						x_premaligk = grow_clonek$xpremalig
						y_premaligk = grow_clonek$ypremalig
						for (p in 1:length(x_premaligk)){
							if (abs(x_premaligk[p])<3/2 && abs(y_premaligk[p])<5/2){
								crypt_intersectk[j]=crypt_intersectk[j]+1
							}
						}
						if (mcryptsk[j]!=0){
							## place center of malignant clone s.t. its entire area is contained in P+M clone
							N=length(x_premaligk)
							## Makig malig clone
							## pick random point in premalig clone
							m_center = sample(1:N,1,prob=rep(1/N,N))
							dist = rep(0,length(x_premaligk))
							for (t in 1:length(x_premaligk)){
								xd = x_premaligk[m_center]-x_premaligk[t]
							    yd = y_premaligk[m_center]-y_premaligk[t]
							    dist[t] = sqrt(xd*xd + yd*yd)
							}
							cutoff_size=sort(dist)[mcryptsk[j]]
							x_maligk= x_premaligk[dist<=cutoff_size]
							y_maligk =y_premaligk[dist<=cutoff_size]
							for (p in 1:length(x_maligk)){
								if (abs(x_maligk[p])<3/2 && abs(y_maligk[p])<5/2){
									mcrypt_intersectk[j]=mcrypt_intersectk[j]+1
								}
							}
						}
					}
					else {
						print(kclonesize[j]-totalcells_quadk)
						positive.biop[,k]= rep(1,length(coverage))
						if (mcryptsk[j]!=0){
							mgrow_clonek =genShape(xy,nbr,size=mcryptsk[j], gamma=gamma, col="lightpink2")
							x_maligk = mgrow_clonek$xpremalig
							y_maligk = mgrow_clonek$ypremalig
							for (p in 1:length(x_maligk)){
								if (abs(x_maligk[p])<3/2 && abs(y_maligk[p])<5/2){
									mcrypt_intersectk[j]=mcrypt_intersectk[j]+1
								}
							}
						}
					}
				}
				total.area = rep(0,(n*4))
				mtotal.area = rep(0,(n*4))
				for (t in 1:(n*4)){
					for (l in 1:length(kclonesize)){
						if (quadrant[l]==t){
							total.area[t] = total.area[t] + crypt_intersectk[l]
							mtotal.area[t] = mtotal.area[t] + mcrypt_intersectk[l]
						}
					}
				}
				for (s in 1:10){
					for (m in 1:length(total.area)){
						if (total.area[m] >= (coverage[s]*biop_crypts)){
							hitcheck=1
							hit=1
							positive.biop[s,k] = positive.biop[s,k] + hit
						}
						else {hitcheck = 0}
						## check if malignant clone is detected in biopsy
						if (mtotal.area[m] >= (mcoverage[s]*biop_crypts) && hitcheck==1){
							hitcheckm=1
							mhit = 1
							mpositive.biop[s,k] = mpositive.biop[s,k]+mhit
						}
						else {hitcheckm=0}	
					}	
					if (sum(mtotal.area)>0){
						if(positive.biop[s,k]>0 && mpositive.biop[s,k]==0){
							missedm[s,k] = missedm[s,k]+1
						}
					}
				}
			}
			else if(max(kmclonesize)>=(totalcells_quadk) && BEdistmm[k]>0){
				positive.biop[,k]= rep(1,length(coverage))
				mpositive.biop[,k]= rep(1,length(mcoverage))
				cat('total coverage', max(kclonesize), 'max size', k, 'person with BE: ', (BEdistmm[k]*kstem*5000), '\n' )
			}
		}
	}
	for (f in 1:totalpop){
		if (length(extramaligpop[extramaligpop==f])>0 && BEdistmm[f]>0){
			maligind = which(extramaligpop==f)
			maligs = extramalig[1,maligind]
			mclone_crypts = round(maligs/kstem)
			n = ceiling(BEdist[f]/2)
			quadrant = rep(0,length(maligs))
			### dimensions of quadrant for appropriate BE length for individual
			dims= rep(0,2)
			dims[1] = -BEdistmm[f]/(2*n)
			dims[2]= BEdistmm[f]/(2*n)
			## area of individual clone intersect with biopsy
			mcrypt_intersectk = rep(0,length(maligs))
			for (j in 1:length(maligs)){
				quadrant[j] = sample(1:(n*4),1,prob=rep(1/(n*4),n*4)) 
				mgrow_clonek =genShape(xy,nbr,size=mclone_crypts[j], gamma=gamma, col="lightpink2")
				x_maligk = mgrow_clonek$xpremalig
				y_maligk = mgrow_clonek$ypremalig
				for (p in 1:length(x_maligk)){
					if (abs(x_maligk[p])<3/2 && abs(y_maligk[p])<5/2){
						mcrypt_intersectk[j]=mcrypt_intersectk[j]+1
					}
				}
			}
			mtotal.area = rep(0,(n*4))
			for (t in 1:(n*4)){
				for (l in 1:length(maligs)){
					if (quadrant[l]==t){
						mtotal.area[t] = mtotal.area[t] + mcrypt_intersectk[l]
					}
				}
			}
			for (s in 1:10){
				for (m in 1:length(mtotal.area)){
					## greater than HGD type coverage threshold 
					if (mtotal.area[m] >= (coverage[s]*biop_crypts)){
						mhit = 1
						mpositive.biop[s,f] = mpositive.biop[s,f]+mhit
						## Also count as HGD
						positive.biop[s,f] = positive.biop[s,f]+mhit
					}					
				}
				if (sum(mtotal.area)>0){
					if(positive.biop[s,f]>0 && mpositive.biop[s,f]==0){
						missedm[s,f] = missedm[s,f]+1
					}
				}	
			}	
		}	
	}
 	return(list(pbiop=positive.biop,mbiop=mpositive.biop, missedmalig=missedm))
}

##########################################################################################################
##      Functions for neoplastic prevalences with a hypothetical high-resolution imaging screen		   ###
## 			requires an input for area resoluton threshold size on surface in mm^2 					   ###
##########################################################################################################

imagingprev = function(k_per_crypt,areathreshold,numbers, Psizes, Msizes, totalpop, extramalig){
	## ASSUMPTION mm^2 surface area threshold for detection via imaging
	min_crypts = areathreshold/((0.0658037)^2*sqrt(12))
	min_cells= min_crypts*k_per_crypt
	sizes = Psizes+Msizes
	positiveHGD= rep(0,totalpop)
	positiveMALIG=rep(0,totalpop)
	missedm = rep(0,totalpop)
	extramaligpop = extramalig[2,]
		if (numbers[1]>0){
		clonesize1 = sizes[1:numbers[1]]
		mclonesize1 = Msizes[1:numbers[1]]
		## usually size 0 for malignant compartment in each P clone
		for (j in 1:length(clonesize1)){
			if (clonesize1[j]>=min_cells){
				positiveHGD[1]=1
			}
			if (mclonesize1[j]>=min_cells){
				positiveMALIG[1]=1
			}
		}
		if (positiveHGD[1]==1 && positiveMALIG[1]==0 && sum(mclonesize1)>0){
			missedm[1]=1
		}
	}	
	for (k in 2:totalpop){
		if (numbers[k]>0){	
			kclonesize= sizes[(sum(numbers[1:(k-1)])+1):(sum(numbers[1:(k-1)])+numbers[k])]
			kmclonesize =Msizes[(sum(numbers[1:(k-1)])+1):(sum(numbers[1:(k-1)])+numbers[k])]
			## usually size 0 for malignant compartment in each P clone
			for (j in 1:length(kclonesize)){
				if (kclonesize[j]>=min_cells){
					positiveHGD[k]=1
				}
				if (kmclonesize[j]>=min_cells){
					positiveMALIG[k]=1
				}
			}
			if (positiveHGD[k]==1 && positiveMALIG[k]==0 && sum(kmclonesize)>0){
				missedm[k]=1
			}
		}		
	}
	for (f in 1:totalpop){
		if (length(extramaligpop[extramaligpop==f])>0){
			maligind = which(extramaligpop==f)
			maligs = extramalig[1,maligind]
			for (j in 1:length(maligs)){
				if (maligs[j]>=min_cells){
					positiveMALIG[f]=1
				}	
			}	
		}
	}		
	return(list(pimage=positiveHGD,mimage=positiveMALIG, missedmalig=missedm))
}

