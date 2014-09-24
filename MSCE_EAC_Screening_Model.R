# Model parameters: can depend on race/gender/GERD symptoms etc.
source('ex_parameter_list.R')
# Density for Barrett's esophagus onset times
source('BE_density_func.R')
# functions for MSCE-EAC model simulation and screening protocol functions
source('MSCE_EAC_screening_functions.R')

params = amparams
# age of index screen of BE patients	
screen_age = 60
# total simulation size
totalpop = 1000
# decide stem cell packing for biopsy protocol / mm^2 : rho = c_sigma * kstem_biop
kstem_biop = 50
#####################################################################################################################
######### 					BE cumulative function to draw using Barry O'Connor cohort        	 		  ###########
###########    mean = 5 cm , range [1,16], 22% SSBE and 78% LSBE, using beta distribution  a= 12/11, b=4  ###########
#####################################################################################################################
BEdist = rbeta(totalpop,16/11,4)*15+1 #in cm
BEdistmm = BEdist*10   #in mm
## Vectorize BE patients with random onset times, Tfinal = screen age - BE onset age : total time for process
## must develop BE before screening ages
BEparams = BEdensity(screen_age,gerd1,gerd2,gerd3,nu0)
nu_cum = BEparams$nu_cum
ages = 1:screen_age
BEstart = rep(100,totalpop)
for (k in 1:totalpop){
	while (BEstart[k] >= screen_age){
		x = runif(1,0,nu_cum[length(nu_cum)])
		if (x<nu_cum[1]){
			x = nu_cum[1]
		}
		temp_BE = approx(nu_cum, ages, xout=x)$y
		BEstart[k] = round(temp_BE)
	}
}
Tfinal = screen_age-BEstart
## KEEP TRACK OF FOLLOWING OUTPUT: ####
## Number of preinitiations
preinit = rep(0,totalpop)
## Number of total initiations 
totalinit = rep(0,totalpop)
## Premalignant clone number and sizes
Pnumbers = rep(0,totalpop)
Psizes = rep(0,500*totalpop)
## Malignant clones : Matrix with size of malignant clone, which individual k it is in, ###
###  and which premalignant clone (characterized by its size) it is in person k       ####
Msizes = matrix(rep(0,totalpop*50*3), nrow=3, ncol=(50*totalpop))
# detection of EA check
EAdetect = rep(0,totalpop)
sizecount = 0
msizecount = 0
#####################################################################################################################
######### 					MSCE-EAC multistage process hybrid simulation				      	 		  ###########
#####################################################################################################################
for ( k in 1:totalpop){
	BE = BEdistmm[k]*kstem*5000
	X = BE
	lambda1 = params[1]*X*Tfinal[k]
	N1 = rpois(1,lambda1)
	preinit[k] = N1
	## times of preinitiations are uniformly distributed:
	tau1 = 0
	tau1 = runif(N1,BEstart[k],screen_age)
	Ptimes = 0
	Ptimes = initiations(tau1,params[2], screen_age)
	if (Ptimes[1]>0){
		totalinit[k] = length(Ptimes)
		tmax = rep(0,length(Ptimes))
		tmax = screen_age-Ptimes
		pclonesizes = rep(0,length(Ptimes))
		mclonesizes = rep(0,length(Ptimes))
		for (i in 1:length(Ptimes)){
			out = 0
			out = pclone(params,n0=1,M=10000,tmax[i])
			pclonesizes[i] = out$finalsize
			mclonesizes[i] = out$maligclone
			if (out$EAdetect==1){ 
				## keep track of EA detections
				EAdetect[k] = 1
				break
			}	    
		}
		# keep track of malig clone sizes from which individual in which clone
		mcloneind = 0
		mtotal = sum(mclonesizes)
		if (mtotal>0){
			mcloneind = which(mclonesizes>0)
			Msizes[1,(msizecount+1):(msizecount+length(mcloneind))] = mclonesizes[mcloneind]
			Msizes[2,(msizecount+1):(msizecount+length(mcloneind))] = k
			Msizes[3,(msizecount+1):(msizecount+length(mcloneind))] = pclonesizes[mcloneind]
			msizecount = msizecount + length(mcloneind)
		}
		## keep track of premalig clone numbers and sizes
		pcloneind = 0
		if (sum(pclonesizes)>0){
			pcloneind = which(pclonesizes>0)
			Pnumbers[k] = length(pcloneind)
			Psizes[(sizecount+1):(sizecount+Pnumbers[k])] = pclonesizes[pcloneind]
			sizecount = sizecount + Pnumbers[k]
		}
	}	
}
## set up matrices with clone numbers and sizes for only those patients who did not develop clinical EAC before screening age
source('pre_biop_setup.R')

#####################################################################################################################
######### 	MSCE-EAC biopsy screening simulation for cancer-free BE patients at index endoscopy			  ###########
#####################################################################################################################
### NEOPLASIA PREVALENCES VIA BIOPSY SAMPLING (change to "biop.sample_abm" with final input for gamma diffusivity parameter for diffusive clone outputs)
pos.biops = biop.sample(nonEAPnumbers, nonEAPsizes, Msizes2, nonEAtotalpop, nonEABEdist,nonEABEdistmm, extramalig,kstem_biop)
## Following vectors have length of totalpop - clinical EAC cases: 
## each component 'i' is number of biopsies for BE patient 'i' that is positive for neoplasia, malignancy, or missed malignancy in neoplastic patient
premaligbiops = pos.biops$pbiop
maligbiops = pos.biops$mbiop
missedmalig = pos.biops$missedmalig
# cell counts adjusted if cells removed by biopsy specimen: used for projections of EAC into future or surveillance screening simulation
nonEAPsizes = pos.biops$Psizes
Msizes2 = pos.biops$Msizes
nonEAPnumbers = pos.biops$Pnumbers
extramalig = pos.biops$extramalig
pos.percentage = number.posbiops = mpos.percentage = mnumber.posbiops= missednumber= missedpercent = rep(0,10)
## For 10 different biopsy sensitivities, compute prevalences of neoplasia, malignancy, missed malignancy in neoplastic, and HGD without detected malignancy
for (k in 1:10){
 	kbiops = premaligbiops[k,]
 	number.posbiops[k] = length(kbiops[kbiops>0])
 	pos.percentage[k] = number.posbiops[k]/nonEAtotalpop
 	mkbiops = maligbiops[k,]
 	mnumber.posbiops[k] = length(mkbiops[mkbiops>0])
 	mpos.percentage[k] = mnumber.posbiops[k]/nonEAtotalpop
 	missed = missedmalig[k,]
 	## to get missed malignancy percentage, don't count patients with one missed malignancy but another detected
 	for (j in 1:nonEAtotalpop){
 		if (missed[j]>0 && mkbiops[j]>0){
 			missed[j] = 0
 		}
 	}
 	missednumber[k] = length(missed[missed>0])
 	missedpercent[k] = missednumber[k]/(number.posbiops[k]-mnumber.posbiops[k])
}
nonEACpop = nonEAtotalpop - mnumber.posbiops
### look at HGD percentage in the absense of prevalent, screen-detected EAC cases upon screen
justhgd = (number.posbiops-mnumber.posbiops)/nonEACpop



#####################################################################################################################
######### 	MSCE-EAC OCT imaging screening simulation for cancer-free BE patients at index endoscopy	  ###########
#####################################################################################################################
### PREVALENCES VIA IMAGING THRESHOLD
pos.imaging = imagingprev(kstem_biop,1,nonEAPnumbers, nonEAPsizes, Msizes2, nonEAtotalpop, extramalig)
premalig_image = pos.imaging$pimage
malig_image = pos.imaging$mimage
missed_image = pos.imaging$missedmalig
mnumber.posimage  = length(malig_image[malig_image>0])
mpos.percentage_i = mnumber.posimage/nonEAtotalpop
number.posimage  = length(premalig_image[premalig_image>0])
pos.percentage_i = number.posimage/nonEAtotalpop
print(pos.percentage_i)
nonEACpop_i = nonEAtotalpop - mnumber.posimage
### look at HGD percentage in the absense of prevalent EAC cases upon screen
justhgd_i = (number.posimage-mnumber.posimage)/nonEACpop_i
missed_number_i = length(missed_image[missed_image>0])
missed_percent_i = missed_number_i/number.posimage
