##########################################################################################
##                  FUTURE EAC INCIDENCE PROJECTIONS AFTER ABLATION                     ##
## PERFORM ABLATION ON HGD PATIENTS: See Supplementary for hazard function derivations  ##
##########################################################################################


## all male parameter example
val = c(allmales$X,allmales$nu, allmales$mu0,allmales$peff,allmales$qeff,allmales$tlag, allmales$alphaP)
## ODE solver package
library('deSolve')
## 4-stage model WITH BE CONVERSION 
phidot = function(s, phi, parms){
  phidot=numeric(6)
  with(as.list(parms),{
  tau = t - s
  RR=5
  GERDpr=1-exp(-gerd1*min(gerd3,tau)-gerd2*max(0,tau-gerd3))
  nu = nu0*((1-GERDpr)+RR*GERDpr)
  phidot[1]=betam -(alpham+betam+rho)*phi[1]+alpham*(phi[1]**2)
  phidot[2]=2*alpham*phi[1]*phi[2]- (alpham+betam+rho)*phi[2]
  phidot[3]=betap+ mu2*phi[1]*phi[3]-(alphap+betap+mu2)*phi[3]+alphap*(phi[3]**2)
  phidot[4]=2*alphap*phi[3]*phi[4]+mu2*(phi[4]*phi[1]+phi[3]*phi[2]) - (alphap+betap+mu2)*phi[4]
  phidot[5]=mu1*phi[5]*(phi[3]-1)
  phidot[6]=mu1*(phi[6]*(phi[3]-1)+phi[5]*phi[4])
  phidot[7]=mu0*phi[7]*(phi[5]-1)
  phidot[8]=mu0*(phi[8]*(phi[5]-1)+phi[7]*phi[6])
  ## A-D transition
  phidot[9]=nu*(phi[7]-phi[9])
  phidot[10]=nu*(phi[8]-phi[10])
  list(c(phidot))
  })
}
#conditional on NOT HAVING BE before screen_age
## 3-stage hazard function ##
h3=function(val,t){
   X=val[1] 
   tlag = val[6]
   nu=val[2]
   mu0=mu1=val[3]
   p=val[4]
   q=val[5]
   alpha = val[7]
   df=q*exp(-p*(t-tlag))-p*exp(-q*(t-tlag))
   ef=mu1/alpha
   return(mu0*X*(1.-((q-p)/df)**ef))
}
## 3-stage survival function ###
s3 = function(val,t) {
   X=val[1]
   nu=val[2]
   mu0=mu1=val[3]
   p=val[4]
   q=val[5]
   tlag = val[6]
   alpha = val[7]
   H3 = numeric(length(t))
   for(i in 1:length(t)) {
		gquad2 = legauss(0,t[i],50)
		xg2 = gquad2$mesh
		wg2 = gquad2$weights
		f1 = q*exp(-p*(xg2-tlag))
		f2 = p*exp(-q*(xg2-tlag))
		f = f1-f2
		ef=mu1/alpha
		y2.tu = ((q-p)/f)**ef
		H3[i]= sum(wg2*(1-y2.tu))
	}
	return(exp(-mu0*X* H3))
}

density5stagecond=function(val,t,gquad,screen_age){
   X= val[1]
   nu=val[2]
   mu0=mu1=val[3]
   p=val[4]
   q=val[5]
   mtlag = val[6]
   alpha=val[7]    ##alpha p for premalignant
   RR=5
   xg = gquad$mesh
   wg = gquad$weights
   f_BE = rep(0,length(xg))
   for ( j in 1:length(xg)){
    GLquad2=legauss(0 ,xg[j],20)
    wg2 = GLquad2$weights
    xg2 = GLquad2$mesh
    int1 = rep(0,20)
    for( k in 1:20){
      GERDpr1=1-exp(-gerd1*min(gerd3,xg2[k])-gerd2*max(0,(xg2[k]-gerd3)))
      amnu_temp =amnu0*((1-GERDpr1)+RR*GERDpr1)
      int1[k] = amnu_temp
    }
    nu_surv = exp(-sum(wg2*int1))
    GERDpr=1-exp(-gerd1*min(gerd3,xg[j])-gerd2*max(0,(xg[j]-gerd3)))
    amnu_temp2 =amnu0*((1-GERDpr)+RR*GERDpr)
    amnu_dens1= amnu_temp2*nu_surv
    f_BE[j] = amnu_dens1
  }
   f.n=f.d = numeric(length(t))
   ts = t-xg      # t-u
   f.n = 1/(1-F_a)*sum(wg*f_BE*h3(val,ts)*s3(val,ts))
   return(f.n)
}
X=10^6
RR = 5
GLquad2=legauss(0 ,screen_age,20)
wg2 = GLquad2$weights
xg2 = GLquad2$mesh
int1 = rep(0,20)
for( k in 1:20){
  GERDpr1=1-exp(-gerd1*min(gerd3,xg2[k])-gerd2*max(0,(xg2[k]-gerd3)))
  nu_temp =nu0*((1-GERDpr1)+RR*GERDpr1)
  int1[k] = nu_temp
}
F_a = 1-exp(-sum(wg2*int1))

#### FOR UNSCREENED BACKGROUND CUMUL:
parms = c(nu=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

ages = 1:80
unscreened_cumul_f= rep(0,80)
unscreened_cumul_haz= rep(0,80)

for (j in 1:length(ages)){
  GLquad1=legauss(0,ages[j],20)
  xg1 = GLquad1$mesh
  wg1 = GLquad1$weights
  temphaz = NULL
  tempf = NULL
  for (i in 1:20){
    parms['t'] = xg1[i]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
    times = c(0,xg1[i])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    # hazard
    temphaz = c(temphaz,-ode$phi10[2]/ode$phi9[2])
    # density
    tempf = c(tempf, -ode$phi10[2])
  } 
  unscreened_cumul_f[j] = sum(wg1*tempf)
  unscreened_cumul_haz[j] = sum(wg1*temphaz)
}

## Make matrix A of cell types for m individuals screened time \sigma before ablation
A = matrix(rep(0,nonEAtotalpop*4),ncol=4,nrow=nonEAtotalpop)
## BE stem cells
A[,1] = nonEABEdistmm*5000*kstem
## P^* preinitiated cells
A[,2] = nonEApreinit
A[1,3] = sum(nonEAPsizes[1:nonEAPnumbers[1]])
A[1,4] = sum(Msizes2[1:nonEAPnumbers[1]])
## find how many malignant clones with no premalignant surviving clones
nonEApop = which(EAdetect==0)
extramalignum = length(extramalig[1,])
for (k in 2:nonEAtotalpop){
  if (nonEAPnumbers[k]>0){
    kclonesize= nonEAPsizes[(sum(nonEAPnumbers[1:(k-1)])+1):(sum(nonEAPnumbers[1:(k-1)])+nonEAPnumbers[k])]
    kmclonesize =Msizes2[(sum(nonEAPnumbers[1:(k-1)])+1):(sum(nonEAPnumbers[1:(k-1)])+nonEAPnumbers[k])]
    ## P initiated clones together
    A[k,3]= sum(kclonesize)
    ## M malignant clones together
    A[k,4] = sum(kmclonesize)  
    for (j in 1:extramalignum){
      	if (extramalig[2,j]==nonEApop[k]){
    	  A[k,4] = A[k,4]+extramalig[2,j]
      	}
    }
  }     
}
##########################################################################################
##   Examples for various ablation proportion vectors omega: Recreates Figure 7.        ##
##########################################################################################

## NEED TO INTEGRATE \int_0^t h(x)dx = F(t) to compare cumulative hazards for partial ablation cases
B=A
## example omega matrix
proportion=matrix(c(.5,.5,.5,.5,.01,.01,.01,.01,1,1,1,0,1,1,0,0,0,0,0,0),nrow=5,byrow=T)
## 40% coverage threshold 
k=5
kbiops = premaligbiops[k,]
popdetect = which(kbiops>0)
mbiops = maligbiops[k,]
popmaligdetect = which(mbiops>0)
combined_prev = c(popdetect,popmaligdetect)
# gives indices of popdetect that are also in popmaligdetect (we only want to treat HGD patients without malig patients)
malig_ind = which(duplicated(combined_prev, fromLast=TRUE))
popdetect = popdetect[-malig_ind]
## nonEApop should be denominator: 1) including screen-detected, old version
#nonEApop = totalpop-sum(EAdetect)
# or 2) without screen-detected
nonEApop = totalpop-sum(EAdetect)-length(popmaligdetect)
for (ablate in 1:4){
  A=B
  ### imperfect ablation, leave, 10% of each cell type
  ## Remove pos. malignant population
  A[popmaligdetect,1]= 0*B[popmaligdetect,1]
  A[popmaligdetect,2]= 0*B[popmaligdetect,2]
  A[popmaligdetect,3]= 0*B[popmaligdetect,3]
  A[popmaligdetect,4]= 0*B[popmaligdetect,4]
  ## HGD without pos. malig
  A[popdetect,1]= proportion[ablate,1]*B[popdetect,1]
  A[popdetect,2]= proportion[ablate,2]*B[popdetect,2]
  A[popdetect,3]= proportion[ablate,3]*B[popdetect,3]
  A[popdetect,4]= proportion[ablate,4]*B[popdetect,4]

  prescreen_cumul_f = unscreened_cumul_f[60]
  prescreen_cumul_haz = unscreened_cumul_haz[60]

  ages = 61:80
  postscreen_cumul_f= rep(0,20)
  postscreen_cumul_haz = rep(0,20)

  for (j in 1:length(ages)){
    guass_points=20
    GLquad1=legauss(60,ages[j],guass_points)
    xg1 = GLquad1$mesh
    wg1 = GLquad1$weights
    temphaz = rep(0,guass_points)
    postscreendensity= rep(0,guass_points)
    screeneddensity= rep(0,guass_points)

    for (i in 1:guass_points){
      ## for one individual i, all values A_i={X_i,P^*_i, P_i,M_i} stored in ith row of matrix A m x 4 matrix. m=totalpop
      screenedsurv = rep(0,nonEApop)
      screenedhaz = rep(0,nonEApop)
      parms = c(nu=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)
      parms['t'] = xg1[i]-60
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
      times = c(0,xg1[i]-60)
      ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
      s4 = ode$phi7[2]
      su3 = ode$phi5[2]
      s2 = ode$phi3[2]
      s1 = ode$phi1[2]
      h4 = -ode$phi8[2]/ode$phi7[2]
      hz3= -ode$phi6[2]/ode$phi5[2]
      h2= -ode$phi4[2]/ode$phi3[2]
      h1 = -ode$phi2[2]/ode$phi1[2]
      for (g in 1:nonEApop){
        cells = A[g,] 
        screenedsurv[g] = ((s4)^cells[1])*(su3^cells[2])*(s2^cells[3])*(s1^cells[4])
        screenedhaz[g] = (h4*cells[1])+(hz3*cells[2])+(h2*cells[3])+(h1*cells[4])
      }
      screeneddensity[i]= (1/(nonEApop))*sum(screenedsurv*screenedhaz)
      if (xg1[i]>val[6]){ 
        GLquad2=legauss(60 ,(xg1[i]-val[6]),20)
        hz=density5stagecond(val,xg1[i],GLquad2)
        postscreendensity[i]= hz 
      }
      else{
        postscreendensity[i]= 0
      }
    } 
    ## f(t): must multiply by fractions of people with BE <60 and BE > 60  
    ## GERD nu 
    temphaz= (F_a)*screeneddensity+(1-F_a)*postscreendensity
    postscreen_cumul_f[j] = sum(wg1*temphaz)
    # post screen haz with cumul sum of density rather than full integral
    postscreen_cumul_haz[j] = -log(1-(postscreen_cumul_f[j]+prescreen_cumul_f))
  }
  if (ablate==1){
    amcumul_inc_50ablate = c(unscreened_cumul_haz[25:60],postscreen_cumul_haz)
  }
  if (ablate == 2){
    amcumul_inc_99ablate2 = c(unscreened_cumul_haz[25:60],postscreen_cumul_haz)
  }
  if (ablate==3){
    amcumul_inc_maligablate2 = c(unscreened_cumul_haz[25:60],postscreen_cumul_haz)  
  }
  if (ablate==4){
    amcumul_inc_hgdmaligablate2 = c(unscreened_cumul_haz[25:60],postscreen_cumul_haz) 
  }
  if (ablate==5){
    amcumul_inc_allablate2 = c(unscreened_cumul_haz[25:60],postscreen_cumul_haz)
  }
}
inc=100000
g_rangehgd2=range(unscreened_cumul_haz*inc, amcumul_inc_hgdmaligablate2*inc)
quartz("Quartz", width=3.27,height=3.4)
par(ps=8,mar=c(2,2,.5,.5), mgp=c(1,.25,0))
colors_treat=c("blue","red" ,"darkmagenta", "dodgerblue"  ,"darkgreen")
plot(2000:2030, unscreened_cumul_haz[50:80]*inc, col=colors_treat[1],type="l", cex=1, lwd=2, lty=1, tck = -.02,ylab = "EAC Cumulative Incidence/100K", xlab= "Calendar Year",cex.lab = 1, cex.axis=1)
lines(2000:2030,amcumul_inc_50ablate[26:56]*inc,type="o", pch=16, cex=1,lty=1,col=colors_treat[2])
lines(2000:2030,amcumul_inc_99ablate2[26:56]*inc,type="o", pch=18, cex=1,lty=1,col=colors_treat[3])
lines(2000:2030,amcumul_inc_maligablate2[26:56]*inc,type="o", pch=17, cex=1,lty=1,col=colors_treat[4])
lines(2000:2030,amcumul_inc_hgdmaligablate2[26:56]*inc,type="o", pch=15, cex=1,lty=1,col=colors_treat[5])
legend('topleft', c("background","ablate 50% of all cell types", "ablate 99% of all cell types", "ablate 100% of malignant cells","ablate 100% of HGD"," & malignant cells"), cex=1, lwd=c(2,1,1,1,1),pch=c(NA,16,18,17,15), lty=c(1),box.lwd = 0,col=c(colors_treat,"white"))
