## 3 morph evolutionary dynamics model for T. cristinae

## function for one generation evolutionary change
evolv<-function(W=NA,N=100,ff=NA){
	## W vector of fitnesses for GG, SS, MM, GS, GM, SM
	## N population size
	## ff lift of G = genotype frequencies: GG, SS, MM, GS, GM, SM
	## and p allele frequencies: G, S, M
	gen<-rmultinom(n=1,size=N,prob=ff$G*W)
	fG<-(gen[1]+.5*gen[4]+.5*gen[5])/N
	fS<-(gen[2]+.5*gen[4]+.5*gen[6])/N
	fM<-(gen[3]+.5*gen[5]+.5*gen[6])/N
	out<-list(G=c(fG^2,fS^2,fM^2,2*fG*fS,2*fG*fM,2*fS*fM),p=c(fG,fS,fM))
	return(out)
}

## function for NFDS for green and stripe
## returns green and stripe fitnesses based on frequencies
## and NFDS fitness function
nfds<-function(muw=0,betaw=0,mubar=0,betabar=0,nG=NA,nS=NA){
	## logit(wbar) = mubar + betabar * pstripe
	## wstripe/wbar = exp(muw + betaw * pstripe)
	## wgreen = 2wbar - wstripe
	## stripe frequency p = pstripe
	p<-nS/(nS+nG)

	## compute wbar
	lwbar<-mubar+betabar*p
	wbar<-1/(1+exp(-lwbar))

	## copmute wstripe
	wstripe<-wbar * exp(muw+ betaw * p)

	## compute wgreen
	wgreen<-2 * wbar - wstripe

	ww<-c(wgreen,wstripe)
	return(ww)
}

## posterior meands from nfdsfit.rdat
## mubar = -.87
## muw = 0.70 
## betabar = -0.26 
## betaw = -0.90
## consistent with replay paper
## nfds(muw = 0.7, betaw = -.9, mubar = -.87, betabar = -.26, nG=10, nS=10)



## main/control function
mainevolv<-function(p0=c(1/3,1/3,1/3),W=NA,Nmuw=0,betaw=0,mubar=0,betabar=0,gen=100){

	## p0 = initial frequency of G, S, and M
	## W = fixed or baseline fitness values
	## loop over generations
	for(i in 1:Ngen){

}


## test simulation
P<-matrix(NA,nrow=3,ncol=1000)
## drift
p<-1/3
gen<-rep(c(p^2,2*p*(1-p)),each=3)
X<-list(G=gen)

ww<-rep(1,6)

X<-evolv(W=ww,N=1000,ff=X)
P[,1]<-X$p

for(i in 2:1000){
	X<-evolv(W=ww,N=1000,ff=X)
	P[,i]<-X$p
}	

## overdominance
p<-1/3
gen<-rep(c(p^2,2*p*(1-p)),each=3)
X<-list(G=gen)

ww<-c(.5,.5,.5,1,1,1)

X<-evolv(W=ww,N=200,ff=X)
P[,1]<-X$p

for(i in 2:1000){
	X<-evolv(W=ww,N=200,ff=X)
	P[,i]<-X$p
}	
