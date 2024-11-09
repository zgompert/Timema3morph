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
## returns vector of w values based on genotype frequencies
## and NFDS fitness function
nfds<-function(W=NA,a=0.5,b=0,dom=1,ff=NA){
	## W is base matrix, use for melanic
	## a is interecept, b is slope
	## both on logit scale
	## dom dominance of stripe relative to green (0 or 1)

	## calcualte number of greens, stripes and melanics

	## compute fitnesses based on phenotype frequencies
	return(W)
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
