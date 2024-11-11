## 3 morph evolutionary dynamics model for T. cristinae

## function for one generation evolutionary change
evolv<-function(W=NA,N=100,ff=NA){
	## W vector of fitnesses for GG, SS, MM, GS, GM, SM
	## N population size
	## ff list of G = genotype frequencies: GG, SS, MM, GS, GM, SM
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
mainevolv<-function(p0=c(1/3,1/3,1/3),Wb=NA,N=200,muw=0,betaw=0,mubar=0,betabar=0,gen=100,pen=1,msd=0.01,odm=0.05,
		    fds=FALSE,fluct=FALSE){

	## p0 = initial frequency of G, S, and M
	## Wb = fixed or baseline fitness values for G, S, M
	## N = population size
	## logit(wbar) = mubar + betabar * pstripe
	## wstripe/wbar = exp(muw + betaw * pstripe)
	## wgreen = 2wbar - wstripe
	## pen = stripe penetrance
	## msd = sd on logit scale for melanic fitness
	## odm = overdominance increase on logit scale
	## loop over generations
	PM<-matrix(NA,nrow=gen+1,ncol=3)
	PM[1,]<-p0
	p<-p0
	for(i in 1:gen){
		## W vector of fitnesses for GG, SS, MM, GS, GM, SM
		W<-c(Wb,(Wb[2]*pen + Wb[1]*(1-pen)),Wb[1],Wb[2])
		Ngen<-N*c(p^2,2*p[1]*p[2],2*p[1]*p[3],2*p[2]*p[3]) 
		## determine stripe and green fitness from NFDS
		if(fds){
			nG<-Ngen[1]+(1-pen)*Ngen[4]+Ngen[5]
			nS<-Ngen[2]+pen*Ngen[4]+Ngen[6]
			wfds<-nfds(muw,betaw,mubar,betabar,nG,nS)
			W[c(1:2,4:6)]<-c(wfds,(wfds[2]*pen + wfds[1]*(1-pen)),wfds)
		}
		## add overdominance
		W[5]<-1/(1+exp(-1*(odm+log(W[5]/(1-W[5])))))	
		W[6]<-1/(1+exp(-1*(odm+log(W[6]/(1-W[6])))))
		## add environmental variation
		if(fluct){
			a<-rnorm(1,0,msd)
			W[3]<-1/(1+exp(-1*(a+log(W[3]/(1-W[3])))))
		}	
		ff<-list(G=Ngen/N,p=p)
		ff<-evolv(W=W,N=N,ff=ff)
		PM[i+1,]<-ff$p
		p<-ff$p
		if(max(p)==1){
			return(PM)
		}
	}
	return(PM)
}

o<- mainevolv(Wb=rep(.5,3),N=200,muw=.7,betaw=-.9,mubar=-.87,betabar=-.26,gen=1000,fds=TRUE)
o<- mainevolv(Wb=rep(.3,3),N=200,muw=.7,betaw=-.9,mubar=-.87,betabar=-.26,gen=1000,odm=.05,fds=TRUE)
o<- mainevolv(Wb=rep(.2,3),N=200,muw=.7,betaw=-.9,mubar=-.87,betabar=-.26,gen=1000,odm=.07,fds=TRUE)

pdf("t3example.pdf",width=5,height=5)
par(mar=c(5,5,.5,.5))
plot(o[,1],type='l',ylim=c(0,1),col="green",xlab="Generation",ylab="Frequency",cex.lab=1.4,lwd=1.5)
lines(o[,2],col="cadetblue",lwd=1.5)
lines(o[,3],col="brown",lwd=1.5)
dev.off()

