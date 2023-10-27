
############# Implementation of the conditional permutation test ##################

generate_X_CPT_gaussian = function(nstep,M,X0,mu,sig2){
	# Runs the conditional permutation test using the distribution X | Z=Z[i] ~ N(mu[i],sig2[i])
	log_lik_mat = -(X0^2)%*%t(1/2/sig2) + X0%*%t(mu/sig2)
	# log_lik_mat[i,j] = density at X=X0[i] when Z=Z[j]
	Pi_mat = generate_X_CPT(nstep,M,log_lik_mat)
	X_mat = X0[Pi_mat]
	dim(X_mat) = c(M,length(X0))
	return(t(X_mat))
}

generate_X_CPT = function(nstep,M,log_lik_mat,Pi_init=NULL){
	# log_lik_mat is the n-by-n matrix with entries log(q(X_i|Z_j))
	# this function produces M exchangeable permutations, initialized with permutation Pi_init
	n = dim(log_lik_mat)[1]
	if(length(Pi_init)==0){
		Pi_init = 1:n
	}
	Pi_ = generate_X_CPT_MC(nstep,log_lik_mat,Pi_init)
	Pi_mat = matrix(0,M,n)
	for(m in 1:M){
		Pi_mat[m,] = generate_X_CPT_MC(nstep,log_lik_mat,Pi_)
	}
	return(Pi_mat)
}

generate_X_CPT_MC = function(nstep,log_lik_mat,Pi_){
	# log_lik_mat is the n-by-n matrix with entries log(q(X_i|Z_j))
	# this function runs the MC sampler, initialized with permutation Pi_
	n = length(Pi_)
	npair = floor(n/2)
	for(istep in 1:nstep){
	  perm = sample(n)
		inds_i = perm[1:npair]
		inds_j = perm[(npair+1):(2*npair)]
		# for each k=1,...,npair, deciding whether to swap Pi_[inds_i[k]] with Pi_[inds_j[k]]
    log_odds = (log_lik_mat[cbind(Pi_[inds_i],inds_j)] + log_lik_mat[cbind(Pi_[inds_j],inds_i)]
          - log_lik_mat[cbind(Pi_[inds_i],inds_i)] - log_lik_mat[cbind(Pi_[inds_j],inds_j)])
    swaps = rbinom(npair,1,1/(1+exp(-pmax(-500,log_odds))))
    Pi_[c(inds_i,inds_j)] = Pi_[c(inds_i,inds_j)] + swaps*(Pi_[c(inds_j,inds_i)]-Pi_[c(inds_i,inds_j)])
	}
	return(Pi_)
}

############# Function to generate X,Y,Z data ###################
generate_XYZ = function(type,param){
	a = rnorm(p)/p; b = rnorm(p)
	Z = matrix(rnorm(n*p),n,p)
	if(param==0 | type=='power'){
		X = Z%*%b + rnorm(n)
		Y = Z%*%a + X*param + rnorm(n)
	}else{
		if(type=='quadratic'){
			X = Z%*%b + param * (Z%*%b)^2 + rnorm(n)
		}
		if(type=='cubic'){
			X = Z%*%b - param * (Z%*%b)^3 + rnorm(n)
		}
		if(type=='tanh'){
			X = tanh(param*(Z%*%b))/param + rnorm(n)
		}
		if(type=='t'){
			X = Z%*%b + rt(n,1/param)*sqrt(1-2*param)
		}
		if(type=='skewnormal'){
			X = Z%*%b + (skewnormal(n,param)-sqrt(2/pi)*param/
			  sqrt(1+param^2))/sqrt(1-2/pi*param^2/(1+param^2))
		}
		if(type=='heteroskedastic'){
			X = Z%*%b +  rnorm(dim(Z)[1]) * abs(Z%*%b)^theta / 
             (sum(b^2))^(theta/2) / sqrt((2^theta)*gamma(theta+0.5)/sqrt(pi))
		}
		Y = Z%*%a + rnorm(n)
	}
	X_CPT = generate_X_CPT_gaussian(nstep,M,X,Z%*%b,rep(1,n))
	X_CRT = (Z%*%b)%*%t(rep(1,M)) + matrix(rnorm(n*M),n,M)
	XYZ = list()
	XYZ$X = X; XYZ$Y = Y; XYZ$Z = Z; XYZ$X_CPT = X_CPT; XYZ$X_CRT = X_CRT
	return(XYZ)
}

# rejection sampling from skew-normal distribution ( density 2*phi(t)*Phi(theta*t) )
skewnormal = function(n,theta){
	samples = NULL
	while(length(samples)<n){
		samples_ = rnorm(n)
		samples = c(samples,samples_[which(runif(n)<=pnorm(theta*samples_))])
	}
	return(samples[1:n])
}

############# Plotting function ###################
plot_mean_se_binary = function(x,y,col,pch){
  n = length(x); nrep = dim(y)[2]
  means = rowMeans(y)
  SEs = sqrt(means * (1-means) / nrep)
  for(i in 1:n){
    points(x[i],means[i],col=col,pch=pch)
    segments(x[i],means[i]-SEs[i],x[i],means[i]+SEs[i],col=col)
  }
}

############# Simulations under the null (testing robustness) ###################

set.seed(12345)
alpha = 0.05; nrep = 1000; n = 50; p = 20; M = 500; nstep = 50

types = c('quadratic','cubic','tanh','t','skewnormal','heteroskedastic')
theta_max = c(0.15,0.025,1,0.3,10,0.2); ntheta = 11
ymax = c(0.35,0.35,0.35,0.1,0.1,0.1)

for(itype in 1:6){
	pvals = array(0,c(ntheta,nrep,2))
	thetas = theta_max[itype] * (0:(ntheta-1))/(ntheta-1)
	for(itheta in 1:ntheta){
		theta = thetas[itheta]
		for(irep in 1:nrep){
			XYZ = generate_XYZ(types[itype],theta)
			X = XYZ$X; Y = XYZ$Y; Z = XYZ$Z; X_CPT = XYZ$X_CPT; X_CRT = XYZ$X_CRT
			T_true = c(abs(cor(X,Y)))
			T_CPT = c(abs(cor(X_CPT,Y)))
			pvals[itheta,irep,1] = (1 + sum(T_CPT >= T_true)) / (1+M)
			T_CRT = c(abs(cor(X_CRT,Y)))
			pvals[itheta,irep,2] = (1 + sum(T_CRT >= T_true)) / (1+M)
		}
	}
	filename = paste0('plot_results_',types[itype],'.pdf')
	pdf(filename,width=3.5,height=4)
	if(itype %% 3 == 1){ylab='Prob. of rejection'}else{ylab=''}
	plot(thetas,thetas,type='n',ylim=c(0,ymax[itype]),xlab=expression(theta),ylab=ylab,las=1)
	jitter = 0.03 * (thetas[2]-thetas[1])
	segments(-max(thetas),alpha,2*max(thetas),alpha,lty=2)
	plot_mean_se_binary(thetas,pvals[,,1]<=alpha,'black',18)
 	plot_mean_se_binary(thetas+jitter,pvals[,,2]<=alpha,'red',20)
	if(itype %% 3 == 1){legend('topleft',legend=c('CPT','CRT'),col=c('black','red'),pch=c(18,20))}
	dev.off()
}

############# Simulations under the alternative (testing power) #################

cs = (0:10)*0.075
pvals = array(0,c(length(cs),nrep,2))
for(ic in 1:length(cs)){
	c = cs[ic]
	for(irep in 1:nrep){
		XYZ = generate_XYZ('power',c)
		X = XYZ$X; Y = XYZ$Y; Z = XYZ$Z; X_CPT = XYZ$X_CPT; X_CRT = XYZ$X_CRT
		T_true = c(abs(cor(X,Y)))
		T_CPT = c(abs(cor(X_CPT,Y)))
		pvals[ic,irep,1] = (1 + sum(T_CPT >= T_true)) / (1+M)
		T_CRT = c(abs(cor(X_CRT,Y)))
		pvals[ic,irep,2] = (1 + sum(T_CRT >= T_true)) / (1+M)
	}
}		
filename = 'plot_results_power.pdf'
pdf(filename,width=3.5,height=4)
plot(cs,cs,type='n',ylim=c(0,1),xlab='c',ylab='Prob. of rejection',las=1)
jitter = 0.03 * (cs[2]-cs[1])
plot_mean_se_binary(cs,pvals[,,1]<=alpha,'black',18)
plot_mean_se_binary(cs+jitter,pvals[,,2]<=alpha,'red',20)
legend('topleft',legend=c('CPT','CRT'),col=c('black','red'),pch=c(18,20))
dev.off()
  
