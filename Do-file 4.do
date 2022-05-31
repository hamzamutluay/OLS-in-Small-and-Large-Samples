

cls
clear 		all
se mo 		off, perm
se se 		1990

se obs 		100
loc 		obs = _N


loc 		alpha = 10
loc 		beta = 1
loc 		theta = -2

sca 		mux1 = 0.5
sca			sigmax1 = 1.5

sca			mux2 = 0.6
sca			sigmax2 = 1.8

loc			replication = 10000



g 			a1 = runiform(0,1)
g 			a2 = runiform(0,1)

g 			z1 = (a1-0.5)/(sqrt(1/12))
g 			z2 = (a2-0.5)/(sqrt(1/12))

g 			x1 = mux1 + sigmax1*z1
g 			x2 = mux2 + sigmax2*z2


mat 		betahat = J(`replication',1,.)
mat 		thetahat = J(`replication',1,.)
g 			sigmasquare = runiform(1,3)


forv 		t = 1 / `replication' {

g 		    epsilon = .

forv 		tt = 1 / `obs' {
	
replace		epsilon = rnormal(0,sigmasquare[`tt']) if `tt' == _n

}

g 			y = `alpha' + `beta'*x1 + `theta'*x2 + epsilon
reg 		y x1 x2 

mat 		betahat[`t',1] = _b[x1]
mat 		thetahat[`t',1] = _b[x2]

drop 		epsilon y
					
}

sca 		summation1 = 0	
sca 		summation2 = 0		
	

forv 		i = 1 / `replication' {
	
sca 		summation1 = 	summation1 + betahat[`i',1]
sca 		summation2 = 	summation2 + thetahat[`i',1]

}

sca 	    mean_betahat = summation1/`replication'
sca 	    mean_thetahat = summation2/`replication'


sca 		bias_betahat = (mean_betahat -`beta')*100
sca  		bias_thetahat = (mean_thetahat -`theta')*100


di 			bias_betahat
di 			bias_thetahat







