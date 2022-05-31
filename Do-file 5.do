


cls
clear 		all
se mo       off, perm
se se       1990


se 		    obs 100
loc 		obs = _N
loc 		beta = 0.5
sca 		mux = 1
sca         sigmax = 2
loc 		rho = 0.8
loc 		replication = 10000

g 		    a = runiform(0,1)
g 			z = (a-0.5)/sqrt(1/12)
g           x = mux + sigmax*z

mat         betahat = J(`replication',1,.)

forv 		t = 1/`replication' {
	
	
g 			u = .	
g 		    epsilon = rnormal(0,1)
replace 	u = epsilon[1]  if _n == 1

forv   	    j = 2/`obs' {
	
replace 	u = `rho'*u[`j'-1] + epsilon[`j']  if _n == `j'
	
}

g 			y = `beta'*x + u

reg 		y x

mat 		betahat[`t',1] = _b[x]

drop		y  u epsilon
	
}

sca 		summation = 0

forv 		i = 1/`replication' {

sca 		summation = summation + betahat[`i',1]	
	
}

sca 	    mean_betahat = summation/`replication'

sca 		bias = (mean_betahat -`beta')*100

di  		bias
