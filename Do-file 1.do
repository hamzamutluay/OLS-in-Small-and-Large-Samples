

cls
clear 		  all
se se 		  1990


se 	 	      obs 100
loc 		  beta = 0.5
loc 		  mux = 1
loc           sigmax = 2
loc 		  replication = 10000

mat 		  betahat  = J(`replication',1,.)


g 			  a = runiform(0,1)
g 			  z = (a-0.5)/( sqrt(1/12))
g 			  x = `mux' + `sigmax'*z




forv 		  t = 1/`replication' {


g 			  e = rnormal(0,1)
g 			  y = `beta'*x + e 
reg 		  y x

mat 		  betahat[`t', 1] = _b[x]

drop 		  e y 

}


sca	          summation = 0


forv 		  tt = 1/`replication' {

sca           summation = summation + betahat[`tt',1]

}

sca 	      mean_betahat = summation/`replication'

sca 		  bias = (mean_betahat -`beta')*100

di  		  bias



