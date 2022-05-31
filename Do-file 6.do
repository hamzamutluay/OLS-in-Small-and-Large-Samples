


cls
clear 		all
se mo 		off, perm
se se 		1990


numlist 	"100(1000)20000"
loc			numlist `r(numlist)'
loc         numlist : subinstr loc numlist " " ",", all
mat			n_all = (`numlist')
loc         dim_col = colsof(n_all)
loc         dim_row = rowsof(n_all)
loc  		col 1

mat 		probability_batehat  = J(`dim_col',2,.)


loc 		beta = 0.5
loc 		rho = 0.5
loc 		epsilon = 0.01
loc 		replication = 2000
loc 		mux = 1
loc 		sigmax = 2


forv 		t = 1/`dim_col' {

se 		    obs `=n_all[`dim_row',`t']'	
loc 		obs = _N

mat 		betahat = J(`replication',1,.)


g 			a = runiform(0,1)
g 			z = (a-0.5)/sqrt(1/12)
g 			x = `mux' + `sigmax'*z


forv 		i = 1/`replication' {
	
g 			e = rnormal(0,1)

g 			u = .
replace 	u = e[1] if _n == 1

forv 		tt = 2/`obs' {
	
replace     u = `rho'*u[`tt'-1] + e[`tt']	if _n == `tt'

}

g 			y = `beta'*x + u

reg 		y x 
mat 		betahat[`i',1] = _b[x]

drop 		y u e

}

sca	        summation = 0

forv        aa= 1/`replication' {
	
sca         summation = summation + ((abs(betahat[`aa',1]-`beta')) < `epsilon')

}

sca 	    mean_betahat = summation/`replication'

mat 		probability_batehat[`t',1] = mean_betahat
mat 		probability_batehat[`t',2] = n_all[`dim_row',`t']		

mat 		drop betahat
drop 		a z x
sca drop    summation mean_betahat
}


svmat       probability_batehat, names(col)
la var 		c1 "Probability"
la var 		c2 "Sample size"

scatter     c1 c2, jitter(7)
