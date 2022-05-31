

cls
clear 		all
se mo 		off,perm
se se 		1990


loc 		beta = 1
loc 		mux  = 1
loc 		sigmax = 2
loc 		replication = 2000
loc 		epsilon = 0.01

numlist 	"100(500)20000"
loc			numlist `r(numlist)'
loc         numlist : subinstr loc numlist " " ",", all
mat			n_all = (`numlist')
loc         dim_col = colsof(n_all)
loc         dim_row = rowsof(n_all)



mat 		probability = J(`dim_col',2,.)

di  		`dim_row' ,`dim_col'



forv 		t = 1/`dim_col' {

set 		obs `=n_all[`dim_row',`t']'	

mat 		betahat = J(`replication',1,.)


g 			a = runiform(0,1)
g 			z = (a-0.5)/( sqrt(1/12))
g 			x = `mux' + `sigmax'*z

forv 		i = 1/`replication' {

g 			e = rnormal(0,1)
g 			y = `beta'*x + e 
reg 		y x
mat 		betahat[`i', 1] = _b[x]


drop 		e y
}

sca	        summation = 0

forv        tt= 1/`replication' {
	
sca         summation = summation + ((abs(betahat[`tt',1]-`beta')) < `epsilon')

}


sca 	    mean_betahat = summation/`replication'

mat 		probability[`t',1] = mean_betahat
mat 		probability[`t',2] = n_all[`dim_row',`t']	


drop 		a z x 
mat         drop betahat


	
}


mat li      probability
svmat       probability, names(col)
la var 		c1 "Probability"
la var 		c2 "Sample size"

scatter     c1 c2, jitter(7)


