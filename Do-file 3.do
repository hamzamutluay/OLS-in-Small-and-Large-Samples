





cls
clear 		all
se 		    mo off,perm
se se 		1990


loc 		beta0 = 10
loc 		beta1 = 1
loc 		beta2 = -2


loc 		mux1  = 0.5
loc 		sigmax1 = 1.5
loc 		mux2  = 0.6
loc 		sigmax2 = 1.8

loc 		replication = 2000
loc 		epsilon = 0.01

numlist 	"100(1000)20000"
loc			numlist `r(numlist)'
loc         numlist : subinstr loc numlist " " ",", all
mat			n_all = (`numlist')
loc         dim_col = colsof(n_all)
loc         dim_row = rowsof(n_all)


mat 		probability_batehat1  = J(`dim_col',2,.)
mat 		probability_batehat2  = J(`dim_col',2,.)



forv 	    t = 1/`dim_col' {

set 		obs `=n_all[`dim_row',`t']'
loc			obs = _N
	

g 			sigmasquare = runiform(1,3)

g 			a1 = runiform(0,1)
g 			a2 = runiform(0,1)

g 			z1 = (a1-0.5)/( sqrt(1/12))
g 			z2 = (a2-0.5)/( sqrt(1/12))

g 			x1 = `mux1' + `sigmax1'*z1
g 			x2 = `mux2' + `sigmax2'*z2

mat 		betahat1  = J(`replication',1,.)
mat 		betahat2  = J(`replication',1,.)


forv 		i = 1/`replication' {
	
g 		    u = .

forv 		tt = 1 / `obs' {
	
replace		u = rnormal(0,sigmasquare[`tt']) if `tt' == _n

}


g 			y = `beta0' + `beta1'*x1 + `beta2'*x2 + u 
reg 		y x1 x2

mat 		betahat1[`i', 1] = _b[x1]
mat 		betahat2[`i', 1] = _b[x2]

drop 		u y
}

sca	        summation1 = 0
sca	        summation2 = 0

forv        ii= 1/`replication' {
	
sca         summation1 = summation1 + ((abs(betahat1[`ii',1]-`beta1')) < `epsilon')
sca         summation2 = summation2 + ((abs(betahat2[`ii',1]-`beta2')) < `epsilon')

}


sca 	    mean_betahat1 = summation1/`replication'
sca 	    mean_betahat2 = summation2/`replication'


mat 		probability_batehat1[`t',1] = mean_betahat1
mat 		probability_batehat1[`t',2] = n_all[`dim_row',`t']	

mat 		probability_batehat2[`t',1] = mean_betahat2
mat 		probability_batehat2[`t',2] = n_all[`dim_row',`t']	


drop 		a1 a2 z1 z2 x1 x2 sigmasquare
mat         drop betahat1 betahat2 
sca drop    summation1 summation2 mean_betahat1 mean_betahat2

	
}


preserve
svmat       probability_batehat1, names(col)
scatter     c1 c2, jitter(7)
gr			save prob1, replace
restore

svmat       probability_batehat2, names(col)
scatter     c1 c2, jitter(7)
gr			save prob2, replace

