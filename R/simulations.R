sim.data=function(sim, distr, n, seed=1) {
    
    set.seed(seed)
    
    sigmasq=.01 
    dil.x=5; dil.y=1
    dil.r=dil.x/dil.y
    theta=c(c=13,d=27560,b=-2.5,f=0.42)
    
    ###############################
    # generate r.x
    
    # generate r.half first and compute r.x from it
    if (distr=="unif1") {
        # the original simulation. non-random uniformly spaced between bounds that are inside c and d
        r.half=(seq(log(15),log(27550),length=n))        
        r.x=four_pl_prc(theta["c"],theta["d"],theta["b"],theta["f"],r.half,dil.r^(-0.5))
    } else if (distr=="rnorm1") {
        r.half=rnorm(n*5, (log(15)+log(27550))/2, sd=(log(27550)-log(14))/4)
        r.half=r.half[r.half>log(15) & r.half<log(27550)][1:n]
        if (length(r.half)!=n) stop("not enough data generated")        
        r.x=four_pl_prc(theta["c"],theta["d"],theta["b"],theta["f"],r.half,dil.r^(-0.5))
        
    # generate r.x directly
    } else if (distr=="unif2") {
        # non-random uniformly spaced between c and d
        r.x=seq(log(theta["c"]),log(theta["d"]),length=n+2)[2:(n+1)]
    } else if (distr=="runif") {
        r.x=sort(runif(n, log(theta["c"]), log(theta["d"])))
    } else if (startsWith(distr, "rbeta_")) {
        # beta distribution
        nu=as.numeric(strsplit(distr,"_")[[1]][-1])
        r.x=sort(rbeta(n,nu[1],nu[2])) * (log(theta["d"])-log(theta["c"])) + log(theta["c"])        
    } else if (startsWith(distr, "mix_")) {
        # in mix_90, the  number is percent responders
        perc=as.numeric(strsplit(distr,"_")[[1]][2])/100
        r.x=sort(c(
            rbeta(n*perc,4,2), # responders distr
            rbeta(n-n*perc,1,1.5) # nonresponders distr, strangely n*(1-perc) may give one less number
        )) * (log(theta["d"])-log(theta["c"])) + log(theta["c"])       
    
    } else {
        stop("distr not supported: "%+%distr)
    }
    
        
    ###############################
    # generate rest
    
    r.y=four_pl_prc(theta["c"],theta["d"],theta["b"],theta["f"],r.x,dil.r)
    xvar=r.x+rnorm(n,sd=sigmasq^.5)
    yvar=r.y+rnorm(n,sd=sigmasq^.5)        
               
    data.frame(xvar, yvar)
}
