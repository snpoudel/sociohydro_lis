SocHydModel_SurgePrecip <- function(
        W,             #annual maxima discharge
        Precip,        #annual maxima rainfall
        Surgethreshold = 2,      #flooding threshold (storm surge)
        Rainthreshold = 0.1, #flooding threshold (precipitation)
        alphad = 1,    #risk taking attitude
        alphaa = 15,   #anxiousness
        alphap = 1,    #activeness
        alphar = 0.25, #effectiveness of preparedness
        Br = 1,        #flood to loss
        mewa = 0.10,   #forgetfullness
        mewp = 0.04,   #decay rate of precautionary measures
        Dmax = 1.0,    #maximum settlement density
        Pmax = 1.0,    #maximum preparedness
        Rmax = 1.0,    #maximum relative loss
        Amax = 1.0,    #maximum awareness
        U_rate = 0.5,  #growth rate
        POT_S_max = 3, #3.5 m
        POT_R_max = 0.15 #0.15 m
){
    #Initialize vectors
    POT = rep(0,length(W))          #Peak over Threshold
    POT_S = max(W) - Surgethreshold
    POT_R = max(Precip) - Rainthreshold
    L = rep(0,length(W))            #L - loss
    D = rep(0.1 * Dmax,length(W))     #D - population density
    R = rep(0,length(W))            #R - relative loss
    U = rep(U_rate,length(W))       #U - growth rate of settlement density
    A = rep(0.2*Amax,length(W))     #A - awareness
    P = rep(0.1*Pmax,length(W))     #P - preparedness
    dDdt = rep(0,length(W))
    dPdt = rep(0,length(W))
    dAdt = rep(0,length(W))
    
    
    for (t in 1:length(W))
    { 
        POT_S[t] = max(0,1*W[t] - Surgethreshold)
        POT_R[t] = max(0,1*Precip[t] - Rainthreshold)
    }
    
    POT_S_Scaled = POT_S/POT_S_max
    POT_R_Scaled = POT_R/POT_R_max
    
    
    for (t in 1:length(W))
    {
        POT[t] = max(POT_S_Scaled[t],POT_R_Scaled[t])
    }
    
    #sociohydro loop
    for (t in 2:length(W))
    {
        #POT[t] = max(POT) #this line was comment before
        #Relative loss
        if( POT[t] > 0)
        {R[t] = max(0,Rmax - Br*exp(-alphar*(Pmax - P[t-1])*POT[t]))}
        else
        {R[t] = 0}
        
        L[t] = R[t]*D[t-1]
        
        dDdt[t] = U[t]*(1 - alphad*A[t-1])*D[t-1]*(1 - D[t-1]/Dmax)
        
        dAdt[t] = alphaa*L[t]*(1 - A[t-1]/Amax) - mewa*A[t-1]
        
        if (R[t] > 0)
        {dPdt[t] = alphap*dAdt[t]*(1 - P[t-1]/Pmax) - mewp*P[t-1]}
        else
        {dPdt[t] = -1*mewp*P[t-1]}
        
        D[t] = max(0,min(D[t-1] + dDdt[t],Dmax))
        A[t] = max(0,min(A[t-1] + dAdt[t],Amax))
        P[t] = max(0,min(P[t-1] + dPdt[t],Pmax))
        
    }
    
    SocHydResults<-data.frame(W,L,D,R,U,A,P) 
    return(SocHydResults)
}

