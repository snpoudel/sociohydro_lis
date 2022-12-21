setwd("C:/SANDEEP/SocioHydro_Model/Data")

library(dataRetrieval)
library(Metrics)
library(corrplot)
library(tidyverse)
library(Hmisc)

source("SHModel_Housing.R")


# Returns numIter length list of entries to be peturbed
probPeturb<-function(x, numIter)
{
    # Input is xBounds & numIter.  
    # Returns numIter entry list with the indices which will be peturbed
    xDims<-nrow(x)
    probabilityVector<-1-log(1:numIter)/log(numIter)
    peturbIdx<-apply(matrix(unlist(lapply(probabilityVector, function(x) as.logical(rbinom(xDims, 1, x)))), byrow=TRUE, ncol=xDims), 1, which)
    return(peturbIdx)
}

#Read in index
fips = read.table("FIPS.txt")
Battery = read.csv("Max_Water_Level_Battery.csv")
precip = read.table("Max_precip_LIS.csv")

#DDS parameters
r= 0.3 #Default value of r(0.2) was finding the local optimum solution( like 2 out of 10 times), so r is increased.
numIter =100000

#running code for all census tracts of interest
n_sim = 1
Coeff_SH_mat <- data.frame(Doubles=double())
Results_SH_mat <- data.frame(Doubles=double())
RMSE_best_mat <- data.frame(Doubles=double())
missing_tract <- c(34017980100,34025990000,36005031900,36047001800,36047015400,36047990100,36061000100,36061000500,36061029700,36061031900,36081091602,36081990100,36085015400,36119982000,36085001800,36085022800,36047016400,36061011700,36119000104,36047035200,36119984000,36081099900,36119005600,9009140200,36005027600,36005050400)
#these are the census tracts with missing values on housing information

for(sim in 1:n_sim)
{
for (ct in 305:305) #length(fips$V1))  #436
{
    #Read in index
    #Census tracts with missing values are skipped first
    if(any(fips$V1[ct] == missing_tract))
    {next}
    
    #Reading of fema claims is divided for HRE(295 census tracts) and LIS as they have different filename and nrows
    else if(ct < 296)
    {fema_claims = read.table(paste("C:/SANDEEP/SocioHydro_Model/Data/FEMA_Claims/HRE_and_LIS/FEMA_Claims_",fips$V1[ct],".csv",sep = ''),sep = ',',skip = 1)
    fema_policies = read.table(paste("C:/SANDEEP/SocioHydro_Model/Data/FEMA_Policies/HRE_and_LIS/HRE_FEMA_Policies_",fips$V1[ct],".csv",sep = ''),sep = ',',skip = 1)}
    
    else
    {fema_claims = read.table(paste("C:/SANDEEP/SocioHydro_Model/Data/FEMA_Claims/HRE_and_LIS/FEMAClaims_",fips$V1[ct],".csv",sep = ''),sep = ',',skip = 2)
    fema_policies = read.table(paste("C:/SANDEEP/SocioHydro_Model/Data/FEMA_Policies/HRE_and_LIS/FEMA_Policies_",fips$V1[ct],".csv",sep = ''),sep = ',',skip = 1)}
    #Policy is missing one census tract(34017980100), a decoy value is added now should be removed later
    
    nhouses <- read.table(paste("C:/SANDEEP/SocioHydro_Model/Data/Housing_Number/nhousing_",fips$V1[ct],".csv",sep = ''),sep=',',skip = 1)
    vhouses <- read.table(paste("C:/SANDEEP/SocioHydro_Model/Data/Housing_Value/vhousing_",fips$V1[ct],".csv",sep = ''),sep=',',skip = 1)
    #missing values of nhouses and vhouses are assigned NA
    nhouses$V3[nhouses$V3==0 | nhouses$V3=="null" | nhouses$V3=="-"] <- NA
    vhouses$V3[vhouses$V3==0 | vhouses$V3=="null" | vhouses$V3=="-"] <- NA
    
    Claims_Rescaled_v = fema_claims$V4/mean(as.numeric(vhouses$V3), na.rm = T)
    Claims_Rescaled_n = fema_claims$V3/mean(as.numeric(nhouses$V3), na.rm = T)
    Claims_Rescaled_d = 0.1 * nhouses$V3/max(nhouses$V3) #Based on literature about 10 percent of total land in ct is used for residential
    Policy_Rescaled_a = fema_policies$V3/mean(nhouses$V3, na.rm = T)
    
    
    #Set calibration parameter boundaries
    xBounds.df = data.frame(matrix(ncol=2,nrow=12))
    colnames(xBounds.df)<-c("min", "max")
    
    #Surge threshold
    xBounds.df$min[1] = 0.001 #what about putting the min value of surge in lower bound(1)?
    xBounds.df$max[1] = 3.5 #wmax 
    
    #Rain threshold
    xBounds.df$min[2] = 0.001 #what about putting the min value of Rain in lower bound(0.03)?
    xBounds.df$max[2] = 0.15 #precipmax
    
    #alphad
    xBounds.df$min[3] = 0.001 
    xBounds.df$max[3] = 25 #constrained based on barplot from 100
    
    #alphaa
    xBounds.df$min[4] = 0.001
    xBounds.df$max[4] = 25
    
    #alphap
    xBounds.df$min[5] = 0.001
    xBounds.df$max[5] = 10 #increased to 10 from 3 based on sh output
    
    #alphar
    xBounds.df$min[6] = 0.001
    xBounds.df$max[6] = 1 #constrained based on barplot from 3
    
    #mewa
    xBounds.df$min[7] = 0.001
    xBounds.df$max[7] = 1 #constrained based on barplot from 3
    
    #mewp
    xBounds.df$min[8] = 0.001
    xBounds.df$max[8] = 2 #constrained based on barplot from 3
    
    #U_rate
    xBounds.df$min[9] = 0.01
    xBounds.df$max[9] = 3 #constrained based on barplot from 5
    
    #POT_S_max
    xBounds.df$min[10] = 0.001 
    xBounds.df$max[10] = 3.5 #max w - min threshold w
    
    #POT_R_max
    xBounds.df$min[11] = 0.001
    xBounds.df$max[11] = 0.15
   
    #phi
    xBounds.df$min[12] = 0.7 
    xBounds.df$max[12] = 1 # Result are good when the phi ranges from 0.7 to 1 (But we can reach to better range if we can calibrate housing price) 
    
    

    # Generate initial first guess
    x_init<-c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,0.1)
    x_best = data.frame(x_init)
    RMSE_best <- 1000000
    
    peturbIdx<-probPeturb(xBounds.df, numIter)
    # Peturb each entry by N(0,1)*r(x_max - x_min) reflecting if beyond boundaries
    sigma<-xBounds.df$max - xBounds.df$min
    
    for (i in 1:numIter)
    {
        # Set up test parameter values as x_test
        x_test<-as.matrix(x_best)
        
        # Get entries we will peturb
        idx<-peturbIdx[[i]]
        if (sum(idx) == 0) {idx = round(runif(1,1,12))}
        
        # Initialize vector of peturbations initially zeros with same length of x so we will add this vector to peturb x
        peturbVec<-rep(0, length(x_test))
        # Generate the required number of random normal variables
        N<-rnorm(length(x_test), mean=0, sd=1)
        
        # Set up vector of peturbations
        peturbVec[idx]<-r*N[idx]*sigma[idx]
        
        # Temporary resulting x value if we peturbed it
        testPeturb<-x_test + peturbVec  
        # Find the values in testPeturb that have boundary violations.  Store the indices in boundaryViolationsIdx
        boundaryViolationIdx<-which(testPeturb<xBounds.df$min | testPeturb > xBounds.df$max)
        
        # Reset those violated indices to the opposite peturbation direction
        peturbVec[boundaryViolationIdx]<-(-1*r*N[boundaryViolationIdx]*sigma[boundaryViolationIdx])
        
        # Find values still at violations of min or max and set them to the minimum or maximum values
        x_test<-x_test + peturbVec
        minViolationIdx<-which(x_test<xBounds.df$min)
        maxViolationIdx<-which(x_test>xBounds.df$max)
        x_test[minViolationIdx]<-xBounds.df$min[minViolationIdx]
        x_test[maxViolationIdx]<-xBounds.df$max[maxViolationIdx]
        
        #Socio-hydrological model    
        Results = SocHydModel_SurgePrecipHousing(W = Battery$Max_Wl_Battery[10:52], Precip = precip$precip_max[10:52], x_test[1], x_test[2], x_test[3], x_test[4], x_test[5], x_test[6], 1, x_test[7], x_test[8], 
                                          1, 1, 1, 1, x_test[9], x_test[10], x_test[11], x_test[12], 1 )
        
        #Objective Function 
        #Rmses are in different scale, so they are normalized by dividing with the mean of observed values.
        #For some census tracts there is no claims and/or policy, so dividing by mean to normalize rmse would cause nan error.
        if(mean(Claims_Rescaled_v) == 0){
            RMSE_Test_V = rmse(Claims_Rescaled_v , Results$L[(44-length(Claims_Rescaled_v)):43] )
        } 
        else {
            
            RMSE_Test_V = rmse(Claims_Rescaled_v , Results$L[(44-length(Claims_Rescaled_v)):43] ) / mean(Claims_Rescaled_v)
        }
        if(mean(Policy_Rescaled_a) == 0){
            
            RMSE_Test_A = rmse(Policy_Rescaled_a , Results$A[31:43] ) / mean(Policy_Rescaled_a)
        }
        else {
            
            RMSE_Test_A = rmse(Policy_Rescaled_a , Results$A[31:43] ) / mean(Policy_Rescaled_a)
        }
        RMSE_Test_D = rmse(Claims_Rescaled_d , Results$D[33:(32 + length(Claims_Rescaled_d))] ) / mean(Claims_Rescaled_d)
        
        
        RMSE_Test <- sqrt( RMSE_Test_V^2 + RMSE_Test_A^2 + RMSE_Test_D^2)
        #read about multiobjective optimization that could help solve the problem
        #Check if this simulation is better
        if (RMSE_Test < RMSE_best) 
        {
            x_best = x_test
            RMSE_best = RMSE_Test 
            
        }
        
        #Print to console
        print_str = paste(ct, "Tract:", fips$V1[ct],"RMSE BEST",i,":",RMSE_best)
        print(print_str)
        
    }
    ###write results
    ##Result of SH coefficients for all census tracts
    x_bestt=t(x_best)
    colnames(x_bestt) <- c("Surge_threshold","Rain_threshold","alpha_d","alpha_a","alpha_p","alpha_r","mew_a","mew_p","U_rate","POT_S_max","POT_R_max", "phi")
    x_bestt <- cbind(census_tract = fips$V1[ct] , x_bestt)
    
    
    ##Result of SH results for all census tracts
    
    Results_SH <- cbind(census_tract = fips$V1[ct], Results)
    
    
    ##Results of RMSE for all census tracts
    RMSE_best <- cbind( census_tract = fips$V1[ct], RMSE_best)
    
    
}
    Coeff_SH_mat <- rbind(Coeff_SH_mat,x_bestt)
    #write.table(Coeff_SH_mat,paste("C:/SANDEEP/SocioHydro_Model/Results/single_tract/SH_Parameters.csv"),sep=",")
    
    Results_SH_mat <-  rbind(Results_SH_mat,Results_SH)
    #write.table(Results_SH_mat,paste("C:/SANDEEP/SocioHydro_Model/Results/single_tract/SH_Results.csv"),sep=",")
    
    RMSE_best_mat <- rbind(RMSE_best_mat,RMSE_best)
    #write.table(RMSE_best_mat,paste("C:/SANDEEP/SocioHydro_Model/Results/single_tract/SH_RMSE.csv"),sep=",")
}
#corrplot(cor(Coeff_SH_mat), method = "number")

# hist.data.frame(Coeff_SH_mat,nclass=5)

##Plot
par(mfrow=c(3,3))

#storm surge
plot(Battery$Year[10:52] ,Battery$Max_Wl_Battery[10:52], type = "o")


#awareness 
plot(Battery$Year[10:52] ,Results$A , type = "o", ylim = c(0,max(Policy_Rescaled_a)))
lines(Battery$Year[40:52],Policy_Rescaled_a, col = "red")


#precipitation
plot(precip$yr[10:52] ,precip$precip_max[10:52], type = "o" )

#preparedness
plot(Battery$Year[10:52] ,Results$P, type = "o",ylim = c(0,1))

#loss
plot(Battery$Year[10:52] ,Results$L, type = "o", ylim = c(0,max(Claims_Rescaled_v)))
lines(fema_claims$V2 , Claims_Rescaled_v, col = "red")

#density
plot(Battery$Year[10:52] ,Results$D, type = "o" , ylim = c(0,max(Claims_Rescaled_d)))
lines(Battery$Year[42:52], Claims_Rescaled_d, col = "red")

#housing price
plot(Battery$Year[10:52] ,Results$HP, type = "o",ylim = c(0,1) )
