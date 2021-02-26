# -------------- Preparation; getting the model implied matrices

# load required R packages
library( "OpenMx" )
library( "psychonetrics" )
library( "qgraph" )
library( "dplyr" )
library( "MASS")

# clear work space, not run
rm( list = ls() )

# load some helper functions
source( "https://raw.githubusercontent.com/KJKan/nwsem/master/nwsem_helper_functions.R" )

# Load WAIS-IV correlation matrices
load( url ( "https://github.com/KJKan/mcfarland/blob/master/WAIS_US.Rdata?raw=true" ) )
load( url ( "https://github.com/KJKan/pame_I/blob/main/WAIS_Germany.Rdata?raw=true" ) )

# WAIS-IV sample Sizes
n_US      <- 1800 
n_Germany <- 1425 

# The WAIS-IV subtests; observed variables 
yvars     <- colnames( WAIS_US )

# Number of observed variables
ny        <- length( yvars ) 

# covariance matrix used in psychonetrics models
cov       <- ( n_Germany - 1 )/n_Germany*WAIS_Germany

# latent constructs to be measured (etas)
lvars     <- c( "P", # Perceptual
                "V", # Verbal
                "W", # Working Memory
                "S"  # Speed 
               )

# Number of latent constructs
ne        <- length( lvars ) 

# theoretical pattern of factor loadings (simple structure)
lambda_g <- matrix( c ( #P  V  W  S   
                        2, 0, 0, 0, # BD
                        0, 2, 0, 0, # SI
                        0, 0, 2, 0, # DS
                        1, 0, 0, 0, # MR
                        0, 1, 0, 0, # VC
                        0, 0, 1, 0, # AR 
                        0, 0, 0, 2, # SS
                        1, 0, 0, 0, # VP
                        0, 1, 0, 0, # IN
                        0, 0, 0, 1, # CD
                        0, 0, 1, 0, # LN
                        0, 0, 1, 0, # FW  
                        0, 1, 0, 0, # CO
                        0, 0, 0, 1, # CA
                        1, 0, 0, 0  # PC 
                        ),
                     ncol = ne, 
                     byrow = TRUE )

lambda_b <- cbind( lambda_g, g = 1 ) 
lambda_b[ ny, ne + 1 ] <- 2


# Extract a network from the US Sample
NWModel_US <- ggm( covs = ( n_US - 1 )/n_US*WAIS_US,
                   omega = "Full",
                   nobs = n_US )

# Prune it
NWModel_US <- NWModel_US %>% prune( alpha = 0.01, recursive = TRUE )

# Aim for further improvement
NWModel_US <- NWModel_US %>% stepup 

# Extract the adjacency matrix and use it as confirmatory network in the German sample
omega      <- 1*( getmatrix( NWModel_US, "omega" ) !=0 )


# ------------------------------ OpenMx models

# --- Higher order g factor model

# The data as OpenMx object
Data <- mxData( WAIS_Germany, type = "cov", numObs = n_Germany )

# Matrix containing the first order factor loadings 
Lambda_g <- mxMatrix( name = 'Lambda',
                      type = 'Full',
                      nrow = ny,
                      ncol = ne,
                      free = lambda_g==1,
                      values = lambda_g/2,
                      labels = label( 'lambda', ny, ne ) )

# Matrix containing second order factor loadings 
Gamma <- mxMatrix( name = 'Gamma',
                   type = 'Full',
                   nrow = ne,
                   ncol = 1,
                   free = TRUE,
                   values = 1,
                   labels = label( 'gamma', ne, 1 ) )

# Matrix containing the variance of the general factor
Phi <- mxMatrix( name = 'Phi',
                 type = 'Symm',
                 nrow = 1,
                 free = FALSE,
                 values = 1,
                 labels = 'phi' )

# Matrix containing the residual variances of the first order factors
Psi_g <- mxMatrix( name = 'Psi',
                   type = 'Diag',
                   nrow = ne,
                   free = TRUE,
                   values = diag(ne),
                   labels = label ( 'psi', ne ) )

# Matrix containing the residual variances of the observed variables
Theta <- mxMatrix( name = 'Theta',
                   type = 'Diag',
                   nrow = ny,
                   free = TRUE,
                   values = 1,
                   labels = label ( 'theta', ny ) )

# The factor model implied variance-covariance matrix
ExpCov_g <- mxAlgebra( name = 'ExpCov',
                       expression = Lambda %*% ( Gamma %*% Phi %*% t( Gamma ) + Psi ) %*% t( Lambda ) + Theta )

# ny by ny identity matrix
I <- mxMatrix( name = 'I',
               type = 'Iden',
               nrow = ny )

# Correlation matrix
ExpCor <- mxAlgebra( name = 'ExpCor',
                     expression = solve( sqrt( I * ExpCov ) ) %*% ExpCov %*% t( solve( sqrt( I * ExpCov ) ) ),
                     dimnames = list( yvars, yvars ) )

# Objective
Obj <- mxExpectationNormal( covariance = 'ExpCov',
                            dimnames = yvars )

# The complete model
gModel <-mxModel( name = 'g Model',
                  Data,
                  Lambda_g,
                  Gamma,
                  Phi,
                  Psi_g,
                  Theta,
                  ExpCov_g,
                  I,
                  ExpCor,
                  Obj,
                  mxFitFunctionML() )

# --- Bifactor model

# Matrix containing the first order factor loadings 
Lambda_b <- mxMatrix( name = 'Lambda',
                      type = 'Full',
                      nrow = ny,
                      ncol = ne + 1,
                      free = lambda_b==1,
                      values = lambda_b/2,
                      labels = label( 'lambda', ny, ne + 1 ) )

# Matrix containing the variances of the factors
Psi_b <- mxMatrix( name = 'Psi',
                   type = 'Diag',
                   nrow = ne + 1,
                   free = TRUE,
                   values = diag( ne + 1 ),
                   labels = label ( 'psi', ne + 1 ) )

# The factor model implied variance-covariance matrix
ExpCov_b <- mxAlgebra( name = 'ExpCov',
                       expression = Lambda %*% Psi %*% t( Lambda ) + Theta )

# The complete model
bModel <-mxModel( name = 'Bifactor model',
                  Data,
                  Lambda_b,
                  Psi_b,
                  Theta,
                  ExpCov_b,
                  I,
                  ExpCor,
                  Obj,
                  mxFitFunctionML() )


# --- Network model

# Matrix containing the scaling parameters
Delta <- mxMatrix( name = 'Delta',
                   type = 'Diag',
                   nrow = ny,
                   ncol = ny,
                   free = TRUE,
                   values = 1,
                   labels = label ('delta', ny ) )

# Matrix containing the partial relations (except the diagonal contains zeroes)
Omega <- mxMatrix( name = 'Omega',
                   type = 'Symm',
                   nrow = ny,
                   ncol = ny,
                   free = omega!=0,
                   values = 0,
                   labels = label( 'omega', ny ) )

# Expected partial correlation Matrix
ExpPcor <- mxAlgebra( name = 'ExpPcor',
                      expression = Omega + I )

# Expected Covariance matrix
ExpCovNW <- mxAlgebra( name = 'ExpCov',
                       expression = Delta %*% solve( I - Omega ) %*% t( Delta ) )

# The complete model
NWModel <- mxModel(name = 'Network Model',
                   Data,
                   I,
                   Delta,
                   Omega,
                   ExpCovNW,
                   ExpPcor,
                   ExpCor,
                   Obj,
                   mxFitFunctionML() )

# ------------------------------ Fit and compare the OpenMx models

# Fit the models
gFit  <- mxRun( gModel )
bFit  <- mxRun( bModel )
NWFit <- mxRun( NWModel )

# Print results
( gRes  <- summary( gFit ) )
( bRes  <- summary( bFit ) )
( NWRes <- summary( NWFit )  )

# Compare
mxCompare( bFit,  gFit )
mxCompare( NWFit, gFit )
mxCompare( NWFit, bFit )


# ------------------------------ Extract the model implied correlation matrices

# standardized sigmas
st_sigma_g  <- matrix( gFit@output$algebras$`g Model.ExpCor`,        ny, dimnames = list( yvars, yvars ) )
st_sigma_b  <- matrix( bFit@output$algebras$`Bifactor model.ExpCor`, ny, dimnames = list( yvars, yvars ) )
st_sigma_nw <- matrix( NWFit@output$algebras$`Network Model.ExpCor`, ny, dimnames = list( yvars, yvars ) )

st_sigmas <- list( g = st_sigma_g, b = st_sigma_b, nw = st_sigma_nw ) 
