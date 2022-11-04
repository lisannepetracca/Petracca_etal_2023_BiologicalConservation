#### COMMENTS SECTION ####
# put some stuff here about who you are, when you wrote this
# and why. what does it need to run

#### LOAD LIBRARIES #####
library(nimble)
library(here)
library(coda)

#### SOURCE FILES ####
# load data, nimble functions, etc

#### MODEL CODE ####
code1 <- nimbleCode({
  # Priors and constraints #####
  
  # use subsections in here for various processes
  
  # END priors and constraints

  
  # Likelihood #####
  
  # use subsections in here for various processes
  
  # END likelihood
  
  
  # Derived quantities #####
  
  # END derived quantities
  
})

#### DATA ####
dat1 <- list() 

#### CONSTANTS ####
const1 <- list()

#### INITIAL VALUES ####
inits1 <- list()

#### PARAMETERS TO MONITOR ####
params1 <- c()

#### MCMC SETTINGS ####
nb <- 5000 #burn-in
ni <- 10000 + nb #total iterations
nt <- 1  #thin
nc <- 3  #chains

#### COMPILE CONFIGURE AND BUILD ####
Rmodel1 <- nimbleModel(code = code1, constants = const1, data = dat1, 
                       check = FALSE, calculate = FALSE, inits = inits1)
conf1 <- configureMCMC(Rmodel2, monitors = params2, thin = nt, 
                       control = list(maxContractions = 1000)) 
# lots of initial model checking you can do by exploring conf1
# if you wanted to change samplers this is where you would do that
Rmcmc1 <- buildMCMC(conf1)  
Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel2)

#### RUN MCMC ####
t.start <- Sys.time()
out2 <- runMCMC(Cmcmc2, niter = ni , nburnin = nb , nchains = nc, inits = inits2,
                setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  
t.end <- Sys.time()
(runTime <- t.end - t.start)

#### MAKE BEAUTIFUL PLOTS AND STUFF ####
summary(out1)
traceplot(out1)