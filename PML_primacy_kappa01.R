library(readxl)
library(mlogit)
library(R.matlab)
library(MASS)
set.seed(1234)

# Effects coding function.
frow2xrow <- function(flevelvalues, nlevels){
  nf <- nrow(nlevels)
  df <- (sum(nlevels) - nf)
  xrow <- matrix(0, nrow=1, ncol=df)
  startidx <- 1
  for (i in 1:nf) {
    nl <- nlevels[i]
    xtmp <-  matrix(0, nrow = 1, ncol = nl-1)
    if(flevelvalues[i] < nl){
      xtmp[flevelvalues[i]] <- 1
    } else {
      for (k in 1:nl-1) {
        xtmp[k] = -1
      }
    }
    xrow[startidx:(startidx+nl-2)] <- xtmp
    startidx <- startidx + nl -1
  }
  xrow
}

## Simulation parameters
ndraws <- 500
cset <- 15
nalts <- 3
nlevels <- c(6, 3,4, 2,3)
nattributes <- 5
mu <- c(0.2141, -0.0135, -0.0067, 0.0395, -0.1169, -0.0180, -0.0247, -0.0308, 0.0174, 0.0389, 0.0417, 0.3, 0)
sigma <- c(0.2264,0.0908,0.0030,0.2226,0.3111,0.1339,0.0180,0.1251,0.1828,0.1914,0.2668,0.1,0.1)
nparam <- length(mu)
nresps <- 150
cov_matrix <- diag(sigma^2)  
beta_m <- matrix(0, nrow = nresps, ncol = length(mu))
for (i in 1:nresps) {
  beta_m[i, ] <- mvrnorm(1, mu = mu, Sigma = cov_matrix)
}


## designs
#naive design
X_data <- readMat(".../X_naive.mat")
X <- X_data$X.naive
cs_id <- rep(1:cset, each = nalts)
ff <- cbind(cs_id, X)
f <- ff[rep(1:nrow(ff), nresps), ]
nruns <- nrow(ff)
nsets <- nruns / nalts
nr <- nresps * nruns
nrcs <- nr / nalts

Xtotal <- matrix(0, nrow = nr, ncol = nparam)
for (i in 1:nr) {
  Xtotal[i,] <- frow2xrow(f[i,2:6], matrix(nlevels, nrow=nattributes, ncol=1))
}

utility <- matrix(0, nrow = nr, ncol = 1)
exputility <- matrix(0, nrow = nr, ncol = 1)
for(i in 1:nresps) {
  X <- Xtotal[((i-1)*nruns+1):(i*nruns),]
  utility[((i-1)*nruns + 1):(i*nruns)] <- X%*%beta_m[i,]
}

exputility <- exp(utility)


beta_PML_naive_primacy <- data.frame()
sd_PML_naive_primacy <- data.frame()
b <- 0

while (b < ndraws) {
  set.seed(2024 + b)
  start_time <- Sys.time()
  estimation_success <- FALSE
  
  while (!estimation_success) {
    # Generate choices
    r <- matrix(0, nrow = nr, ncol = 1)
    for (i in 1:nrcs) {
      u1 <- exputility[3 * i - 2]
      u2 <- exputility[3 * i - 1]
      u3 <- exputility[3 * i]
      p1 <- u1 / (u1 + u2 + u3)
      p2 <- u2 / (u1 + u2 + u3)
      sump <- p1 + p2
      rnd <- runif(1)
      if (p1 >= rnd) {
        r[3 * i - 2, 1] <- 1
      } else if (sump >= rnd) {
        r[3 * i - 1, 1] <- 1
      } else {
        r[3 * i, 1] <- 1
      }
    }
    
    # Create data frame 'm' with choices
    m <- data.frame(
      subject_id = rep(1:nresps, each = nruns),
      chosen = r,
      choice_set_id = rep(ff[,1], nresps),
      X1 = f[, 2],
      X2 = f[, 3],
      X3 = f[, 4],
      X4 = f[, 5],
      X5 = f[, 6]
    )
    
    # Effect coding
    choices_effect_coded <- matrix(0, nrow = nrow(m), ncol = nparam)
    m_matrix <- as.matrix(m[, 4:8])
    for (i in 1:nrow(m)) {
      choices_effect_coded[i, ] <- frow2xrow(m_matrix[i, ], matrix(nlevels, nrow = nattributes, ncol = 1))
    }
    
    # Prepare choices for mlogit
    choices <- matrix(0, nrow = nsets * nresps, ncol = nparam * nalts)
    for (i in 1:nrow(choices)) {
      choices[i, 1:13] <- choices_effect_coded[i * 3 - 2, ]
      choices[i, 14:26] <- choices_effect_coded[i * 3 - 1, ]
      choices[i, 27:39] <- choices_effect_coded[i * 3, ]
    }
    
    # Indicate chosen profile
    profile_choice <- rep(NA, nsets * nresps)
    for (i in 1:length(profile_choice)) {
      if (m[i * 3 - 2, "chosen"] == 1) {
        profile_choice[i] <- "A"
      } else if (m[i * 3 - 1, "chosen"] == 1) {
        profile_choice[i] <- "B"
      } else if (m[i * 3, "chosen"] == 1) {
        profile_choice[i] <- "C"
      }
    }
    
    # Create subject IDs and choice set numbers
    id_and_choice_set <- data.frame(
      subject_id = rep(1:nresps, each = nsets),
      choice_set_id = rep(1:nsets, times = nresps)
    )
    
    info <- cbind(id_and_choice_set, profile_choice)
    choices_df <- cbind(info, as.data.frame(choices))
    colnames(choices_df) <- c(
      "subject_id", "choice_set_id", "chosen",
      paste0("V", 1:13, "_A"), paste0("V", 1:13, "_B"), paste0("V", 1:13, "_C")
    )
    
    # Convert to mlogit.data
    data_z <- mlogit.data(choices_df, choice = "chosen", shape = "wide", varying = 4:42, sep = "_", id.var = "subject_id")
    
    # Try model estimation
    tryCatch({
      fit_PML_with_order_naive_primacy <- mlogit(chosen ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13|0, 
                                                 data = data_z,
                                                 rpar = c(V1 = "n", V2 = "n", V3 = "n", V4 = "n", V5 = "n", V6 = "n", V7 = "n", V8 = "n", V9 = "n", V10 = "n", V11 = "n", V12 = "n", V13 = "n"), 
                                                 panel = TRUE, R = 1000,halton = list(dim = 13))
      
      beta_PML_naive_primacy <- rbind(beta_PML_naive_primacy, as.data.frame(t(fit_PML_with_order_naive_primacy[["coefficients"]][1:13])))
      sd_PML_naive_primacy <- rbind(sd_PML_naive_primacy, as.data.frame(t(fit_PML_with_order_naive_primacy[["coefficients"]][14:26])))
      estimation_success <- TRUE
      duration <- Sys.time() - start_time
      cat('draw:', b+1, 'duration:', duration, "\n")
    }, error = function(e) {
      
    })
  }
  
  b <- b + 1
}


# Random order
b <- 0
# Generate random order permutations
generate_WSPrmutations <- function() {
  unlist(replicate(cset, sample(1:nalts), simplify = FALSE))
}
random_matrix <- matrix(unlist(replicate(nresps, generate_WSPrmutations(), simplify = FALSE)), nrow = cset * nalts, ncol = nresps)

f_random <- cbind(ff[, 1:5], order = random_matrix[, 1])
for (i in 2:nresps) {
  sf <- cbind(ff[, 1:5], order = random_matrix[, i])
  f_random <- rbind(f_random, sf)
}

Xtotal_random <- matrix(0, nrow = nr, ncol = nparam)
for (i in 1:nr) {
  Xtotal_random[i, ] <- frow2xrow(f_random[i, 2:6], matrix(nlevels, nrow = nattributes, ncol = 1))
}

utility_random <- matrix(0, nrow = nr, ncol = 1)
exputility_random <- matrix(0, nrow = nr, ncol = 1)

for(i in 1:nresps) {
  X <- Xtotal_random[((i-1)*nruns+1):(i*nruns),]
  utility_random[((i-1)*nruns + 1):(i*nruns)] <- X%*%beta_m[i,]
}
exputility_random <- exp(utility_random)

beta_PML_random_primacy <- data.frame()
sd_PML_random_primacy <- data.frame()

while (b < ndraws) {
  set.seed(2024 + b)
  start_time <- Sys.time()
  estimation_success <- FALSE
  
  while (!estimation_success) {
    
    # Generate choices
    r <- matrix(0, nrow = nr, ncol = 1)
    for (i in 1:nrcs) {
      u1 <- exputility_random[3 * i - 2]
      u2 <- exputility_random[3 * i - 1]
      u3 <- exputility_random[3 * i]
      p1 <- u1 / (u1 + u2 + u3)
      p2 <- u2 / (u1 + u2 + u3)
      sump <- p1 + p2
      rnd <- runif(1)
      if (p1 >= rnd) {
        r[3 * i - 2, 1] <- 1
      } else if (sump >= rnd) {
        r[3 * i - 1, 1] <- 1
      } else {
        r[3 * i, 1] <- 1
      }
    }
    
    # Create data frame 'm' with choices
    m <- data.frame(
      subject_id = rep(1:nresps, each = nruns),
      chosen = r,
      choice_set_id = rep(ff[,1], nresps),
      X1 = f_random[, 2],
      X2 = f_random[, 3],
      X3 = f_random[, 4],
      X4 = f_random[, 5],
      X5 = f_random[, 6]
    )
    
    # Effect coding
    choices_effect_coded <- matrix(0, nrow = nrow(m), ncol = nparam)
    m_matrix <- as.matrix(m[, 4:8])
    for (i in 1:nrow(m)) {
      choices_effect_coded[i, ] <- frow2xrow(m_matrix[i, ], matrix(nlevels, nrow = nattributes, ncol = 1))
    }
    
    # Prepare choices for mlogit
    choices <- matrix(0, nrow = nsets * nresps, ncol = nparam * nalts)
    for (i in 1:nrow(choices)) {
      choices[i, 1:13] <- choices_effect_coded[i * 3 - 2, ]
      choices[i, 14:26] <- choices_effect_coded[i * 3 - 1, ]
      choices[i, 27:39] <- choices_effect_coded[i * 3, ]
    }
    
    # Indicate chosen profile
    profile_choice <- rep(NA, nsets * nresps)
    for (i in 1:length(profile_choice)) {
      if (m[i * 3 - 2, "chosen"] == 1) {
        profile_choice[i] <- "A"
      } else if (m[i * 3 - 1, "chosen"] == 1) {
        profile_choice[i] <- "B"
      } else if (m[i * 3, "chosen"] == 1) {
        profile_choice[i] <- "C"
      }
    }
    
    # Create subject IDs and choice set numbers
    id_and_choice_set <- data.frame(
      subject_id = rep(1:nresps, each = nsets),
      choice_set_id = rep(1:nsets, times = nresps)
    )
    
    info <- cbind(id_and_choice_set, profile_choice)
    choices_df <- cbind(info, as.data.frame(choices))
    colnames(choices_df) <- c(
      "subject_id", "choice_set_id", "chosen",
      paste0("V", 1:13, "_A"), paste0("V", 1:13, "_B"), paste0("V", 1:13, "_C")
    )
    
    # Convert to mlogit.data
    data_z <- mlogit.data(choices_df, choice = "chosen", shape = "wide", varying = 4:42, sep = "_", id.var = "subject_id")
    
    # Try model estimation
    tryCatch({
      fit_PML_with_order_random_primacy <- mlogit(chosen ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13|0, 
                                                  data = data_z,
                                                  rpar = c(V1 = "n", V2 = "n", V3 = "n", V4 = "n", V5 = "n", V6 = "n", V7 = "n", V8 = "n", V9 = "n", V10 = "n", V11 = "n", V12 = "n", V13 = "n"), 
                                                  panel = TRUE, R = 1000,halton = list(dim = 13))
      
      beta_PML_random_primacy <- rbind(beta_PML_random_primacy, as.data.frame(t(fit_PML_with_order_random_primacy[["coefficients"]][1:13])))
      sd_PML_random_primacy <- rbind(sd_PML_random_primacy, as.data.frame(t(fit_PML_with_order_random_primacy[["coefficients"]][14:26])))
      estimation_success <- TRUE
      duration <- Sys.time() - start_time
      cat('draw:', b+1, 'duration:', duration, "\n")
    }, error = function(e) {
      # Do nothing on error, retry
    })
  }
  
  b <- b + 1
}



write.csv(beta_PML_random_primacy, ".../beta_PML_random_primacy.csv", row.names = FALSE)
write.csv(sd_PML_random_primacy, ".../sd_PML_random_primacy.csv", row.names = FALSE)

write.csv(beta_PML_naive_primacy, ".../beta_PML_naive_primacy.csv", row.names = FALSE)
write.csv(sd_PML_naive_primacy, ".../sd_PML_naive_primacy.csv", row.names = FALSE)








