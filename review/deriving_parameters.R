# Simulate data using JAGS for the puposes of deriving parameters, demonstration
# Josh Nowak
# 12/2022
# library(jagsUI)
# library(tidybayes)
# library(dplyr)

# The model doesn't matter much, so make it simple
model_file <- tempfile()

writeLines("
model{
  for(p in 1:n_pop) {
    N[p,1] ~ dnorm(p*10, 0.001)
  }

  lambda ~ dnorm(1, 1)

  for(p in 1:n_pop) {
    for(y in 2:n_year){
      N[p,y] <- N[p,y-1] * lambda
    }
  }

  # Derive total population
  for(y in 1:n_year) {
    Ntot[y] <- sum(N[,y])
  }

  # Derive mean of each population
  for(p in 1:n_pop) {
    Nmu[p] <- mean(N[p,])
  }

}
", con = model_file)



# We will define data for 3 populations over 10 years
jdat <- list(
  n_pop = 3,
  n_year = 10
)

out <- jagsUI::jags(
  data = jdat,
  inits = NULL,
  parameters.to.save = c("Ntot", "Nmu", "N"),
  model.file = model_file,
  n.iter = 2000,
  n.chains = 3,
  DIC = FALSE
)


# out is now our reference for evaluating order of operations on various
#  parameter derivations

# First let's look at the posteriors
str(out$samples)

# They are also stored in...
str(out$sims.list)

# Let's start easy and get the mean of a each herd in each year
# Since the original script used tidybayes let's keep that convention
post <- out |>
  tidybayes::gather_draws(N[i,j])

# Using tidybayes
post |>
  tidybayes::mean_qi()

out$mean$N
out$q2.5$N
out$q97.5$N

# Looks good, but no idea what happened
# Two step version, first get the mean per iteration across all chains
iter_mu <- post |>
  dplyr::group_by(i, j, .iteration) |>
  dplyr::summarize(
    mean = mean(.value),
    .groups = "drop"
  )

# Now get mean over all iterations
mu_herd_year <- iter_mu |>
  dplyr::group_by(i, j) |>
  dplyr::summarize(
    mu = mean(mean)
  )

mu_herd_year |>
  dplyr::filter(i == 1)

# Compare to summary done by jagsUI
out$mean$N[1,]

# That worked, but only because the mean is forgiving
# Now let's see if we can get the sum correct
# What if we start with the mean over the iterations
iter_mu |>
  dplyr::group_by(j, .iteration) |>
  dplyr::summarize(
    totN = sum(mean),
    .groups = "drop"
  ) |>
  dplyr::group_by(j) |>
  dplyr::summarize(
    totN_mu = mean(totN),
    lcl = quantile(totN, 0.025),
    ucl = quantile(totN, 0.975)
  )

# Compare it to model quantity, the point estimate is correct, but the quantiles
#  are different
out$mean$Ntot
out$q2.5$Ntot
out$q97.5$Ntot

# This time let's start from the beginning, sum per draw and then summarize
# A draw is different because it is unique among chains whereas the iteration is
# unique within a chain, so there are 3 iterations with value 1 if 3 chains
tot_draw <- post |>
  dplyr::group_by(j, .draw) |>
  dplyr::summarize(
    mean = sum(.value),
    .groups = "drop"
  )

# Now that we have the sum per draw, summarize that posterior
tot_draw |>
  dplyr::group_by(j) |>
  dplyr::summarize(
    lcl = quantile(mean, 0.025),
    tot = mean(mean),
    ucl = quantile(mean, 0.975)
  )

# Perfect, that worked well
out$q2.5$Ntot
out$mean$Ntot
out$q97.5$Ntot

# One more time to get the mean per population, the other derived quantity in
#  the model, one step
post |>
  dplyr::group_by(i, .draw) |>
  dplyr::summarize(
    mean = sum(.value),
    .groups = "drop"
  ) |>
  dplyr::group_by(i) |>
  dplyr::summarize(
    lcl = quantile(mean, 0.025),
    tot = mean(mean),
    ucl = quantile(mean, 0.975)
  )

out$q2.5$Nmu
out$mean$Nmu
out$q97.5$Nmu

# Oldschool, just use apply and the sims.list
# Total N
dd <- apply(out$sims.list$N, c(1, 3), sum)
apply(dd, 2, mean)
out$mean$Ntot

apply(dd, 2, quantile, 0.025)
out$q2.5$Ntot

################################################################################
# End