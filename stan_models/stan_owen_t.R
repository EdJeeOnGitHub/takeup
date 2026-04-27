library(OwenQ)

# OwenT(h, a) requires scalar a; wrap to match the stan_owen_t(vector, vector) interface
exposed_stan_func <- list(
  stan_owen_t = function(h, a) mapply(OwenT, h, a)
)
