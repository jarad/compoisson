
com.fit = function(x)
{
	xbar = (x[,1] %*% x[,2]) / sum(x[,2]);
	options(warn = -1);
	result = optim(c(xbar,1), function(p) {return (-com.loglikelihood(x, p[1], p[2]));},
		method="L-BFGS-B", lower=c(1e-10,0));
	options(warn = 0);
	
	lambda = result$par[1];
	nu = result$par[2];
	fit = list( lambda = lambda,
	            nu = nu,
	            z = com.compute.z(lambda, nu),
	            fitted.values = sum(x[,2]) * dcom(x[,1], lambda, nu),
				log.likelihood = com.loglikelihood(x, lambda, nu) );

	return (fit);
}



com.compute.log.z = function(lambda, nu, log.error = 0.001)
{
  # Perform argument checking
  if (lambda <= 0 || nu < 0)
    stop("Invalid arguments, only defined for lambda > 0, nu >= 0");
	
  # Initialize values
  j = 0;
  llambda = log(lambda)                                 # precalculate for speed
  inclfact = 0                                          # log(factorial(0))
  z = j * llambda - nu * inclfact;                      # first term in sum
  z.last = -Inf                                         # to ensure entering the loop

  # Continue until we have reached specified precision
  while (abs(z - z.last) > log.error) {
    z.last = z;                                         # For comparison in while statement
    z = com.log.sum(z, j * llambda - nu * inclfact );   # Log of current sum
    j = j + 1;                                          # Next term in sum
    inclfact = inclfact+log(j)                          # Calculate increment for log factorial
  }
  return (z);
}



com.compute.z = function(lambda, nu, log.error = 0.001)
{
  return (exp(com.compute.log.z(lambda,nu,log.error)));
}

dcom = function(x, lambda, nu, z = NULL)
{
	# Perform argument checking
	if (lambda < 0 || nu < 0)
		stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
	if (x < 0 || x != floor(x))
		return (0);
	if (is.null(z) || z <= 0)
		z = com.compute.z(lambda, nu);
	
	# Return pmf
	return ((lambda ^ x) * ((factorial(x)) ^ -nu) / z);
}

com.log.density = function(x, lambda, nu, log.z = NULL)
{
	# Perform argument checking
	if (lambda < 0 || nu < 0)
		stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
	if (x < 0 || x != floor(x))
		return (0);
	if (is.null(log.z)) { log.z = com.compute.log.z(lambda, nu); }
	
	# Return log pmf
	return ((x * log(lambda) - nu * lfactorial(x)) - log.z);
}

com.loglikelihood = function(x, lambda, nu)
{
	# Perform argument checking
	if (lambda < 0 || nu < 0)
		return (-Inf);

	log.z = com.compute.log.z(lambda, nu);
	return (x[,2] %*% ( x[,1] * log(lambda) - nu * lfactorial(x[,1]) - log.z ));
}

com.expectation = function(f, lambda, nu, log.error = 0.001)
{
	log.z = com.compute.log.z(lambda, nu);

	# Initialize variables
	ex = -Inf;
	ex.last = 0;
	j = 0;

	# Continue until we have reached specified precision
	while ((ex == -Inf && ex.last == -Inf) || abs(ex - ex.last) > log.error)
	{
		ex.last = ex;
		ex = com.log.sum(ex, log(f(j)) + com.log.density(j, lambda, nu, log.z));

		j = j + 1;
	}
	return (exp(ex));
}

com.mean = function(lambda, nu)
{
	return ( com.expectation(function (x) {x}, lambda, nu) );
}

com.var = function(lambda, nu)
{
	return ( com.expectation(function(x) {x^2}, lambda, nu) - (com.mean(lambda,nu))^2 );
}

rcom = function(n, lambda, nu, log.z = NULL) {
  # Check arguments
  if (lambda < 0 || nu < 0)
    stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  if (nu == 1) return(rpois(n, lambda))

  if (is.null(log.z))
    log.z = com.compute.log.z(lambda, nu);

  r = NULL;	# Vector of random values

  for (i in 1:n) {
    # Get a uniform random variable and find the smallest value x such that
    # the cdf of x is greater than the random value
    log.prob = log(runif(1));
    j = 0;
    log.dens = com.log.density(j, lambda, nu, log.z)
    while (log.prob>log.dens) {
      j = j+1
      log.prob = log(exp(log.prob)-exp(log.dens))     # Decrement the remaining probability
      log.dens = com.log.density(j,lambda,nu,log.z)
    }

    r = c(r, j);
  }

  return (r);
}


com.confint = function(data, level=0.95, B=1000, n=1000)
{
	# B = Number of repetitions of bootstrap
	# n = number of obs. in each sample used to find bootstrap mean
	# data = the original data
	# level - confidence interval level


	# Check arguments
	if (level <= 0 || level >= 1)
		stop("Invalid arguments, 0 < level < 1");
	if (n <= 0)
		stop("Invalid arguments, n > 0");
	if (B <= 0)
		stop("Invalid arguments, B > 0");


	boot.lambda = matrix(0,B,1)
	boot.nu = matrix(0,B,1)

	# Sample statistic
	COMobject = com.fit(data)
	lambda.mle = COMobject$lambda 
	nu.mle = COMobject$nu

	# Changing data into form we can sample from
	fulldata = matrix(0,sum(data[,2]),1)
	index = 0
	for (i in 1:length(data[,1]))
	{
		index2 = index + data[i,2]
		fulldata[(index+1):index2] = data[i,1]
	}

	# Creates a vector of means (the bootstrap vector)
	for (i in 1:B)
	{
		samplewR = sample(fulldata, n, replace=TRUE);
		sample = data.matrix( as.data.frame(table(samplewR)) );
		sample = cbind( sample[,1] - 1, sample[,2] );
		COMsampleobject = com.fit(sample);
		boot.lambda[i] = COMsampleobject$lambda;
		boot.nu[i] = COMsampleobject$nu;
		remove(samplewR);
	}

	# Pivotal method calculation
	boot.lambda = sort(boot.lambda)
	boot.nu = sort(boot.nu)

	lower.bound = (1 - level) / 2;
	upper.bound = 1 - lower.bound;

	lower.index = floor(lower.bound * B)+1;
	upper.index = floor(upper.bound * B);

	pivotal.ci = matrix(0,2,2);
	pivotal.ci[1,1] = max(0, 2*lambda.mle - boot.lambda[upper.index]);
	pivotal.ci[1,2] = 2*lambda.mle - boot.lambda[lower.index];

	pivotal.ci[2,1] = max(0, 2*nu.mle - boot.nu[upper.index]);
	pivotal.ci[2,2] = 2*nu.mle - boot.nu[lower.index];

	rownames(pivotal.ci) = c("lambda", "nu");
	colnames(pivotal.ci) = c(lower.bound, upper.bound);

	return (pivotal.ci);
}

