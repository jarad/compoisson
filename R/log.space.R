com.log.sum = function(x,y)		# log.sum(x,y) = log( exp(x) + exp(y) )
{
	if (x == -Inf)
		{ return (y); }
	else if (y == -Inf)
		{ return (x); }
	else if (x > y)
		{ return (x + log( 1 + exp(y - x) ) ); }
	else
		{ return (y + log( 1 + exp(x - y) ) ); }
}

com.log.difference = function(x,y)	# log.difference(x,y) = log( exp(x) - exp(y) )
{
	if (x == -Inf)
		{ return (NaN); }
	else if (y == -Inf)
		{ return (x); }
	else if (x > y)
		{ return (x + log( 1 - exp(y - x) ) ); }
	else
		{ return (NaN); } # y<x returns NaN?
}

