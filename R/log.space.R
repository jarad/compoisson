com.log.sum = function(x,y)		# log.sum(x,y) = log( exp(x) + exp(y) )
{
	if (x == -Inf)
		{ return (y); }
	else if (y == -Inf)
		{ return (x); }
	else if (x > y)
		{ return (x + log1p( exp(y - x) ) ); }
	else
		{ return (y + log1p( exp(x - y) ) ); }
}


