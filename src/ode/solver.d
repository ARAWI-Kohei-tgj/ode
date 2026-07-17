module ode.solver;

import ode.plot;
/*************************************************************
 * Classical Runge-Kutta method
 *
 * Params:
 *  R= type of floating point
 *  Order= order of the ode
 *
 *  f= function ddot(x) = f(t, x, dot(x)) 
 *  tInit= initial time
 *  tMax= terminal time
 *  xInit= initial value of x
 *  vInit= initial value of dot(x)
 *  h= time step
 *************************************************************/
ScatterPlot!(R, Order) odeRungeKutta4(R, int Order)(in R function(R, R, R) f,
	in R tInit,
	in R tMax,
	in R xInit,
	in R vInit,
	in R h)
if(Order == 2){
	import std.numeric: dotProduct;
	enum CurveOdr= 4;
	enum R[CurveOdr] weight= [1.0, 2.0, 2.0, 1.0];

	R x= xInit;
	R v= vInit;
	R[CurveOdr] k, l;
	ScatterPlot!(R, Order) result;

	result.append(tInit, xInit, vInit);

	for(R t= tInit+h; t <= tMax; t += h){
		k[0]= v;
		l[0]= f(t, x, v);

		k[1]= (v+l[0]*h/2);
		l[1]= f(t+h/2, x+k[0]*h/2, v+l[0]*h/2);

		k[2]= (v+l[1]*h/2);
		l[2]= f(t+h/2, x+k[1]*h/2, v+l[1]*h/2);

		k[3]= (v+l[2]*h);
		l[3]= f(t+h, x+k[2]*h, v+l[2]*h);

		x += h/6*dotProduct(k, weight);
		v += h/6*dotProduct(l, weight);

		result.append(t, x, v);
	}

	return result;
}

/// ditto
ScatterPlot!(R, Order) odeRungeKutta4(R, int Order)(in R function(R, R) f,
	in R tInit,
	in R tMax,
	in R xInit,
	in R h)
if(Order == 1){
	import std.numeric: dotProduct;
	enum CurveOdr= 4;
	enum R[CurveOdr] weight= [1.0, 2.0, 2.0, 1.0];

	R x= xInit;
	R[CurveOdr] k;
	auto result= ScatterPlot!(R, Order)();

	result.append(tInit, xInit);

	for(R t= tInit+h; t <= tMax; t += h){
		k[0]= f(t, x);
		k[1]= f(t+h/2, x+h*k[0]/2);
		k[2]= f(t+h/2, x+h*k[1]/2);
		k[3]= f(t+h, x+h*k[2]);

		x += h/6.0*dotProduct(k, weight);

		result.append(t, x);
	}

	return result;
}
 
