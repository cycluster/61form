// Code for the case of the K3 Fibration (projecting w first) 

// What is l61 for again? 

l61:=[1, -7, -3, 25, 3, 21, -9, -63, 6, -21, -4, -75, -3, 63, -9, 169, 37, -42, -75, 75, 27, 28, 10, 189, -76, 21, -90, -225, 212, 63, -6, -623, 12, -259, -27, 150, -88, 525, 9, -189, -3, -189, 547, -100, 18, -70, -147, -507, 25];

getEC := function(W,p);
	K<Z> := FunctionField(GF(p));
	R<X> := PolynomialRing(K);
	return(EllipticCurve(X^3 + (-1/48*Z^8*W^8 - 1/12*Z^8*W^7 - 1/12*Z^7*W^\
	8 - 1/8*Z^8*W^6 - 1/3*Z^7*W^7 - 1/8*Z^6*W^8 - 1/12*Z^8*W^5 - 5/12*Z^7*W^6 - 1/\
	3*Z^6*W^7 - 1/12*Z^5*W^8 - 1/48*Z^8*W^4 - 1/12*Z^7*W^5 + 1/6*Z^6*W^6 + 1/6*Z^5\
	*W^7 - 1/48*Z^4*W^8 + 1/6*Z^7*W^4 + 5/4*Z^6*W^5 + 7/4*Z^5*W^6 + 5/12*Z^4*W^7 +\
	 1/12*Z^7*W^3 + 7/6*Z^6*W^4 + 49/12*Z^5*W^5 + 13/8*Z^4*W^6 + 1/6*Z^3*W^7 + 1/6\
	*Z^6*W^3 + 53/12*Z^5*W^4 + 19/6*Z^4*W^5 - 1/8*Z^6*W^2 + 9/4*Z^5*W^3 + 197/48*Z\
	^4*W^4 - Z^3*W^5 - 1/3*Z^2*W^6 + 1/2*Z^5*W^2 + 19/6*Z^4*W^3 - 4/3*Z^3*W^4 - 4/\
	3*Z^2*W^5 + 1/12*Z^5*W + 31/24*Z^4*W^2 - 1/3*Z^3*W^3 - 2*Z^2*W^4 + 1/4*Z^4*W +\
	 1/3*Z^3*W^2 - 4/3*Z^2*W^3 - 1/48*Z^4 + 1/6*Z^3*W - 1/3*Z^2*W^2)*X + (1/864*Z^\
	12*W^12 + 1/144*Z^12*W^11 + 1/144*Z^11*W^12 + 5/288*Z^12*W^10 + 1/24*Z^11*W^11\
	 + 5/288*Z^10*W^12 + 5/216*Z^12*W^9 + 7/72*Z^11*W^10 + 13/144*Z^10*W^11 + 5/21\
	6*Z^9*W^12 + 5/288*Z^12*W^8 + 5/48*Z^11*W^9 + 41/288*Z^10*W^10 + 5/72*Z^9*W^11\
	 + 5/288*Z^8*W^12 + 1/144*Z^12*W^7 + 5/144*Z^11*W^8 - 5/144*Z^10*W^9 - 5/36*Z^\
	9*W^10 - 5/144*Z^8*W^11 + 1/144*Z^7*W^12 + 1/864*Z^12*W^6 - 1/36*Z^11*W^7 - 25\
	/72*Z^10*W^8 - 103/108*Z^9*W^9 - 169/288*Z^8*W^10 - 7/72*Z^7*W^11 + 1/864*Z^6*\
	W^12 - 1/36*Z^11*W^6 - 55/144*Z^10*W^7 - 65/36*Z^9*W^8 - 67/36*Z^8*W^9 - 13/24\
	*Z^7*W^10 - 1/16*Z^6*W^11 - 1/144*Z^11*W^5 - 13/96*Z^10*W^6 - 115/72*Z^9*W^7 -\
	 815/288*Z^8*W^8 - 137/144*Z^7*W^9 - 7/96*Z^6*W^10 - 1/72*Z^5*W^11 + 1/48*Z^10\
	*W^5 - 5/8*Z^9*W^6 - 313/144*Z^8*W^7 - 11/48*Z^7*W^8 + 329/432*Z^6*W^9 + 5/36*\
	Z^5*W^10 + 5/288*Z^10*W^4 - 5/72*Z^9*W^5 - 53/96*Z^8*W^6 + 73/36*Z^7*W^7 + 467\
	/144*Z^6*W^8 + 7/8*Z^5*W^9 + 1/18*Z^4*W^10 - 1/36*Z^9*W^4 + 5/18*Z^8*W^5 + 305\
	/72*Z^7*W^6 + 1007/144*Z^6*W^7 + 9/4*Z^5*W^8 + 1/9*Z^4*W^9 - 5/216*Z^9*W^3 + 5\
	5/288*Z^8*W^4 + 589/144*Z^7*W^5 + 937/96*Z^6*W^6 + 97/24*Z^5*W^7 - 5/18*Z^4*W^\
	8 - 2/27*Z^3*W^9 + 7/144*Z^8*W^3 + 95/48*Z^7*W^4 + 1237/144*Z^6*W^5 + 203/36*Z\
	^5*W^6 - 10/9*Z^4*W^7 - 4/9*Z^3*W^8 + 5/288*Z^8*W^2 + 7/18*Z^7*W^3 + 611/144*Z\
	^6*W^4 + 391/72*Z^5*W^5 - 4/3*Z^4*W^6 - 10/9*Z^3*W^7 + 43/48*Z^6*W^3 + 28/9*Z^\
	5*W^4 - 5/9*Z^4*W^5 - 40/27*Z^3*W^6 - 1/144*Z^7*W - 5/288*Z^6*W^2 + 7/8*Z^5*W^\
	3 + 1/6*Z^4*W^4 - 10/9*Z^3*W^5 - 1/48*Z^6*W + 1/18*Z^5*W^2 + 2/9*Z^4*W^3 - 4/9\
	*Z^3*W^4 + 1/864*Z^6 - 1/72*Z^5*W + 1/18*Z^4*W^2 - 2/27*Z^3*W^3)));
 end function;
// Type(getEC(3)); is CrvEll

// Equation from monodromy of K3 PF equation
// Z := IntegerRing();
// R<t> := PolynomialRing(Z);
// MonPol := t*(t+1)*(t^5 + 28*t^4 + 156*t^3 + 240*t^2 + 108*t -4);

GetMonPolRoots := function(p)
	R<t>:=PolynomialRing(GF(p));
	pol := t*(t+1)*(t^5 + 28*t^4 + 156*t^3 + 240*t^2 + 108*t -4);
	roots := Roots(pol);
	v := [];
	for i in roots do
		Append(~v,i[1]);
	end for;
	return(v);
end function;

GoodFibreVals := function(p)
	F := GF(p);
	fp := [i : i in [0..p-1]];
	vec := [x : x in fp | not x in GetMonPolRoots(p)];
	return vec;
end function;

LFunFacEC := function(w,p) 
	E := getEC(w,p);
	return(Factorisation(LFunction(E)));
end function;

quarticTerm := function(w,p)
 print w;
 pols:=LFunFacEC(w,p);
 for x in pols do
 if (Degree(x[1]) eq 4) then
 return(x[1]);
 end if;
 end for;
 return(pols);
 end function;

printquartic := function(p)
	for i in GoodFibreVals(p) do
		quarticTerm(i,p);
	end for;
	return(0);
end function;

SumOverFibres := function(p)
print "Generating sum over fibres for", p;
print "Bad fibres are:", GetMonPolRoots(p);
s:= [];
for i in GoodFibreVals(p) do
print "Checking", i;
pols:=LFunFacEC(i,p);
for x in pols do
if (Degree(x[1]) eq 4) then
Append(~s,-1*Coefficient(x[1],1)*p^4);
end if;
end for;
end for;
print "Sum over fibres:";
return(&+s);
end function;

// Example: p = 29
SumOverFibres(29);