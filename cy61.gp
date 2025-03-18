ns_xyz(x,y,z,  p) = {
 A0 = (1+x)*(1+y)*(1+x+z)*(1+z+y*z);
 A1 = (1+x)*(1+y)*(1+x+z)*(2+2*z+y*z) - x*y*z;
 A2 = (1+x)*(1+y)*(1+x+z)*(1+z);
 if(A2%p,
   if(A0%p, 1+kronecker(A1^2 - 4*A0*A2,p), A1%p != 0),
   if(A1%p, A0%p!=0, if(A0%p,0,p-1));
 );
}

num_sols(p) = {
 su=0;
 forvec(X=[[1,p-1],[1,p-1],[1,p-1]],
   su+=ns_xyz(X[1],X[2],X[3], p);
 );
 su;
}

tr_frob(p) = {
 p^3 - 8*p^2 + 25*p - 33 - num_sols(p);
}

forprime(p=2,61,print("a("p") = "tr_frob(p)));

\q