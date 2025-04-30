M:=CuspidalSubspace(ModularSymbols(61,2,-1));
f := NewformDecomposition(M)[2];
K:=HeckeEigenvalueField(f);
L:=Compositum(K,CyclotomicField(9));
N:=ModularSymbols(61,2);
G<a>:=DirichletGroup(19,L);
A:=HeckeOperator(f,2);
// ClassNumber(K);
B:=Matrix(K,A);
V, U, F:=JordanForm(B);
R:=Matrix(L,V);
S:=Matrix(L,U);
C:=0;
v43:=N!<1,[Cusps()|-1/43,0]>;
v20:=N!<1,[Cusps()|-1/20,0]>;
v38:=N!<1,[Cusps()|-1/38,0]>;
v26:=N!<1,[Cusps()|-1/26,0]>;

print("M : "); M;
print("");
print("f : "); f;
print("");
print("A : "); A;
print("");
print("K : "); K;
print("");
print("Class Number of K: "); ClassNumber(K);
print("");
print("L : "); L;
print("");
// print("F: "); F;
// print("");
// print("U : "); U;
// print("");
// print("V: "); V;
// print("");

for i in {1..18} do
D:=N!<1,[Cusps()|0,i/19]>;
C:=C+Evaluate(a^2,i)*((S[1][1])*(IntersectionPairing(D,v43)+IntersectionPairing(D,v20))
+(S[1][2])*(IntersectionPairing(D,v38)+(1/2)*IntersectionPairing(D,v20))
+(S[1][3])*(IntersectionPairing(D,v26)-(1/2)*IntersectionPairing(D,v20)));
end for;

print("C: ");C;
print("");
Z:=IntegerRing();
Factorization(Z!(Norm(C)));
for j in {1..12} do
P:=Decomposition(L,19)[j];
Evaluate(C, P[1]), Evaluate(Evaluate(a,2), P[1]);
end for;