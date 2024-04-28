using Test
using ToricAtiyahBott
using Oscar

v = projective_space(NormalToricVariety, 3);
line = cohomology_class(toric_divisor(v, [1,0,0,0]))^2;
P = ev(1, a_point(v))*ev(2, a_point(v));
t = IntegrateAB(v, line, 2, P)[1];

@test t == 1

x = domain(blow_up(v, [1,1,1]; coordinate_name="Ex1"));
X = domain(blow_up(x, [-1,0,0]; coordinate_name="Ex2"));
mg = moment_graph(X);
(H, E1, E2) = (mg[7,8], mg[4,5], mg[1,2]);
(d, e1, e2) = (2, -2, -2);
beta = d*H + e1*E1 + e2*E2;
P = class_one();
s = IntegrateAB(X, beta, 0, P)[1];

@test s == 1//8