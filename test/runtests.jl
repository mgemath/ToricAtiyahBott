using Test
using ToricAtiyahBott
using Oscar

v = projective_space(NormalToricVariety, 3);
line = cohomology_class(toric_divisor(v, [1,0,0,0]))^2;
P = ev(1, a_point(v))*ev(2, a_point(v));
t = IntegrateAB(v, line, 2, P)[1];

@test t == 1