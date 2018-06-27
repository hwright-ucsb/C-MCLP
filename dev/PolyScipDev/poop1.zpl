set I := {1..3};
param c1[I] := <1> 3, <2> 2, <3> -4;
param c2[I] := <3> 2 default 1;
param low[I] := <1> 1, <2> 0.4, <3> 0;
var x[I] integer >= 0;
maximize Obj1: sum <i> in I: c1[i]*x[i];
Obj2: sum <i> in I: c2[i]*x[i];
subto Eqn: sum <i> in I: x[i] == 2;
subto Lower: sum <i> in I: low[i]*x[i] <= 1.5;