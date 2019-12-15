
#
# Read("~/Workspace/groupsSB/epi/G2/generic.UU.gi");
#

type:="G";
rank:=2;
nr_pos_roots:=6;

Read("~/Workspace/groupsSB/epi/group.gi");



handleUaUbUc:=function()
	local rels, vals,i,v;
	# UaUb=Uc = U(a1+b1,a2+b2,a3+b3+a2*b1)
	rels:=Set(Concatenation(Ua*Ub-Uc));
	
	vals:=[];
	Append(vals,[[[vars[22]],[vars[2]+vars[12]]]]);
	Append(vals,[[[vars[21]],[vars[1]+vars[11]]]]);
	Append(vals,[[[vars[23]],[vars[3]+vars[13]+vars[2]*vars[11]]]]);
#a_{2}^2*b_{1}+2*a_{2}*b_{1}*b_{2}-2*a_{3}*b_{2}-a_{4}-b_{4}
#-a_{2}*b_{1}^2-2*a_{3}*b_{1}-a_{4}-b_{4}+c_{4},
	Append(vals,[[[vars[24]],[vars[2]*vars[11]^2+2*vars[3]*vars[11]+vars[4]+vars[14]]]]);
#a_{2}*b_{1}^3+3*a_{3}*b_{1}^2+3*a_{4}*b_{1}+a_{5}+b_{5}-c_{5}
	Append(vals,[[[vars[25]],[vars[2]*vars[11]^3+3*vars[3]*vars[11]^2+3*vars[4]*vars[11]+vars[5]+vars[15]]]]);
#a_{2}^2*b_{1}^3
#-a_{2}*b_{1}^3*b_{2}
#+3*a_{2}*a_{3}*b_{1}^2
#+3*a_{2}*b_{1}^2*b_{3}
#-3*a_{3}*b_{1}^2*b_{2}

#+3*a_{3}^2*b_{1}
#+6*a_{3}*b_{1}*b_{3}
#-3*a_{4}*b_{1}*b_{2}
#+3*a_{4}*b_{3}
#-a_{5}*b_{2}
#-a_{6}-b_{6}+c_{6}
	Append(vals,[[[vars[26]],[
-vars[2]^2*vars[11]^3
+vars[2]*vars[11]^3*vars[12]
-3*vars[2]*vars[3]*vars[11]^2
-3*vars[2]*vars[11]^2*vars[13]
+3*vars[3]*vars[11]^2*vars[12]
-3*vars[3]^2*vars[11]
-6*vars[3]*vars[11]*vars[13]
+3*vars[4]*vars[11]*vars[12]
-3*vars[4]*vars[13]
+vars[5]*vars[12]
+vars[6]+vars[16]]]]);


	for i in [1..Length(rels)] do
		for v in vals do
			Print(rels[i],"\n");
			rels[i]:=One(APR)*Value(rels[i],v[1],v[2]);
		od;
	od;
	return [vals,rels];
end;