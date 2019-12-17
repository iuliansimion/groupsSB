
#
# Read("~/Workspace/groupsSB/epi/A2/handle.ext.gi");
#

type:="A";
rank:=2;
nr_pos_roots:=3;

Read("~/Workspace/groupsSB/epi/group.gi");

highest_weight:=[1,1];# corresponding highest weight module is the Lie algebra

Read("~/Workspace/groupsSB/epi/ext.gi");


ext_e:=function(k,e)
	local W,extb,result,v;
    W:=extA[k];
    extb:=extA_basis[k];

	result:=[];
	for v in extb do
		Append(result,[Coefficients(extb,e^v)]);
	od;
	result:=TransposedMat(result);
	return result;
end;


ext_root_group:=function(k,index,t)
	local ee,tmp,result,i;
	ee:=ext_e(k,cb0[index]);
	tmp:=ee;
	result:=tmp^0;
	i:=1;
	while Length(Set(Concatenation(tmp)))<>1 do
		result:=result+t^i*tmp/Factorial(i);
		i:=i+1;
		tmp:=tmp*ee;
	od;
	return result;
end;

handleReg:=function()
    local k,ua,ub,uab,u,e,h,h1,h2;
    k:=3; # power of exterior algebra
    ua:=ext_root_group(k,1,1);
    ub:=ext_root_group(k,2,1);
    uab:=ext_root_group(k,3,1/2); #p\neq 2
    u:=ua*ub*uab;
    #
    #uua:=ext_root_group(1,2);
    #uub:=ext_root_group(2,2);
    #uuab:=ext_root_group(3,2);
    #uu:=uua*uub*uuab;
    #
    id_mat:=u^0;
    e:=u-id_mat;

    ese:=Eigenspaces(Rationals,e)[1]; #eigenvalues 0, centralizer
    h1:=ext_e(k,cb0[7]);
    h2:=ext_e(k,cb0[8]);;
    h:=h1*h2;   
    esh1:=Eigenspaces(Rationals,h);
    evh1:=Eigenvalues(Rationals,h);

    esh10:=[];
    index_ev0:=Position(evh1,0);
    esh10:=esh1[index_ev0];

    #
    # sum and intersection
    #
    si:=SumIntersectionMat(Basis(ese),Basis(esh10));
    inter:=si[2];
    return inter;
end;

#
# c_2=0
#
#handle2:=function()
#   local k,ua,ub,uab,u,e;
k:=3;
ua:=ext_root_group(k,1,1);
uab:=ext_root_group(k,3,1);
u:=ua*uab;
e:=u-id_mat;

ese:=Eigenspaces(Rationals,e)[1]; #eigenvalues 0, centralizer
h1:=ext_e(k,cb0[7]);;           
esh1:=Eigenspaces(Rationals,h1);
evh1:=Eigenvalues(Rationals,h1);

esh10:=[];
index_ev0:=Position(evh1,0);
esh10:=esh1[index_ev0];

#
# sum and intersection
#
si:=SumIntersectionMat(Basis(ese),Basis(esh10));
inter:=si[2];
#    return inter;
#end;