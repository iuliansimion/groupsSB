
#
# Read("~/Workspace/groupsSB/epi/A2/case.1.gi");
#

type:="A";
rank:=2;
nr_pos_roots:=3;

Read("~/Workspace/groupsSB/epi/group.gi");

highest_weight:=[1,1];# corresponding highest weight module is the Lie algebra

Read("~/Workspace/groupsSB/epi/ext.gi");

#
# p>2
#

getRep:=function(k,t)
    local u,ua,ub,uab;
    ua:=ext_root_group(k,1,t);
    ub:=ext_root_group(k,2,t);
    uab:=ext_root_group(k,3,t^2/2); 
    u:=ua*ub*uab;

    return u;   
end;

#
# k:=3; # power of exterior algebra
#
handle:=function(k)
    local id_mat,u,e,h,h1,h2,ese,evh,esh,esh0,index_ev0,si,inter;
    u:=getRep(k,1);
    id_mat:=u^0;
    e:=u-id_mat;
    ese:=Eigenspaces(Rationals,e)[1]; #eigenvalues 0, centralizer
    
    # torus
    h1:=ext_e(k,cb0[7]);
    h2:=ext_e(k,cb0[8]);;
    h:=h1+h2;   
    esh:=Eigenspaces(Rationals,h);
    evh:=Eigenvalues(Rationals,h);

    # centralizer of torus
    esh0:=[];
    index_ev0:=Position(evh,0);
    esh0:=esh[index_ev0];

    #
    # sum and intersection
    #
    si:=SumIntersectionMat(Basis(ese),Basis(esh0));
    inter:=si[2];
    return inter;
end;