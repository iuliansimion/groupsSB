
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
    #ub:=ext_root_group(k,2,t);
    uab:=ext_root_group(k,3,t); 
    u:=ua*uab;

    return u;   
end;

#
# k:=3; # power of exterior algebra
#
handle:=function(k)
    local id_mat,u,e,h,h1,h2,ese,evh1,esh1,evh2,esh2,esh01,index_ev01,esh02,index_ev02,si,si1,si2,inter;
    u:=getRep(k,1);
    id_mat:=u^0;
    e:=u-id_mat;
    ese:=Eigenspaces(Rationals,e)[1]; #eigenvalues 0, centralizer
    
    # torus
    # centralizer of T1
    h1:=ext_e(k,cb0[7]);
    esh1:=Eigenspaces(Rationals,h1);
    evh1:=Eigenvalues(Rationals,h1);
    esh01:=[];
    index_ev01:=Position(evh1,0);
    esh01:=esh1[index_ev01];

    # centralizer of T1
    h2:=ext_e(k,cb0[8]);
    esh2:=Eigenspaces(Rationals,h2);
    evh2:=Eigenvalues(Rationals,h2);
    esh02:=[];
    index_ev02:=Position(evh2,0);
    esh02:=esh2[index_ev02];

    #
    # sum and intersection
    #
    #return [ese,esh01,esh02];
    si1:=SumIntersectionMat(Basis(ese),Basis(esh01));
    si2:=SumIntersectionMat(Basis(ese),Basis(esh02));
    si:=SumIntersectionMat(si1[2],si2[2]);
    inter:=si[2];
    return inter;
    #return [inter,esh01,esh02];
end;