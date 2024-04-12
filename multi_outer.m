%%% note that the first two arguments to this function should be 
%%% mptPolytopes otherwise the intersection operation (i.e. and) is not
%%% exact in CORA 2021
%%% this function builds a box-aligned overapproximation to polytopes by
%%% using Algorithm 6 in the paper Inner and outer approximations of polytopes 
%%% using boxes by Alberto Bemporad, Carlo Filippi, Fabio D. Torrisi

function out = multi_outer(P, originalmpt, epsilon)

B = interval( P );

if( volume( B ) < epsilon || in( originalmpt, B ) )
    out = B;
else
    u = max( vertices( B ), [], 2 );
    l = min( vertices( B ), [], 2 );
    [ ~, k ] = max( u - l );
    g  = (u( k ) + l( k ))/2;
    
    u1 = u; u2 = u;
    l1 = l; l2 = l;
    
    l1( k ) = g;
    u2( k ) = g;
    
    Q1 = and( P, interval( l1, u1 ) ); %% intersection
    Q2 = and( P, interval( l2, u2 ) ); %% intersection
    
    out = [ multi_outer( Q1, originalmpt, epsilon ), ...
        multi_outer( Q2, originalmpt, epsilon ) ]; %% note we do not union the output because that would throw away the interval info
end
    