#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 18:58:41 2024

@author: legras

Matlab text of oscul.m

function [yp,ym] = oscul(chem,xx,p,off,ip)
%test that NaNs are not in chem because LLTd silently fails in this case
if (any(isnan(chem)))
    disp('NaN in data, osculation fails');
end
if (nargin==3)
    off=0; ip=0;
end
if length(xx) ~= length(chem)
    disp('chem and xx not same length');
end
[Sp,Cp] = LLTd(xx,0.5 * xx .* xx - p * (chem-off),xx);
[Sm,Cm] = LLTd(xx,0.5 * xx .* xx + p * (chem+off),xx);
yp = (- 0.5 * xx .* xx + Cp)/p;
ym = (  0.5 * xx .* xx - Cm)/p;
if (ip==1)
   figure;
   plot(xx,chem,'-r',xx,yp,'-g',xx,ym,'-g');
end
return


"""
import numpy as np
from numba import jit

def oscul(yy, xx, p, off, function = '+-'):
     
    if len(yy) != len(xx):
        print('Vectors of unequal length')
        return -1

    if '+' in function:
        [Sp,Cp] = LLTd(xx,0.5 * xx**2 - p * (yy-off),xx)
        yp = (- 0.5 * xx**2 + Cp)/p
    if '-' in function:
        [Sm,Cm] = LLTd(xx,0.5 * xx**2 + p * (yy+off),xx)
        ym = (  0.5 * xx**2 - Cm)/p
    
    match function:
        case '+': return yp
        case '-': return ym
        case '+-': return [yp,ym]
        case '-+': return [ym,yp]
        case _: return -1

def LLTd(X,Y,S):
    """
    % LLTd		Compute numerically the discrete Legendre transform
    %		of a set of planar points (X(i),Y(i)) at 
    %		slopes S, i.e., Conj(j)=max_{i}[S(j)*X(i)-Y(i)].
    %		It uses the Linear-time Legendre Transform algorithm
    %		to speed up computations.
    %		
    %		All slopes not needed to build the affine interpolation may be
    %		removed by using fusionca.m instead of fusionma.m in
    %		fusion.m but computation time will take longer (SS may
    %		then be smaller than S). 
    %
    % Usage:
    %  [SS,Conj]=LLTd(X,Y,S)
    %
    % Input parameters:
    %  X,Y - Column vectors. Usually Y(i)=f(X(i))
    %  S - Column vector. We want to know the conjugate at S(j)
    %
    %  
    % Output parameter:
    %  SS - Slopes needed to define the conjugate. Removing affine parts
    %	(by using fusionca.m instead of fusionma.m in fusion.m) may
    %	result in needing less slopes. This is a column vector. 
    %  Conj - numerical values of the conjugate at slopes SS. This is
    %	a column vector.
    %
    % Functions called:
    %  bb - Compute the planar convex hull.
    %  fusion - sort two increasing sequences of slopes. This is the core
    %  of the LLT algorithm.
    %
    % Being called by:
    %  LFt - Compute the conjugate of a function.
    %
    % Example:
    %  X=[-5:0.5:5]'
    %  Y=X.^2
    %  S=(Y(2:size(Y,1))-Y(1:size(Y,1)-1))./(X(2:size(X,1))-X(1:size(X,1)-1))
    %  [SS Conj]=LLTd(X,Y,S)
    """
   # Compute the convex hull (input sensitive algorithm).
    [XX,YY] = bb(X,Y) 
    
    #h = len(XX) #m h=size(XX,1)
    # Compute the slopes associated with the primal points
    C = (YY[1:] - YY[:-1]) / (XX[1:] - XX[:-1])  			   
    # Merge sequences C and S
    [SS,H] = fusionma(C,S) #!! H is an index used to slice XX and YY and must 
    # be used here withe a shift -1 
    # speed (use fusionma.m) or to emphasize compression of output data (use
    # fusionca.m). 
    Conj = SS * XX[H-1] - YY[H-1]

    return [SS,Conj]

from numba import jit
@jit(nopython=True,cache=True)
def bb (X,Y):
    """
    % bb		Compute the planar convex hull of the set (X(i),Y(i)).
    %		By requesting X to be sorted we obtain a linear-time
    %		algorithm. 
    %		
    % Usage:
    %  [bbx,bby] = bb (X,Y)
    %
    % Input parameters:
    %  X,Y - Column vectors of same size. (X(i),Y(i)) is the planar set
    %	whose convex hull we compute. X must be sorted increasingly
    %	with distinct points: X(i)<X(i+1).
    %
    %  
    % Output parameter:
    %  bbx,bby - Verteces of the convex hull. Obviously bbx is a subset of
    %	X and bby a subset of Y. Column vectors.
    %
    %
    % Being called by:
    %  LLTd - Compute the discrete Legendre transform.
    %
    % Example:
    %  X=[-2:0.25:2]'
    %  Y=(X.^2-ones(size(X))).^2
    %  [bbx,bby]=bb(X,Y)
    %  plot(X,Y,'x');hold on;plot(bbx,bby,'o');
    %
    % bb uses the Beneat-Beyond algorithm [Preparata and Shamos,
    % Computational Geometry, Springer Verlag, 1990]
    """
    
    # Initialisation
    n = len(X) 
    v = 1 
    CX = np.empty(n) 
    CY = np.empty(n) 
    CX[:2] = X[:2] 
    CY[:2] = Y[:2] 
    
    # We look at each point only once, this gives the linear-time
    # algorithm.
    for i in range(2,n): 
        y=((CY[v]-CY[v-1])/(CX[v]-CX[v-1]))*(X[i]-CX[v])+CY[v]
        while (v > 1) & (Y[i] <= y): 
            # We erase points which are not vertices of the convex hull
            v -= 1
            y=((CY[v]-CY[v-1])/(CX[v]-CX[v-1]))*(X[i]-CX[v])+CY[v]
        if v>1: 
            if Y[i] == y: 
                CX[v] = X[i]
                CY[v] = Y[i]
            else:
                CX[v+1] = X[i]
                CY[v+1] = Y[i]
                v += 1
        else: 
            if Y[i] > y: # Trivial convex hull
                CX[2] = X[i]
                CY[2] = Y[i]
                v = 2
            else: 
                CX[1] = X[i]
                CY[1] = Y[i] 
    return   [CX[:v+1],CY[:v+1]]

def fusionma(C,S):
    """
    % fusionma	merge two increasing sequences.
    %		fusionma(C,S) gives the resulting slopes fS and
    %		indices fH where each slope support the epigraph. 
    %		All in all, it amounts to finding the first indice i
    %		such that C(i-1)<=S(j)<C(i)
    %
    % 		This function uses matlab syntax to speed up
    %		computation, see fusionca.m for a more classical
    %		programming of the same function
    % 
    %		All slopes not needed to build the affine interpolation may be
    %		removed by using fusionca.m instead of fusionma.m but
    %		computation time will take longer (SS may then be
    %		smaller than S).  
    %
    % Usage:
    %  [fS,fH]=fusionma(C,S)
    %
    % Input parameters:
    %  C,S - Increasing column vectors. We want to know the conjugate at
    %	S(j). C gives the slope of planar input points.
    %
    %  
    % Output parameter:
    %  fS - Slopes supporting the epigraph. This is a subset of S. Column vector.
    %  fH - Index at which the slope support the epigraph. Column vector.
    %
    %
    % Being called by:
    %  fusion - Choose either of fusionca.m or fusionma.m
    %
    % Example:
    %  X=[-5:0.5:5]'
    %  Y=X.^2
    %  C=(Y(2:size(Y,1))-Y(1:size(Y,1)-1))./(X(2:size(X,1))-X(1:size(X,1)-1))
    %  S=C-0.25*ones(size(C))
    %  [fS,fH]=fusionma(C,S)
    """
    
    # The BUILT-IN sort function is faster than
    #looking at each slope with a for loop 
    I = np.argsort(np.concatenate((C,S))) + 1
    I = np.array(np.where(I>(len(C)))[0])
    J = I[1:] - I[:-1] -1
    L = np.cumsum(J)
    # This (may be) obscure computation ONLY aims at using Matlab syntax
    # to do exactly the same as fusionca.m whose programming sticks to the
    # LLT as written in ["Faster than the Fast Legendre Transform, the
    # Linear-time Legendre Transform", Y. Lucet]
    return [S, np.insert(L+I[0]+1,0,I[0]+1)]
    # This return matches exactly fusionca when there is no elimination

