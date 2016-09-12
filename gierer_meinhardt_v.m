%{
Reaction Diffusion solver (Forward time, center space)
Matt Bovyn
September 6 2016
%}

function out=gierer_meinhardt_v(u,v,v_bar,a,beta)
%Evaluates the function part of the Gierer-Meinhardt equation for v

out=a*u.^2+v_bar-beta*v;

end