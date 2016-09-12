%{
Reaction Diffusion solver (Forward time, center space)
Matt Bovyn
September 6 2016
%}

function out=gierer_meinhardt_u(u,v,u_bar,a,alpha)
%Evaluates the function part of the Gierer-Meinhardt equation for v

out=a*(u.^2./v)+u_bar-alpha*u;

end