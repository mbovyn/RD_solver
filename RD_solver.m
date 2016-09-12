%{
Reaction Diffusion solver (Forward time, center space)
Matt Bovyn
September 6 2016

This code implements a FTCS scheme to solve the Gierer-Meinhardt
equations. It calls three functions. The first creates the FTCS matrix and
is called FTCS_matrix.m. The second and third are gierer_meinhardt_u.m and 
gierer_meinhardt_v.m which return values of the functions associated with
the system.
%}

%% Set Parameters

%--------------------------------------------------------------------------
% user inputs

%number of points in space (stay <=100 or the matrix creation takes too
%long. It's not coded efficiantly.)
N=100;
%number of points in time
n_t=20000;

%choose whether or not we want to see each time step as the code is
%running. Makes things very slow if we do
animate=0;

%choose which parameters we want to use, as it was useful to keep several
%sets handy
%set 2 makes spots
params=1;

%--------------------------------------------------------------------------

%forms a labarynth pattern
if params==1;
    %diffusion
    nu_u=.05;
    nu_v=.5;

    %function
    a=1000;
    u_bar=1;
    v_bar=0;
    alpha=27;
    beta=40;
end

%forms spots
if params==2
    %diffusion
    nu_u=.05;
    nu_v=.5;

    %function
    a=5;
    u_bar=.5;
    v_bar=0;
    alpha=45;
    beta=50;
end

%doesn't form patterns
if params==3
    %diffusion
    nu_u=.5;
    nu_v=.5;

    %function
    a=5;
    u_bar=1;
    v_bar=0;
    alpha=1;
    beta=50;
end

%forms few, small spots
if params==4
    %diffusion
    nu_u=.05;
    nu_v=.5;

    %function
    a=-5;
    u_bar=.5;
    v_bar=0;
    alpha=55;
    beta=20;
end

%% set inital condition
%choose the intial conditions we want to use as well as if we'd like to add
%the reaction into the solver

%diff_only==1 used for debugging purposes, doesn't include the reaction
%term
diff_only=0;

%random perterbation of steady state if we are using the reaction
if diff_only==0
    %steady state values for Gierer Meinhardt system
    u_ss=(beta+u_bar)/alpha;
    v_ss=(a/beta)*((beta+u_bar)/alpha)^2;

    %set them as the initial conditions with 10% initial perterbation
    ic_u=u_ss*ones(N^2,1);
    ic_u=abs(ic_u+.1*randn(size(ic_u)));
    ic_v=v_ss*ones(N^2,1);
    ic_v=abs(ic_v+.1*.1*randn(size(ic_v)));
end

%to test diffusion only
if diff_only==1
    ic_u=zeros(N^2,1);
    ic_u(5050)=100; %high at one point in the center
    ic_v=zeros(N^2,1);
    ic_v(5050)=1; %high at one point in the center
end


%% set time and space grids
%solver only implemented for periodic boundary conditions

%spatial grid
side_length=2*pi;
space_step=side_length/N;

%calculate dt
%mu's are forced to be less than 1/2 for stability
% dt_u=space_step^2/(2.1);
% dt_v=space_step^2/(2.1);
% dt=min([dt_u;dt_v]);

%mu=1/2 turned out to be too big with reaction term
dt=space_step^2/5;

%calculate mu's for u and v (seperate grid spacings not implemented)
mu_u_x=nu_u*dt/space_step^2;
mu_u_y=nu_u*dt/space_step^2;

mu_v_x=nu_v*dt/space_step^2;
mu_v_y=nu_v*dt/space_step^2;

%intialize matricies
u=ic_u; %set inital conditions
v=ic_v; 

clear ic_u ic_v %save memory by clearing unnessicary variables

%% create FTCS matricies using FTCS_matrix function written for this


%for both u and v, create a sparse matrix which encodes the scheme
%sparsity is very important here as it makes the computation very fast
M_u=FTCS_matrix(N,mu_u_x,mu_u_y);
M_v=FTCS_matrix(N,mu_v_x,mu_v_y);


%% Run solver forward through time

%create a figure if we are animating
if animate==1
    figure(1)
end

%for each time step
for n=1:n_t
    
    %if we are using the full scheme, we multiply the matrix we created
    %with the vector of values at the previous time point and add the
    %function evaluation
    if diff_only==0
        u=M_u*u...
            +dt*gierer_meinhardt_u(u,v,u_bar,a,alpha);
        v=M_v*v...
            +dt*gierer_meinhardt_v(u,v,v_bar,a,beta);
    end
    
    %if we are using only diffusion, we simply perform the matrix
    %multiplication
    if diff_only==1
        u=M_u*u;
        v=M_v*v;
    end
    
    %if we want to see the animation while the simulation is running, draw
    %the patterns of both species
    if animate==1 && mod(n,500)==0
        subplot(1,2,1)
        imagesc(reshape(u,[N,N]))
        title1=sprintf('Concentration of Species 1\nt=%g',n);
        title(title1)
        xlabel('Space dimension 1')
        ylabel('Space dimension 2')
        colorbar
        subplot(1,2,2)
        imagesc(reshape(v,[N,N]))
        title2=sprintf('Concentration of Species 2\nt=%g',n);
        title(title2)
        xlabel('Space dimension 1')
        ylabel('Space dimension 2')
        colorbar
        drawnow
    end
    
    %display counters to make sure we are running
    if mod(n,1000)==0
        disp(n)
    end
end

%% Display results

%set if we want to display the last frame for a quick check
last_frame=0;

%if we want the last frame displayed, do so for both species
if last_frame==1
    
    figure(4)
    
    subplot(1,2,1)
        imagesc(reshape(u,[N,N]))
        title1=sprintf('Concentration of Species 1\nt=%g',n);
        title(title1)
        xlabel('Space dimension 1')
        ylabel('Space dimension 2')
        colorbar
        subplot(1,2,2)
        imagesc(reshape(v,[N,N]))
        title2=sprintf('Concentration of Species 2\nt=%g',n);
        title(title2)
        xlabel('Space dimension 1')
        ylabel('Space dimension 2')
        colorbar
end

%% surf animation - now defunct

%I had implemented animation before, but changed the scheme to be more
%efficiant. If you want to have a movie of the concentrations changing, you
%could save u and v to u_record and v_record every 100 time steps or
%something like that, then this code would animate it for you.

%{
%decide if we want the time evolution of either species animated on a
%surface plot
surf_animate_u=0;
surf_animate_v=0;

%if we do, animate the moving surf plot
if surf_animate_u==1
    
    figure(2)

    %set initial frame
    sh=surf(reshape(u_record(:,1),[N,N]));
    set(gca,'zlim',[min(min(u_record(:,50:end))) ...
        max(max(u_record(:,50:end)))])

    %update for subsequent time steps
    for n=2:n_t
        set(sh,'zdata',reshape(u_record(:,n),[N,N]))
        pause(.05) %slow things down or it goes too fast
    end
end

%the same for the other species
if surf_animate_v==1
    
    figure(3)

    sh=surf(reshape(v_record(:,1),[N,N]));
    set(gca,'zlim',[min(min(v_record(:,50:end))) ...
        max(max(v_record(:,50:end)))])

    for n=2:n_t
        set(sh,'zdata',reshape(v_record(:,n),[N,N]))
        pause(.05)
    end
end

%}

%% Display the steady state results on a surface plot

figure(5)
label=linspace(0,side_length,N);
surf(label,label,reshape(u,[N,N]),'linestyle','none')
xlabel('Space dimension 1')
ylabel('Space dimension 2')
zlabel('Concentration of Species 1')
axis tight
title('Final state of Species 1')

%% Display the steady state results on a contour plot

figure(6)
contourf(label,label,reshape(u,[N,N]),'linestyle','none')
h=colorbar;
xlabel('Space dimension 1')
ylabel('Space dimension 2')
ylabel(h,'Concentration of Species 1')
title('Final state of Species 1')