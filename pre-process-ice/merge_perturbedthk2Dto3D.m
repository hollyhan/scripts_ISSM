% This script prescribes perturbed ice thickness on local 2D-projected domain
% onto the global 3D domain. For the script to work correctly, one needs to
% check

addpath('/Users/kyhan/Desktop/LIA-Project/Pre-process/') % where the code 'perturb_thk_GRF.m' is saved 
addpath('/Users/kyhan/Desktop/ISSM_Tutorial_Exercise/LIA-thickness/');
addpath('/Users/kyhan/Desktop/ISSM_Tutorial_Exercise/SolidEarthModule-Tutorial/SLC_tutorials/GIA_simple_run/');	
addpath('/Users/kyhan/Desktop/ISSM_Tutorial_Exercise/LIA-thickness/mesher_and_coastlines');
cmap_b2r = load('/Users/kyhan/Desktop/Data/others/blue_white_red_colormap.mat').cmap;

% load model files
% md_global_ready_for_solve.mat already has the initial ice thickness set
% up. This script simply replaces the spcthickness field in the file with
% the new (perturbed) ice thickness data
md_global = loadmodel('/Users/kyhan/Desktop/LIA-Project/Pre-process/md_global_ready_for_solve.mat');
md_local = loadmodel('/Users/kyhan/Desktop/ISSM_Tutorial_Exercise/LIA-thickness/750CE_Initial_GrISmodel.mat'); % output 750CE is on the 3D (curved) local field. Need to project onto 2D field because I'll be using `transfer` values from the 2D field back onto 3D global field. 
md_world = loadmodel('/Users/kyhan/Desktop/LIA-Project/Pre-process/md_world.mat');

ice_thickness = load('Thickness_1050CE_2014.mat').thickness; % original file is struct with fields 
time = load('/Users/kyhan/Desktop/ISSM_Tutorial_Exercise/LIA-thickness/Time_1050CE_2014.mat').time; 
load('/Users/kyhan/Desktop/LIA-Project/Pre-process/transfer.mat');

% perturbed the ice thickness field from ISSM
disp('Perturbing the given ice thickness field by generating a Gaussian Random Field')
perturbed_ice_thk = perturb_thk_GRF(ice_thickness);

% transfer data from the local mesh onto the global merged mesh	
disp('Projecting initial model file onto 2D')
projected_md_local_thk = project2d(md_local, md_local.geometry.thickness,1); %this is needed only when inputs are not already projected on 2D field


%Define loading time series
disp('Assigning new ice thickness time series')
md = md_global;
    
n_time=size(time,2);
for t=1:n_time
	for i=1:size(projected_md_local_thk) 
		md.masstransport.spcthickness(transfer(i+md_world.mesh.numberofvertices,2),t) = perturbed_ice_thk(i,t);
	end
end

for i=1:size(projected_md_local_thk)
    md.geometry.thickness(transfer(i+md_world.mesh.numberofvertices,2))= perturbed_ice_thk(i,1);
end

% Note that 'md.masstransport.spcthickness(end,:)' has time tags already
% from pre-processing of md_world file. Just need to make sure the
% timestamps in the field are exactly the same as `time` variable in this
% script. 
% i.e. sum(md_global.masstransport.spcthickness(end,:)-time)% this should be zero

figure;
plotmodel(md,'data',md.masstransport.spcthickness(1:end-1,1) - md.geometry.thickness); % initial ice thickness diff should be zero

figure;
plotmodel(md,'data',md.masstransport.spcthickness(1:end-1,1),'coastlines','on') % initial thickness

figure;
deltaIce = md.masstransport.spcthickness(1:end-1,150) - md.masstransport.spcthickness(1:end-1,1);
plotmodel(md,'data',deltaIce, 'coastlines','on','colormap', flipud(cmap_b2r) ,'caxis',[-200 200])% total ice thickness change

figure; % plot the perturbation
t = 150;
thk_perturbation = md_global.masstransport.spcthickness(1:end-1,t)-md.masstransport.spcthickness(1:end-1,t);
plotmodel(md,'data', thk_perturbation, 'coastlines','on','colormap', flipud(cmap_b2r) ,'caxis',[-200 200])

% save model
save md_global_ready_for_solve_perturbed_ice.mat md;