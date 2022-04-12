%--------------------------------------------------------------------------
% Human branching cholangiocyte organoids recapitulate functional bile duct 
% formation. 
% Roos et al. Cell Stem Cell (2022)
%-------------------------------------------------------------------------- 
% Simulation of the three-dimensional (3d) branching-annihilating random walk
% considering: tip-duct repulsion, a refractory period for branching and
% terminations, and a typical tip size at branching. For details, see
% Methods in paper.
% 
% INPUTS: see MODEL PARAMETERS below, and Methods. The default parameters 
% provided, correspond to the parameters used to simulate the BRCO condition. 
%
% OUTPUTS:
% edge_list: Array of size number_of_links*3, containing the source node id (column 1), 
%            target node id (column 2), and distance from source to target (column 3).
%            For all nodes in the simulation.
% node_positions : Array of size number_of_nodes*4. Column 1: node id. 
%            Columns 2-4: coordinates (x,y,z) for each node, respectively.  
% edge_list_tree: Reduced edge_list, where each node corresponds to either
%                 a branching point or a termination point or the root node
% node_positions_tree : node positions for the reduced branching tree.
%
% For questions: ib443 (at) cam.ac.uk
%-------------------------------------------------------------------------- 
%% Working directory 
% define a working directory to save output files (only used if save_files == 1)
% in the form 'D:\Documents\BARW_model\sim_output'
work_path = 'D:\Documents\BARW_model\sim_output';
save_files = 0; % set to 1 if want to save output files in work_path folder
plot_over_time = 1; % set to 1 if want to plot simulation (set n_reps to 1)
%--------------------------------------------------------------------------
%% SIMULATION PARAMETERS
n_reps = 1; % number of realisations
avoid = 1; % set to 1 to enable tip-duct repulsion (0 for no repulsion)
% MODEL PARAMETERS
sigma = 0.04; % branching to elongation rate ratio
noise_amplitude = 0.37; % noise aplitude for the active tips
annihil_radius = 16; % R_a : annihilation radius
t_max = 192;% maximum simuation time 
max_refractory_time = 10; % tip refractory period
repulsion_strength = 0.2; % tip-duct repulson strength
repulsion_scale = 1.2; % distance for tip-duct repulsion = repulsion_scale*annihil_radius
% Branching angle distribution: truncated Gaussian
pretruncMean = 90*pi/180;
pretruncSD = 50*pi/180;
untruncated = makedist('Normal',pretruncMean,pretruncSD);
branch_angle_dist = truncate(untruncated,80*pi/180,100*pi/180); 
% distribution of typical tip sizes at branching:
tip_size_dist = makedist('Uniform',1,8);
%--------------------------------------------------------------------------
%% STORE PARAMETERS if save_files == 1
repulsion_thresh = annihil_radius.*repulsion_scale;
if save_files == 1
    sim_params = containers.Map;
    sim_params('noise_amplitude') = noise_amplitude;
    sim_params('sigma') = sigma;
    sim_params('annihil_radius') = annihil_radius;
    sim_params('n_reps') = n_reps;
    sim_params('T_max') = t_max;
    sim_params('max_refractory_time') = max_refractory_time;
    sim_params('repulsion_thresh') = repulsion_thresh;
    sim_params('repulsion_strength') = repulsion_strength;
    [exi,save_path] = create_save_folder(work_path,sim_params);
else
    exi = 0;
end
%--------------------------------------------------------------------------
%% RUN simulation if the parameter set has not been simulated yet
if exi == 0
    if save_files == 1
        save([save_path,'sim_params.mat'],'sim_params')
    end
    % Run n_reps of the process.
    % NOTE: change this "for" by "parfor" to run
    % realisations in parallel.
    parfor n_rep = 1:n_reps
        % initialise the system with a "rod" with two active tips
        [edge_list,node_positions,edge_list_clean,node_positions_clean,tot_part_number,n_active_part,active_part_number,active_part_pos,active_part_angle,refractory_time] = initial_condition();
        % initial time
        t = 0;
        
        if plot_over_time == 1
            figure;
            plot_3d_network(edge_list,node_positions);
            axis equal; set(gca,'Visible','off'); set(gcf,'color','w');
            hold on;
        end
        % run while t < t_max and there are active particles in the system
        while t < t_max && n_active_part>0
            % update time and refractory times
            t = t + 1;
            refractory_time(:,2) = refractory_time(:,2)+1;
            % random ordering of tip particles to perform action
            part_order = randperm(n_active_part);
            % draw random number to choose which action to perform
            ran = rand(n_active_part,1);
            % reorder particles in list
            active_part_number = active_part_number(part_order);
            active_part_pos = active_part_pos(part_order,:);
            active_part_angle = active_part_angle(part_order,:);
            % iterate through all active particles
            for part_num = n_active_part:-1:1
                % choose a tip cell to elongate or branch
                rt = refractory_time(part_num,2);
                % tip elongates if random number is < 1/(1+sigma) and
                % refractory period has no elapsed for the given particle 
                if ran(part_num) <= 1/(1+sigma) || rt < max_refractory_time % tip elongation
                    part_pos = active_part_pos(part_num,2:end);
                    part_angle = active_part_angle(part_num,2:end);
                    % update particle position
                    if avoid == 0
                        [new_part_pos, new_part_angle] = update_particle_position(part_pos,part_angle,dr,noise_amplitude);
                        % check ditance to trace if refractory period has elapsed
                        if rt < max_refractory_time
                            h = 0;
                        else
                            h = check_distance_to_trace(part_pos,new_part_pos,annihil_radius,node_positions);
                        end
                    elseif avoid == 1
                        [new_part_pos, new_part_angle,h] = update_particle_position_avoid(node_positions,part_pos,part_angle,noise_amplitude,annihil_radius,repulsion_thresh,repulsion_strength);
                        if rt < max_refractory_time
                            h = 0;
                        end
                    end
                    % add particle and update data
                    [edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,tot_part_number,refractory_time] = add_particle(h,edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,tot_part_number,part_num,new_part_pos,new_part_angle,refractory_time);
                    % plot if plot_over_time == 1
                    if plot_over_time == 1
                        pause(0.01)
                        if h
                            plot([part_pos(1),new_part_pos(1)],[part_pos(2),new_part_pos(2)],'-k','linewidth',2);
                            scatter(new_part_pos(1),new_part_pos(2),20,'sr','filled');
                        else
                            plot([part_pos(1),new_part_pos(1)],[part_pos(2),new_part_pos(2)],'-k','linewidth',2);
                        end
                    end
                    
                % tip branches otherwise    
                else % branch
                    part_pos = active_part_pos(part_num,2:end);
                    part_angle = active_part_angle(part_num,2:end);
                    % draw a random number from truncated gaussian
                    ang = random(branch_angle_dist);
                    % branch at equal opposite angles:
                    branch_angles = [ang/2 -ang/2];
                    % random rotation angle
                    rotate_angle = pi*rand();
                    % draw a tip size
                    tip_size = random(tip_size_dist);
                    % find offspring tip position
                    [new_part_pos1,new_part_angle1,new_part_pos2,new_part_angle2] = branch_tip(part_pos,part_angle,rotate_angle,branch_angles,tip_size);
                    
                    [edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,tot_part_number,refractory_time] = add_two_particles(edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,tot_part_number,part_num,new_part_pos1,new_part_angle1,new_part_pos2,new_part_angle2,refractory_time);
                    
                    if plot_over_time == 1
                        pause(0.01)
                        plot([part_pos(1),new_part_pos1(1)],[part_pos(2),new_part_pos1(2)],'-k','linewidth',2);
                        plot([part_pos(1),new_part_pos2(1)],[part_pos(2),new_part_pos2(2)],'-k','linewidth',2);
                    end
                end
            end
        end
        % save results
        if save_files == 1
            [edge_list_large,node_positions_large] = clean_edge_list(edge_list_clean,node_positions_clean);
            parsave([save_path,'edge_list_',num2str(n_rep.','%04d'),'.mat'],edge_list)
            parsave([save_path,'node_positions_',num2str(n_rep.','%04d'),'.mat'],node_positions)
            parsave([save_path,'edge_list_tree_',num2str(n_rep.','%04d'),'.mat'],edge_list_large)
            parsave([save_path,'node_positions_tree_',num2str(n_rep.','%04d'),'.mat'],node_positions_large)
        end
    end
end

%% FUNCTION DEFINITIONS

function [edge_list,node_positions,edge_list_clean,node_positions_clean,tot_part_number,n_active_part,active_part_number,active_part_pos,active_part_angle,refractory_time] = initial_condition()
% construct the initial condition
init_pos = zeros(1,3); init_pos(1) = 0;
n_active_part = 2;
% active particle positions:
active_part_pos = [2,1/2,0,0];
active_part_pos(end+1,:) = [3,-1/2,0,0];
active_part_angle = [2,pi/2,0];
active_part_angle(end+1,:) = [3,pi/2,pi];
% total number of particles
tot_part_number = 3;
% idx of ative tips:
active_part_number = [2,3];
% network info:
edge_list = [1 2 1/2;1 3 1/2];
node_positions = [1,init_pos];
node_positions(end+1,:) = [2,1/2,0,0];
node_positions(end+1,:) = [3,-1/2,0,0];
edge_list_clean = edge_list;
node_positions_clean = node_positions;
refractory_time = [2,0;3,0];
end

function [new_pos,new_angle] = update_particle_position(part_pos,part_angle,noise_amplitude)
% update particle position, considering a random noise (on the sphere) with
% amplitude noise_amplitude
part_angle(1) = part_angle(1) + 2*(rand()-0.5)*noise_amplitude;
part_angle(2) = part_angle(2) + acos(1 - 2*rand())*noise_amplitude;
new_pos = part_pos + [sin(part_angle(1))*cos(part_angle(2)) sin(part_angle(1))*sin(part_angle(2)) cos(part_angle(1))];
new_angle = part_angle;
end

function [new_part_pos,new_part_angle,h] = update_particle_position_avoid(node_positions,part_pos,part_angle,noise_amplitude,annihil_radius,repulsion_thresh,repulsion_strength)
% update particle position, considering a random noise (on the sphere) with
% amplitude noise_amplitude, and considering tip-duct repulsion, which acts
% over a radius annihil_radius*repulsion_thresh with "strength" repulsion_strength
new_part_angle(1) = part_angle(1) + 2*(rand()-0.5)*noise_amplitude;
new_part_angle(2) = part_angle(2) + 2*(rand()-0.5)*noise_amplitude;
new_part_pos = part_pos + [sin(new_part_angle(1))*cos(new_part_angle(2)) sin(new_part_angle(1))*sin(new_part_angle(2)) cos(new_part_angle(1))];
% find tracers that act on the moving tip
trace_pos = node_positions(:,2:end);
dists = pdist2(new_part_pos,trace_pos);
v = trace_pos - new_part_pos;
w = new_part_pos - part_pos;
w = repmat(w,size(v,1),1);
angles = atan2(vecnorm(cross(v,w,2),2,2),dot(v,w,2))';
ind = find(dists < repulsion_thresh & abs(angles) < pi/4);
h = any(dists < annihil_radius & abs(angles) < pi/4,'all');
% lets move away from the centre of mas of nearby tracers
numel(ind)
if h == 0 && numel(ind) > 0
    in_trace = trace_pos(ind,:);
    cm = mean(in_trace,1);
    vect = cm - part_pos; vect = vect/norm(vect);% vector connecting CM and part_pos
    % now we push the new part a quantity fr*(-vect)
    new_part_pos = new_part_pos - repulsion_strength.*vect;
    v1 = new_part_pos - part_pos;
    % compute new particle angle
    new_part_angle = [atan2(sqrt(v1(1)^2 + v1(2)^2),v1(3)), atan2(v1(2),v1(1))];
end
end
                      
function [new_pos1,new_angle1,new_pos2,new_angle2] = branch_tip(part_pos,part_angle,rotate_angle,branch_angles,tip_size)
% create two new branches that form and angle branch_angles between them
% and they are revolved and angle rotate_angle around the axis formed by
% the direction of the parent tip.
% normalised reference axis:
ref_ax = [sin(part_angle(1))*cos(part_angle(2)) sin(part_angle(1))*sin(part_angle(2)) cos(part_angle(1))];
% random rotation angle
%     rotate_angle = pi*rand();
% branch angles
theta1 = part_angle(1)+branch_angles(1);
theta2 = part_angle(1)+branch_angles(2);
% generate new branches by rotating the branches around the axis formed by
% the parent branch.
pos = tip_size*[sin(theta1)*cos(part_angle(2)) sin(theta1)*sin(part_angle(2)) cos(theta1)];
v1 = rotate_3D(pos', 'any', rotate_angle, ref_ax')';
pos = tip_size*[sin(theta2)*cos(part_angle(2)) sin(theta2)*sin(part_angle(2)) cos(theta2)];
v2 = rotate_3D(pos', 'any', rotate_angle, ref_ax')';

new_angle1 = [atan2(sqrt(v1(1)^2 + v1(2)^2),v1(3)), atan2(v1(2),v1(1))];
new_angle2 = [atan2(sqrt(v2(1)^2 + v2(2)^2),v2(3)), atan2(v2(2),v2(1))];
new_pos1 = v1 + part_pos;
new_pos2 = v2 + part_pos;
end

function [edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,tot_part_number,refractory_time] = add_particle(h,edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,tot_part_number,part_num,new_part_pos,new_part_angle,refractory_time)
% add particle to the list, and update all data
if h
    active_part_pos(part_num,:) = [];
    active_part_angle(part_num,:) = [];
    active_part_number(part_num) = [];
    refractory_time(part_num,:) = [];
    n_active_part = n_active_part - 1;
else
    tot_part_number = tot_part_number+1;
    edge_list(end+1,:) = [active_part_pos(part_num,1) tot_part_number 1];
    
    [r,c] = find(edge_list_clean == active_part_pos(part_num,1));
    
    d = pdist2(active_part_pos(part_num,2:4),new_part_pos);
    if numel(r) == 0
        edge_list_clean(end+1,:) = [active_part_pos(part_num,1) tot_part_number d];
        node_positions_clean(end+1,:) = [tot_part_number new_part_pos];
    else
        edge_list_clean(r,c) = tot_part_number;
        edge_list_clean(r,3) = edge_list_clean(r,3) + d;
        
        node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos];
    end
    
    active_part_pos(part_num,:) = [tot_part_number,new_part_pos];
    active_part_angle(part_num,:) = [tot_part_number,new_part_angle];
    node_positions(end+1,:) = [tot_part_number new_part_pos];
    active_part_number(part_num) = tot_part_number;
    
    refractory_time(part_num,1) = tot_part_number;
end
end

function [edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,tot_part_number,refractory_time] = add_two_particles(edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,tot_part_number,part_num,new_part_pos1,new_part_angle1,new_part_pos2,new_part_angle2,refractory_time)
% add two particles (at branching) to the list, and update all data
d1 = pdist2(active_part_pos(part_num,2:4),new_part_pos1);
d2 = pdist2(active_part_pos(part_num,2:4),new_part_pos2);

edge_list(end+1,:) = [active_part_pos(part_num,1) tot_part_number+1 d1];
edge_list(end+1,:) = [active_part_pos(part_num,1) tot_part_number+2 d2];
edge_list_clean(end+1,:) = [active_part_pos(part_num,1) tot_part_number+1 d1];
edge_list_clean(end+1,:) = [active_part_pos(part_num,1) tot_part_number+2 d2];

tot_part_number = tot_part_number+1;
active_part_pos(part_num,:) = [tot_part_number,new_part_pos1];
active_part_angle(part_num,:) = [tot_part_number,new_part_angle1];
node_positions(tot_part_number,:) = [tot_part_number new_part_pos1];
node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos1];
active_part_number(part_num) = tot_part_number;
refractory_time(part_num,:) = [tot_part_number,0];

tot_part_number = tot_part_number+1;
active_part_pos(end+1,:) = [tot_part_number,new_part_pos2];
active_part_angle(end+1,:) = [tot_part_number,new_part_angle2];
node_positions(tot_part_number,:) = [tot_part_number new_part_pos2];
node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos2];
active_part_number(end+1) = tot_part_number;
refractory_time(end+1,:) = [tot_part_number,0];

n_active_part = n_active_part + 1;

% add tracer particles in case the tip particles are at a distance >1
% from their parent tip (this is for the tip-duct sensing to work properly) 
r = 0;
while ~isempty(r)
    r = find(edge_list(:,3) > 1);
    n1 = edge_list(r,1);
    n2 = edge_list(r,2);
    d = edge_list(r,3);
    
    new_parts_num = [(tot_part_number + 1):(tot_part_number+length(r))]';
    tot_part_number = tot_part_number + length(r);
    edge_list(end+1:end+length(r),:) = [n1 new_parts_num d/2];
    edge_list(end+1:end+length(r),:) = [new_parts_num n2 d/2];
    pos_new = (node_positions(n1,2:4) + node_positions(n2,2:4))/2;
    node_positions(new_parts_num,:) = [new_parts_num pos_new];
    edge_list(r,:) = [];
end

end

function [edge_list_large,node_positions_large] = clean_edge_list(edge_list_large,node_positions_large)
% this is to put the nodes in edge_list_clean and node_positions_clean in
% sequencial order.
un_nodes = unique(edge_list_large(:,1:2));
complete_list = 1:max(un_nodes);
rm_nodes = setdiff(complete_list,un_nodes); % these are the nodes that must be deleted
for i = 1:numel(un_nodes)
    node_positions_large(un_nodes(i),1) = i;
    ind = find(edge_list_large(:,1) == un_nodes(i));
    edge_list_large(ind,1) = i;
    ind = find(edge_list_large(:,2) == un_nodes(i));
    edge_list_large(ind,2) = i;
end
node_positions_large(rm_nodes,:) = [];
end

function parsave(fname, x)
% for saving in parfor (parallel computing)
save(fname, 'x')
end

function h = check_distance_to_trace(part_pos,new_part_pos,annihil_radius,node_positions)
% check any tracer (duct) particle is at and angle pi/4 or lower with
% respect the the direction of motion of the tip and if their distance is
% lower than annihil_radius.
trace_pos = node_positions(:,2:end);
dists = pdist2(new_part_pos,trace_pos);
v = trace_pos - new_part_pos;
w = new_part_pos - part_pos;
w = repmat(w,size(v,1),1);
angles = atan2(vecnorm(cross(v,w,2),2,2),dot(v,w,2))';
h = any(dists < annihil_radius & abs(angles) < pi/4, 'all');
end

function [exi,save_path] = create_save_folder(work_path,newparams)
% this function check if a simulation with parameters given in "newparams"
% exists in the working directory work_path. If it does not, then it
% creates a folder to store the data.
exi = 0;
save_path =[];
dinfo = dir(work_path);
dinfo(ismember( {dinfo.name}, {'.', '..','.DS_Store'})) = [];
if ~ispc; slsh = '/'; else; slsh = '\'; end
if length(dinfo) == 0
    n_strPadded = sprintf( '%05d', 1);
    save_path = [work_path,slsh,n_strPadded,slsh];
    mkdir(save_path)
else
    for i = 1:length(dinfo)
        fold = dinfo(i).name;
        load([work_path,slsh,fold,slsh,'sim_params.mat']);
        ks = keys(sim_params);
        check_if_exists = [];
        for j = 1:length(ks)
            check_if_exists(j) = all(sim_params(ks{j}) == newparams(ks{j}));
        end
        if all(check_if_exists)
            disp(['Simulation data exists for these parameters: ',dinfo(i).name]);
            exi = 1;
            break
        end
    end
    n = 0;
    if exi == 0
        for i = 1:length(dinfo)
            n(i) = str2num(dinfo(i).name);
        end
        n_strPadded = sprintf( '%05d', max(n) + 1 );
        save_path = [work_path,slsh,n_strPadded,slsh];
        mkdir(save_path)
    end
end
end

function [R, Rm] = rotate_3D(V, mode, theta, u, angle_unit)
% rotate_3D : function to compute the rotation of a vector or an array of vectors in 2D or 3D space.
% Source:
% Nicosahedron (2022). Any 3D rotation (https://github.com/NicolasDouillet/rotate_3D/releases/tag/v2.6), GitHub. Retrieved April 12, 2022.
%
% Syntax
% R = rotate_3D(V, mode, theta);
% R = rotate_3D(V, mode, theta, u);
% R = rotate_3D(V, mode, theta, u, angle_unit);
% [R,Rm] = rotate_3D(V, mode, theta, u, angle_unit);
%
%
% Description
% R = rotate_3D(V, mode, theta) computes the vector R, which results
% from the rotation of V vector around one of the the basis vectors, which
% is choosen in the mode : 'x', 'y', or 'z'.
%
% R = rotate_3D(V, mode, theta, u) computes the vector R, which results
% from the rotation of V vector around u vector and of theta angle in radian.
%
% R = rotate_3D(V, mode, theta, u, angle_unit) uses angle_unit for theta
% unit (radian or degree).
%
% [R,Rm] = rotate_3D(V, mode, theta, u, angle_unit) also returns the
% rotation matrix.
%
% Important NB : in 2D -(xOy) plan- mandatory rotation axis is 'z'. It will
% be set as so by default if input is different. Also in 2D, in case u is missing it
% is automatically set to the origin [0,0]' by default.
%
% Input parsing
Ndim = size(V,1);
assert(nargin > 2, 'Not enough input arguments.');
assert(nargin < 6, 'Too many input arguments.');
assert(Ndim > 1 && Ndim < 4, 'Input argument V must have between one and three rows : 1 < size(V,1) <= 3.');
assert(strcmpi(mode,'x') || strcmpi(mode,'y') || strcmpi(mode,'z') || strcmpi(mode,'any'),...
       'Bad mode argument : mode must be a string in the set {''x'',''X'',''y'',''Y'',''z'',''Z'',''any'',''ANY''}.');
if nargin < 5    
    angle_unit = 'radian';    
    if nargin < 4        
        if Ndim == 2           
            u = [0,0]';    
        elseif Ndim == 3    
            switch mode   
                case {'x', 'X'}
                    u = [1 0 0]';  
                case {'y', 'Y'}
                    u = [0 1 0]';
                case {'z', 'Z'}
                    u = [0 0 1]';     
            end 
        end
    else
        assert(Ndim < 3 || ~strcmpi(mode,'any') || norm(u) > 0,'3D rotation axis u must not equal null vector.');  
    end
else
    assert(strcmpi(angle_unit,'radian') || strcmpi(angle_unit,'degree'),'angle_unit value must be either ''radian'' or ''degree''.');
    if strcmpi(angle_unit,'degree')
        theta = pi * theta / 180;
    end
end
% Body
% Rotation matrix construction and resulting rotated vector computation
switch Ndim
    case 2 % rotate around a point (2D vector) in (xOy) plan -> mandatory rotation axis is 'z' 
        Rm = [cos(theta) -sin(theta);
              sin(theta)  cos(theta)];
        W = V - u;
        R = Rm * W;
        R = R + u;                
    case 3
        switch mode
            case {'x', 'X'} % X axis rotation matrix ; u = i = [1 0 0]'
                Rm = [1          0           0;
                      0 cos(theta) -sin(theta);
                      0 sin(theta)  cos(theta)];
            case {'y', 'Y'} % Y axis rotation matrix ; u = j = [0 1 0]'
                Rm = [cos(theta)   0  -sin(theta);
                      0            1           0;
                      sin(theta)  0  cos(theta)];
            case {'z', 'Z'} % Z axis rotation matrix ; u = k = [0 0 1]'
                Rm = [cos(theta) -sin(theta) 0;
                      sin(theta)  cos(theta) 0;
                      0           0          1];
            case {'any', 'ANY'} % Any u axis rotation matrix
                u = u/norm(u);
                Rm = [u(1,1)^2+cos(theta)*(1-u(1,1)^2) (1-cos(theta))*u(1,1)*u(2,1)-u(3,1)*sin(theta) (1-cos(theta))*u(1,1)*u(3,1)+u(2,1)*sin(theta);
                      (1-cos(theta))*u(1,1)*u(2,1)+u(3,1)*sin(theta) u(2,1)^2+cos(theta)*(1-u(2,1)^2) (1-cos(theta))*u(2,1)*u(3,1)-u(1,1)*sin(theta);
                      (1-cos(theta))*u(1,1)*u(3,1)-u(2,1)*sin(theta) (1-cos(theta))*u(2,1)*u(3,1)+u(1,1)*sin(theta) u(3,1)^2+cos(theta)*(1-u(3,1)^2)];
            otherwise
                error('Bad mode argument : mode must be a string in the set {''x'',''X'',''y'',''Y'',''z'',''Z'',''any'',''ANY''}.');   
        end
        R = Rm * V;            
end
end % rotate_3D

function plot_3d_network(edge_list,node_positions)
% Plot the 3d network given by edge_list and node positions
hold on
for i=1:length(edge_list(:,1))
    n1 = edge_list(i,1);
    n2 = edge_list(i,2);
    r1 = node_positions(n1,2:4);
    r2 = node_positions(n2,2:4);
    plot3([r1(1) r2(1)],[r1(2) r2(2)],[r1(3) r2(3)],'-k','linewidth',1.5)
end
hold off
end

function [paths,dists,D,G] = node_level(edge_list,node_positions,source_node)
% Extract branching tree information:
% paths : path between source_node and all other nodes
% dists : distance from the source_node to every other node (in the case of
%         the branching tree, this equates to the generation of "level" in 
%         the tree.
% D     : degree of every node
% G     : graph object
G = graph(edge_list(:,1),edge_list(:,2));
paths = {};
dists = [];
for i = 1:size(node_positions,1)
    nod = node_positions(i,1);
    if nod~=i
        disp('Warning: incompatible node_positions file')
        break;
    end
    if isfinite(nod)
        [P,d] = shortestpath(G,source_node,nod);
        paths{i} = P;
        dists(i) = d;
    else
        paths{i} = nan;
        dists(i) = nan;
    end
end
D = degree(G);
end

function [qis,undists] = termination_prob_vs_level(dists,D)
% Extract the termination probabilities as fucntion of level, using the 
% distances (dists) and degrees (D) of each node, as produced by the
% function "node_level".
undists = unique(dists);
undists = undists(isfinite(undists));
qis = [];
for i = 1:numel(undists)
    di = undists(i);
    if di == 0
    else
        indi = find(dists == di);
        degs = D(indi);
        nd = numel(degs);
        d1 = sum(degs == 1);
        d3 = sum(degs == 3);
        qis(di) = d1/nd;
    end
end
end

function add_generations_to_tree(max_level)
% Add level lines to plot of the branching tree
grid on
yticks(1:(max_level+1))
yTickLabels = arrayfun(@num2str,max_level:-1:0,'uni',false);
yticklabels(yTickLabels)
box off
ax1 = gca;
ax1.XAxis.Visible = 'off';
ylabel('Tree level (i)')
end