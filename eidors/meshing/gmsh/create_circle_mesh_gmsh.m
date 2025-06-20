function mdl = create_circle_mesh_gmsh(basename, center, rad, elem_size )
% Create a 2D Circular FEM using GMSH
% mdl= CREATE_GMSH_2D_CIRCLE(rad, n_elec)
%
% mdl - EIDORS forward model
% rad - model radius

% (C) 2009 Bartosz Sawicki. License: GPL version 2 or version 3
% $Id: create_circle_mesh_gmsh.m 6929 2024-05-31 15:36:50Z bgrychtol $

if ischar(basename) && strcmp(basename,'UNIT_TEST'),do_unit_test,return,end

geo_filename = sprintf('%s.geo', basename);
geo_fid= fopen(geo_filename,'w');

point_no = 1;
% Points to define circle
center_no = point_no;
fprintf(geo_fid,'Point(%d) = {%f, %f, 0, %f};\n',center_no, ...
    center(1), center(2), 1);
point_no = point_no + 1;
inpoint1_no = point_no;
fprintf(geo_fid,'Point(%d) = {%f, %f, 0, %f};\n',inpoint1_no, ...
    center(1)+rad, center(2), elem_size);
point_no = point_no + 1;
inpoint2_no = point_no;
fprintf(geo_fid,'Point(%d) = {%f, %f, 0, %f};\n',inpoint2_no, ...
    center(1)-rad, center(2), elem_size);
point_no = point_no + 1;


line_no = 1;
% Internal circle and loop
circle1_no = line_no;
line_no = line_no + 1;
fprintf(geo_fid,'Circle(%d) = {%d, %d, %d};\n', circle1_no, inpoint1_no, ...
    center_no, inpoint2_no);
circle2_no = line_no;
line_no = line_no + 1;
fprintf(geo_fid,'Circle(%d) = {%d, %d, %d};\n', circle2_no, inpoint2_no, ...
    center_no, inpoint1_no);

loop_no = line_no;
line_no = line_no + 1;
fprintf(geo_fid,'Line Loop(%d) = {%d,%d};\n', loop_no, circle1_no, ...
    circle2_no);

fprintf(geo_fid, 'Plane Surface(%d) = {%d};\n', line_no, loop_no);
line_no = line_no + 1;

fprintf(geo_fid, 'Mesh 2;');

fclose(geo_fid);

% Call Gmsh 
disp('Calling Gmsh. Please wait ...');
call_gmsh( geo_filename);

msh_filename = sprintf('%s.msh', basename);

disp(['Now reading back data from file: ' msh_filename])

[srf,vtx,fc,bc,simp,edg,mat_ind] = gmsh_read_mesh('circle.msh');
mdl.type     = 'fwd_model';
mdl.name = 'GMSH Circle';
mdl.nodes    = vtx;
mdl.elems    = simp;
mdl= eidors_obj('fwd_model', mdl);


function do_unit_test
object_center = [0.0, 0.1];
object_radius = 0.6;

mdl = create_circle_mesh_gmsh('circle', object_center, object_radius, 0.04 );
show_fem(mdl);
