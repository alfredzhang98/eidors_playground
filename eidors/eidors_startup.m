function eidors_startup( path_array )
% Script to start EIDORS
% Set path and variables correctly
% USAGE:
%   startup - setup basic eidors usage functions
%   startup( { dev directory paths })

% NOTE: this is a function, so that we don't put variables into the
% workspace

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: eidors_startup.m 6828 2024-04-28 09:02:21Z bgrychtol $

if nargin == 0
    path_array = {};
end

eidors_msg('=== SETTING UP EIDORS version %s ===',  ...
     eidors_obj('eidors_version'), 0);

HOMEDIR=strrep(pwd,'\','/');
ver = version_check;
archdir = set_paths(HOMEDIR, ver, path_array);
eidors_cache('init');
% make EIDORS less verbose by default
eidors_msg('log_level', 2); 
set_defaults(HOMEDIR);
print_welcome(HOMEDIR,archdir, ver);

function set_defaults(HOMEDIR)
    % default functions
    eidors_default('set','fwd_solve','fwd_solve_1st_order');
    eidors_default('set','calc_system_mat','system_mat_1st_order');
    eidors_default('set','calc_jacobian','jacobian_adjoint');
    eidors_default('set','inv_solve','inv_solve_diff_GN_one_step');
    eidors_default('set','calc_RtR_prior','prior_laplace');
    eidors_default('set','calc_R_prior','prior_laplace');
    eidors_default('set','fwd_solve','fwd_solve_1st_order');
    eidors_default('set','GREIT_desired_img','GREIT_desired_img_sigmoid');
    eidors_default('set','mk_coarse_fine_mapping','mk_analytic_c2f');
    %models are NOT normalized by default
    eidors_default('set','mdl_normalize',@(x) 0); 
    %netgen finelevel is coarse
    eidors_default('set','call_netgen_finelevel','');

    calc_colours('defaults'); % default calc_colours

    %  Set max cache size. Not completely sure about this
    %  but 1GB should be available in most modern machines
    eidors_cache('cache_size', 1024^3 );
    eidors_cache('boost_priority', 0 ); % set default priority

    % Set default model cache location
    mk_library_model('LIBRARY_PATH',[HOMEDIR, '/models/cache']);
    eidors_cache('cache_path',[HOMEDIR, '/models/cache']);

    eidors_cache('eidors_path',HOMEDIR);


function ver = version_check
    ver= eidors_obj('interpreter_version');

    if ver.isoctave
        if ver.ver < 3.008
            warning(['EIDORS REQUIRES AT LEAST OCTAVE V3.8.0\n' ...
                'Several functions may not work with your version']);
        end
        % EIDORS does this a lot, so we remove the warning
        warning('off','Octave:load-file-in-path') 
    else
        if ver.ver < 7.006
            warning(['EIDORS REQUIRES AT LEAST MATLAB V7.6.\n' ...
                'Several functions may not work with your version']);
        end
    end

function archdir = set_paths(HOMEDIR, ver,path_array)
    % we will be overloading built-in functions. Disable the warning.
    warning off MATLAB:dispatcher:nameConflict
    addpath( HOMEDIR );
    addpath([HOMEDIR, '/solvers']);
    addpath([HOMEDIR, '/solvers/inverse']);
    addpath([HOMEDIR, '/solvers/forward']);
    addpath([HOMEDIR, '/solvers/forward/tools-beta']); % Tools for MC's solvers which are not yet general
    addpath([HOMEDIR, '/algorithms']);
    addpath([HOMEDIR, '/interface']);
    addpath([HOMEDIR, '/models']);
    addpath([HOMEDIR, '/meshing']);
    addpath([HOMEDIR, '/meshing/netgen']);
    addpath([HOMEDIR, '/meshing/distmesh']);
    addpath([HOMEDIR, '/meshing/gmsh']);
    addpath([HOMEDIR, '/meshing/stl']);
    addpath([HOMEDIR, '/sample_data']);
    addpath([HOMEDIR, '/examples']);
    addpath([HOMEDIR, '/tools']);
    addpath([HOMEDIR, '/graphics/matlab']);
    addpath([HOMEDIR, '/graphics/vtk']);
    addpath(genpath([HOMEDIR, '/external'])); %add subdirectories
    addpath([HOMEDIR, '/deprecated']);
    addpath([HOMEDIR, '/overloads']);
    if ~ver.isoctave && ver.ver < 9.008
        addpath([HOMEDIR, '/overloads/sparse']);
    end
% This was not an overload, but a function missing until 2017b 
%     if ~ver.isoctave && ver.ver < 9.0073
%         addpath([HOMEDIR, '/overloads/vecnorm']);
%     end
    if ver.isoctave 
        warning off Octave:shadowed-function
        addpath([HOMEDIR, '/overloads/octave']);
        warning on Octave:shadowed-function
    end
    warning on MATLAB:dispatcher:nameConflict

    %addpath([HOMEDIR, '/tests']);

    DEVDIR = [HOMEDIR(1:find(HOMEDIR == '/',1,'last')-1) '/dev'];
    for i = 1:length(path_array)
        p = genpath([DEVDIR, '/', path_array{i}]);
        addpath(p);
    end

    % We need to add an architecture specific directory for mex files
    if ver.isoctave
        archdir= strcat('/arch/octave/',computer);
% Make sure that the archdir exists
        [status] = mkdir([HOMEDIR,archdir]);
    else
        % problem was that MSVC compiles depended on versions of MSVCRT1xx.dll.
	% which was not available with all windows/matlab versions.
	% New concept is to try building with MINGW compiler
        archdir= '/arch/matlab';
    end
    addpath([HOMEDIR, archdir]);
    fname = [HOMEDIR, archdir, '/eidors_var_id.', mexext];
    
    if ~exist(fname, 'file')
       eidors_msg('STARTUP: missing a required, pre-compiled mex file: eidors_var_id', 1);
       compile_mex(HOMEDIR,archdir,ver);
    end

    % check that the compiled mex file is newer than the source file
    srcf = strcat(HOMEDIR,'/arch/eidors_var_id.cpp');
    mexf = strcat(fname);
    if exist(srcf) == 2 && exist(mexf) == 3
        srcd=dir(srcf);
        mexd=dir(mexf);


        % We thank MATLAB for their version issues
        newer_src = false;
        try newer_src = datenum(srcd.date) > datenum(mexd.date);
        catch
           newer_src = srcd.datenum > mexd.datenum;
        end

        if newer_src
           warning('eidors_var_id.mex file is older than source file and should be recompiled.');
        end

        ok = eidors_var_id_ok;
        if newer_src || ~ok
           while 1
              if ~ok
                 resp = input('Would you like to compile now? [Y/n]: ','s');
              else
                 resp = input('Would you like to compile now? [y/N]: ','s');
                 if isempty(resp)
                    resp = 'n';
                 end
              end
              
              switch lower(resp)
                 case {'n', 'no'}
                    if ver.isoctave
                       eidors_msg([...
                          '  Please compile it using:\n'...
                          '    cd ',HOMEDIR,'/arch\n'...
                          '    mkoctfile -v --mex eidors_var_id.cpp\n'...
                          '    mkdir -p ..',archdir,'\n'...
                          '    mv *.mex ..',archdir ...
                          ],1);
                    else
                       eidors_msg([ ...
                          '  Please compile it using:\n'...
                          '     cd ',HOMEDIR,'/arch\n'...
                          '     mex "',HOMEDIR,'/arch/eidors_var_id.cpp"\n'...
                          '     mv *.mex* ',HOMEDIR,'/arch/matlab\n' ...
                          'If you have a 64 bit machine, please use "mex -largeArrayDims ..."' ...
                          ],1);
                    end
                    break;
                 case {'','y','yes'}
                   compile_mex(HOMEDIR,archdir,ver);
                   break;
              end
           end
        end
    end

function compile_mex(HOMEDIR,archdir, ver)
    eidors_msg('Attempting to compile eidors_var_id',2);
    c = computer;
    flags = [];

    if ver.isoctave
         curdir = cd;
         cd(sprintf('%s/arch',HOMEDIR));
         mex eidors_var_id.cpp
% We create the directory earlier
%        system(sprintf('mkdir -p ..%s',archdir));
%    Has to be absolute paths because Matlab coders are so stupid!!
         movefile(sprintf('%s/arch/*.mex',HOMEDIR), ...
                  sprintf('%s%s/',HOMEDIR,archdir));
         cd(curdir)
         return
    end
    
    if strcmp(c(end-1:end),'64')
       flags = '-largeArrayDims';
    end  
    cmd = sprintf('mex %s "%s/arch/eidors_var_id.cpp"', flags, HOMEDIR);
% it seems to be better to use matlabs mex, especially since
% there is a latex derivative called mex to interfere with us
    eval(cmd);

% the assholes at matlab don"t respect the 'f' flag in their own
% documentation. this means we need to rewrite the whole file move.
% after 60 years of pcs you would think that copying files is 
% understood technology!
    targ = sprintf('%s%s/eidors_var_id.%s',HOMEDIR,archdir,mexext);
    try
    delete( targ );
    end
    movefile(sprintf('%s/eidors_var_id.%s',HOMEDIR, mexext), targ)

    ok = eidors_var_id_ok; % test it
    if ~ok
       fprintf([ ...
    'After compilation, eidors_var_id does not seem to be working.' ...
    'Sorry, you will need to debug this yourself. Some ideas are:\n\n' ...
    'On windows, try "mex -setup". You may need to install a compiler.' ...
    'For your matlab version (ie R2013a), see:' ...
    'http://www.mathworks.com/support/compilers/R2013a\n\n' ...
    'On linux, you may need to install older compilers like gcc-4.4.' ...
    'These can be used by writing\n' ...
    '   mex CC=gcc-4.4 CXX=g++-4.4 -largeArrayDims eidors_var_id.cpp\n']);
    end

function print_welcome(HOMEDIR,archdir,ver)
    eidors_ver = eidors_obj('eidors_version');
    if eidors_ver(end) == '+' % post release version
       % THIS IS HORRIBLE, HORRIBLE CRAP IN SVN. LOTS OF USERS WANT GlobalRev
       % BUT THE ARROGANT SVN AUTHORS REFUSE TO PROVIDE IT!!!!
       % If available, we could write it into each checkin without having
       % to rely on a broken system call.
       [status, result] = system('svnversion');
       if status==0;
          eidors_ver = [eidors_ver, ' SVN_ID=', result(1:end-1)];
       else
          eidors_ver = [eidors_ver, ' SVN_ID=unknown'];
       end
    end
    eidors_msg('Installed EIDORS (Ver: %s)', eidors_ver,1);

    eidors_msg('Parameter: cache_size=%.0f MB',eidors_cache('cache_size')/(1024*1024),1);
    eidors_msg('Parameter: mapped_colour=%d',calc_colours('mapped_colour'),1);
    if calc_colours('greylev')>=0
        eidors_msg('Default background colour: black',1);
    else
        eidors_msg('Default background colour: white',1);
    end
    eidors_msg('EIDORS mex folder: %s%s',HOMEDIR,archdir,1);
    eidors_msg('EIDORS cache folder: %s (must be writable)', ...
         eidors_cache('cache_path'),1);
    eidors_msg('EIDORS model cache: %s', mk_library_model('LIBRARY_PATH'),1);


    % helpful messages
    % TODO: test if desktop is available
    if ver.isoctave
        canwritehtml=0;
    else try
            mf = com.mathworks.mde.desk.MLDesktop.getInstance.getMainFrame;
            if isempty(mf)
                canwritehtml=0;
            else
                canwritehtml=1;
            end
        catch
            canwritehtml=0;
        end
    end
    if canwritehtml
        tutorials = '<a href="http://eidors3d.sf.net/tutorial/tutorial.shtml">Tutorials</a>';
    else
        tutorials = 'Tutorials';
    end
    eidors_msg(['New to EIDORS? Have a look at the ',tutorials,'.'],1);

function ok = eidors_var_id_ok;
    id0 = '';
    try id0 = eidors_var_id([]); end
    if strcmp(id0, ...
      'id_DA39A3EE5E6B4B0D3255BFEF95601890AFD80709')  % SHA1 of []
       ok = 1;
    elseif strcmp(id0, ...
      'id_99AA06D3014798D86001C324468D497F8C010DD2')  % xxHash of []
       ok = 1;
    else
       ok = 0;
    end
    if ok==0
       warning('caching (function eidors_var_id) is not working');
    else
       eidors_msg('tested function eidors_var_id: OK',1);
    end
