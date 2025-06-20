function status= call_netgen(geo_file, vol_file, msz_file, finelevel)
% CALL_NETGEN: call netgen to create a vol_file from a geo_file
% status= call_netgen(geo_file, vol_file, msz_file, finelevel)
%  staus = 0 -> success , negative -> failure
%
% geo_file = geometry file (input)
% vol_file = FEM mesh file (output)
% msz_file = (deprecated) Meshsize file in netgen format
%
% Finelevel controls the fineness of the mesh
%   default is '' -> coarse
%   valid values are 'fine' or 'veryfine'
% 
% Finelevel can also be set using eidors_default
%   eidors_default call_netgen_finelevel 'fine'
%
% To debug the shapes and view the netgen interface, use
%   eidors_debug on call_netgen

% $Id: call_netgen.m 7112 2024-12-28 17:17:14Z aadler $
% (C) 2006 Andy Adler. Licensed under GPL V2

if nargin<3
   msz_file= '';
end

if ~isempty(msz_file)
   eidors_msg('call_netgen: Warning. Using an *.msz file. netgen rarely uses this file without care. Use the MSZPOINTS option in ng_write_opt instead.');
end

if nargin<4
   %  finelevel= '-veryfine';
   %  finelevel= '-fine';
   finelevel= eidors_default('get','call_netgen_finelevel');
end

% Netgen executable filename
cache_path = eidors_cache('cache_path');
eidors_path= eidors_cache('eidors_path');
if  isunix
   % If the file ng.sh has been created, call it instead
   % If we're calling netgen.exe, then use 
   % NETGEN=/path/to/netgen-5.3_x64/bin NETGENDIR=/path/to/netgen-5.3_x64/bin wine "/path/to/Netgen-5.3_x64/bin/netgen.exe" $1 $2 $3
   if isfile([cache_path,'/ng.sh'])
      ng_name = [cache_path,'/ng.sh'];
   elseif isdir([eidors_path,'/../Netgen-5.3_x64/bin'])
      % Check for a version of netgen shipped with eidors
      % Fixme if we ship a different version of netgen
      netgen_path = [eidors_path,'/../Netgen-5.3_x64/bin'];
      ng_name = [cache_path, '/ng.sh'];

      fprintf('Creating batch file to access Netgen. Ensure wine is installed\n')
      [status] = system('wine --version'); % test wine works
      if status~=0;
         error('<wine> must be installed to use provided netgen.exe');
      end
      fid= fopen(ng_name,'w');
      if fid<0; error('Unable to write to %s',cache_path); end
      fprintf(fid,['NETGEN="',netgen_path,'" NETGENDIR="',netgen_path, ...
                   '" wine "',netgen_path,'/netgen.exe" $1 $2 $3\n']);
      fclose(fid);
      system(['chmod a+x "',ng_name,'"']); % make executable
   else
      ng_name = 'netgen';
   end
else % windows
   % If the file ng.bat has been created, call it instead
   % If we're calling netgen.exe, then use 
   % NETGEN=/path/to/netgen-5.3_x64/bin NETGENDIR=/path/to/netgen-5.3_x64/bin wine "/path/to/Netgen-5.3_x64/bin/netgen.exe" $1 $2 $3
   if isfile([cache_path,'/ng.bat'])
      ng_name = [cache_path,'/ng.bat'];
   elseif isdir([eidors_path,'/../Netgen-5.3_x64/bin'])
      % Check for a version of netgen shipped with eidors
      % Fixme if we ship a different version of netgen
      netgen_path = [eidors_path,'/../Netgen-5.3_x64/bin'];
      ng_name = [cache_path, '/ng.bat'];
      fprintf('Creating batch file to access Netgen. Ensure wine is installed\n')
      fid= fopen(ng_name,'w');
      if fid<0; error('Unable to write to %s',cache_path); end
      fprintf(fid,'set NETGEN=%s\n', netgen_path);
      fprintf(fid,'set NETGENDIR=%s\n', netgen_path);
      fprintf(fid,'"%s/netgen.exe" %%*\n', netgen_path);
      fclose(fid);
   else
      ng_name = 'netgen.exe';
   end
end
    
if ~isempty(msz_file)
    warning('EIDORS:Deprecated', 'Use ng_write_opt instead of specifying a meshsize file');
    warning('Overwriting ng.opt file');
    fid= fopen('ng.opt','w'); %create ng.opt file in local dir
    if fid==-1
        error(['Netgen requires writing files in the current directory(%s). ', ...
            'Unfortunately, you don''t have permission. ' ...
            'Your options are: 1) change your working directory to one in which you have write permission, or ' ...
            '2) change the permissions on the current working directory.'], pwd);
    end
    %     fprintf(fid,'options.segmentsperedge 5\n'); % Another
    %                                                   potentially useful parameter
    %                                                   except netgen ignores it
    fprintf(fid,'options.meshsizefilename= %s\n',msz_file);
    fclose(fid);
end

while( 1 )
   if ~isunix
      % on Linux, Netgen runs in the current directory
      % enforce this behaviour in Windows
      oldpath = getenv('NETGEN_USER_DIR');
      setenv('NETGEN_USER_DIR', cd);
   end

   if eidors_debug('query','call_netgen')
      sys_cmd = sprintf('"%s" %s  -geofile=%s  -meshfile=%s ', ...
         ng_name, finelevel,geo_file,vol_file);
   else
      sys_cmd = sprintf('"%s" %s -batchmode -geofile="%s"  -meshfile="%s" ', ...
         ng_name, finelevel,geo_file,vol_file);
   end
   status= system_cmd( sys_cmd );

   if status==0; break; end % 
   try
      if isunix
         global eidors_objects;
         ignore_errors = 0;
         try
            ignore_errors = eidors_objects.ng_ignore_errors;
         end
         fprintf(['It seems that netgen has not worked (on linux). ' ...
            'If you downloaded a eidors-v??-ng, check <wine> is installed. ' ...
            'Otherwise, Check that netgen is installed and on the path. ' ...
            'Check the environment variable NETGENDIR is set. ' ...
            ' This can be set via the following commands:\n' ...
            '  setenv(''NETGENDIR'',''/path/to/netgen/bin'');\n' ...
            '  setenv(''PATH'',[''/path/to/netgen/bin:'',getenv(''PATH'')]);\n' ...
            'Please enter a new netgen file name\n' ]);
         if ignore_errors == 1;
            eidors_msg('Ignoring Netgen Error (as requested)',1);
            break;
         end
         ng_name = input( ...
           'netgen file name (with path)? [or i=ignore, e=error, a=always ignore]','s');
         if strcmp(ng_name,'i'); break;end
         if strcmp(ng_name,'a');
             eidors_objects.ng_ignore_errors = 1;
             break;
         end
         if strcmp(ng_name,'e'); error('user requested'),end;
      else % Windows
         % Check for a version of netgen shipped with eidors
         % Fixme if we ship a different version of netgen
         netgen_path = [eidors_path,'/../Netgen-5.3_x64/bin'];
         if exist(netgen_path)
            fprintf('Creating batch file to access Netgen.\n')
            fid= fopen([cache_path, '/ng.bat'],'w');
            if fid<0; error('Unable to write to %s',cache_path); end
            fprintf(fid,'set TCL_LIBRARY=%s/lib/tcl8.3\n', netgen_path); % REQ for ng <= 4.4
            fprintf(fid,'set TIX_LIBRARY=%s/lib/tix8.1\n', netgen_path); % REQ for ng <= 4.4
            fprintf(fid,'set NETGENDIR=%s\n', netgen_path); % REQ for ng >= 4.9
            fprintf(fid,'"%s/netgen.exe" %%*\n', netgen_path);
            fclose(fid);
         else
            fprintf([ ...
               'Netgen call failed. Is netgen installed and on the search path?\n' ...
               'If you are running under windows, I can attempt to create\n' ...
               'a batch file to access netgen.\n' ...
               'Please enter the directory in which to find netgen.exe.\n' ...
               'A typical path is "C:\\Program Files (x86)\\Netgen-5.0_x64\\bin"\n' ...
               'If you don''t have a copy, download it from' ...
               'http://sourceforge.net/projects/netgen-mesher/ \n\n']);
            netgen_path = input('netgen_path? [or i=ignore, e=error] ','s');
            if strcmp(netgen_path,'i'); break;end
            if strcmp(netgen_path,'e'); error('user requested'),end;
            if exist( sprintf('%s/netgen.exe',netgen_path) , 'file' ) || ...
                  exist( sprintf('%s/bin/netgen.exe',netgen_path) , 'file' )
               disp('Found netgen version 4.4 or higher');
               netgen_exe = netgen_path;
               if exist( sprintf('%s/bin/netgen.exe',netgen_path) , 'file' )
                  netgen_exe = [netgen_path '/bin'];
               end
               
               fid= fopen([cache_path, '/ng.bat'],'w');
               if fid<0; error('Unable to write to %s',cache_path); end
               fprintf(fid,'set TCL_LIBRARY=%s/lib/tcl8.3\n', netgen_path); % REQ for ng <= 4.4
               fprintf(fid,'set TIX_LIBRARY=%s/lib/tix8.1\n', netgen_path); % REQ for ng <= 4.4
               fprintf(fid,'set NETGENDIR=%s\n', netgen_path); % REQ for ng >= 4.9
               fprintf(fid,'"%s/netgen.exe" %%*\n', netgen_exe);
               fclose(fid);
            elseif exist( sprintf('%s/ng431.exe',netgen_path) , 'file' )
               disp('Found netgen version 4.3.1');
               fid= fopen([cache_path, '/ng.bat'],'w');
               if fid<0; error('Unable to write to %s',cache_path); end
               fprintf(fid,'set TCL_LIBRARY=%s/lib/tcl8.3\n', netgen_path);
               fprintf(fid,'set TIX_LIBRARY=%s/lib/tcl8.2\n', netgen_path);
               fprintf(fid,'"%s/ng431.exe" %%*\n', netgen_path);
               fclose(fid);
            else
               warning(['cannot find a version of netgen that I know about\n' ...
                  'Install netgen or check the path\n']);
            end
         end
      end
   catch e
      if strncmp(computer,'PC',2)
         % restore Netgen settings on Windows
         setenv('NETGEN_USER_DIR', oldpath);
      end
      rethrow(e)
   end
end
