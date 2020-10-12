function xASL_system(command, source, destination, bForce)
% xASL function to copy or move files.
%
% FORMAT: RES = xASL_adm_GetFsList([strDirectory, strRegEx, bGetDirNames, bExcludeHidden, bIgnoreCase, nRequired])
%
% INPUT:
%   command        'move':   Move file from source to destination.
%                  'copy':   Copy file from source to destination.
%   source         Source location of the corresponding file.
%                  Please use forward slashes for your paths (should work on windows and unix.)
%   destination    destination of the corresponding file.
%                  Please use forward slashes for your paths (should work on windows and unix.)
%   bForce         Overwrite file if it already exists in the corresponding destination.
%   
%
% OUTPUT:          n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Wrapper script for the Matlab movefile and copyfile functions with built in overwrite option.
%
% EXAMPLE:     command = 'move';
%              source = '/home/folder_A/my_file.txt';
%              destination = '/home/folder_B/my_file.txt';
%              bForce = true;
%              xASL_system(command, source, destination, bForce);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2020 ExploreASL

    % Check if source exists
    if ~isfile(source) % xASL_exist does not seem to be able to do this with spaces
        error('Source file does not exist...');
    end

    % Check if destination exists
    destination_exists=isfile(destination);

    % Check if file does not exist already or if it should be overwritten
    if ~destination_exists || bForce        
        switch command
           case 'move'
              movefile(source,destination);
           case 'copy'
              copyfile(source,destination);
           otherwise
              warning('No known command...')
        end
    else
        warning('Trying to move/copy with existing destination and without force...')
    end
    