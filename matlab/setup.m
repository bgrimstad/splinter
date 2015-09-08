% This file is part of the SPLINTER library.
% Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.

function setup()
    % Change this to the directory where the MatLab interface of SPLINTER
    % is installed.
    %splinter_path = '/home/bjarne/Code/C++/splinter2/splinter/build/splinter-matlab';
    currentDirectory = pwd;
    [upperPath, deepestFolder, ~] = fileparts(currentDirectory) ;
    splinter_path = upperPath;
    
    % Read version file. Name of library file depends on the version.
    versionFile = fullfile(splinter_path, 'version');
    versionFileId = fopen(versionFile, 'r');
    version = fscanf(versionFileId, '%d-%d');
    fclose(versionFileId);
    
    majorVersion = version(1);
    minorVersion = version(2);

    % Add the directory containing the MatLab interface of SPLINTER to the
    % search path that MatLab searches through to find .m files.
    addpath(fullfile(splinter_path, 'matlab'));

    windows = ispc();
    mac = ismac();
    linux = isunix() && ~mac;

    % Detect architecture. Linux and MAC does not have x86 builds.
    arch = 'x86-64';
    if(strcmp('PCWIN', computer()))
        arch = 'x86';
    end

    % Header file is at the same location no matter the OS
    headerFile = fullfile(splinter_path, 'include', 'matlab.h');

    libBaseName = strcat('splinter-matlab-', int2str(majorVersion));
    libBaseName = strcat(libBaseName, '-');
    libBaseName = strcat(libBaseName, int2str(minorVersion));
    
    if(windows)
        libFileDir = fullfile(splinter_path, 'lib', 'windows', arch);
        libFile = fullfile(libFileDir, strcat(libBaseName, '.dll'));
    elseif(linux)
        libFileDir = fullfile(splinter_path, 'lib', 'linux', arch);
        libFile = fullfile(libFileDir, strcat('lib', strcat(libBaseName), '.so'));
    elseif(mac)
        libFileDir = fullfile(splinter_path, 'lib', 'osx', arch);
        libFile = fullfile(libFileDir, strcat('lib', strcat(libBaseName), '.so'));
    else
        libFileDir = fullfile(splinter_path, 'lib', 'linux', arch);
        libFile = fullfile(libFileDir, strcat('lib', strcat(libBaseName), '.so'));
    end

    % The Splinter class is implemented as a Singleton
    s = Splinter.getInstance();
    s.load(libFile, headerFile);
end
