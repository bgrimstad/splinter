% This file is part of the SPLINTER library.
% Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.

function setup()
    % Change this to the directory where the MatLab interface of SPLINTER
    % is installed.
    splinter_path = '/home/anders/C++/SPLINTER/SPLINTER/build/splinter-matlab';
    %splinter_path = 'C:/Users/Anders/Documents/Github/SPLINTER/SPLINTER/build/splinter-matlab';

    % Add the directory containing the MatLab interface of SPLINTER to the
    % search path that MatLab searches through to find .m files.
    addpath(fullfile(splinter_path, 'matlab'));

    windows = ispc();
    mac = ismac();
    linux = isunix() && ~mac;

    % Detect bitness. Linux and MAC does not have 32 bit builds.
    wordSize = 64;
    if(strcmp('PCWIN', computer()))
        wordSize = 32;
    end

    % Header file is at the same location no matter the OS
    headerFile = fullfile(splinter_path, 'include', 'matlab.h');

    % Library file location is dependent on arch and OS
    wordSizeString = strcat(int2str(wordSize), 'bit');
    if(windows)
        libFileDir = fullfile(splinter_path, 'lib', 'Windows', wordSizeString);
        libFile = fullfile(libFileDir, 'splinter-matlab-1-4.dll');
    elseif(linux)
        libFileDir = fullfile(splinter_path, 'lib', 'Linux', wordSizeString);
        libFile = fullfile(libFileDir, 'libsplinter-matlab-1-4.so');
    elseif(mac)
        libFileDir = fullfile(splinter_path, 'lib', 'OSX', wordSizeString);
        libFile = fullfile(libFileDir, 'libsplinter-matlab-1-4.so');
    else
        libFileDir = fullfile(splinter_path, 'lib', 'Linux', wordSizeString);
        libFile = fullfile(libFileDir, 'libsplinter-matlab-1-4.so');
    end

    % The Splinter class is implemented as a Singleton
    s = Splinter.getInstance();
    s.load(libFile, headerFile);
end