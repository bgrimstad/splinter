function setup()
    % Change this to the directory where the MatLab interface of SPLINTER
    % is installed.
    splinter_path = '/home/anders/splinter-matlab';

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
        libFile = fullfile(libFileDir, 'splinter-matlab-1-2.dll');
    elseif(linux)
        libFileDir = fullfile(splinter_path, 'lib', 'Linux', wordSizeString);
        libFile = fullfile(libFileDir, 'libsplinter-matlab-1-2.so');
    elseif(mac)
        libFileDir = fullfile(splinter_path, 'lib', 'OSX', wordSizeString);
        libFile = fullfile(libFileDir, 'libsplinter-matlab-1-2.so');
    else
        libFileDir = fullfile(splinter_path, 'lib', 'Linux', wordSizeString);
        libFile = fullfile(libFileDir, 'libsplinter-matlab-1-2.so');
    end

    % The Splinter class is implemented as a Singleton
    s = Splinter.getInstance();
    s.load(libFile, headerFile);
end