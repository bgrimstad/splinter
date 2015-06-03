Splinter_alias = '';
loaded = false;

function obj = Splinter(objectFile, header, alias)
    obj.Splinter_alias = alias;

    loadlibrary(objectFile, header, 'alias', alias)
end

