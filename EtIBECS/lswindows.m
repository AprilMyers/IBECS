function varargout=lswindows(varargin)
%LS List folder contents.
%   LS displays the results of the 'ls' command on UNIX. On UNIX, 
%   LS returns a character row vector of filenames separated 
%   by tab and space characters. On Windows, LS returns an m-by-n 
%   character array of filenames, where m is the number of filenames 
%   and n is the number of characters in the longest filename found. 
%   Filenames shorter than n characters are padded with space characters.
%
%   You can pass any flags to LS as well that your operating system supports.
%
%   See also DIR, MKDIR, RMDIR, FILEATTRIB, COPYFILE, MOVEFILE, DELETE.

%   Copyright 1984-2019 The MathWorks, Inc.
%=============================================================================
% validate input parameters
if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
    
    if ~iscellstr(varargin)
        error(message('MATLAB:ls:InputsMustBeStrings'));
    end

    try
        isURL = any(matlab.virtualfileio.internal.validators.isIRI(varargin)); 
    catch
        isURL = false;
    end
    
    if isURL
        error(message('MATLAB:ls:URLNotAllowed'));
    end
end

% check output arguments
if nargout > 1
    error(message('MATLAB:ls:TooManyOutputArguments'));
end


    if nargin == 0
        % Display output of dir in wide format.  
        % dir; prints out info.  
        % d = dir; does not.
        if nargout == 0
            dir;
        else
            d = dir;
            listing = char(d.name);
        end
    elseif nargin == 1
        if nargout == 0
            dir(varargin{1});
        else
            d = dir(varargin{1});
            listing = char(d.name);
        end
    else
        error(message('MATLAB:ls:TooManyInputArguments'));
    end
%end

% determine output mode, depending on presence of output arguments
i%f nargout == 0 && isunix
    disp(listing)
%else
if nargout > 0
    varargout{1} = listing;
end

%---------------------------------------------------------------------------
function quotedArgs = quoteUnixCmdArg(tildeArgs)
% Algorithm: Start and end each argument with a single quote (squote).
%            Within each argument:
%            1. squote -> squote '\' squote squote
%            2. '!'    -> squote '\' '!' squotex
%            3. '*'    -> squote '*' squote	(MATLAB globbing character)
%
if isempty(tildeArgs)
    quotedArgs = '';
    return
end
% Do any tilde expansion first 
ix = find(strncmp(tildeArgs,'~',1)); 
if ~isempty(ix) 
  tildeArgs(ix) = unix_tilde_expansion(tildeArgs(ix)); 
end 

% Special cases to maintain as literal: single quote or ! with '\thing_I_found'
quotedArgs= regexprep(tildeArgs,'[''!]','''\\$&''');

% Special cases to maintain as NOT literal: Replace * with 'thing_I_found'
quotedArgs= regexprep(quotedArgs,'[*]','''$&''');

quotedArgs = strcat(' ''', quotedArgs, '''');
quotedArgs = [quotedArgs{:}];
