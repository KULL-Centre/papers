function [data, seq] = fastaread(filename,varargin)
%FASTAREAD reads FASTA format file.
%
%   S = FASTAREAD(FILENAME) reads a FASTA format file FILENAME, returning
%   the data in the file as a structure. FILENAME can also be a URL or
%   MATLAB character array that contains the text of a FASTA format file.
%   S.Header is the header information. S.Sequence is the sequence stored
%   as a string of characters.
%
%   [HEADER, SEQ] = FASTAREAD(FILENAME) reads the file into separate
%   variables HEADER and SEQ. If the file contains more than one sequence,
%   then HEADER and SEQ are cell arrays of header and sequence information.
%
%   FASTAREAD(...,'IGNOREGAPS',TF) removes any gap symbol ('-' or '.')
%   from the sequence(s) when TF is true. Default is false.
%
%   FASTAREAD(...,'BLOCKREAD', M) allows you to read in a single entry or
%   block of entries from a file containing multiple sequences. If M is a
%   scalar then the M'th entry in the file is read. If M is a two element
%   vector then the block of entries starting at entry M(1) and ending at
%   entry M(2) will be read.  Use Inf for M(2) to read all entries in the
%   file starting at position M(1). 
%
%   FASTAREAD(...,'TRIMHEADERS',TF) trims the header after the first
%   whitespace when TF is true. White space characters include a horizontal
%   tab (char(9)) and a space (char(32)). Default is false. 
%
%   FASTA format specified here:
%   http://www.ncbi.nlm.nih.gov/BLAST/fasta.shtml
%
%   Examples:
%
%       % Read the sequence for the human p53 tumor gene.
%       p53nt = fastaread('p53nt.txt')
%
%       % Read the sequence for the human p53 tumor protein.
%       p53aa = fastaread('p53aa.txt')
%
%       % Read a block of entries from a file
%       pf2_5_10 = fastaread('pf00002.fa','blockread',[ 5 10], ...
%                            'ignoregaps',true)
%
%   See also EMBLREAD, FASTAINFO, FASTAWRITE, FASTQINFO, FASTQREAD,
%   FASTQWRITE, GENBANKREAD, GENPEPTREAD, HMMPROFDEMO, MULTIALIGNREAD,
%   SEQALIGNVIEWER, SEQPROFILE, SEQVIEWER, SFFINFO, SFFREAD.

%   Copyright 2002-2012 The MathWorks, Inc.


% check input is char
% in a future version we may accept also cells
if ~ischar(filename)
    error(message('bioinfo:fastaread:InvalidInput'))
end

% default
ignoreGaps = false;
trimh = false;
blockRead = false;
% get input arguments
if  nargin > 1
    if rem(nargin,2) ~= 1
        error(message('bioinfo:fastaread:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'ignoregaps','blockread','trimheaders'};
    for j=1:2:nargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:fastaread:UnknownParameterName', pname));
            %elseif length(k)>1
            %    error('bioinfo:fastaread:AmbiguousParameterName',...
            %        'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % ignore gaps
                    ignoreGaps = bioinfoprivate.opttf(pval);
                    if isempty(ignoreGaps)
                        error(message('bioinfo:fastaread:IgnoreGapsNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 2  % range
                    range = pval;
                    if ~isnumeric(range) || numel(range)> 2 || isempty(range)
                        error(message('bioinfo:fastaread:BadBlockRange'))
                    end
                    blockRead = true;
                    range = sort(range);
                case 3 % trimheaders
			        trimh = bioinfoprivate.opttf(pval, okargs{k}, mfilename);    
            end
        end
    end
end


if size(filename,1)>1  % is padded string
    if blockRead
        warning(message('bioinfo:fastaread:IgnoredRange'))
    end
    ftext = cell(size(filename,1),1);
    for i=1:size(filename,1)
        ftext(i,1)=strtrim(strread(filename(i,:),'%s','whitespace','','delimiter','\n'));
    end
    % try then if it is an url
elseif (strfind(filename(1:min(10,end)), '://'))
    if blockRead
        warning(message('bioinfo:fastaread:IgnoredRange'))
    end
    if (~usejava('jvm'))
        error(message('bioinfo:fastaread:NoJava'))
    end
    try
        ftext = urlread(filename);
    catch allExceptions
        error(message('bioinfo:fastaread:CannotReadURL', filename));
    end
    ftext = strread(ftext,'%s','delimiter','\n');

    % try then if it is a valid filename
elseif  (exist(filename,'file') || exist(fullfile(pwd,filename),'file'))
    if blockRead
        blockText = getfileblock(filename,range,'>');
        try
           ftext = strread(blockText,'%s','delimiter','\n');
        catch theErr
           if strcmpi(theErr.identifier,'MATLAB:nomem')
                error(message('bioinfo:fastaread:BlockTooBig'));
           else
                 rethrow(theErr);
           end
        end
    else %=== read entire file
        fid = fopen(filename);
		c = onCleanup(@()fclose(fid));
        try
            ftext = textscan(fid,'%s','delimiter','\n');
            ftext = ftext{:}; 
        catch theErr
            if strcmpi(theErr.identifier,'MATLAB:nomem')
                error(message('bioinfo:fastaread:FileTooBig'));
            else 
                rethrow(theErr);
            end
        end
    end
else  % must be a string with '\n', convert to cell
    if blockRead
        warning(message('bioinfo:fastaread:IgnoredRange'))
    end
    ftext = strread(filename,'%s','delimiter','\n');
end

% it is possible that there will be multiple sequences
commentLines = strncmp(ftext,'>',1);

if ~any(commentLines)
    error(message('bioinfo:fastaread:FastaNotValid'))
end

numSeqs = sum(commentLines);
seqStarts = [find(commentLines); size(ftext,1)+1];
data(numSeqs,1).Header = '';

try
    for theSeq = 1:numSeqs
        % Check for > symbol ?
        data(theSeq).Header = ftext{seqStarts(theSeq)}(2:end);
        % convert 1x0 empty char array to ''; 
        if isempty(data(theSeq).Header) 
           data(theSeq).Header = ''; 
        end 

        firstRow = seqStarts(theSeq)+1;
        lastRow = seqStarts(theSeq+1)-1;
        numChars = cellfun('length',ftext(firstRow:lastRow));
        numSymbols = sum(numChars);
        data(theSeq).Sequence = repmat(' ',1,numSymbols);
        pos = 1;
        for i=firstRow:lastRow,
            str = strtrim(ftext{i});
            len =  length(str);
            if len == 0
                break
            end
            data(theSeq).Sequence(pos:pos+len-1) = str;
            pos = pos+len;
        end
        data(theSeq).Sequence = strtrim(data(theSeq).Sequence);
        if ignoreGaps
            data(theSeq).Sequence = strrep(data(theSeq).Sequence,'-','');
            data(theSeq).Sequence = strrep(data(theSeq).Sequence,'.','');
        end
    end

    % trim headers
    if trimh
       for i = 1:numSeqs
          data(i).Header = sscanf(data(i).Header,'%s',1);
       end
    end
    
    % in case of two outputs
    if nargout == 2
        if numSeqs == 1
            seq = data.Sequence;
            data = data.Header;
        else
            seq = {data(:).Sequence};
            data = {data(:).Header};
        end
    end

catch allExceptions
    error(message('bioinfo:fastaread:IncorrectDataFormat'))
end

